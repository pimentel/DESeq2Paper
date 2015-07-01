library("GenomicRanges")
library("Biobase")
library("BiocParallel")
load("../../../data/bottomly_sumexp.RData")
randomSubsets <- read.table("random_subsets.txt",strings=FALSE)
se <- bottomly[,match(randomSubsets[1,],colnames(bottomly))]
colData(se)$run <- colnames(se)
eset <- ExpressionSet(assay(se),
                      AnnotatedDataFrame(as.data.frame(colData(se))))
pData(eset)$condition <- pData(eset)$strain
levels(pData(eset)$condition) <- c("A","B")

library("DESeq")
library("DESeq2")
library("edgeR")
library("limma")
library("samr")
library("DSS")
library("EBSeq")
source("../runScripts.R")

algos <- list("DESeq"=runDESeq,"DESeq2"=runDESeq2,"edgeR"=runEdgeR,"edgeR-robust"=runEdgeRRobust,
              "DSS"=runDSS,"voom"=runVoom,"SAMseq"=runSAMseq,"EBSeq"=runEBSeq)

namesAlgos <- names(algos)
names(namesAlgos) <- namesAlgos

resTest <- list()
resHeldout <- list()
lfcTest <- list()
lfcHeldout <- list()

nreps <- 30

register(MulticoreParam(workers=5,verbose=TRUE))

res <- bplapply(1:nreps, function(i) {   
    cat(i," ")
    testSet <- as.character(randomSubsets[i,1:6])
    heldOutSet <- as.character(randomSubsets[i,-(1:6)])
    eTest <- eset[,testSet]
    eHeldout <- eset[,heldOutSet]
    resTest0 <- lapply(namesAlgos, function(n) algos[[n]](eTest))
    resHeldout0 <- lapply(namesAlgos, function(n) algos[[n]](eHeldout))
    cuffDir <- "/g/huber/projects/bottomly/cuffdiff_results/"
    cuffResTest <- read.table(paste0(cuffDir,"cuffdiff_random_testset_",i,"/gene_exp.diff"),header=TRUE)
    rownames(cuffResTest) <- cuffResTest$gene_id
    cuffResHeldout <- read.table(paste0(cuffDir,"cuffdiff_random_heldoutset_",i,"/gene_exp.diff"),header=TRUE)
    rownames(cuffResHeldout) <- cuffResHeldout$gene_id
    common <- intersect(rownames(eset),rownames(cuffResTest))
    resIdx <- match(common,rownames(eset))
    # cuffdiff2 q-values are equal to BH adjusted p-values over genes with status = "OK"
    cuffResTestCommon <- cuffResTest[match(common,rownames(cuffResTest)),]
    cuffResTestCommon$padj <- cuffResTestCommon$q_value
    cuffResHeldoutCommon <- cuffResHeldout[match(common,rownames(cuffResHeldout)),]
    cuffResHeldoutCommon$padj <- cuffResHeldoutCommon$q_value
    stopifnot(all(cuffResTest$gene_id == cuffResHeldout$gene_id))
    resTest <- as.data.frame(c(lapply(resTest0, function(z) z$padj[resIdx]),
                               list("cuffdiff2" = cuffResTestCommon$padj)))
    resHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$padj[resIdx]),
                                  list("cuffdiff2" = cuffResHeldoutCommon$padj)))
    lfcTest <- as.data.frame(c(lapply(resTest0, function(z) z$beta[resIdx]),
                               list("cuffdiff2" = cuffResTestCommon$log2.fold_change.)))
    lfcHeldout <- as.data.frame(c(lapply(resHeldout0, function(z) z$beta[resIdx]),
                                  list("cuffdiff2" = cuffResHeldoutCommon$log2.fold_change.)))
    rownames(resTest) <- common
    rownames(resHeldout) <- common
    rownames(lfcTest) <- common
    rownames(lfcHeldout) <- common
    list(resTest=resTest,resHeldout=resHeldout,lfcTest=lfcTest,lfcHeldout=lfcHeldout)
})

resTest <- lapply(res, "[[", "resTest")
resHeldout <- lapply(res, "[[", "resHeldout")
lfcTest <- lapply(res, "[[", "lfcTest")
lfcHeldout <- lapply(res, "[[", "lfcHeldout")

namesAlgos <- c(names(algos),"cuffdiff2")
names(namesAlgos) <- namesAlgos

save(resTest,resHeldout,lfcTest,lfcHeldout,namesAlgos,file="sensitivityPrecision.RData")

