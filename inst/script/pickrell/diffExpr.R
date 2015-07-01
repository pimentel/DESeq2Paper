library("GenomicRanges")
library("Biobase")
library("BiocParallel")
load("../../../data/pickrell_sumexp.RData")
randomSubsets <- read.table("random_subsets.txt",strings=FALSE)
se <- pickrell
colData(se)$run <- colnames(se)
eset <- ExpressionSet(assay(se),
                      AnnotatedDataFrame(as.data.frame(colData(se))))

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

resList <- list()

nreps <- 30

register(MulticoreParam(workers=8,verbose=TRUE))

resList <- bplapply(1:nreps, function(i) {
  cat(i," ")
  testSet <- as.character(randomSubsets[i,])
  eTest <- eset[,testSet]
  pData(eTest)$condition <- factor(rep(c("A","B"),each=ncol(eTest)/2))
  resTest0 <- lapply(namesAlgos, function(n) algos[[n]](eTest))
  cuffDir <- "/g/huber/projects/pickrell/cuffdiff_results/"
  cuffResTest <- read.table(paste0(cuffDir,"cuffdiff_random_",i,"/gene_exp.diff"),header=TRUE)
  rownames(cuffResTest) <- cuffResTest$gene_id
  common <- intersect(rownames(eset),rownames(cuffResTest))
  cuffResTestCommon <- cuffResTest[match(common,rownames(cuffResTest)),]
  resTest <- as.data.frame(c(lapply(resTest0, function(z) z$pvals[match(common,rownames(eset))]),
                             list("cuffdiff2" = cuffResTestCommon$p_value)))
  rownames(resTest) <- common
  resTest
})

namesAlgos <- c(names(algos),"cuffdiff2")
names(namesAlgos) <- namesAlgos

save(resList, namesAlgos, file="specificity.RData")

