setwd("..")
library("GenomicFeatures")
genesGTF <- "/tophatindices/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf"
hse <- makeTranscriptDbFromGFF(genesGTF,format="gtf",dataSource="Illumina iGenomes: Ensembl, NCBIM37, release 66, downloaded March 2013")
saveDb(hse,file="hseDb.RData")
exonsByGene <- exonsBy(hse, by="gene")

options(mc.cores=8)
library("Rsamtools")
goodSet <- list.files("tophat_out")
fls <- paste0("tophat_out/",goodSet,"/accepted_hits.bam")
names(fls) <- goodSet
bamLst <- BamFileList(fls,yieldSize=100000)
bottomly <- summarizeOverlaps(exonsByGene, bamLst,
                              mode="Union",
                              singleEnd=TRUE,
                              ignore.strand=TRUE,
                              fragments=TRUE)
load("samples_merged.RData")
idx <- match(colnames(bottomly),rownames(samplesMerged))
colData(bottomly) <- DataFrame(samplesMerged[idx,])
library("annotate")
exptData(bottomly) <- list(pmid2MIAME("21455293"))
save(bottomly, file="bottomly_sumexp.RData")


