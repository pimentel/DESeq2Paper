setwd("..")
library("GenomicFeatures")
genesGTF <- "/tophatindices/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
hse <- makeTranscriptDbFromGFF(genesGTF,format="gtf",dataSource="Illumina iGenomes: Ensembl, GRCh37, release 70, downloaded March 2013")
saveDb(hse,file="hseDb.RData")
exonsByGene <- exonsBy(hse, by="gene")

options(mc.cores=8)
library("Rsamtools")
goodSet <- read.table("good_set.txt")
fls <- paste0("tophat_out/",goodSet$run,"/accepted_hits.bam")
names(fls) <- goodSet$run
bamLst <- BamFileList(fls,yieldSize=100000)
pickrell <- summarizeOverlaps(exonsByGene, bamLst,
                              mode="Union",
                              singleEnd=TRUE,
                              ignore.strand=TRUE,
                              fragments=TRUE)

library("annotate")
exptData(pickrell) <- list(pmid2MIAME("20220758"))
save(pickrell, file="pickrell_sumexp.RData")
