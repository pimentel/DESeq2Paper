setwd("..")
samples <- read.table("bottomly_phenodata.txt",header=TRUE)
rownames(samples) <- samples$sample.id
samples$experiment <- samples$sample.id

library("SRAdb")
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
conversion <- sraConvert(in_acc = samples$sample.id, out_type =
                         c("sra","submission","study","sample","experiment","run"),
                         sra_con = sra_con)
write.table(conversion,file="sra_conversion.txt")

conversion <- read.table("sra_conversion.txt")
samplesMerged <- merge(samples, conversion)
samplesMerged <- samplesMerged[order(samplesMerged$strain),]
samplesMerged$run <- as.character(samplesMerged$run)
rownames(samplesMerged) <- samplesMerged$run
  
library(caret)
c57Idx <- which(samplesMerged$strain == "C57BL/6J")
dbaIdx <- which(samplesMerged$strain == "DBA/2J")

printIdx <- function() {
  c57Sub <- createDataPartition(samplesMerged$experiment.number[c57Idx], p=2/10)[[1]]
  dbaSub <- createDataPartition(samplesMerged$experiment.number[dbaIdx], p=2/11)[[1]]
  idx <- c(c57Idx[c57Sub], dbaIdx[dbaSub])
  c(samplesMerged$run[idx], samplesMerged$run[-idx])
}

numIdx <- function() {
    c57Sub <- createDataPartition(samplesMerged$experiment.number[c57Idx], p=2/10)[[1]]
    dbaSub <- createDataPartition(samplesMerged$experiment.number[dbaIdx], p=2/11)[[1]]
    c(c57Idx[c57Sub], dbaIdx[dbaSub])
}

set.seed(5)
randomSubsets <- t(replicate(30, printIdx()))

set.seed(5)
numSubsets <- t(replicate(30, numIdx()))
all(dist(numSubsets) > 0)

write.table(randomSubsets,file="random_subsets.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

as.character(samplesMerged[randomSubsets[1,],"strain"])
write.table(as.character(samplesMerged[randomSubsets[1,],"strain"]),
            file="random_strains.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
