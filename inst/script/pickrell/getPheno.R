setwd("..")
library("GEOquery")
gse <- getGEO("GSE19480")[[1]]
pDataGSE <- pData(gse)[,c("title","supplementary_file_1")]
gseSamples <- data.frame(id=sub("(.*)","\\1",pDataGSE$title),
                      experiment=sub(".*\\/(SRX.*)","\\1",pDataGSE$supplementary_file_1))
gseSamples$hapmap_id <- sub("(NA[0-9]*)_.*","\\1",gseSamples$id)
hapmap <- read.csv("HapMap_samples.csv",strings=FALSE)
names(hapmap)[1] <- "hapmap_id"
samples <- merge(gseSamples,hapmap,by="hapmap_id")
write.table(samples,"geo_samples.txt")

library("SRAdb")
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
conversion <- sraConvert(in_acc = samples$sra, out_type =
                         c("sra","submission","study","sample","experiment","run"),
                         sra_con = sra_con)
write.table(conversion,file="sra_conversion.txt")

samples <- read.table("geo_samples.txt")
conversion <- read.table("sra_conversion.txt")
readlength <- read.csv("readlength.txt",header=FALSE)
names(readlength) <- c("run","readlength")
samplesMerged <- merge(merge(samples, conversion),readlength)

goodSet <- samplesMerged[samplesMerged$readlength == 46 & samplesMerged$sex == 1,]
goodSet <- goodSet[-grep("_2_argonne",as.character(goodSet$id)),]
goodSet <- goodSet[order(goodSet$run),]
write.table(goodSet,"good_set.txt")

goodSet <- read.table("good_set.txt")
goodSet$run <- as.character(goodSet$run)
set.seed(5)
randomSubsets <- t(replicate(30,goodSet$run[sample(nrow(goodSet),10,replace=FALSE)]))
write.table(randomSubsets,file="random_subsets.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
