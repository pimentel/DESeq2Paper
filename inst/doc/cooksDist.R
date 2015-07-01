
## ----lib-----------------------------------------------------------------
library("DESeq2")
library("DESeq2paper")


## ----cooksSetup, cache=TRUE----------------------------------------------
data("bottomly_sumexp")
dds <- DESeqDataSetFromMatrix(assay(bottomly), DataFrame(colData(bottomly)), ~ strain)
dds <- dds[,c(8:11,15:17,12:14,18:21)]
as.data.frame(colData(dds))
dds <- DESeq(dds)


## ----cooksNoReplace, cache=TRUE------------------------------------------
ddsNoReplace <- DESeq(dds, minReplicatesForReplace=Inf)


## ----cooksFiltering------------------------------------------------------
names(assays(dds))
maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
idx <- which(rownames(dds) == "ENSMUSG00000076609")
unname(counts(dds)[idx,])


## ----cooksDist, dev="pdf", fig.align="center", fig.width=5.5, fig.height=2, fig.cap="Demonstration using the Bottomly et al dataset, of detection of outlier counts (A) using Cook's distances (B), and refitting after outliers have been replaced by the trimmed median over all samples (C). The dotted line indicates the fitted mean on the common scale."----
makeColors <- function(y=c(-1e6,1e6)) {
  polygon(c(-1,-1,7.5,7.5),c(y,y[2:1]),
          col=rgb(0,1,0,.1),border=NA)
  polygon(c(7.5,7.5,50,50),c(y,y[2:1]),
          col=rgb(0,0,1,.1),border=NA)
}
line <- 0.6
adj <- -.5
cex <- 1

par(mfrow=c(1,3),mar=c(4.3,4.3,3,1))
out <- assays(dds)[["cooks"]][idx,] > qf(.99, 2, ncol(dds) - 2)
plot(counts(dds,normalized=TRUE)[idx,],main="With outlier",
     ylab="normalized counts",xlab="samples",
     pch=as.integer(colData(dds)$strain) + 1,
     ylim=c(0,max(counts(dds,normalized=TRUE)[idx,])),
     col=ifelse(out,"red","black"))
makeColors()
q0 <- 2^(mcols(ddsNoReplace)$Intercept[idx] + mcols(ddsNoReplace)$strainC57BL.6J[idx])
q1 <- 2^(mcols(ddsNoReplace)$Intercept[idx] + mcols(ddsNoReplace)$strainDBA.2J[idx])
segments(1,q0,7,q0,lty=3)
segments(8,q1,14,q1,lty=3)
mtext("A",side=3,line=line,adj=adj,cex=cex)

plot(assays(dds)[["cooks"]][idx,],
     main="Cook's distances",
     ylab="",xlab="samples",
     log="y",
     pch=as.integer(colData(dds)$strain) + 1,
     col=ifelse(out,"red","black"))
makeColors(y=c(1e-5,1e5))
abline(h=qf(.99, 2, ncol(dds) - 2))
mtext("B",side=3,line=line,adj=adj,cex=cex)

plot(assays(dds)[["replaceCounts"]][idx,]/sizeFactors(dds),
     main="Outlier replaced",
     ylab="normalized counts",xlab="samples",
     ylim=c(0,max(assays(dds)[["replaceCounts"]][idx,]/sizeFactors(dds))),
     pch=as.integer(colData(dds)$strain) + 1,
     col=ifelse(out,"red","black"))
makeColors()
q0 <- 2^(mcols(dds)$Intercept[idx] + mcols(dds)$strainC57BL.6J[idx])
q1 <- 2^(mcols(dds)$Intercept[idx] + mcols(dds)$strainDBA.2J[idx])
segments(1,q0,7,q0,lty=3)
segments(8,q1,14,q1,lty=3)
mtext("C",side=3,line=line,adj=adj,cex=cex)



## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


