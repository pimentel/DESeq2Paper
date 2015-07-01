
## ----lib-----------------------------------------------------------------
library("DESeq2")
library("DESeq2paper")


## ----lfcThresholdSetup, cache=TRUE---------------------------------------
data("bottomly_sumexp")
dds <- DESeqDataSetFromMatrix(assay(bottomly), DataFrame(colData(bottomly)), ~ strain)
dds <- DESeq(dds)
ddsNoPrior <- DESeq(dds, betaPrior=FALSE)
theta <- 1
res <- results(dds,lfcThreshold=theta,altHypothesis="greaterAbs")
padj <- res$padj
padj[is.na(padj)] <- 1
resNoPrior <- results(ddsNoPrior,lfcThreshold=theta,altHypothesis="lessAbs")
padjAB <- resNoPrior$padj
padjAB[is.na(padjAB)] <- 1


## ----lfcThreshold, dev="png", fig.align="center", fig.width=7, fig.height=3.5, dpi=150, fig.cap="Testing for logarithmic fold changes above (A) and below (B) a positive threshold, using the Bottomly et al dataset. Testing for logarithmic fold changes below a threshold (B) requires that a prior on logarithmic fold changes not be used."----

line <- 0.75
adj <- -.35
cex <- 1.5

par(mfrow=c(1,2),mar=c(4.5,4.5,3,1))
ymax <- 3

plotMA(res, ylim=c(-ymax, ymax),
       colNonSig=rgb(0,0,0,.5),
       colSig=rgb(1,0,0,.5),
       main=expression(H[A]:~abs(beta) > 1),
       colLine=NULL, ylab=expression(log[2]~fold~change))
abline(h=c(-1,1)*theta, col="dodgerblue",lty=3,lwd=3)
legend("bottomright","adj. p < .1",pch=16,col="red",bg="white",cex=.8)
mtext("A",side=3,line=line,adj=adj,cex=cex)

plotMA(resNoPrior, ylim=c(-ymax,ymax),
       colNonSig=rgb(0,0,0,.5),
       colSig=rgb(1,0,0,.5),
       main=expression(H[A]:~abs(beta) < 1),
       colLine=NULL, ylab=expression(log[2]~fold~change))
abline(h=c(-1,1)*theta, col="dodgerblue",lty=3,lwd=3)
legend("bottomright","adj. p < .1",pch=16,col="red",bg="white",cex=.8)
mtext("B",side=3,line=line,adj=adj,cex=cex)


## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


