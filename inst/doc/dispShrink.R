
## ----lib-----------------------------------------------------------------
library("DESeq2")
library("DESeq2paper")


## ----setupDispShrink, cache=TRUE-----------------------------------------
data("pickrell_sumexp")
idx <- 1:5
ddsP <- DESeqDataSetFromMatrix(assay(pickrell)[,idx], 
                               DataFrame(colData(pickrell)[idx,]), ~ 1)
ddsP <- estimateSizeFactors(ddsP)
ddsP <- estimateDispersions(ddsP)

data("bottomly_sumexp")
idx <- c(1:3,5:7)
ddsB <- DESeqDataSetFromMatrix(assay(bottomly)[,idx], 
                               DataFrame(colData(bottomly)[idx,]), ~ strain)
ddsB <- estimateSizeFactors(ddsB)
ddsB <- estimateDispersions(ddsB)



## ----pick----------------------------------------------------------------
plotDispShrink <- function(dds) {
  # pick 40 equally spaced genes along the base mean
  bins <- 10^seq(from=0, to=5, length=20)
  pickone <- function(x) {
    if (sum(x) == 0) return(NULL)
    if (sum(x) == 1) return(which(x))
    sample(which(x), 1)
  }
  up <- sapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst > 1e-4 & !mcols(dds)$dispOutlier & mcols(dds)$dispGeneEst > mcols(dds)$dispFit & mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean < bins[i+1]))
  down <- sapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst > 1e-4 & !mcols(dds)$dispOutlier & mcols(dds)$dispGeneEst < mcols(dds)$dispFit & mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean < bins[i+1]))
  # pick 5 outliers
  bins <- 10^seq(from=1, to=4, length=6)
  outliers <- do.call(c,lapply(seq_along(bins[-1]), function(i) pickone(mcols(dds)$dispGeneEst / mcols(dds)$dispMAP > 2 & mcols(dds)$dispOutlier & mcols(dds)$baseMean > bins[i] & mcols(dds)$baseMean < bins[i+1])))
  s <- c(up, down, outliers)
  s <- s[!is.na(s)]
  with(mcols(dds[s,]),
       plot(baseMean, dispGeneEst,log="xy",
            pch=16,xlab="mean of normalized counts",
            ylab="dispersion estimate",yaxt="n",ylim=c(.001,100)))
  axis(2,at=10^(-3:2),label=10^(-3:2))
  xs <- 10^(-20:50/10)
  lines(xs,dispersionFunction(dds)(xs),col="red",lwd=2)
  with(mcols(dds[s,][!mcols(dds[s,])$dispOutlier,]),
       arrows(baseMean, dispGeneEst, baseMean, dispersion,
              length=.075, col="dodgerblue",lwd=2))
  with(mcols(dds[s,][mcols(dds[s,])$dispOutlier,]),
       segments(baseMean, dispGeneEst, baseMean, dispMAP,
                col="dodgerblue",lwd=2, lty=3))
  with(mcols(dds[s,][mcols(dds[s,])$dispOutlier,]),
       points(baseMean, dispersion,
              cex=2,col="dodgerblue",lwd=2))
  legend("topright",c("MLE","prior mean","MAP"),pch=20,
         col=c("black","red","dodgerblue"),bg="white")
}



## ----dispShrink, dev="pdf", fig.align="center", fig.width=8, fig.height=4, fig.cap="Plot of dispersion estimates over the base mean for a subset of the (A) Bottomly et al and (B) Pickrell et al dataset. Dispersion outliers are circled in blue with dotted lines indicating the effect shrinkage would have had on the estimate. Genes were selected for ease of visualization, including an enrichment of dispersion outliers."----
line <- -.1
adj <- -.3
cex <- 1.5
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5))
set.seed(1)
plotDispShrink(ddsB)
mtext("A",side=3,line=line,adj=adj,cex=cex)
set.seed(1)
plotDispShrink(ddsP)
mtext("B",side=3,line=line,adj=adj,cex=cex)


## ----dispEsts, dev="png", fig.align="center", fig.width=8, fig.height=4, dpi=150, fig.cap="Plot of dispersion estimates over the base mean for a subset of the (A) Bottomly et al and (B) Pickrell et al dataset. Dispersion outliers are circled in blue with dotted lines indicating the effect shrinkage would have had on the estimate. Vertical lines indicate the reciprocal of dispersion on the scale of the samples with the smallest and largest size factor; the estimation of dispersion to the left of this line is difficult as described in the Methods section."----
line <- -.1
adj <- -.3
cex <- 1.5
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5))
set.seed(1)
plotDispEsts(ddsB,cex=.1,ylim=c(1e-8,1e1))
asymptDispB <- attr(dispersionFunction(ddsB),"coefficients")["asymptDisp"]
abline(v=1/(asymptDispB * range(sizeFactors(ddsB))))
mtext("A",side=3,line=line,adj=adj,cex=cex)
set.seed(1)
plotDispEsts(ddsP,cex=.1,ylim=c(1e-8,1e1))
asymptDispP <- attr(dispersionFunction(ddsP),"coefficients")["asymptDisp"]
abline(v=1/(asymptDispP * range(sizeFactors(ddsP))))
mtext("B",side=3,line=line,adj=adj,cex=cex)


## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


