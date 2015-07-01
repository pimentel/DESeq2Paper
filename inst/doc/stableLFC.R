
## ----lib-----------------------------------------------------------------
library("DESeq2")
library("DESeq2paper")
library("caret")


## ----stableLFCsetup, cache=TRUE------------------------------------------
data("bottomly_sumexp")
dds <- DESeqDataSetFromMatrix(assay(bottomly), colData(bottomly), ~ experiment.number + strain)
dds <- dds[rowSums(counts(dds)) > 0,]
cond1 <- which(colData(dds)$strain == "C57BL/6J")
cond2 <- which(colData(dds)$strain == "DBA/2J")

set.seed(1)
idx1 <- cond1[createDataPartition(colData(dds)$experiment.number[cond1],p=.5)[[1]]][1:5]
idx2 <- cond2[createDataPartition(colData(dds)$experiment.number[cond2],p=.5)[[1]]][1:5]
idx <- c(idx1,idx2)
table(1:ncol(dds) %in% idx, colData(dds)$strain)
table(1:ncol(dds) %in% idx, colData(dds)$experiment.number)

ddsList <- list(dds[,idx], dds[,-idx])
ddsListNoPrior <- list()
ddsListPC <- list()

for (i in 1:2) {
  ddsList[[i]] <- DESeq(ddsList[[i]])
  ddsListNoPrior[[i]] <- nbinomWaldTest(ddsList[[i]],betaPrior=FALSE)
  ddsListPC[[i]] <- ddsList[[i]]
  counts(ddsListPC[[i]]) <- counts(ddsListPC[[i]]) + 1L
  ddsListPC[[i]] <- DESeq(ddsListPC[[i]],betaPrior=FALSE)
}


## ----rmse----------------------------------------------------------------
rmseNoPriorNoPrior <- sqrt(mean((results(ddsListNoPrior[[1]])$log2FoldChange - 
                                 results(ddsListNoPrior[[2]])$log2FoldChange)^2, na.rm=TRUE))
rmsePriorPrior <- sqrt(mean((results(ddsList[[1]])$log2FoldChange - 
                             results(ddsList[[2]])$log2FoldChange)^2, na.rm=TRUE))
rmseNoPriorPrior <- sqrt(mean((results(ddsListNoPrior[[1]])$log2FoldChange - 
                               results(ddsList[[2]])$log2FoldChange)^2, na.rm=TRUE))
rmseNoPriorPC <- sqrt(mean((results(ddsListNoPrior[[1]])$log2FoldChange - 
                            results(ddsListPC[[2]])$log2FoldChange)^2, na.rm=TRUE))


## ----stableLFC, dev="png", fig.align="center", fig.width=6, fig.height=3.25, dpi=150, fig.cap="Comparing stability of logarithmic fold changes across two balanced, random subsets of the Bottomly et al dataset, (A) without a prior on logarithmic fold changes, (B) with a zero-centered Normal prior on logarithmic fold changes"----

line <- 0.4
adj <- -.3
cex <- 1.5

plotMM <- function(x,y,s,l=4,...) {
  idx <- abs(x) < 10 & abs(y) < 10
  x <- x[idx]
  y <- y[idx]
  s <- s[idx]
  plot(x,y,xlim=c(-l,l),ylim=c(-l,l),type="n",...)
  abline(0,1,col=rgb(0,0,0,1))
  abline(v=0,h=0,col=rgb(0,0,0,1))
  cols <- ifelse(s,rgb(1,0,0,.5),rgb(0,0,0,.2))
  points(x,y,cex=.6,col=cols,pch=20)
}

rmselg <- function(rmse) {
  legend("bottomright",legend=paste("RMSE:",round(rmse,2)),
         bg="white",adj=c(.2,.5),cex=.9)  
}


par(mfrow=c(1,2),mar=c(4.5,4.5,2,1))
plotMM(results(ddsListNoPrior[[1]])$log2FoldChange,
       results(ddsListNoPrior[[2]])$log2FoldChange,
       (results(ddsListNoPrior[[1]],ind=FALSE)$padj < .1 &
        results(ddsListNoPrior[[2]],ind=FALSE)$padj < .1),
       xlab=expression(MLE~log[2]~fold~change),
       ylab=expression(MLE~log[2]~fold~change))
rmselg(rmseNoPriorNoPrior)
mtext("A",side=3,line=line,adj=adj,cex=cex)

plotMM(results(ddsList[[1]])$log2FoldChange,
       results(ddsList[[2]])$log2FoldChange,
       (results(ddsList[[1]],ind=FALSE)$padj < .1 &
        results(ddsList[[2]],ind=FALSE)$padj < .1),
       xlab=expression(MAP~log[2]~fold~change),
       ylab=expression(MAP~log[2]~fold~change))
rmselg(rmsePriorPrior)
mtext("B",side=3,line=line,adj=adj,cex=cex)


## ----catPlot, dev="pdf", fig.align="center", fig.width=5, fig.height=5, fig.cap="Concordance at the top (CAT) plot showing that the shrunken LFCs and and pseudocount-based estimators for logarithmic fold change provided more stable rankings compared to ranking based on the unhsrunken LFCs."----

catFn <- function(n, l) {
  top1 <- head(order(-abs(results(l[[1]])$log2FoldChange)),n)
  top2 <- head(order(-abs(results(l[[2]])$log2FoldChange)),n)
  length(intersect(top1,top2))
}

ns <- c(50,100,200,1000,2000)
priorCat <- sapply(ns, catFn, ddsList)/ns
noPriorCat <- sapply(ns, catFn, ddsListNoPrior)/ns
pcCat <- sapply(ns, catFn, ddsListPC)/ns

data.frame(n=ns, MAP=priorCat*ns, MLE=noPriorCat*ns, PC=pcCat*ns)

plot(ns, priorCat, type="b", col="blue",
     ylim=c(0,1),lwd=2,log="x",
     xlab="size of list",
     ylab="proportion in common")
points(ns, noPriorCat, type="b", col="purple",lwd=2)
points(ns, pcCat, type="b", col="forestgreen",lwd=2)
legend("topright",legend=c("MAP","pseudocount","MLE"),
       col=c("blue","forestgreen","purple"),
       lwd=2, pch=1)


## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


