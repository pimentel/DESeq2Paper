
## ----setup---------------------------------------------------------------
library("DESeq2")
ms <- rep(c(6,8,16),c(2,4,4))
ps <- rep(c(2,3,2,3),c(4,2,2,2))
alphas <- rep(c(.05,.2), 5)

set.seed(1)
d <- data.frame()
for (i in seq_along(ms)) {
  m <- ms[i]
  p <- ps[i]
  alpha <- alphas[i]
  theorvar <- trigamma((m-p)/2)
  dds <- makeExampleDESeqDataSet(n=4000,m=m,interceptMean=10,
                                 interceptSD=0,dispMeanRel=function(x) alpha)
  colData(dds)$group <- factor(rep(c("X","Y"),times=m/2))
  design(dds) <- if (p == 2) { ~ condition } else { ~ group + condition }
  sizeFactors(dds) <- rep(1,ncol(dds))
  dds <- estimateDispersionsGeneEst(dds)
  disp <- mcols(dds)$dispGeneEst
  # exclude the dispersions which head to -Infinity
  samplevar <- var(log(disp[disp > 1e-7]))
  d <- rbind(d, data.frame(m=m,p=p,alpha=alpha,theorvar=theorvar,samplevar=samplevar))
}
 


## ----out, results="asis"-------------------------------------------------
library("xtable")
names(d) <- c("m","p","disp.","theor. var.","sample var.")
print(xtable(d,digits=c(0,0,0,2,3,3)),include.rownames=FALSE)


## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


