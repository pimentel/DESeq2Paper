
## ----loadData------------------------------------------------------------
library("DESeq2paper")
data("specificity")
alpha <- .01
# for cuffdiff2: set the genes with zero row sum to have p-value = NA
# as we do not count these to the denominator of the FPR for other algos
for (i in seq_len(length(resList))) {
  resList[[i]]$cuffdiff2[is.na(resList[[i]]$DESeq2)] <- NA 
}
resMat <- t(sapply(resList, function(z) colMeans(z < alpha, na.rm=TRUE)))
colnames(resMat) <- namesAlgos
resMat <- subset(resMat,select=-EBSeq) # EBSeq does not produce p-values
library("ggplot2")
d <- data.frame(fpr=as.vector(resMat),
                algorithm=factor(rep(colnames(resMat),each=nrow(resMat)),
                levels=colnames(resMat)))
# these points are outliers
d[d$fpr >= .05,]


## ----renameAtoB----------------------------------------------------------
renameAtoB <- function(f,a,b) {
  levels(f)[levels(f) == a] <- b
  f
}


## ----rename--------------------------------------------------------------
d$algorithm <- renameAtoB(d$algorithm, "DESeq", "DESeq (old)")
d$algorithm <- renameAtoB(d$algorithm, "cuffdiff2", "Cuffdiff 2")


## ----specificityBoxplot, dev="pdf", fig.align="center", fig.width=4.5, fig.height=3, fig.cap="Estimate of E($p$-value $< .01$), constructed by counting the number of $p$-values less than $.01$ and dividing by the total number of tests. Genes with all zero counts not included. Plot cropped to remove outliers printed in code above."----
p <- ggplot(d, aes(x=reorder(algorithm, fpr, median), y=fpr, color=algorithm))
p + geom_boxplot(outlier.colour=rgb(0,0,0,0)) + theme_bw() +
  geom_point(position = position_jitter(w = 0.1, h = 0), color="grey50", size=1) +
  geom_hline(aes(yintercept=alpha)) + 
  ylab("1 - specificity (FPR)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("") + scale_colour_discrete(guide="none") + 
  coord_cartesian(ylim=c(0,.05))



## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


