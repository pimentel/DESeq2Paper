
## ----lib-----------------------------------------------------------------
library("DESeq2paper")


## ----loadDE--------------------------------------------------------------
data("results_simulateDE")
# SAMseq has no calls for the 3 vs 3 comparison
res <- res[!(res$m == 6 & res$algorithm == "SAMseq"),]
res$m <- factor(res$m)
levels(res$m) <- paste0("m=",levels(res$m))
res$effSize <- factor(res$effSize)
levels(res$effSize) <- c("fold change 2","fold change 3","fold change 4")
res$algorithm <- factor(res$algorithm)
levels(res$algorithm)[levels(res$algorithm) == "DESeq"] <- "DESeq (old)"
resMinusEBSeq <- res[res$algorithm != "EBSeq",] # EBSeq does not produce p-values
resMinusEBSeq$algorithm <- droplevels(resMinusEBSeq$algorithm)


## ----simulateDE, dev="pdf", fig.width=7, fig.height=5, fig.align="center", fig.cap="Use of simulation to assess the sensitivity and specificity of algorithms across combinations of sample size and effect size. The sensitivity was calculated as the fraction of genes with adjusted p-value less than 0.1 among the genes with true differences between group means. The specificity was calculated as the fraction of genes with p-value greater than 0.01 among the genes with no true differences between group means.  The p-value was chosen instead of the adjusted p-value, as this allows for comparison against the expected fraction of p-values less than a critical value given the uniformity of p-values under the null hypothesis. DESeq2 often had the highest sensitivity of those algorithms which control the false positive rate, i.e., those algorithms which fall on or to the left of the vertical black line (1\\% p-values less than 0.01 for the non-DE genes)."----
library("ggplot2")
p <- ggplot(resMinusEBSeq, aes(y=sensitivity, x=oneminusspecpvals, color=algorithm, shape=algorithm))
p + geom_point() + theme_bw() + facet_grid(effSize ~ m) +
  scale_shape_manual(values=1:9) +
  xlab("1 - specificity (false positive rate)") + 
  coord_cartesian(xlim=c(-.003,.035)) + 
  geom_vline(xintercept=.01) + 
  scale_x_continuous(breaks=c(0,.02))



## ----simulateDEPrec, dev="pdf", fig.width=7, fig.height=5, fig.align="center", fig.cap="Sensitivity and precision of algorithms across combinations of sample size and effect size. The sensitivity was calculated as the fraction of genes with adjusted p-value less than 0.1 among the genes with true differences between group means.  The precision was calculated as the fraction of genes with true differences between group means among those with adjusted p-value less than 0.1. DESeq2 often had the highest sensitivity of those algorithms which controlled the false discovery rate, i.,e., those algorithms which fall on or to the left of the vertical black line."----
library("ggplot2")
p <- ggplot(res, aes(y=sensitivity, x=oneminusprec, color=algorithm, shape=algorithm))
p + geom_point() + theme_bw() + facet_grid(effSize ~ m) +
  scale_shape_manual(values=1:9) +
  xlab("1 - precision (false discovery rate)") + 
  coord_cartesian(xlim=c(-.03, .3)) + 
  geom_vline(xintercept=.1)



## ----sensStratify--------------------------------------------------------
library("reshape")
id.vars <- c("algorithm","effSize","m")
measure.vars <- c("sens0to20","sens20to100","sens100to300","sensmore300")
melted <- melt(res[,c(id.vars,measure.vars)], id.vars=id.vars, measure.vars=measure.vars)
names(melted) <- c(id.vars, "aveexp", "sensitivity")
levels(melted$aveexp) <- c("<20","20-100","100-300",">300")


## ----simulateDESensStratify, dev="pdf", fig.width=7, fig.height=5, fig.align="center", fig.cap="The sensitivity of algorithms across combinations of sample size and effect size, and further stratified by the mean of counts of the differentially expressed genes in the simulation data. Points indicate the average over 6 replicates. Algorithms all show a similar dependence of sensitivity on the mean of counts. The height of the sensitivity curve should be compared with the previous plot indicating the total sensitivity and specificity of each algorithm."----
p <- ggplot(melted, aes(y=sensitivity, x=aveexp, group=algorithm, color=algorithm, shape=algorithm))
p + stat_summary(fun.y="mean", geom="line") +
  stat_summary(fun.y="mean", geom="point") +
  theme_bw() + facet_grid(effSize ~ m) +
  scale_shape_manual(values=1:9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("mean counts")



## ----loadOut-------------------------------------------------------------
data("results_simulateOutliers")
# when < 7 replicates DESeq does not replace
res <- res[!(res$algorithm == "DESeq2-noRepl" & res$m < 14),]
# when >= 7 replicates DESeq does not filter
res <- res[!(res$algorithm == "DESeq2-noFilt" & res$m >= 14),]
res$m <- factor(res$m)
levels(res$m) <- paste0("m=",levels(res$m))
res$percentOutlier <- 100 * res$percentOutlier
res$percentOutlier <- factor(res$percentOutlier)
levels(res$percentOutlier) <- paste0(levels(res$percentOutlier),"% outlier")


## ----outlierPadj---------------------------------------------------------
resSensPadj <- res[res$senspadj < .1,]
resSensPadj <- resSensPadj[nrow(resSensPadj):1,]
resSensPadj <- resSensPadj[!duplicated(with(resSensPadj, paste(algorithm, m, percentOutlier))),]
summary(resSensPadj$senspadj)


## ----simulateOutliersSens, dev="pdf", fig.width=8, fig.height=4, fig.align="center", fig.cap="Sensitivity-specificity curves for detecting true differences in the presence of outliers. Negative Binomial counts were simulated for 4000 genes and total sample sizes (m) of 10 and 20, for a two-group comparison.  80\\% of the simulated genes had no true differential expression, while for 20\\% of the genes true logarithmic (base 2) fold changes of -1 or 1. The number of genes with simulated outliers was increased from 0\\% to 15\\%. The outliers were constructed for a gene by multiplying the count of a single sample by 100.  Sensitivity and specificity were calculated by thresholding on p-values. Points indicate an adjusted p-value of 0.1. DESeq2 with the default settings and edgeR with the robust setting had higher area under the curve compared to running edgeR without the robust option, turning off DESeq2 gene filtering, and turning off DESeq2 outlier replacement. DESeq2 filters genes with potential outliers for samples with 3 to 6 replicates, while replacing outliers for samples with 7 or more replicates, hence the filtering can be turned off for the top row (m=10) and the replacement can be turned off for the bottom row (m=20)."----
library("ggplot2")
p <- ggplot(res, aes(x=oneminusspec, y=sensitivity, color=algorithm))
p + scale_x_continuous(breaks=c(0,.1,.2)) + 
  scale_y_continuous(breaks=c(0,.2,.4,.6,.8)) + 
  geom_line() + theme_bw() +
  facet_grid(m ~ percentOutlier) + xlab("1 - specificity") + 
  coord_cartesian(xlim=c(-.03, .25), ylim=c(-.05, .9)) +
  geom_point(aes(x=oneminusspec, y=sensitivity, shape=algorithm),
           data=resSensPadj)



## ----simulateOutliersPrec, dev="pdf", fig.width=8, fig.height=4, fig.align="center",fig.cap="Outlier handling: One minus the precision (false discovery rate) plotted over various thresholds of adjusted p-value. Shown is the results for the same simulation with outliers described in the previous figure. Points indicate an adjusted p-value of 0.1. edgeR run with the robust setting had false discovery rate generally above the nominal value from the adjusted p-value threshold (black diagonal line). DESeq2 run with default settings was generally at or below the line, which indicates control of the false discovery rate."----
p <- ggplot(res, aes(x=precpadj, y=oneminusprec, color=algorithm))
p + scale_x_continuous(breaks=c(0,.1,.2)) + 
  scale_y_continuous(breaks=c(0,.1,.2)) + 
  geom_line() + theme_bw() + 
  facet_grid(m ~ percentOutlier) + 
  geom_abline(intercept=0,slope=1) +
  xlab("adjusted p-value") + ylab("1 - precision (FDR)") + 
  coord_cartesian(xlim=c(-.03, .25), ylim=c(-.05, .25)) + 
  geom_point(aes(x=precpadj, y=oneminusprec, shape=algorithm),
             data=res[res$precpadj == .1,])




## ----lfcAccuracyHist, fig.height=2.5, fig.width=8, fig.align="center", fig.cap="Benchmarking LFC estimation: Models for simulating logarithmic (base 2) fold changes. For the bell model, true logarithmic fold changes were drawn from a Normal with mean 0 and variance 1. For the slab bell model, true logarithmic fold changes were drawn for 80\\% of genes from a Normal with mean 0 and variance 1 and for 20\\% of genes from a Uniform distribution with range from -4 to 4. For the slab spike model, true logarithmic fold changes were drawn similarly to the slab bell model except the Normal is replaced with a spike of logarithmic fold changes at 0. For the spike spike model, true logarithmic fold changes were drawn according to a spike of logarithmic fold changes at 0 (80\\%) and a spike randomly sampled from -2 or 2 (20\\%). These spikes represent fold changes of 1/4 and 4, respectively."----
par(mfrow=c(1,4),mar=c(3,3,3,1))
n <- 1000
brks <- seq(from=-4,to=4,length.out=20)
trimit <- function(x) x[x > -4 & x < 4] # for visualization only
hist(trimit(rnorm(n)), breaks=brks, col="black", main="bell", xlab="", ylab="")
hist(trimit(c(rnorm(n * 8/10), runif(n * 2/10, -4, 4))), breaks=brks, col="black", main="slab bell", xlab="", ylab="")
hist(c(rep(0, n * 8/10), runif(n * 2/10, -4, 4)), breaks=brks, col="black", main="slab spike",xlab="", ylab="")
hist(c(rep(0, n * 8/10), sample(c(-2, 2), n * 2/10, TRUE)), breaks=brks, col="black", main="spike spike",xlab="", ylab="")



## ----lfcAccuracy, fig.height=5, fig.width=7, fig.align="center",fig.cap="Root mean squared error (RMSE) for estimating logarithmic fold changes under the four models of logarithmic fold changes and varying total sample size m. Simulated Negative Binomial counts were generated for two groups and for 1000 genes. Points and error bars are drawn for the mean and 95\\% confidence interval over 10 replications. DESeq2 and GFOLD, which both implement posterior logarithmic fold change estimates, had lower root mean squared error to the true logarithmic fold changes over all genes, compared to predictive logarithmic fold changes from edgeR, either using the default value of 0.125 for the edgeR argument prior.count, or after increasing prior.count to 10 (edgeR predFC10)."----
data("results_simulateLFCAccuracy")
library("ggplot2")
library("Hmisc")
p <- ggplot(data=res, aes(x=m, y=RMSE, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point") + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  theme_bw() + xlab("total sample size") + facet_wrap(~ type) + 
  scale_x_continuous(breaks=unique(res$m))




## ----lfcAccuracyDE, fig.height=2.5, fig.width=6, fig.align="center", fig.cap="Root mean squared error (RMSE) of logarithmic fold change estimates, only considering genes with non-zero true logarithmic fold change. For the same simulation, shown here is the error only for the 20\\% of genes with non-zero true logarithmic fold changes (for bell and slab bell all genes have non-zero logarithmic fold change).  DESeq2 had generally lower root mean squared error, compared to GFOLD which had higher error for large sample size and edgeR which had higher error for low sample size."----
p <- ggplot(data=res[grepl("spike",res$type),], aes(x=m, y=DiffRMSE, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point") + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  theme_bw() + xlab("total sample size") +  ylab("RMSE only of DE genes") + 
  facet_wrap(~ type) + 
  scale_x_continuous(breaks=unique(res$m))



## ----lfcAccuracyMAE, fig.height=5, fig.width=7, fig.align="center", fig.cap="Mean absolute error (MAE) of logarithmic fold change estimates. Results for the same simulation, however here using median absolute error in place of root mean squared error. Mean absolute error places less weight on the largest errors. For the bell and slab bell models, DESeq2 and GFOLD had the lowest mean absolute error, while for the slab spike and spike spike models, GFOLD and edgeR with a prior.count of 10 had lowest mean absolute error."----
p <- ggplot(data=res, aes(x=m, y=MAE, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point") + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  theme_bw() + xlab("total sample size") +  ylab("MAE") + 
  facet_wrap(~ type) + 
  scale_x_continuous(breaks=unique(res$m))



## ----lfcAccuracyDiffMAE, fig.height=2.5, fig.width=6,fig.align="center",fig.cap="Mean absolute error (MAE) of logarithmic fold change estimates, only considering those genes with non-zero true logarithmic fold change. While in the previous figure, considering all genes for the slab spike and spike spike models, GFOLD and edgeR with a prior.count of 10 had lowest median absolute error, the median absolute error for these methods was relatively large for large sample size, when considering only the 20\\% of genes with true differentially expression. DESeq2 and edgeR generally had the lowest mean absolute error."----
p <- ggplot(data=res[grepl("spike",res$type),], aes(x=m, y=DiffMAE, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point") + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  theme_bw() + xlab("total sample size") +  ylab("MAE only of DE genes") + 
  facet_wrap(~ type) + 
  scale_x_continuous(breaks=unique(res$m))



## ----simulateCluster, fig.width=8, fig.height=5, fig.align="center", fig.cap="Adjusted Rand Index from clusters using various transformation and distances compared to the true clusters from simulation. The methods assessed were Euclidean distance on counts normalized by size factor, log2 of normalized counts plus a pseudocount of 1, and after applying the rlog and variance stabilizing transformation. Additionally, the Poisson Distance from the PoiClaClu package and the Biological Coefficient of Variation (BCV) distance from the plotMDS function of the edgeR package were used for hierarchical clustering. The points indicate the mean from 20 simulations and the bars are 95 percent confidence intervals. In the equal size factor simulations, the Poisson Distance, variance stablizing transformation (VST), and the rlog transformation had the highest accuracy in recovering true clusters. In the unequal size factor simulations, the size factors for the 4 samples of each group were set to [1, 1, 1/3, 3]. Here, the rlog outperformed the Poisson distance and the VST."----
data("results_simulateCluster")
library("ggplot2")
library("Hmisc")
res$sizeFactor <- factor(res$sizeFactor)
levels(res$sizeFactor) <- paste("size factors", levels(res$sizeFactor))
res$dispScale <- factor(res$dispScale)
levels(res$dispScale) <- paste(levels(res$dispScale),"x dispersion")
p <- ggplot(res, aes(x=rnormsd, y=ARI, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point", aes(shape=method)) + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  facet_grid(sizeFactor ~ dispScale, scale="free") + theme_bw() + 
  ylab("adjusted Rand Index") + xlab("SD of group differences")



## ----simulateCluster1xDisp, fig.width=6, fig.height=3,fig.align="center",fig.cap="Adjusted Rand Index from clusters using various transformation and distances compared to the true clusters from simulation. The same results as the previous figure, only showing the panels with dispersion equal to the estimates from the Pickrell et al dataset (1 x dispersion)."----
data("results_simulateCluster")
library("ggplot2")
library("Hmisc")
res$sizeFactor <- factor(res$sizeFactor)
levels(res$sizeFactor) <- paste("size factors", levels(res$sizeFactor))
res <- res[res$dispScale == 1,]
p <- ggplot(res, aes(x=rnormsd, y=ARI, color=method, shape=method))
p + stat_summary(fun.y=mean, geom="point", aes(shape=method)) + 
  stat_summary(fun.y=mean, geom="line") + 
  stat_summary(fun.data=mean_cl_normal, geom="errorbar") + 
  facet_wrap(~ sizeFactor) + theme_bw() + 
  ylab("adjusted Rand Index") + xlab("SD of group differences")



## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


