
## ----lib-----------------------------------------------------------------
library("DESeq2")
library("DESeq2paper")


## ----simulate------------------------------------------------------------
makeSimScript <- system.file("script/makeSim.R", package="DESeq2paper", mustWork=TRUE)
source(makeSimScript)
data("meanDispPairs")
n <- 4000
m <- 10
condition <- factor(rep(c("A","B"), each = m/2))
x <- model.matrix(~ condition)
beta <- rep(0,n)
sim <- makeSim(n,m,x,beta,meanDispPairs)
mat <- sim$mat
dds <- DESeqDataSetFromMatrix(mat, DataFrame(condition), ~ condition)
dds <- DESeq(dds)
data <- log10(cbind("genewise"=mcols(dds)$dispGeneEst,
                    fit=mcols(dds)$dispFit,
                    maximum=pmax(mcols(dds)$dispGeneEst, mcols(dds)$dispFit),
                    MAP=dispersions(dds),
                    true=sim$disp))
data <- data[data[,"genewise"] > -7,]
data <- data[rowSums(is.na(data)) == 0,]


## ----dispAccuracyPairs, dev="png", fig.align="center", fig.width=10, fig.height=10,fig.cap="Scatterplot of various estimates of dispersion using DESeq2, against the true dispersion in the logarithmic scale (base 10) from simulated counts. The blue, red, and yellow colors indicate regions of increasing density of points. Counts for 4000 genes and for 10 samples were simulated for two groups with no true difference in means. The Negative Binomial counts had mean and dispersion drawn from the joint distribution of the mean and gene-wise dispersion estimates from the Pickrell et al dataset.  The estimates shown are genewise, the CR-adjusted maximum likelihood estimate; fit the value from the fitted curve; maximum, the maximum of the two previous values (the estimate used in the older version of DESeq); and MAP, the maximum a posteriori estimate used in DESeq2. The correlations shown in the bottom panels do not include the very low gene-wise estimates of dispersion which can result in potential false positives. The MAP, shrunken estimates used in DESeq2 were closer to the diagonal, while the maximum estimate was typically above the true value of dispersion, which can lead to overly-conservative inference of differential expression."----
library("LSD")
heatpairs(data, xlim=c(-2,2),ylim=c(-2,2), main="")


## ----sessInfo, echo=FALSE, results="asis"--------------------------------
toLatex(sessionInfo())


