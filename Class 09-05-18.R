#Class 09-05-18

library(SpikeInSubset)

data(spikein95)

image(spikein95)

ids <- geneNames(spikein95)

ids[1:10]

mas5.eset <- mas5(spikein95)

mas5.e <- log2(exprs(mas5.eset))

boxplot(spikein95)

x11()

boxplot(mas5.e, col = 2:5)

density1 <- density(mas5.e[, 1])

plot(density1, main = "MAS5 expression measure distributions")

density2 <- density(mas5.e[, 2])

lines(density2, col = "red")

density3 <- density(mas5.e[, 3])

lines(density3, col = "blue")



# Making MA plots 

# M: difference in average log intensities

# A: average log intensities



d <- rowMeans(mas5.e[,1:3]) - rowMeans(mas5.e[,4:6])

a <- rowMeans(mas5.e)

plot(a, d, ylim = c(-5, 5), main = "MAS 5.0 MA plot", xlab = "A", ylab = "M", pch = ".")

abline(h = c(-1, 1))



#Finding specific probesets in the plot
#() function pdata- phenotype/meta data - probes and ids
#colnames- columns are probe names and rows are samples
#to find specific gene, look up affy probe id ESR1 <- "1298_at", "1682_at", "1893_s_at" #map on plot with green dot
#use gene annot website genecards.weizmann

genes_of_interest <- c("100_g_at", "1000_at", "1002_f_at")

spikedin <- colnames(pData(spikein95))

spikedIndex <-match(genes_of_interest,featureNames(mas5.eset))

points(a[spikedIndex], d[spikedIndex],pch = 19, col = "green")

genes_of_interest <- c("669_s_at", "670_s_at", "671_at", "672_at", "673_at")

spikedIndex <-match(genes_of_interest,featureNames(mas5.eset))

points(a[spikedIndex], d[spikedIndex],pch = 19, col = "orange")


spikedin <- colnames(pData(spikein95))

spikedIndex <- match(spikedin,featureNames(mas5.eset))

points(a[spikedIndex], d[spikedIndex],pch = 19, col = "red")

genes_of_interest
[1] "669_s_at" "670_s_at" "671_at"   "672_at"   "673_at"  
> spikedIndex <- match(genes_of_interest, featureNames(mas5.eset))
> points(a[spikedIndex], d[spikedIndex], pch=19, col="green")


# Trying RMA normalization



rma.eset <- rma(spikein95)

rma.e <- exprs(rma.eset)



# Making MA plots 

x11()

d <- rowMeans(rma.e[,1:3]) - rowMeans(rma.e[,4:6])

a <- rowMeans(rma.e)

plot(a, d, ylim = c(-2, 2), main = "RMA MA plot", xlab = "A", ylab = "M", pch = ".")

abline(h = c(-1, 1))


abline(h = c(-.5, .5), col="red")
abline(h = c(-.3, .3), col="orange")

#go to bioconductor and DL 
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")

#Load

library("genefilter")

#Plot







pData(rma.eset) <- pData(mas5.eset)

tt <- rowMeans(rma.e)

lod <- tt

plot(d, lod, cex = 0.25, main = "Volcano plot for MA", xlim = c(-2, 2), xlab = "M", ylab = "A", yaxt = "n")

axis(2, at = seq(0, 3, by = 1), labels = 10^(-seq(0, 3, by = 1)))

points(d[spikedIndex], lod[spikedIndex], pch = 19, col = "red")

abline(h = 2, v = c(-1, 1))

