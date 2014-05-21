## Make some graphics to show rgr fitting results
dat <- read.csv("~/work/data/data/hazards/hazards-bc-firs.csv")

## Table of different fits used to calculate the rGr
library(plotrix)
degrees <- as.data.frame(table(dat$degree))
names(degrees) <- c("degree", "number")

## rGR residuals from hazard fitting
pdf("~/work/hazards/pics/hazard-rgr-vs-size.pdf")
plot(dat$priorbv, dat$rgr, xlab = "size", ylab = "rGR", main = "Hazard: rGR vs. Size")
abline(h = 1, lty = 2, col = "red", lwd = 2)
addtable2plot(15, 4, degrees, cex=1.5)
dev.off()
