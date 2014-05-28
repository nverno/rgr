################################################################################
##
##   Look at the effect of adding more degrees to bspline poly fits
##
################################################################################

## data: show for install 1 plot 10 time 76
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")

## Test data, only look at ALIVE firs that grew
bclong <- subset(bclong1, spec == "FD" & stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees
tst <- bclong[bclong$install == 1 & bclong$plot == 10 & bclong$time == 85,]
pind <- tst[,"priorbv"]
ind <- tst[,"bv"]
dep <- tst[,"bvgrowth"]
degree <- 9

################################################################################
##
##        Make Plot Showing B-splines fits of different degrees
##
################################################################################
plot(ind, dep, pch=1, main="Fitting B-spline Polynomials with\nVarying Degrees",
     ylab = "BV growth", xlab = "BV")
dd <- data.frame(ind=ind, dep=dep)
dd <- dd[order(dd$ind),]
for (i in 1:degree) {
    bmat <- bs(dd$ind, degree = i)
    fit <- lm(dd$dep ~ bs(dd$ind, degree = i) - 1)
    print(paste("df: ", i))
    lines(dd$ind, predict(lm(dd$dep ~ bs(dd$ind, degree = i) - 1)), col = i)
}
legend("bottomright", legend = c(1:degree), col = 1:degree, lty = 1)
