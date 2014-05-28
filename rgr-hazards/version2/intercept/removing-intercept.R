################################################################################
##
##   Test different methods of removing intercept term in regression
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

## Centered data
indCent <- ind - mean(ind)
depCent <- dep - mean(dep)

################################################################################
##
##                                   fits
##
################################################################################
## Poly fits
fit1 <- lm(depCent ~ poly(indCent, degree = 3))
fit2 <- lm(depCent ~ poly(indCent, degree = 3) - 1)
fit3 <- lm(dep ~ poly(ind, degree = 3))
fit4 <- lm(dep ~ poly(ind, degree = 3) - 1)

## bs fits
fit5 <- lm(dep ~ bs(ind, degree = 3))
fit6 <- lm(dep ~ bs(ind, degree = 3) - 1)
fit7 <- lm(depCent ~ bs(indCent, degree = 3))
fit8 <- lm(depCent ~ bs(indCent, degree = 3) - 1)

## linear
fit9 <- lm(dep ~ ind)
fit10 <- lm(depCent ~ indCent)
fit11 <- lm(dep ~ ind - 1)
fit12 <- lm(depCent ~ indCent - 1)

## Best is fit6, B-spline & remove intercept from linear model

