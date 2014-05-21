## Messing around trying to figure out how to fit piecewise nonlinear to
## bvgrowth ~ priorbv
library(splines)
source("~/work/functions/functions.R")
source("~/work/functions/functions-hazards.R")

## data - use long-bc-derived, just doug first
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")
## Remove outlier from install 10, plot 2, year 97
bclong <- subset(bclong1, spec == "FD" & stat == "ALIVE" & bvgrowth>=0)

## some test data
tst <- subset(bclong, install == 1 & plot == 10 & time == 79)

## B-spline basis for polynomial splines
basis <- bs(tst$priorbv, degree = 9)
fm1 <- lm(bvgrowth ~ bs(priorbv, degree = 9), data = tst)

## Examine
plot(tst$priorbv, tst$bvgrowth)
points(tst$priorbv, predict(fm1), col = "blue", pch = 2)
summary(fm1)
cor.test(tst$bvgrowth/predict(fm1), tst$bv)
plot(tst$priorbv, tst$bvgrowth/predict(fm1))

## Get a correlated plot
## run fit-rgr3.R to get list of correlated
tst <- subset(bclong, install==72 & plot==14 & time==97)
basis <- bs(tst$priorbv, degree = 9)
fm1 <- lm(bvgrowth ~ bs(priorbv, degree = 3), data = tst)
fm2 <- lm(bvgrowth ~ bs(priorbv, degree = 2), data = tst)
AIC(fm1)
AIC(fm2)
anova(fm1, fm2)

## Examine
plot(tst$priorbv, tst$bvgrowth)
points(tst$priorbv, predict(fm1), col = "blue", pch = 2)
summary(fm1)
cor.test(tst$bvgrowth/predict(fm1), tst$bv)
plot(tst$priorbv, tst$bvgrowth/predict(fm1))

