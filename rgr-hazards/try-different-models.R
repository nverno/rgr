## Looking for a decent model to predict bvgrowth from priorbv
## Needs to be good fit (lot MSE, distributed residuals) and to be
##  uncorrelated from bole volume at the end of the prediction period
library(ggplot2)
library(segmented)
source("~/work/simulation/hazard-model.R")
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")
bc <- subset(bclong1, stat == "ALIVE" & bvgrowth>=0)
bc$sdpclass <- factor(as.numeric(bc$sdpclass))

tst <- subset(bc, install == 1 & plot == 16 & time == 76)
pb <- tst$priorbv
bv <- tst$bvgrowth
bounds <- cont2class(tst, "priorbv", 3)

## Try log fit compared to seg vs lin vs pow
tst <- tst[tst$bvgrowth>0,]
ggplot(tst, aes(priorbv, bvgrowth)) + geom_point() + geom_smooth()
fit1 <- lm(bvgrowth ~ priorbv, data = tst)
fit2 <- nls(bvgrowth ~ a*priorbv^b, start = list(a = 0.5, b=0.5), data = tst)
fit3 <- segmented(fit1, seg.Z = ~priorbv, psi = list(priorbv = 25))
                  control = seg.control(stop.if.error=FALSE))

fit4 <- lm(log(bvgrowth) ~ log(priorbv), data = tst)
fit5 <- nls(bvgrowth ~ a1*priorbv + a2*priorbv^(1/2) + a3*priorbv^2,
            start = list(a1 = .5, a2 = 0.5, a3 = 0.5), data = tst)
fit6 <- lm(bvgrowth ~ poly(priorbv, 3), data = tst)
fit7 <- lm(bvgrowth ~ poly(priorbv, 7), data = tst)

## lets see
plot(tst$priorbv, tst$bvgrowth)
points(tst$priorbv, fit6$fitted, col = "red")
points(tst$priorbv, fit7$fitted, col = "blue")
abline(fit1, col = "blue")
plot(fit3, col = "green", add = TRUE)
curve(coef(fit2)[[1]]*x^coef(fit2)[[2]], add = TRUE, col = "red")
curve(exp(coef(fit4)[[1]]) * x ^ coef(fit4)[[2]], add = TRUE, col = "purple")
curve(coef(fit5)[[1]] * x + coef(fit5)[[2]] * x^(1/2) + coef(fit5)[[3]]*x^2,
      col = "orange", add = TRUE)
