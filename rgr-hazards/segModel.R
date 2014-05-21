library(segmented)
source("~/Allometry/simulation/hazard-model.R")
bclong1 <- read.csv("~/Allometry/data/long-bc-derived.csv")
bc <- subset(bclong1, stat == "ALIVE" & bvgrowth>=0)
bc$sdpclass <- factor(as.numeric(bc$sdpclass))

tst <- subset(bc, install == 1 & plot == 16 & time == 76)
pb <- tst$priorbv
bv <- tst$bvgrowth
bounds <- cont2class(tst, "priorbv", 3)

## try out segmented model
fit <- lm(bvgrowth ~ priorbv, data = tst)
##tst[,c("priorbv","bvgrowth")] <- tst[,c("priorbv","bvgrowth")] * 1000
seg.fit <- segmented(fit, seg.Z = ~priorbv, psi = list(priorbv = c(.1)),
                     seg.control(K = 5))
par(mfrow = c(1,1))
plot(tst$priorbv, tst$bvgrowth)
plot(seg.fit, add = TRUE)
lines(seg.fit, col = "blue")
abline(fit, col = "red")
slope(seg.fit)
plot(seg.fit)

## res
plot(tst$bvgrowth, residuals(seg.fit))
points(tst$bvgrowth, residuals(fit), col = "blue")
rmse1 <- sqrt(sum(residuals(fit)^2))
rmse2 <- sqrt(sum(residuals(seg.fit)^2))

## lrs, rgr
lrs1 <- (predict(fit)+tst$priorbv)/max(tst$priorbv + predict(fit))
lrs2 <- (predict(seg.fit)+tst$priorbv)/max(tst$priorbv + predict(seg.fit))
rgr1 <- predict(fit)/tst$bvgrowth
rgr2 <- predict(seg.fit)/tst$bvgrowth

# cor
par(mfrow = c(1,3))
plot(rgr1, tst$bv, main="Single fit")
plot(rgr2, tst$bv, main="Segmented")
rgr3 <- rgr2; rgr3[rgr3<0] <- 0;
plot(rgr3, tst$bv, main="Segmented, Below zero rGR set to 0")
cor1 <- cor.test(rgr1, tst$bv)
cor2 <- cor.test(rgr2, tst$bv)
cor3 <- cor.test(rgr3, tst$bv)

points(tst$priorbv, tst$bvgrowth)
intercept(seg.fit)
vcov(seg.fit)
draw.history(seg.fit)

library(ggplot2)
ggplot(tst, aes(priorbv, bvgrowth))
## How to fit segmented models to the data?
