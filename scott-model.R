# Reworking scott's model: local relative size (LRS)
#  and residual (ratio) Growth Rate (rGR) as two independent
#  measures of size and growth
#
# All sampling years: 73, 76, 79, 82, 85, 91, 97 (years: 0, 3, 6, 9, 12, 18, 24)
#  Three 6-yr periods: Period 1 (yrs 3- 9), Period 2 (12-18) and Period 3 (18-24)
#
# Bole volumes: calculated from height, dbh, taper equations (Kozak 1988)
#  Note: I think those are the simplified bv equations, need to look into it.
#
# LRS = absolute size of individual/maximum individual size in plot
#
# rGR = obs. growth (in prior three years) / predicted growth (in prior three years)
#  - calculated separately for each 6 yr interval
#  - for period 3, prior size estimated from obs. data by assuming 6 periods
#     of annually compounded growth
#
# Predicted growth = predicted growth based on size per plot, 3 yr growth ~ prior size
#  Models used: chapman-richards growth equation, modified weibull, segmented linear
#   model, gompertz (rejected for poor fit to smaller trees)
#  Model selection: MSE, visual inspection of residuals
#
# Sequence:
#  1. Calculate LRS by year/plot
#  2. Predict growth for each prior yr period based on individual starting size for that
#      period.  Predict separately for each plot using power model.
#  3. rGR: calculate normally for periods prior to Period 1 and 2: (73-76, 82-85). For
#      the third prior period, estimate basal area in 88 by assuming annual
#      compounding of growth.
# #
# Note: using BA and HT separately instead of BV for now, BA should be pretty much
#  the same as BV, will use the same 6 yr sampling intervals with only two 3 yr
#  prior periods.  For the last 6 yr period (Period 3), rGR will be calculated from
#  previous 6 yr period (Period 2).
#
source("functions.R")
dflong <- read.csv("long-df.csv")
dflong <- dflong[order(dflong$pplot,dflong$time),]
dflong <- subset(dflong, stat!="DEAD")

# *LRS*
# Max Tree Ht/plot/year and LRH/plot/yr
#  LRH = local relative height
dflong$pplot <- as.factor(dflong$pplot)
maxhts <- ddply(dflong, .(pplot,time), summarise, max(ht, na.rm = TRUE))
times <- ddply(dflong, .(pplot,time), function(x) nrow(x))
dflong$maxplotht <- rep(maxhts[,3], times = times[,3])
dflong$relht <- dflong$ht/dflong$maxplotht

# Max Tree/BA/plot/yr and LRB/plot/yr
#  LRB = local relative basal area
maxbas <- ddply(dflong, .(pplot,time), summarise, max(ba, na.rm = TRUE))
dflong$maxplotba <- rep(maxbas[,3], times = times[,3])
dflong$relba <- dflong$ba/dflong$maxplotba

# Neighbor density in radius of 4 meters
ind.var = "priorba"
fit.MLE.models(dflong, sr=4, spec="FD", ind.var = "priorba", dep.var = "bagrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)
dflong <- merge(dflong, targets, all=TRUE)

# **Predicted Growth**
#   Model: BA ~ a*BA^b
p1 <- ggplot(data=dflong, aes(x=priorba, y=bagrowth)) + geom_point() +
    stat_smooth(method="lm")
p2 <- ggplot(data=dflong, aes(x=priorba, y=bagrowth)) + geom_point() +
    stat_smooth(lwd = 2, color = "red")
mod <- nls(bagrowth ~ a*priorba^b, start = list(a = 0.01, b = 1),
        data = dflong[!is.na(dflong$priorba),])
plot(dflong$priorba, dflong$bagrowth)
curve(coef(mod)[["a"]] * x ^ coef(mod)[["b"]], add = TRUE, col="red", lwd = 2)

# Calculate power model parameters/plot for periods prior to Periods 1-3
#  (0-3, 9-12, 15-18*)
#  * needs to be estimated by assuming annual compound growth
#
# Find compound growth rate to estimate Period 3 priors, 6 yrs of compound growth
#  y = y0(1 + r)^t
#  y^(1/6) = y0(1+r) # 6 compounded growth periods
#  ln(y) = A + tB + err; where A = ln(y0), B = ln(1 + r)
#  r = exp(B) - 1
get.compr <- function(size, prior, obs.yrs, est.yrs = NULL) {
    ifelse(length(size)>1,
       { comp.fit <- lm(log(size^(1/obs.yrs)) ~ log(prior))
         r <- exp(coef(comp.fit)[[1]]) - 1 },
           r <- (size/prior)^(1/obs.yrs) - 1)
    ifelse(!missing(est.yrs),
           return(prior*(1 + r)^est.yrs),
           return(r))
}

# Estimate 88 data assuming annually compounded growth,
#  fit power models/pplot/yr, and
#  use ba - priorba for trees with time = 76, 85, or 91

priors <- c(76, 85, 91)
dfprior <- subset(dflong, time %in% priors & !is.na(priorba))
rgr <- rbind.fill(ddply(dfprior, .(pplot, time), .fun = function(x) {
    ifelse(x[1,"time"] == 91, {
        x[,"est.priorba"] <- apply(x, 1, function(z) {
            get.compr(as.numeric(z["ba"]), as.numeric(z["priorba"]), 6, 3)
        })
        fit <- nls((ba-est.priorba) ~ a*est.priorba^b, start=list(a=.1,b=1), data=x)
        rgr <- (x[,"ba"]-x[,"est.priorba"])/predict(fit)#, newdata = x[,"est.priorba"])
    }, {
        fit <- nls((ba-priorba) ~ a*priorba^b, start=list(a=.1,b=1), data=x)
        rgr <- (x[,"ba"]-x[,"priorba"])/predict(fit) #, newdata = x[,"priorba"])
    })
                                        # interactively show graphs
    if(x[1,"time"] == 91) {
        x[,"priorba"] <- x[,"est.priorba"]
    }
    x$bagrowth <- x$ba - x$priorba
    plot <- x$pplot[1]
    year <- x$time[1]
    par(ask = TRUE, mfrow = c(2,2) )
    scatter.smooth(fitted(fit), residuals(fit), lpars = list(col="red"),
                   main = paste("Residuals vs Fitted",plot,year))
    abline(h = 0, lty = 2)
    scatter.smooth(x$relba, rgr, main = "rGR vs LRS", lpars = list(col="red"))
    abline(h = 1, lty = 2)
    plot(x$priorba, x$bagrowth, main = "BA growth vs. Prior BA with Predicted")
    curve(coef(fit)[[1]] * x ^ coef(fit)[[2]], add=TRUE, col="blue", lwd=2)
    if(length(x$numnebs[!is.na(x$numnebs)])>1) {
        scatter.smooth(x$numnebs, rgr, main = "rGR vs Neighbor Density",
                       lpars=list(col="red"))
        abline(h=1, lty=2);
        abline(lm(rgr ~ x$numnebs), col="blue")
    }
    data.frame(rgr = rgr)
}))

## tst <- data.frame(x = 1:4, y=c(11,2,NA,5))
## fi <- lm(y ~ x, data = tst)
## plot(tst$x, tst$y)
## abline(fi)
## length(tst$y[!is.na(tst$y)])

## nrow(dfprior[!is.na(dfprior$numnebs),])

dfprior <- cbind(dfprior, rgr = rgr$rgr)

## plot(dfprior$relba, dfprior$numnebs)
## plot(dfprior$numnebs, dfprior$rgr)

## mean(tst$rgr)
## median(tst$rgr)

## ggplot(tst, aes(rgr)) + geom_histogram()
## par(mfrow = c(2,1))
## plot(tst$rgr, tst$bagrowth)
## plot(tst$priorba, tst$bagrowth)
## plot(tst$priorba, tst$rgr)

## plot(tst$relht, tst$rgr); abline(h=1, col="red")
## names(tst)

# Background: priorht bad predictor of ht growth
p1 <- ggplot(dfprior, aes(priorht, htgrowth)) + geom_point() + geom_smooth()
p1
p2 <- ggplot(dfprior, aes(priorba, bagrowth)) + geom_point() + geom_smooth()
p2
p3 <- ggplot(dfprior, aes(priordbh, dbhgrowth)) + geom_point() + geom_smooth()

# as predictors of rgr
p4 <- ggplot(dfprior, aes(priorht, rgr)) + geom_point() + geom_smooth()
p5 <- ggplot(dfprior, aes(priorba, rgr)) + geom_point() + geom_smooth()
p6 <- ggplot(dfprior, aes(priordbh, rgr)) + geom_point() + geom_smooth()
p4
p5
p6

p7 <- ggplot(dfprior, aes(numnebs, rgr)) + geom_point() + geom_smooth()
p7
p8 <- ggplot(dfprior, aes(bagrowth, rgr)) + geom_point() + geom_smooth()
p8
p9

# neighbor model with rgr
# Simple neighbor model as predictor of rgr, no size
simple = function(ps, ind.var = "priorba")
{
    PG = ps[["PG"]]
    alpha = ps[["alpha"]]
    beta = ps[["beta"]]
    C = ps[["C"]]
    D = ps[["D"]]
    nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)
    competition.effect <- exp(-(C) * nci^D)
    PG * competition.effect
}

dat = dfprior
sr = 4
spec = "FD"
ind.var = "rgr"
dep.var = "bagrowth"
models = "simple"
method = "Nelder-Mead"
maxit = 1000
realdist = TRUE
ps <- get.params(sr,spec,ind.var,dep.var,currentmodel = models)
tokeep <- c("alpha","beta","C","D","PG","sd")
ps <- ps[tokeep]

fit <- fit.MLE.models(dat,sr,spec,ind.var,dep.var,models=models,method="SANN",
                      realdist = realdist, maxit = 10000)

# get fit parameters and predict nci
newps <- get.params(sr, spec, ind.var, dep.var, models)
targets$NCI <- simple(newps, ind.var = ind.var)/ps[["PG"]]

plot(targets$NCI, targets$bagrowth)
abline(b= 0, a = ps[["PG"]])

plot(targets$rgr, targets$bagrowth)

ggplot(targets, aes(NCI, bagrowth)) + geom_point(alpha = 0.4) +
    geom_smooth()
ggplot(targets, aes(NCI, rgr)) + geom_point(alpha = 0.4) +
    geom_smooth()



