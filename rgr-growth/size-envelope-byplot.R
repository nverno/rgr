## Looking at top performers by plot
##  calculate top performers per plot/yr
source("functions.R")
source("neighborhood-models.R")

## data (start with just DF data)
dflong <- read.csv("long-df.csv")
dflong <- subset(dflong, stat!="DEAD" & bagrowth>=0)

## sample size by plot/year
samp <- ddply(dflong, .(pplot, time), function(x) nrow(x))
range(samp$V1)

## Neighbor density in radius of 4 meters
ind.var = "priorba"
fit.MLE.models(dflong, sr=6, spec="FD", ind.var = "priorba", dep.var = "bagrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)
dflong <- merge(dflong, targets, all=TRUE)

## rgr: predicting top performers per plot
##  1. break each plot into three quantile based size classes, plot 40 has only 8
##  2. fit power function to top 25% from each size class
##  rgr = obs/predicted
yrs <- c(76,79,82,85,91,97)
dflong$bacl <- rep(NA, nrow(dflong))
dflong$bagrcl <- rep(NA, nrow(dflong))

# ignore plot 40
dflong <- dflong[dflong$pplot != 40.1,]
rgr <- rbind.fill(ddply(dflong, .(pplot, time), .fun = function(x) {
    tryCatch(expr={
                                        # fitting power curve to all points
        fit.reg <- nls((bagrowth) ~ a*priorba^b, start=list(a=.05,b=.5), data=x)
        rgr.reg <- x[,"bagrowth"]/predict(fit.reg)
                                        # create size classes
        x$bacl <- quantclass(x$priorba,3)
                                        # top 25% (by bagrowth) from each size class
        x <- x[order(x$bacl),]
        top <- rbind.fill(ddply(x, .(bacl), .fun = function(d) {
            bagrcl <- as.numeric(quantclass(d$bagrowth, 4, smallest = 0))
            bagrcl = data.frame(bagrcl = bagrcl)
        }))
        x$bagrcl <- as.factor(top$bagrcl)
                                        # fit to top performers
        fit.top <- nls(bagrowth ~ a*priorba^b, start=list(a=.1,b=.8),
                       data=x[!is.na(x$bagrcl) & x$bagrcl == 4,])
        pred.top <- predict(fit.top, newdata = x)
        rgr.top <- x[,"bagrowth"]/pred.top
                                        # calculate residuals from top fit to all tree
        res.top <- pred.top-x$bagrowth
                                        # ** Calculations done **
                                        # ** interactively show graphs **
        plot <- x$pplot[1]
        year <- x$time[1]
        par(ask = TRUE, mfrow = c(2,2) )
        scatter.smooth(x$priorba, res.top, lpars = list(col="red"),
                       main = paste("Residuals vs Fitted",plot,year))
        abline(h = 0, lty = 2)
        plot(x$priorba, x$bagrowth, main = "BA growth vs. Prior BA with Predicted")
        curve(coef(fit.reg)[[1]] * x ^ coef(fit.reg)[[2]], add=TRUE, col="blue", lwd=2)
        curve(coef(fit.top)[[1]] * x ^ coef(fit.top)[[2]], add=TRUE, col="green", lwd=2)
        hist(res.top, main = "Distribution of the residuals")
        abline(v = 1, lty = 2)
        if(length(x$numnebs[!is.na(x$numnebs)])>2) {
            plot(x$numnebs, rgr.top, main = "rGR vs Neighbor Density")
            abline(h=1, lty=2);
            tryCatch(abline(lm(rgr.top ~ x$numnebs), col="blue"),
                     error=function(e) NULL )
        }
        data.frame(pplot=x$pplot, time=x$time, id=x$id, rgr.top = rgr.top)},
             error=function(e) print(x$pplot))
}))

## combine new rgr
dflong <- dflong[order(dflong$pplot,dflong$time,dflong$id),]
rgr <- rgr[order(rgr$pplot,rgr$time,rgr$id),]

dflong$rgrplot <- rgr$rgr.top
plot(dflong$priorba, dflong$rgrplot)

## run MLE using rgrplot
## create new starting parameters
## pars <- read.csv("parameters.csv")
## rowcopy = pars[83:86,]
## rowcopy$model <- "simplest"
## rowcopy$mdep.var <- "rgrplot"
## pars <- rbind(pars, rowcopy)
## write.csv(pars, "parameters.csv", row.names =FALSE)

# specify parameters for MLE fitting, using 'simplest' model
dat = dflong
sr = c(6)
spec = "FD"
ind.var = "priorba" ## will be used in neighborhood calculations
dep.var = "rgrplot"
currentmodel = "simplest"
method = "Nelder-Mead"
maxit = 10000
realdist = TRUE
pars <- read.csv("parameters.csv")
ps <- get.params(sr,spec,ind.var,dep.var,currentmodel = currentmodel)

fit <- fit.MLE.models(dat=dat,sr=sr,spec=spec,ind.var=ind.var,dep.var=dep.var,
                      models=currentmodel,method="SANN",
                      realdist = realdist, maxit = 100000)

# get fit parameters and predict nci
#  nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)
newps <- get.params(sr = sr, spec = spec, ind.var = ind.var, dep.var = dep.var,
                    currentmodel = currentmodel)
alpha = newps[["alpha"]]; beta = newps[["beta"]]
targets$nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)

# graph rgr vs nci, bagrowth vs nci
par(mfrow = c(1,2))
plot(targets$nci, targets$rgr)
plot(targets$nci, targets$bagrowth)

# graph predicted over obs. on graph of rgr vs priorba
targets$predicted <- do.call(currentmodel, list(newps))
plot(targets$nci, targets$rgr)
##points(targets$priorba, targets$predicted, col = "red")
##points(targets$nci, targets$predicted, col = "blue")

pdf('rgrbyplot.pdf')
plot(targets$priorba, targets$rgr, main = "RGR by plot vs. Prior BA",
       xlab = "Prior BA", ylab = "RGR by Plot")
points(targets$priorba, targets$predicted, col = "red")
dev.off()
##points(targets$nci, targets$predicted, col = "blue")

write.csv(dflong, "dflong-rgr.csv", row.names = FALSE)

###################################################################################
## testing
tst <- subset(dflong, pplot == 1.1)
tst$bacl <- quantclass(tst$priorba,3)
tst <- tst[order(tst$bacl),]
top <- rbind.fill(ddply(tst, .(bacl), .fun = function(d) {
    bagrcl <- as.numeric(quantclass(d$bagrowth, 4, smallest = 0))
    bagrcl = data.frame(bagrcl = bagrcl)
}))
table(top)
tst$bagrcl <- as.factor(top$bagrcl)

ggplot(tst, aes(priorba, bagrowth, color=bagrcl)) + geom_point(alpha = 0.4, pwd =2)
                                        # fit to top performers
ggplot(tst, aes(priorba, bagrowth, color=bagrcl)) + geom_point(alpha=.5) +
    facet_wrap(~bagrcl) + geom_smooth()

fit.top <- nls(bagrowth ~ a*priorba^b, start=list(a=.1,b=.8),
               data=tst[!is.na(tst$bagrcl) & tst$bagrcl == 4,])
pred.top <- predict(fit.top, newdata = tst)
rgr.top <- tst[,"bagrowth"]/pred.top
                                        # calculate residuals from top fit to all tree
res.top <- pred.top-tst$priorba
hist(res.top)

# test the neighbor graph
if(length(tst$numnebs[!is.na(tst$numnebs)])>1) {
    scatter.smooth(tst$numnebs, rgr.top, main = "rGR vs Neighbor Density",
                   lpars=list(col="red"))
    abline(h=1, lty=2);
    abline(lm(rgr.top ~ tst$numnebs), col="blue")
}

