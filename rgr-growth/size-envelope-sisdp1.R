## Looking at top performers by a combination of SI and SDP (as a class)
##  calculate top performers per Install / Plot / Period
source("functions.R")

## data - use long-bc-derived, all species included
bclong <- read.csv("long-bc-derived.csv")
bclong <- subset(bclong, stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees
## not marked with status
bclong$pplot <- factor(bclong$pplot)
bclong$sdpclass <- factor(as.numeric(bclong$sdpclass))

## sample size by install / plot / period
samp <- ddply(bclong, .(install, plot, time), function(x) nrow(x))
range(samp$V1)

## Neighbor density in radius of 4 meters
ind.var = "priorbv"
fit.MLE.models(bclong, sr=6, spec="FD", ind.var = ind.var, dep.var = "bvgrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)
bclong <- merge(bclong, targets, all=TRUE)

## rgr: predicting top performers per plot
##  1. break each plot into three quantile based size classes, plot 40 has only 8
##  2. fit power function to top 25% from each size class
##  rgr = obs/predicted
yrs <- c(76,79,82,85,91,97)
bclong$bvcl <- rep(NA, nrow(bclong))
bclong$bvgrcl <- rep(NA, nrow(bclong))
errors <- c()
errors.top <- c()

## test run
samps <- ddply(bclong, .(sdpclass, si), function(x) nrow(x) )
tst <- subset(bclong, sdpclass == "2" & si == 38)
fit.reg <- nls(bvgrowth ~ a*priorbv^b, start=list(a=.05,b=.5), data=tst)

rgr <- ddply(bclong, .(sdpclass, si), .fun = function(x) {
    x <- droplevels(x)
    fit.reg <- NULL ## Try NLS fit
    try(fit.reg <- nls((bvgrowth) ~ a*priorbv^b, start=list(a=.05,b=.5), data=x,
                       control = nls.control(warnOnly = TRUE)))
    if(is.null(fit.reg)) { ## NLS failed, do something
        print(paste(levels(x$sdpclass),mean(x$si)))
    }
    else { ## NLS fit successful
        ## if(length(x[,"bvgrowth"])!=length(predict(fit.reg)))
        ##     print(paste(levels(x$sdpclass), mean(x$si)))
        rgr.reg <- x[,"bvgrowth"]/predict(fit.reg)
        x$bvcl <- quantclass(x$priorbv,3) ## create size classes
        x <- x[order(x$bvcl),] ## top 25% (by bvgrowth) from each size class
        top <- rbind.fill(ddply(x, .(bvcl), .fun = function(d) {
            bvgrcl <- as.numeric(quantclass(d$bvgrowth, 4, smallest = 0))
            bvgrcl = data.frame(bvgrcl = bvgrcl)
        }))
        x$bvgrcl <- as.factor(top$bvgrcl)
        fit.top <- NULL ## try NLS fit on top performers
        try(fit.top <- nls(bvgrowth ~ a*priorbv^b, start=list(a=.1,b=.8),
                           data=x[!is.na(x$bvgrcl) & x$bvgrcl == 4,],
                           control = nls.control(warnOnly = TRUE)), silent = TRUE)
        if(is.null(fit.top)) { ## NLS failed, do something
            print("error")
            errors.top <<- c(errors.top, levels(x$sdpclass))
        }
        pred.top <- predict(fit.top, newdata = x)
        rgr.top <- x[,"bvgrowth"]/pred.top
        res.top <- pred.top-x$bvgrowth
    }
    data.frame(install = mean(x$install), plot = mean(x$plot), id = x$id,
    rgr.top = rgr.top)
})


    ##                                     # ** Calculations done **
    ##                                     # ** interactively show graphs **
    ## plot <- x$pplot[1]
    ## year <- x$time[1]
    ## par(ask = TRUE, mfrow = c(2,2) )
    ## scatter.smooth(x$priorbv, res.top, lpars = list(col="red"),
    ##                main = paste("Residuals vs Fitted",plot,year))
    ## abline(h = 0, lty = 2)
    ## plot(x$priorbv, x$bvgrowth, main = "BV growth vs. Prior BV with Predicted")
    ## curve(coef(fit.reg)[[1]] * x ^ coef(fit.reg)[[2]], add=TRUE, col="blue", lwd=2)
    ## curve(coef(fit.top)[[1]] * x ^ coef(fit.top)[[2]], add=TRUE, col="green", lwd=2)
    ## hist(res.top, main = "Distribution of the residuals")
    ## abline(v = 1, lty = 2)
    ## if(length(x$numnebs[!is.na(x$numnebs)])>2) {
    ##     plot(x$numnebs, rgr.top, main = "rGR vs Neighbor Density")
    ##     abline(h=1, lty=2);
    ##     tryCatch(abline(lm(rgr.top ~ x$numnebs), col="blue"),
    ##              error=function(e) NULL )
    ## }
    ##         data.frame(pplot=x$pplot, time=x$time, id=x$id, rgr.top = rgr.top)}
    ##                     error=function(e) print(x$pplot))
    ##           }))

## combine new rgr
rgr <- rgr[,c("install","plot","id","rgr.top")]
rgr <- rgr[order(rgr$install, rgr$plot, rgr$id),]
bclong <- bclong[order(bclong$install, bclong$plot, bclong$id),]
bclong$rgrsisdp <- rgr$rgr.top
## tst <- merge(bclong, rgr, by = c("install","plot","id"), all.x = TRUE)
## any rgr >1 call 1 for now
bclong[bclong$rgrsisdp>1.5,"rgrsisdp"] <- 1.5
plot(bclong$priorbv, bclong$rgrsisdp)

## run MLE using rgrplot
## create new starting parameters
## pars <- read.csv("parameters.csv")
## rowcopy = pars[83:86,]
## rowcopy$model <- "simplest"
## rowcopy$mdep.var <- "rgrplot"
## pars <- rbind(pars, rowcopy)
## write.csv(pars, "parameters.csv", row.names =FALSE)

# specify parameters for MLE fitting, using 'simplest' model
source("neighborhood-models.R")
dat = bclong
sr = c(6)
spec = "FD"
ind.var = "priorbv" ## will be used in neighborhood calculations
dep.var = "rgrsisdp"
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

# graph rgr vs nci, bvgrowth vs nci
par(mfrow = c(1,2))
plot(targets$nci, targets[,dep.var])
plot(targets$nci, targets$bvgrowth)

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

write.csv(bclong, "bclong-rgr.csv", row.names = FALSE)

###################################################################################
## testing
tst <- subset(bclong, pplot == 1.1)
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

