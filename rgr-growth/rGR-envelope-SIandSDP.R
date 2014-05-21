## Looking at top performers by a combination of SI and SDP (as a class)
##  calculate top performers per Install / Plot / Period
source("functions.R")

## data - use long-bc-derived, all species included
bclong1 <- read.csv("~/Allometry/data/long-bc-derived.csv")
bclong <- subset(bclong1, stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees
## not marked with status
bclong$pplot <- factor(bclong$pplot)
bclong$sdpclass <- factor(as.numeric(bclong$sdpclass))

## sample size by install / plot / period
samp <- ddply(bclong, .(install, plot, time), function(x) nrow(x))
range(samp$V1)

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
    data.frame(install = x$install, plot = x$plot, id = x$id,
    rgr.top = rgr.top, res.top = res.top, pred.top = pred.top, reg.a = coef(fit.reg)[[1]],
               reg.b = coef(fit.reg)[[2]], top.a = coef(fit.top)[[1]],
               top.b = coef(fit.top)[[2]])
})

## write rgr data.frame for graphical analysis
write.csv(rgr, "data/rgr-parameters-sisdp.csv", row.names = FALSE)

## combine new rgr
rgr <- rgr[,c("install","plot","id","rgr.top")]
rgr <- rgr[order(rgr$install, rgr$plot, rgr$id),]
bclong <- bclong[order(bclong$install, bclong$plot, bclong$id),]
bclong$rgrsisdp <- rgr$rgr.top

## tst <- merge(bclong, rgr, by = c("install","plot","id"), all.x = TRUE)
## any rgr >1 call 1 for now
bclong[bclong$rgrsisdp>1.5,"rgrsisdp"] <- 1.5
plot(bclong$priorbv, bclong$rgrsisdp)

## Add new rgr values to original bclong1 (after testing)
bclong1 <- bclong1[order(bclong1$install,bclong1$plot,bclong1$id),]
bclong <- bclong[order(bclong$install,bclong$plot,bclong$id),]
bclong1$rgrsisdp <- rep(NA, nrow(bclong1))
toadd <- c("rgrsisdp")
bclong1[bclong1$stat=="ALIVE" & !is.na(bclong1$bvgrowth) & bclong1$bvgrowth>=0,toadd] <-
    bclong[,toadd]

## Write data to disk?
write.csv(bclong1, "/data/long-bc-derived.csv", row.names = FALSE)
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

