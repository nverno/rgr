## Looking at top performers by a combination of SI and SDP (as a class)
##  calculate top performers per Install / Plot / Period
source("~/work/functions/functions.R")
source("~/work/functions/functions-hazards.R")

## data - use long-bc-derived, all species included
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")
bclong <- subset(bclong1, stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees
## not marked with status
bclong$pplot <- factor(bclong$pplot)
bclong$sdpclass <- factor(as.numeric(bclong$sdpclass))

## sample size by install / plot / period
samp <- ddply(bclong, .(install, plot, time), function(x) nrow(x))
range(samp$V1)

##  fit using power function
##  rgr = obs growth/predicted growth
yrs <- c(76,79,82,85,91,97)
bclong$bvcl <- rep(NA, nrow(bclong))
bclong$bvgrcl <- rep(NA, nrow(bclong))
errors <- c()
correlated <- c()

## test run, make sure starting parameters are decent for nls fit
samps <- ddply(bclong, .(install, plot, time), function(x) nrow(x) )
tst <- subset(bclong, install == 1 & plot == 10 & time == 79)
fit.reg <- nls(bvgrowth ~ a*priorbv^b, start=list(a=.05,b=.5), data=tst)

rgr <- ddply(bclong, .(install, plot, time), .fun = function(x) {
    x <- droplevels(x)
    fit.nls <- NULL ## Try NLS fit
    try(fit.nls <- nls((bvgrowth) ~ a*priorbv^b, start=list(a=.05,b=.5), data=x,
                       control = nls.control(warnOnly = TRUE)))
    if(is.null(fit.nls)) { ## NLS failed, report problem data
        print(paste(unique(tst$install), unique(tst$plot)))
    }
    else { ## NLS fit successful
        ## if(length(x[,"bvgrowth"])!=length(predict(fit.nls)))
        ##     print(paste(levels(x$sdpclass), mean(x$si)))
        rgr.nls <- x[,"bvgrowth"]/predict(fit.nls)
        rgr.nls[rgr.nls>10] <- 10
    }
    ## Linear Fit and correlations
    fit.lin <- lm(bvgrowth ~ priorbv, data = x)
    rgr.lin <- x[,"bvgrowth"]/predict(fit.lin); rgr.lin[rgr.lin<0] <- 0
    rgr.lin[rgr.lin>10] <- 10
    lin.cor <- cor.test(rgr.lin, x$bv)
    nonlin.cor <- cor.test(x[,"bvgrowth"]/predict(fit.nls), x$bv)
    ## Choose fit type by least significant correlation
    ifelse (lin.cor$p.value > nonlin.cor$p.value,
           { rgr <- rgr.lin; fit.type = "linear" },
           { rgr <- rgr.nls; fit.type = "power" })
    ## if there is significant correlation under best fit, report
    if (cor.test(rgr, x[, "bv"])$p.value < 0.05)
        correlated <<- c(correlated, paste("Install ", unique(x$install),
                                          ", Plot ", unique(x$plot),
                                          ", Year ", unique(x$time)))
    ## Results => data.frame
    data.frame(install = x$install, plot = x$plot, id = x$id,
               nls.a = rep(coef(fit.nls)[[1]], nrow(x)),
               nls.b = rep(coef(fit.nls)[[2]], nrow(x)),
               rgr = rgr,
               method = fit.type,
               rgr.power = rgr.nls, rgr.lin = rgr.lin,
               lm.int = rep(coef(fit.lin)[[1]], nrow(x)),
               lm.slope = rep(coef(fit.lin)[[2]], nrow(x)),
               lin.cor = rep(lin.cor$p.value, nrow(x)),
               nonlin.cor = rep(nonlin.cor$p.value, nrow(x)))
})

## write rgr data.frame for graphical analysis
write.csv(rgr, "~/work/data/data/rgr-parameters-forhazard.csv", row.names = FALSE)

## combine new rgr
rgr <- rgr[,c("install","plot","id","rgr","method","rgr.power","rgr.lin")]
rgr <- rgr[order(rgr$install, rgr$plot, rgr$id),]
bclong <- bclong[order(bclong$install, bclong$plot, bclong$id),]
bclong$hazRgr.power <- rgr$rgr.power
bclong$hazRgr.lin <- rgr$rgr.lin
bclong$hazRgr <- rgr$rgr

## If set max rgr value to 10, plot rgr vs. priorbv
bclong[bclong$hazRgr.power > 10,"hazRgr.power"] <- 10
pdf(file = "~/Allometry/hazards/rgrVSpriorbv-power.pdf")
plot(bclong$priorbv, bclong$hazRgr.power, main = "rGR (power fit) vs Prior Bole Volume, rGR
was calculated by install / plot / year.  Max rGR set at 10.",
     xlab = "Prior Bole Volume", ylab = "rGR (power Fit)")
abline(h = 1, lty = 2, col = "cyan", lwd = 2)
abline(h = median(bclong$hazRgr.power), col = "red", lwd = 2, lty = 2)
abline(h = mean(bclong$hazRgr.power), col = "green", lwd = 2, lty = 2)
text(x = 15, y = c(6,7), c(paste0("Median = ", round(median(bclong$hazRgr.power), 2)),
             paste("Mean = ",round(mean(bclong$hazRgr.power), 2))), col = c("red","green"))
dev.off()

## Linear fits, negative rgrs => 0
bclong[bclong$hazRgr.lin < 0,"hazRgr.lin"] <- 0
pdf(file = "~/Allometry/hazards/rgrVSpriorbv-lin.pdf")
plot(bclong$priorbv, bclong$hazRgr.lin, main = "rGR (lin fit) vs Prior Bole Volume, rGR
was calculated by install / plot / year.  Max rGR set at 10.",
     xlab = "Prior Bole Volume", ylab = "rGR (lin Fit)")
abline(h = 1, lty = 2, col = "cyan", lwd = 2)
abline(h = median(bclong$hazRgr.lin), col = "red", lwd = 2, lty = 2)
abline(h = mean(bclong$hazRgr.lin), col = "green", lwd = 2, lty = 2)
text(x = 15, y = c(6,7), c(paste0("Median = ", round(median(bclong$hazRgr.lin), 2)),
             paste("Mean = ",round(mean(bclong$hazRgr.lin), 2))), col = c("red","green"))
dev.off()

## rgr final fits
bclong[bclong$hazRgr < 0,"hazRgr"] <- 0
pdf(file = "~/Allometry/hazards/rgrVSpriorbv-best.pdf")
plot(bclong$priorbv, bclong$hazRgr, main = "rGR (best fit) vs Prior Bole Volume, rGR
was calculated by install / plot / year.  Max rGR set at 10.",
     xlab = "Prior Bole Volume", ylab = "rGR (best Fit)")
abline(h = 1, lty = 2, col = "cyan", lwd = 2)
abline(h = median(bclong$hazRgr), col = "red", lwd = 2, lty = 2)
abline(h = mean(bclong$hazRgr), col = "green", lwd = 2, lty = 2)
text(x = 15, y = c(6,7), c(paste0("Median = ", round(median(bclong$hazRgr), 2)),
             paste("Mean = ",round(mean(bclong$hazRgr), 2))), col = c("red","green"))
dev.off()

## Add new rgr values to original bclong1 (after testing)
bclong1 <- bclong1[order(bclong1$install,bclong1$plot,bclong1$id),]
bclong <- bclong[order(bclong$install,bclong$plot,bclong$id),]
bclong1$hazRgr.power <- rep(NA, nrow(bclong1))
bclong1$hazRgr.lin <- rep(NA, nrow(bclong1))
bclong1$hazRgr <- rep(NA, nrow(bclong1))
toadd <- c("hazRgr.power", "hazRgr.lin", "hazRgr")
bclong1[bclong1$stat=="ALIVE" & !is.na(bclong1$bvgrowth) & bclong1$bvgrowth>=0,toadd] <-
    bclong[,toadd]

## Write data to disk?
write.csv(bclong1, "~/Allometry/data/hazards-bc.csv", row.names = FALSE)

