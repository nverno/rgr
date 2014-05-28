################################################################################
##
##           Interactive Graphics to see rGR fitting for hazards
##
## - To get the parameters, just run fitSplines function with appropriate
##   variables.  The degree of the fit used is stored in "hazard-firs-bc.csv"
##
################################################################################
## Currently all the fits are B-spline polynomials
source("~/work/rgr/rgr-hazards/version2/bsplines/functions.R")

## data
dat <- read.csv("~/work/data/data/hazards/hazards-bc-firs.csv")
combs <- unique(dat[c("install","plot","time")]) # unique combos of install/plot/time

## Example: returns fit object
#tst <- dat[dat$install == 1 & dat$plot == 10 & dat$time == 85,]
fitSplines(ind=tst$priorbv, ind2=tst$bv, dep=tst$bvgrowth)

printout <- FALSE
null.subs <- c() # track empty subsets
for (i in 1:nrow(combs)) {
    par(ask = TRUE, mfrow = c(2,2))
    x <- droplevels(subset(dat, install == combs[i,"install"] &
                           plot == combs[i,"plot"] & time == combs[i,"time"]))
    print(paste("Install: ", unique(x$install), ", Plot: ", unique(x$plot),
                ", Year: ",unique(x$time)))
    if(nrow(x) == 0) { ## skip NULL subsets, store which were null and print info
        print("Null subset");
        null.subs <- c(null.subs, paste(combs[i,"sdpclass"],combs[i,"si"],sep = "."))
        next;
    }
    ## Get the model
    mod <- fitSplines(ind=x$priorbv, ind2=x$bv, dep=x$bvgrowth)[[1]]
    deg <- length(coef(mod))
    preds <- predict(mod)
    res <- residuals(mod)
    stopifnot(deg == unique(x$degree))

    ## Get the correlation between predicted and bv
    corr <- cor(x$bv, res)
    pval <- cor.test(x$bv, res)$p.value

    ## Print to file if desired
    ## if (printout == TRUE) {
    ##     pdf(file = paste0("~/Allometry/hazards/rgrPics/",unique(x$install),".",
    ##         unique(x$plot),".",unique(x$time),".pdf"))
    ##     par(mfrow = c(2,2))
    ## }

    ## *** Start plotting stuff  ***
    ## plot 1: points + predicted line/points
    predsort <- data.frame(priorbv=x$priorbv, preds=preds)
    predsort <- predsort[order(predsort$priorbv),]
    plot(x$priorbv, x$bvgrowth, pch = 2, main = "BV growth vs Prior BV with predicted points/line",
         xlab = "Prior BV", ylab = "BV growth")
    lines(predsort$priorbv, predsort$preds, col = "dark red", type = "o", pch = 3)

    ## plot 2: residuals scatter
    plot(x$priorbv, res, pch = 2)
    abline(h = 0, col = "dark blue", lty = 2)

    ## plot 3: residuals hist
    hist(res)
    abline(v = 0, lty = 2, col ="dark red", lwd = 3)

    ## plot 4: correlation info
    plot(1:10, 1:10, type = "n")
    degr <- paste("Degree: ", deg)
    corinfo <- paste("Correlation: ", format(corr, digits = 3))
    corp <- paste("P-value: ", format(pval, scientific = TRUE, digits = 3))
    text(c(5,5,5), c(8,4,1), labels = c(corinfo, corp, degr), cex = 1.2)

    ## if (printout == TRUE)
    ##     dev.off()
}
