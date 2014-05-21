## Interactive Graphics to see rGR fitting by SI and SDP
## rgr parameters dataset has fitted values, coefficients,
## bclong dataset is used to get numbers of neighbors
source("~/Allometry/functions.R")
bclong <- read.csv("~/Allometry/data/hazards-bc.csv")
bclong <- subset(bclong, stat=="ALIVE" & bvgrowth >= 0)
rgr <- read.csv("~/Allometry/data/rgr-parameters-forhazard.csv")
rgr <- rgr[order(rgr$install, rgr$plot, rgr$id),]

## Neighbor density in radius of 6 meters
ind.var = "priorbv"
fit.MLE.models(bclong, sr=6, spec="FD", ind.var = ind.var, dep.var = "bvgrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)

## add numnebs, bv stuff to rgr for graphics
bclong <- merge(bclong, targets, all=TRUE)
bclong <- bclong[order(bclong$install, bclong$plot, bclong$id),]
toadd <- c("priorbv","bvgrowth","sdpclass","numnebs")
rgr[,toadd] <- bclong[,toadd]
rgr$sdpclass <- factor(rgr$sdpclass)

## Subset by combinations install / plot /year, interactively produce graphics
##  showing residuals rgr vs number of neighbors,
##  fitted line on untransformed data,
bclong$sdpclass <- factor(bclong$sdpclass)
combs <- unique(bclong[c("install","plot","time")])
null.subs <- c()
printout = FALSE

for(i in 1:nrow(combs)) {
    par(ask = TRUE, mfrow = c(2,2))
    x <- droplevels(subset(rgr, install == combs[i,"install"] &
                           plot == combs[i,"plot"] & time == combs[i,"time"]))
    print(paste("Install: ", unique(x$install), ", Plot: ", unique(x$plot),
                ", Year: ",unique(x$time)))
    if(nrow(x) == 0) { ## skip NULL subsets, store which were null and print info
        print("Null subset");
        null.subs <- c(null.subs, paste(combs[i,"sdpclass"],combs[i,"si"],sep = "."))
        next;
    }
    ## get the coefficients
    nls.a = unique(x$nls.a); nls.b = unique(x$nls.b);
    lm.int = unique(x$lm.int); lm.slope = unique(x$lm.slope);
    if (length(nls.a) > 1 | length(nls.b) > 1) stop("Wrong number of parameters (a,b)")
    lin.cor <- unique(x$lin.cor)
    nonlin.cor <- unique(x$nonlin.cor)
    ## Start plotting stuff
    nls.pred <- nls.a * x$priorbv ^ nls.b
    nls.res <- nls.pred - x$bvgrowth
    lin.pred <- lm.int + lm.slope*x$priorbv
    lin.res <- lin.pred - x$bvgrowth
    ## Choose fit with lowest correlation between rgr and priorbv
    ifelse(abs(cor(x$rgr.lin, x$priorbv)) < abs(cor(x$rgr.power, x$priorbv)),
           { low.res <- lin.res; low.rgr <- x$rgr.lin; used <- "linear fit"},
           { low.res <- nls.res; low.rgr <- x$rgr.power; used <- "power fit"} )
    if (printout == TRUE) {
        pdf(file = paste0("~/Allometry/hazards/rgrPics/",unique(x$install),".",
            unique(x$plot),".",unique(x$time),".pdf"))
        par(mfrow = c(2,2))
    }
    scatter.smooth(x$priorbv, low.res, lpars = list(col="red"),
                   main = paste("Residuals vs Fitted, ", used,
                   "\nInstall: ", unique(x$install),
                   ", Plot: ", unique(x$plot), ", Year: ", unique(x$time)))
    abline(h = 0, lty = 2)
    plot(x$priorbv, x$bvgrowth, main = "BV growth vs. Prior BV with Predicted")
    curve(lm.int + x*lm.slope, add=TRUE, col="blue", lwd=2)
    curve(nls.a * x ^ nls.b, add=TRUE, col="green", lwd=2)
    ## histogram
    hist(low.res, main = paste("Distribution of the residuals, ", used))
    abline(v = 0, lty = 2, lwd = 3, col = "red")
    if(length(x$numnebs[!is.na(x$numnebs)])>2) {
        plot(x$numnebs, low.rgr, main = "rGR vs Neighbor Density")
    abline(h=1, lty=2);
    tryCatch(abline(lm(low.rgr ~ x$numnebs), col="blue"),
             error=function(e) NULL )
    }
    if (printout == TRUE)
        dev.off()
}

tst <- lm(bvgrowth ~ priorbv, data = x)
