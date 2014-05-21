## Interactive Graphics to see rGR fitting by SI and SDP
## rgr parameters dataset has fitted values, coefficients,
## bclong dataset is used to get numbers of neighbors
source("~/work/functions/functions.R")
source("~/work/functions/functions-growth.R")

bclong <- read.csv("~/work/data/data/long-bc-derived.csv")
bclong <- subset(bclong, stat=="ALIVE" & bvgrowth >= 0)
rgr <- read.csv("~/work/data/data/rgr-parameters-sisdp.csv")
rgr <- rgr[order(rgr$install, rgr$plot, rgr$id),]

## Neighbor density in radius of 4 meters
ind.var = "priorbv"
fit.MLE.models(bclong, sr=6, spec="FD", ind.var = ind.var, dep.var = "bvgrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)

## add numnebs, bv stuff to rgr for graphics
bclong <- merge(bclong, targets, all=TRUE)
bclong <- bclong[order(bclong$install, bclong$plot, bclong$id),]
toadd <- c("priorbv","bvgrowth","numnebs")
rgr[,toadd] <- bclong[,toadd]
rgr$sdpclass <- factor(rgr$sdpclass)

## Subset by combinations of sdpclass and si, interactively produce graphics
##  showing residuals to rgr.top, rgr.top vs number of neighbors,
##  fitted line on untransformed data,
bclong$sdpclass <- factor(bclong$sdpclass)

tst <- ddply(bclong, .(sdpclass), .fun = function(x) {
    x <- droplevels(x)
    data.frame =c(
    sdps <- nrow(x)
    )})
bclong$sdpclass <- rep(c(1:3), times = tst[,2])

sdps <- as.numeric(levels(factor(bclong$sdpclass)))
sis <- as.numeric(levels(factor(bclong$si)))
combs <- expand.grid(sdpclass = sdps, si = sis)
null.subs <- c()

## Print the results to file
printout = TRUE;

for(i in 1:nrow(combs)) {
    par(ask = TRUE, mfrow = c(2,2))
    x <- droplevels(subset(rgr, sdpclass == combs[i,"sdpclass"] &
                           si == combs[i,"si"]))
    print(paste(levels(x$sdpclass), mean(x$si)))
    if(nrow(x) == 0) { ## skip NULL subsets, store which were null and print info
        print("Null subset");
        null.subs <- c(null.subs, paste(combs[i,"sdpclass"],combs[i,"si"],sep = "."))
        next;
    }
    ## Print out the image to disk
    if (printout == TRUE) {
        pdf(file = paste0("~/work/rgr/rgr-growth/visuals/rgrPics/",unique(x$install),".",
            unique(x$plot),".",unique(x$time),".pdf"))
        par(mfrow = c(2,2))
    }
    ## get the coefficients
    reg.a <- mean(x$reg.a)
    reg.b <- mean(x$reg.b)
    top.a <- mean(x$top.a)
    top.b <- mean(x$top.b)
    ## Start plotting stuff
    scatter.smooth(x$priorbv, x$res.top, lpars = list(col="red"),
                   main = paste("Residuals vs Fitted",levels(x$sdpclass),mean(x$si)))
    abline(h = 0, lty = 2)
    plot(x$priorbv, x$bvgrowth, main = "BV growth vs. Prior BV with Predicted")
 ##   curve(reg.a * x ^ reg.b, add=TRUE, col="blue", lwd=2)
    curve(top.a * x ^ top.b, add=TRUE, col="green", lwd=2)
    hist(x$res.top, main = "Distribution of the residuals")
    abline(v = 1, lty = 2)
    if(length(x$numnebs[!is.na(x$numnebs)])>2) {
        plot(x$numnebs, x$rgr.top, main = "rGR vs Neighbor Density")
    abline(h=1, lty=2);
    tryCatch(abline(lm(x$rgr.top ~ x$numnebs), col="blue"),
             error=function(e) NULL )
    }
    if (printout == TRUE)
        dev.off()
}


