## Looking at top performers by a combination of SI and SDP (as a class)
##  calculate top performers per Install / Plot / Period
## Note: the LRS value calculated here is WRONG, it should be calculated as
##  individuals size relative to the largets tree in the plot including ALL species.
##  It is correctly calculated in the subsequent script that preps this output data
##  for hazard analysis (~/work/hazards/addColumn/lrs.R)
## source("~/work/functions/functions.R")
source("~/work/rgr/rgr-hazards/version2/functions.R")
require(plyr)

## data - use long-bc-derived, all species included
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")
## Remove outlier from install 10, plot 2, year 97
bclong <- subset(bclong1, spec == "FD" & stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees

## not marked with status
bclong$pplot <- factor(bclong$pplot)
bclong$sdpclass <- factor(as.numeric(bclong$sdpclass))

## sample size by install / plot / period
## samp <- ddply(bclong, .(install, plot, time), function(x) nrow(x))
## range(samp$V1)

## Track correlated plots
corrplots <- data.frame(install=NULL,plot=NULL,time=NULL)

## calculate rGRs
rgr <- ddply(bclong, .(install, plot, time), .fun = function(x) {
    x <- droplevels(x)
    maxSize <- max(x$bv)
    print(paste(unique(x$install),unique(x$plot), unique(x$time)))
    rgrs <- removeCorr(x, degree = 3)
    if (unique(rgrs$corr < 0.05)) # report correlated plot
        corrplots <<- rbind(data.frame(install = unique(x$install), plot = unique(x$plot),
                                       time = unique(x$time)),corrplots)
    ## Results => data.frame
    data.frame(install = x$install, plot = x$plot, id = x$id,
               rgr = rgrs$rgr,
               method = rgrs$model,
               corr = rgrs$corr,
               degree = rgrs$degree,
               LRS = x$bv/maxSize)
})

## Check how many plots are significantly correlated
## nrow(subset(rgr, corr < 0.05)) ## 87, 2 plots

## Store correlated plots: currently none
## write.csv(corrplots, "~/work/data/data/hazards/corrplots.csv", row.names=FALSE)

## add rgr column to bclong and make a dataset to run hazards on
rgr <- rgr[order(rgr$install, rgr$plot, rgr$time, rgr$id),]
bc <- bclong[order(bclong$install, bclong$plot, bclong$time, bclong$id),]
bc[,c("rgr","method","degree")] <- rgr[,c("rgr","method","degree")]

## store hazards dataset
write.csv(bc, "~/work/data/data/hazards/hazards-bc-firs.csv", row.names=FALSE)
