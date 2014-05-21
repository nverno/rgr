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
    model <- NULL
    deg <- 11
    while (is.null(model) & deg > 0) {
        model <- try(removeCorr(x, degree = deg), silent = TRUE)
        if (is.null(model)) deg <- deg - 1
    }
    if (!is.null(model)) {
        switch(model$model,
               "lin" = {
                   deg = 0
                   fit <- lm(x$bvgrowth ~ x$priorbv)
                   rmse <- sqrt(sum(residuals(fit)^2))
                   rgr <- x$bvgrowth/predict(fit)
                   corr <- cor.test(rgr, x$bv)$p.value
               },
               "pow" = {
                   deg = 0
                   fit <- nls(x$bvgrowth ~ a*x$priorbv^b, start = list(a=0.5,b=0.5))
                   rmse <- sqrt(sum(residuals(fit)^2))
                   rgr <- x$bvgrowth/predict(fit)
                   corr <- cor.test(rgr, x$bv)$p.value
               },
               "poly" = {
                   deg = model$degree
                   fit <- lm(x$bvgrowth ~ poly(x$priorbv, model$degree))
                   rmse <- sqrt(sum(residuals(fit)^2))
                   rgr <- x$bvgrowth/predict(fit)
                   corr <- cor.test(rgr, x$bv)$p.value
               })
    }
    ## Results => data.frame
    data.frame(install = x$install, plot = x$plot, id = x$id,
               rgr = rgr,
               method = rep(model$model, nrow(x)),
               corr = rep(corr, nrow(x)),
               degree = rep(deg, nrow(x)))
})

