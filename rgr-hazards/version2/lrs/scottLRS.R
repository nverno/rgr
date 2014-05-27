## Compute Scotts's version of LRS
## Questions:
## - Is the predicted fit just a linear fit?

## Returns LRS given prior and current size values (i.e. priorbv and bv)
## Options for prediction are: linear, power
## LRS = (priorbv + predicted bvgrowth) / (max bv at end of period)
scottLRS <- function(prior, current, model = "linear") {
    fit <- lm((current-prior) ~ prior)
    if (model == "power")
        fit <- nls((current-prior) ~ a*prior^b, start = list(a=1, b=1))
    LRS = (prior + predict(fit, data.frame(prior=prior, current=current))) /
        max(current, na.rm=T)
    return( LRS )
}


## Read data
library(plyr)
bclong <- read.csv("~/work/data/data/long-bc-derived.csv")

## Fit lrs values
bclong$lrs <- rep(NA, nrow(bclong))
lrsvals <- ddply(bclong[bclong$time != 73, ], .(install, plot, time), function(x) {
    data.frame(id = x$id, lrs = scottLRS(prior = x$priorbv, current = x$bv, model = "linear"))
})

## add lrs values
lrs <- lrsvals[order(lrsvals$install, lrsvals$plot, lrsvals$time, lrsvals$id),]
bclong <- bclong[order(bclong$install, bclong$plot, bclong$time, bclong$id),]
bclong[bclong$time != 73,]$lrs <- lrs$lrs

## write data
write.csv(bclong, "~/work/data/data/long-bc-derived.csv", row.names=FALSE)
