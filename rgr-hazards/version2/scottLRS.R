## Compute Scotts's version of LRS
## Questions:
## - Is the predicted fit just a linear fit?

## Returns LRS given prior and current size values (i.e. priorbv and bv)
## Here a linear fit is used to predict the estimated growth between periods
## LRS = (priorbv + predicted bvgrowth) / (max bv at end of period)
scottLRS <- function(prior, current) {
    predicted <- lm(current ~ prior)
    LRS = (prior + predicted)/max(current)
    return(LRS)
}

## Test
tst <- bclong[bclong$install == 1 & bclong$plot == 10,]
lrs <- scottLRS(prior=tst$prio)
