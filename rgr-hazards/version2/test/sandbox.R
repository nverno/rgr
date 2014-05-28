################################################################################
##
##                                  Set-up
##
################################################################################

## data
bclong1 <- read.csv("~/work/data/data/long-bc-derived.csv")

## Test data, only look at ALIVE firs that grew
bclong <- subset(bclong1, spec == "FD" & stat == "ALIVE" & bvgrowth>=0) ## not sure how to treat trees
tst <- bclong[bclong$install == 1 & bclong$plot == 10,]
ind <- tst[,"priorbv"]
ind2 <- tst[,"bv"]
dep <- tst[,"bvgrowth"]
degree <- 9
## dat <- tst

## params
debug <- TRUE
corr <- 0.05
intrcpt <- FALSE
info <- TRUE

##########################################################################
##
##              Polynomial regression with B-splines
##
## Parameters:
## - ind: independent variable from prior period (i.e. priorbv)
## - ind2: independent variable at start of period (i.e. bv)
## - dep: dependent variable (i.e bvgrowth)
## - degree: range of degrees of freedom allowed (ie. if degree=3 then
##           polynomial fits of degree 1:3 will be tried and the best chosen)
## - corr: minimum correlation test pvalue allowed (ie. alpha = 0.05)
## - intrcpt: determines if polynomials should include intercept
## - info: If true, returns the following fitting information:
##     * best: degree of best fitting polynomial
##     * smallest:
## - debug: print debugging info
##########################################################################
##
fitSplines <- function(ind, ind2, dep, degree=9, corr=0.05, intrcpt=FALSE,
                       info = FALSE, debug = FALSE) {
    fits <- list() ## store all fits to do analysis at end
    nonSig <- c() ## store number of insignificant coefs in each fit
    for (i in 1:degree) {
        if (debug==TRUE)
            print(paste("Fitting degree",i))
        fit <- NULL
        if (intrcpt) { # fit without intercept
            try(fit <- lm(dep ~ bs(ind, degree = i) + 0), silent=TRUE)
        } else {
            try(fit <- lm(dep ~ bs(ind, degree = i)), silent=TRUE)
        }
        if (!is.null(fit)) { ## Successful fit, check significance of coefs
            summ <- summary(fit)$coefficients[,4]
            nonSig <- length(summ[summ > corr])
            rmse <- sqrt(sum(residuals(fit)^2))
            if (debug==TRUE)
                print(paste("rmse of", rmse,", AIC:",AIC(fit)))
            fits[i] <- list(fit)
        }
        if (is.null(fit)) break; ## stop fitting if one fails
    }
    if (info==TRUE)
        return (bestSpline(fits, ind2, dep, corr))
    else return (fits[bestSpline(fits, ind2, dep, corr)[["best"]]])
}


## Helper function for fitSplines
## Takes a list of fits and returns smallest degree that is uncorrelated to bv
## Also returns the AIC and rmse of the best one
bestSpline <- function(fits, ind2, dep, corr=0.05) {
    aics <- sapply(fits, AIC)
    logLiks <- sapply(fits, logLik)
    rmses <- sapply(fits, function(x) sqrt(sum(residuals(x)^2)/length(x$fitted)))
    corrs <- sapply(fits, function(x) cor.test(dep/predict(x), ind2)$p.value)
    ifelse (length(which(corrs > corr)) > 0,
            { smallest <- min(which(corrs > corr)) },
            { smallest <- which(corrs == max(corrs)) })
    best <- smallest # currently the best
    ## Test smallest against larger for significantly better models
    if (smallest < length(fits)) {
        others <- fits[(smallest+1):length(fits)]
        pvals <- sapply(others, function(x) anova(fits[[smallest]], x)$Pr[2])
        if (length(which(pvals < corr)) > 0) {
            possible <- min(which(pvals < corr)) + smallest
            if (corrs[possible] > corr) ## better model must also be uncorrelated
                best <- possible  ## update best fit
        }
    }
    return(list(best=best, smallest=smallest, degrees=c(1:length(fits)),
                corrs=corrs, rmses=rmses, logLiks=logLiks, aics=aics))
}


################################################################################
##
##                                  Sandbox
##
################################################################################

## with info
fitSplines(ind, ind2, dep, degree = 3, info=TRUE, debug = TRUE, corr = 0.05)

## w/o info
fitSplines(ind, ind2, dep, degree = 2, info=FALSE, debug = TRUE, corr = 0.05)


## correlation between residuals and priorbv
fit <- lm((tst$bv - tst$priorbv) ~ tst$priorbv)
plot()
preds <- predict(fit, data.frame())
