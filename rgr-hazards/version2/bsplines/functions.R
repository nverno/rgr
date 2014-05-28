##########################################################################
##
##              Polynomial regression with B-splines
##
## Returns list of fits from 1:degree
##
## Parameters:
## - ind: independent variable at start of period (i.e. priorbv)
## - ind2: independent variable at end of period (i.e. bv)
## - dep: dependent variable (i.e bvgrowth across periods)
## - degree: range of degrees of freedom allowed (ie. if degree=3 then
##           polynomial fits of degree 1:3 will be tried and the best chosen)
## - corr: minimum correlation test pvalue allowed (ie. alpha = 0.05)
## - intrcpt: determines if polynomials should include intercept
## - info: If true, returns information about the fits, described
##         in the helper function, bestSpline
## - debug: print debugging info during fitting
##
##########################################################################
##
require(splines)
fitSplines <- function(ind, ind2, dep, degree=9, corr=0.05, intrcpt=FALSE,
                       info = FALSE, debug = FALSE) {
    fits <- list() ## store all fits to do analysis at end
    nonSig <- c() ## store number of insignificant coefs in each fit
    for (i in 1:degree) {
        if (debug==TRUE)
            print(paste("Fitting degree",i))
        fit <- NULL
        if (intrcpt) {              # fit with intercept
            try(fit <- lm(dep ~ bs(ind, degree = i)), silent=TRUE)
        } else {                    # fit without intercept
            try(fit <- lm(dep ~ bs(ind, degree = i) - 1), silent=TRUE)
        }
        if (!is.null(fit)) {        # Successful fit, check significance of coefs
            summ <- summary(fit)$coefficients[,4]
            nonSig <- length(summ[summ > corr])
            rmse <- sqrt(sum(residuals(fit)^2))
            if (debug==TRUE)
                print(paste("rmse of", rmse,", AIC:",AIC(fit)))
            fits[i] <- list(fit)
        }
        if (is.null(fit)) break;    # stop fitting if one fails
    }
    if (info==TRUE)
        return (bestSpline(fits, ind2, dep, corr))
    else return (fits[bestSpline(fits, ind2, dep, corr)[["best"]]])
}


################################################################################
##
##                      Helper function for fitSplines
##
## Takes a list of polynomial fits and returns information about them, ranking
##  them by rmse, aic, loglikelihood and checking to see if the residuals
##  after fitting are correlated with the independent variable.
##
## Parameters:
## - ind2: independent variable at end of period (i.e. bv)
## - dep: dependent variable (i.e bvgrowth b/w priorbv and bv)
## - corr: minimum correlation test pvalue allowed (ie. alpha = 0.05)
## - fits: list of polynomial fits to compare
##
## Output:
## - best: degree of best fitting polynomial that is also uncorrelated with ind2
## - smallest: smallest degree with non-significant correlation
## - degrees: degrees of fits
## - corrs: pvalues of correlation tests between predicted values and the residuals
## - rmses: root mean square errors
## - logLiks: log likelihoods
## - aics: AIC of fits
##
################################################################################
bestSpline <- function(fits, ind2, dep, corr=0.05) {
    aics <- sapply(fits, AIC)
    logLiks <- sapply(fits, logLik)
    rmses <- sapply(fits, function(x) sqrt(sum(residuals(x)^2)/length(x$fitted)))
    corrs <- sapply(fits, function(x) cor.test(dep - predict(x), ind2)$p.value)
    ifelse (length(which(corrs > corr)) > 0,
            { smallest <- min(which(corrs > corr)) },
            { smallest <- which(corrs == max(corrs)) })
    current <- best <- smallest                           # currently the best
    ## Test smallest against larger for significantly better models
    candobetter <- TRUE
    while (candobetter) {
        if (current < length(fits)) {
            others <- fits[(current+1):length(fits)]
            pvals <- sapply(others, function(i) anova(fits[[current]], i)["Pr(>F)"][2,])
            if (length(which(pvals < corr)) > 0) {
                possible <- min(which(pvals < corr)) + current
                ifelse (corrs[possible] > corr,
                    { current <- best <- possible },     # better model must also be uncorrelated
                    { candobetter <- FALSE })             # update best fit
            } else { candobetter <- FALSE }
        } else { candobetter <- FALSE }
    }
    return(list(best=best, smallest=smallest, degrees=c(1:length(fits)),
                corrs=corrs, rmses=rmses, logLiks=logLiks, aics=aics))
}

