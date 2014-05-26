## These functions are used to fit rGR for hazard analysis:
## The primary goal was to have final rGR values uncorrelated with LRS for use in
##  hazard modelling.

## Library dependencies
library(splines)
library(segmented)

######################################################################################
##
## rGR functions
##

## Function to select model to calculate rgr.  Goal is to find a well suited
##  model that also predicts bole volume growth that is not correlated to observed
##  bole volume at the end of the prediction period
## Model selection:
## - Choose best model by RMSE
## - If best model significantly correlated to bv at end of period, try next best
## - Possible models: linear, power, polynomials, segmented linear
removeCorr <- function(dat, models = c("lin","pow","poly","seg","bs"),
                       depen = "bvgrowth", indep = "priorbv", indep2 = "bv",
                       degree = 9, debug = FALSE, scottLRS = FALSE,
                       rgrLimits = c(0,10)) {
    if (length(models) < 1) stop("Must specify models: lin, pow, poly, seg, etc.")
    dep <- dat[,depen]
    ind <- dat[,indep]
    ind2 <- dat[,indep2]
    pval <- 0
    numTries = 20
    removed <- c()
    while (pval < 0.05 & numTries > 0) { # While the data is correlated
        if (!"lin" %in% removed) {
            fit.lin <- lm(dep ~ ind)
            rmse.lin <- sqrt(sum(residuals(fit.lin)^2))
            rgr.lin <- dep/predict(fit.lin)
            cor.lin <- cor.test(rgr.lin, ind2)$p.value
            if (scottLRS == TRUE)
                cor.lin <- cor.test(predictLRS(fit.lin, ind, ind2), rgr.lin)$p.value
        }
        if (!"bs" %in% removed) {
            fit.bs <- NULL
            try(fit.bs <- fitSplines(ind, ind2, dep, degree=degree)[[1]], silent = TRUE)
            try(fit.info <- fitSplines(ind, ind2, dep, degree=degree, info=TRUE),
                silent = TRUE)
            if (!is.null(fit.bs)) {
                rmse.bs <- fit.info[["rmses"]][fit.info[["best"]]]
                rgr.bs <- dep/predict(fit.bs)
                cor.bs <- cor.test(rgr.bs, ind2)$p.value
                if (scottLRS == TRUE)
                    cor.bs <- cor.test(predictLRS(fit.bs, ind, ind2), rgr.bs)$p.value
            }
        }
        if (!"seg" %in% removed) {
            fit.seg <- lm(dep ~ ind)
            try(fit.seg <- segmented(fit.seg, seg.Z = ~ind, psi = 2),silent = TRUE)
            if (!all(predict(fit.seg) == predict(fit.lin))) { # seg successful
                rmse.seg <- sqrt(sum(residuals(fit.seg)^2))
                rgr.seg <- dep/predict(fit.seg)
                cor.seg <- cor.test(rgr.seg, ind2)$p.value
                if (scottLRS == TRUE)
                    cor.seg <- cor.test(predictLRS(fit.seg, ind, ind2), rgr.seg)$p.value
            }
            else { # remove seg model
                if (debug == TRUE)
                    print(paste("Removing segmented model: fit failed"))
                removed <- c(removed, "seg")
            }
        }
        if(!"pow" %in% removed) {
            fit.pow <- NULL # control = nls.control(warnOnly = TRUE)
            try(fit.pow <- nls(dep ~ a*ind^b, start = list(a=0.5,b=0.5)), silent = TRUE)
            if (is.null(fit.pow)) fit.pow <- findPowerFit(dep, ind)
            if (!is.null(fit.pow)) {
                rmse.pow <- sqrt(sum(residuals(fit.pow)^2))
                rgr.pow <- dep/predict(fit.pow)
                cor.pow <- cor.test(rgr.pow, ind2)$p.value
            }
            if (scottLRS == TRUE)
                cor.pow <- cor.test(predictLRS(fit.pow, ind, ind2), rgr.pow)$p.value
        }
        if(!"poly" %in% removed) {
            polyn <- bestPoly(dat, polys = degree)
            degree <- polyn[["degree"]]
            fit.poly <- lm(as.formula(polyNoInt(degree,indep,depen)), data = dat)
            rmse.poly <- polyn[["rmse"]]
            rgr.poly <- dep/predict(fit.poly)
            cor.poly <- cor.test(rgr.poly, ind2)$p.value
            if (scottLRS == TRUE)
                cor.poly <- cor.test(predictLRS(fit.poly, ind, ind2), rgr.poly)$p.value
        }
        ## get best RMSE, check p-value, repeat without best if necessary
        stillhere <- models[!models %in% removed]
        if (length(stillhere) < 1) break;
        possibles<- paste0("rmse.",stillhere)
        rmses <- sapply(possibles, function(x) if (exists(x)) get(x))
        best <- rmses[which(rmses==min(rmses))]
        bestMOD <- gsub(".*[.]", "", names(best))
        pval <- get(paste0("cor",".",bestMOD))
        ## debug
        if (debug == TRUE) {
            print(c("Remaining models:",paste(stillhere)))
            print(paste("Poly degree:", degree))
            print(pval)
        }
        ## If p-value < 0.05 repeat without model (or in case of poly, with restricted
        ##  degree untill reaching degree 1
        if (pval < 0.05) {
            numTries = numTries - 1
            print(paste("removing:",bestMOD))
            if (bestMOD == "poly" & degree > 2)
                degree <- degree-1
            else
                removed <- c(removed, bestMOD)
        }
    }
    if (bestMOD == "bs") degree <- fit.info[["best"]]
    rgr = get(paste("rgr", bestMOD, sep = "."))
    ## Set rgr to conform to rgrLimits
    rgr[rgr > rgrLimits[2]] <- rgrLimits[2]
    rgr[rgr < rgrLimits[1]] <- rgrLimits[1]
    data.frame(
        rgr = rgr,
        rmse = rep(get(paste0("rmse.", bestMOD)), length(rgr)),
        corr = rep(get(paste0("cor.", bestMOD)), length(rgr)),
        model = rep(bestMOD, length(rgr)),
        degree = rep(degree, length(rgr)))
}

## helper function to predict LRS
## takes a model, prior and current (i.e. priorbv and bv)
## LRS = priorbv + predicted bvgrowth / max bv
predictLRS <- function(mod, prior, current) {
    LRS = (prior + predict(mod))/max(current)
    return(LRS)
}

## helper function to output results to data.frame for removeCorr
## Returns rgrs, rmse, corrs, model type, degree (0 if not polynomial)
returnDat <- function(bestMOD) {
    rgr = get(paste("rgr", bestMOD, sep = "."))
    data.frame(
        rgr = rgr,
        rmse = rep(get(paste0("rmse.", bestMOD)), length(rgr)),
        corr = rep(get(paste0("cor.", bestMOD)), length(rgr)),
        model = rep(bestMOD, length(rgr)),
        degree = rep(degree, length(rgr)))
}

## Function takes a range of polynomials and returns the polynomial best fitted
## NOTE: allows for the removal of intercept (set the intercept to be 0)
## To determine goodness of fit:
## - Iterate from lowest to highest
## - If coefficients lose significance stop and compare most recent to previous
##   by AIC
## - Choose the best of those and return that
bestPoly <- function(dat, ind = "priorbv", dep = "bvgrowth", polys = 9, corr = 0.05,
                     debug = FALSE) {
    bestRMSE <- Inf
    best <- NULL
    rmse <- NULL
    nonSig <- 0
    for (i in 2:polys) {
        if (debug==TRUE)
            print(paste("Fitting degree",i))
        fit <- NULL
        form <- polyNoInt(i,ind,dep)
        try(fit <- lm(as.formula(form), data = dat), silent=TRUE)
        if (!is.null(fit)) { ## Successful poly fit, check significance of coefs
            summ <- summary(fit)$coefficients[,4]
            nonSig <- length(summ[summ > corr])
            rmse <- sqrt(sum(residuals(fit)^2))
            if (debug==TRUE)
                print(paste("rmse of", rmse))
            if (rmse < bestRMSE) {
                best <- i
                bestRMSE <- rmse
            }
        }
        if (is.null(fit)) { ## Fit failed, return the last polynomial as the best
            i = i-1
            break
        }
        if (nonSig > 0) { ## non-significant coefs, compare to last and return best
            last <- lm(as.formula(polyNoInt(i-1,ind,dep)),data = dat)
            ifelse(AIC(last) < AIC(fit),
               { i = i-1; nonSig <- 0; break }, { break })
        }
    }
    results <- c(degree = i, rmse = rmse, numInsigCoefs = nonSig, aic = AIC(fit))
    return (results)
}

## Testing
## ind <- "priorbv"
## dep <- "bvgrowth"
## tst1 <- lm(as.formula(polyNoInt(2, ind, dep)), data = tst)
## tst2 <- lm(as.formula(polyNoInt(3, ind, dep)), data = tst)
## summary(tst1)
## summary(tst2)
## plot(tst$priorbv, tst$bvgrowth)
## points(tst$priorbv, predict(tst1), col = "blue")
## points(tst$priorbv, predict(tst2), col = "green")


## Helper function for bestPoly to make a formula for a polynomial without an intercept
## Takes integer argument defining the degree of the polynomial, ind. var, and dep. var
## Return form: dep ~ -1 + ind + I(ind^2) + ... + I(ind^n)
polyNoInt <- function(degs, ind, dep) {
    form <- paste(dep,"~","-1")
    for (i in 1:degs) {
        form <- paste0(form," + ","I(",ind,"^",i,")")
    }
    return(form)
}

## If power fitting fails in removeCorr this function will search more extensively
##  for a possible fit, and if it fails it will return NULL
findPowerFit <- function(depen, indep) {
    arange <- seq(.1,1,.1)
    brange <- seq(0.01,.5,length.out = 10)
    possibles <- expand.grid(a = arange,b = brange)
    fit <- NULL
    try(fit <- nls(depen ~ a * indep^b, start = list(a=0.5, b=0.01)),silent=TRUE)
    if (is.null(fit)) {
        for (row in 1:nrow(possibles)) {
            start = list(a = possibles[row,"a"], b = possibles[row, "b"])
            try(fit <- nls(depen ~ a * indep^b, start = start),silent = TRUE)
            if (!is.null(fit)) break;
        }
    }
    return(fit)
}

## Find best fit using splines bs to get basis
## Testing
## ind <- tst[,"priorbv"]
## ind2 <- tst[,"bv"]
## dep <- tst[,"bvgrowth"]
## degree <- 9
## dat <- tst
fitSplines <- function(ind, ind2, dep, degree=9, info = FALSE, debug = FALSE, corr=0.05) {
    fits <- list() ## store all fits to do analysis at end
    nonSig <- c() ## store number of insignificant coefs in each fit
    for (i in 1:degree) {
        if (debug==TRUE)
            print(paste("Fitting degree",i))
        fit <- NULL
        try(fit <- lm(dep ~ bs(ind, degree = i)), silent=TRUE)
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

## dat <- tst
## indep <- "priorbv"
## depen <- "bvgrowth"

## ## Testing
## x <- lm(bvgrowth ~ poly(priorbv, 11), data = tst)
## plot(tst$priorbv, residuals(x))
## rgr <- tst$bvgrowth/predict(x)
## cor.test(tst$bvgrowth, rgr)


## library(polynom)
## plot(poly.calc(1:13))
## plot(tst$priorbv, tst$bvgrowth)
## curve(tst$bvgrowth~ poly(x, 13), add=TRUE)
