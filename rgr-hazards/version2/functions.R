## Functions to fit rgr for hazard analysis
## Goal: to have as little correlation as possible between rgr and LRS for
##  hazard models
source("~/work/rgr/rgr-hazards/version2/bsplines/functions.R")
require(splines)
require(segmented)

################################################################################
##
##                               rGR functions
##
################################################################################

## Function to compute rGR
## Currently it is just the residuals of fitting
## Formerly: observed / predicted
rgrfunc <- function(dep, mod) {
    dep - predict(mod)
}

## Function to select model to calculate rgr.  Goal is to find a well suited
##  model that also predicts bole volume growth that is not correlated to LRS
## Model selection:
## - Choose best model by RMSE
## - If best model significantly correlated to bv at end of period, try next best
## - Possible models: linear, power, polynomials, segmented linear
removeCorr <- function(dat, models = c("lin","pow","poly","seg","bs"),
                       depen = "bvgrowth", indep = "priorbv", indep2 = "bv",
                       degree = 9, debug = FALSE, scottLRS = FALSE,
                       rgrLimits = c(-100,100), intrcpt=FALSE) {
    if (length(models) < 1) stop("Must specify models: lin, pow, poly, seg, etc.")
    dep <- dat[,depen]
    ind <- dat[,indep]
    ind2 <- dat[,indep2]
    pval <- 0
    numTries = 20
    removed <- c()
    while (pval < 0.05 & numTries > 0) { # While the data is correlated
        if (!"lin" %in% removed) {
            ifelse(intrcpt,
               { fit.lin <- lm(dep ~ ind) },
               { fit.lin <- lm(dep ~ ind - 1) })
            rmse.lin <- sqrt(sum(residuals(fit.lin)^2))
            rgr.lin <- rgrfunc(dep, fit.lin)
            cor.lin <- cor.test(rgr.lin, ind2)$p.value
            if (scottLRS == TRUE)
                cor.lin <- cor.test(predictLRS(fit.lin, ind, ind2), rgr.lin)$p.value
        }
        if (!"bs" %in% removed) {
            fit.bs <- NULL
            try(fit.bs <- fitSplines(ind, ind2, dep, intrcpt=intrcpt, degree=degree)[[1]], silent = TRUE)
            try(fit.info <- fitSplines(ind, ind2, dep, degree=degree, info=TRUE, intrcpt=intrcpt),
                silent = TRUE)
            if (!is.null(fit.bs)) {
                rmse.bs <- fit.info[["rmses"]][fit.info[["best"]]]
                rgr.bs <- rgrfunc(dep, fit.bs)
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
                rgr.seg <- rgrfunc(dep, fit.seg)
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
                rgr.pow <- rgrfunc(dep, fit.pow)
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
            rgr.poly <- rgrfunc(dep, fit.poly)
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
