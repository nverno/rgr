## Data: long-bc-derived.csv
## rGR: Envelope modeled to top performers in subsets by SI and SDP
## Neighbor model: "simplest"
## SR: 6
## Target species: Doug Firs
## Neighbor variable: Prior Bole Volume
## Fitting Algorithm: Nelder-Mead
source("functions.R")
source("neighbor-model/neighborhood-models.R")

## Data
dat <- read.csv("long-bc-derived.csv")

## define arguments
rgrcol = "rgrsisdp"
sr = 6
spec = "FD"
ind.var = "priorbv"
dep.var = "rgrsisdp"
currentmodel = "simplest"

## Make neighbor matrices
fit.MLE.models(dat, sr=sr, spec=spec, ind.var = ind.var, dep.var = dep.var,
               realdist = TRUE)

## Retrieve parameters of most recent fit
pars <- read.csv("parameters.csv")
ps <- get.params(sr,spec,ind.var,dep.var,currentmodel = currentmodel)
ps

## Predict NCI
alpha = ps[["alpha"]]; beta = ps[["beta"]]
targets$nci2 <- rowSums(bas > 0, na.rm = TRUE)

targets$nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)

## Graph rGR vs NCI, bvgrowth vs NCI
par(mfrow = c(1,2))
plot(targets$nci, targets[,dep.var])

plot(targets$nci, targets[,dep.var])
points(targets$nci, targets$predicted, col = "red")

## Graph predicted over obs. on graph of rGR vs priorBV
targets$predicted <- do.call(currentmodel, list(ps))
plot(targets$nci, targets[,rgrcol], main = "RGR vs Number of Neighbors",
     xlab = "Number of neighbors", ylab = "RGR")
##points(targets$priorba, targets$predicted, col = "red")
##points(targets$nci, targets$predicted, col = "blue")

plot(targets$nci2, targets[,dep.var])

##pdf('rgrbyplot.pdf')
plot(targets$priorbv, targets[,rgrcol], main = "RGR by plot vs. Prior BV",
       xlab = "Prior BV", ylab = "RGR")
points(targets$priorbv, targets$predicted, col = "red")
##dev.off()
##points(targets$nci, targets$predicted, col = "blue")
