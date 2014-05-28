################################################################################
##
##                    Store rGR Coefficients for Plotting
##
## - To get the parameters, just run fitSplines function with appropriate
##   variables.  The degree of the fit used is stored in "hazard-firs-bc.csv"
##
################################################################################
## Currently all the fits are B-spline polynomials
source("~/work/rgr/rgr-hazards/version2/bsplines/functions.R")

## data
dat <- read.csv("~/work/data/data/hazards/hazards-bc-firs.csv")

tst <- dat[dat$install == 1 & dat$plot == 10 & dat$time == 85,]
fitSplines(ind=tst$priorbv, ind2=tst$bv, dep=tst$bvgrowth)

