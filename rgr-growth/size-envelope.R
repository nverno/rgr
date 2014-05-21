# Create a fit for individuals that grew the most
#  break data into size classe and fit growth vs size for individuals in the top
#  quantile
source("functions.R")
source("neighborhood-models.R")

# data (start with just DF data)
dflong <- read.csv("long-df.csv")
dflong <- subset(dflong, stat!="DEAD" & bagrowth>=0)

# size classes, breaking priorba into 16 quantiles (~1150 per class)
dflong$bacl <- (quantclass(dflong$priorba, 16))
table(dflong$bacl)
ggplot(dflong, aes(priorba, bagrowth, color=bacl)) + geom_point()

# top performers from each size class (grew the most)
dflong <- dflong[order(dflong$bacl),]
top <- rbind.fill(ddply(dflong, .(bacl), .fun = function(x) {
   bagrcl <- as.numeric(quantclass(x$bagrowth, 16, smallest = 0))
   bagrcl = data.frame(bagrcl = bagrcl)
}))
table(top)
dflong$bagrcl <- as.factor(top$bagrcl)
dflong$bacl <- as.numeric(dflong$bacl)

# graph it, add vertical lines to show bacls??
nonas <- subset(dflong, !is.na(priorba) & !is.na(bagrowth) & !is.na(bagrcl))
ggplot(nonas, aes(priorba, bagrowth, color=bagrcl)) + geom_point(alpha = 0.4)

# power fit to top performing class (bagrcl == 16)
topperf <- subset(dflong, bagrcl == 16)
ggplot(topperf, aes(priorba, bagrowth)) + geom_point() + geom_smooth()
top.fit <- nls(formula = bagrowth ~ a*priorba^b, start = list(a=1, b=1),
               data=topperf)
reg.fit <- nls(formula = bagrowth ~ a*priorba^b, start = list(a=1, b=1),
               data=dflong)
# graph it
plot(dflong$priorba, dflong$bagrowth, main = "Top performers and the rest")
points(topperf$priorba, topperf$bagrowth, col ="blue")
curve(coef(top.fit)[1] * x  ^ coef(top.fit)[2], add=TRUE, col = "red", lwd = 2)
curve(coef(reg.fit)[1] * x  ^ coef(reg.fit)[2], add=TRUE, col = "green", lwd = 2)

# residuals from top performers and regular fit
plot(top.fit)
plot(reg.fit)
scatter.smooth(topperf$priorba, residuals(top.fit)); abline(h=0, lty=2)

# calulate residuals from all plants to top performers fits
#  rgr: obs/predicted
dflong$rgr <- dflong[,"bagrowth"]/(coef(top.fit)[1]*dflong["priorba"]^coef(top.fit)[2])

# Neighbor density in radius of 4 meters
ind.var = "priorba"
fit.MLE.models(dflong, sr=4, spec="FD", ind.var = "priorba", dep.var = "bagrowth",
               realdist = TRUE)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)
dflong <- merge(dflong, targets, all=TRUE)

# graph rgr vs numnebs
plot(dflong$numnebs, dflong$rgr)

# MLE fit
# neighbor model with rgr
# Simple neighbor model as predictor of rgr, no size
simplest = function(ps, ind.var = "priorba")
{
    PG = ps[["PG"]]
    alpha = ps[["alpha"]]
    beta = ps[["beta"]]
    C = ps[["C"]]
    nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)
    competition.effect <- exp(-C*nci)
    PG * competition.effect
}

# create new starting parameters
## pars <- read.csv("parameters.csv")
## rowcopy = pars[82,]
## rowcopy$model <- "simplest"
## rowcopy$PG <- NA
## pars <- rbind(pars, rowcopy)
## write.csv(pars, "parameters.csv", row.names =FALSE)

# specify parameters for MLE fitting
dat = dflong
sr = c(3:6)
spec = "FD"
ind.var = "priorba" ## will be used in neighborhood calculations
dep.var = "rgr"
currentmodel = "simplest"
method = "Nelder-Mead"
maxit = 10000
realdist = TRUE
pars <- read.csv("parameters.csv")
ps <- get.params(sr,spec,ind.var,dep.var,currentmodel = currentmodel)

fit <- fit.MLE.models(dat=dat,sr=sr,spec=spec,ind.var=ind.var,dep.var=dep.var,
                      models=currentmodel,method="SANN",
                      realdist = realdist, maxit = 100000)

# get fit parameters and predict nci
#  nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)
sr <-4
newps <- get.params(sr = sr, spec = spec, ind.var = ind.var, dep.var = dep.var,
                    currentmodel = currentmodel)
alpha = newps[["alpha"]]; beta = newps[["beta"]]
targets$nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)

# graph rgr vs nci, bagrowth vs nci
par(mfrow = c(1,2))
plot(targets$nci, targets$rgr)
plot(targets$nci, targets$bagrowth)

# graph predicted over obs. on graph of rgr vs priorba
##targets$predicted <- do.call(currentmodel, list(newps))
plot(targets$nci, targets$rgr)
##points(targets$priorba, targets$predicted, col = "red")
##points(targets$nci, targets$predicted, col = "blue")

plot(targets$priorba, targets$rgr)
points(targets$priorba, targets$predicted, col = "red")
##points(targets$nci, targets$predicted, col = "blue")

# write file
write.csv(dflong, "dflong-rgr.csv", row.names=FALSE)
