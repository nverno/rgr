## compare envelope models: by plot vs by whole
dflong <- read.csv("dflong-rgr.csv")

## 6 meter radius
# specify parameters for MLE fitting, using 'simplest' model
dat = dflong
sr = c(6)
spec = "FD"
ind.var = "priorba" ## will be used in neighborhood calculations
dep.var = "rgrplot"
currentmodel = "simplest"
method = "Nelder-Mead"
maxit = 10000
realdist = TRUE
pars <- read.csv("parameters.csv")

# parameters for rgr by plot
psbyplot <- get.params(sr,spec,ind.var,dep.var,currentmodel = currentmodel)

# parameters for rgr as whole
dep.var = "rgr"
pswhole <- get.params(sr,spec,ind.var,dep.var,currentmodel = currentmodel)

# make neighbor matrices to get predicted values
## Neighbor density in radius of 4 meters
fit.MLE.models(dflong, sr=sr, spec=spec, ind.var =ind.var, dep.var =dep.var,
               realdist = realdist)
targets$numnebs <- rowSums(bas > 0, na.rm = TRUE)
dflong <- merge(dflong, targets, all=TRUE)

# get predicted values using parameters from both fits...
#  nci <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)
alpha = psbyplot[["alpha"]]; beta = psbyplot[["beta"]]
targets$ncibyplot <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)

alpha = pswhole[["alpha"]]; beta = pswhole[["beta"]]
targets$nciwhole <- rowSums(((bas ^ alpha)/(distances ^ beta)), na.rm=TRUE)

# graph rgr vs nci
pdf("rgrVSnci.pdf")
par(mfrow = c(1,2))
plot(targets$ncibyplot, targets$rgrplot, main = "RGR vs NCI, by Plot")
plot(targets$nciwhole, targets$rgrwhole, main = "RGR vs NCI")
dev.off()

## graph bagrowth vs nci
pdf("rgrVSnci.pdf")
par(mfrow = c(1,2))
plot(targets$ncibyplot, targets$bagrowth, main = "RGR vs BA growth, by Plot")
plot(targets$nciwhole, targets$bagrowth, main = "RGR vs NCI")
dev.off()

## graph bagrowth vs nci

plot(targets$nci, targets$bagrowth)
