## Check correlation between rgr and LRS
dat <- read.csv("~/work/data/data/hazards/hazards-bc-firs.csv")

library(plyr)
corrs <- ddply(dat, .(install, plot, time), function(x) {
    cc <- cor.test(x$rgr, x$LRS, alternative = "two.sided", method="pearson")
    data.frame(pval = cc$p.value,
               statistic = cc$statistic,
               corr = cc$estimate,
               degree = unique(x$degree),
               method = unique(x$method),
               num_trees = nrow(x))
})

## Plots still correlated
corrs <- corrs[order(corrs$pval),]
corrs[corrs$pval < 0.05, ]


cc <- cor.test(dat$rgr, dat$LRS)
