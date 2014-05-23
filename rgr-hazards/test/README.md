## Check correlation between rgr and LRS

```r
dat <- read.csv("~/work/data/data/hazards/hazards-bc-firs.csv")

library(plyr)
corrs <- ddply(dat, .(install, plot, time), function(x) {
    cc <- cor.test(x$rgr, x$LRS, alternative = "two.sided", method = "pearson")
    data.frame(pval = cc$p.value, statistic = cc$statistic, corr = cc$estimate, 
        degree = unique(x$degree), method = unique(x$method), num_trees = nrow(x))
})
```


## Plots still correlated

```r
corrs <- corrs[order(corrs$pval), ]
corrs[corrs$pval < 0.05, ]
```

```
##     install plot time      pval statistic   corr degree method num_trees
## 132      19    5   97 8.039e-05     4.233 0.4795      2    lin        62
## 260      72   14   79 2.730e-04     4.182 0.6270      2     bs        29
## 256      72    2   85 4.766e-04     4.255 0.7081      1     bs        20
## 33        3   15   82 5.410e-04     3.889 0.5854      4     bs        31
## 8         1   16   79 6.141e-04     3.764 0.5368      2     bs        37
## 224      65    6   79 9.804e-04     3.382 0.2996      1     bs       118
## 178      51    3   85 1.674e-03     3.558 0.5958      2   poly        25
## 191      52    5   91 4.975e-03     2.899 0.3254      2     bs        73
## 46        4   17   85 5.187e-03     2.875 0.3078      3     bs        81
## 219      65    5   82 6.433e-03     2.780 0.2607      3     bs       108
## 77        8    1   91 2.341e-02     2.421 0.4430      1     bs        26
## 83        8   14   91 3.070e-02     2.247 0.3465      2     bs        39
## 170      51    1   79 4.128e-02     2.098 0.2927      7     bs        49
## 250      71   14   85 4.533e-02     2.023 0.1816      1     bs       122
```

