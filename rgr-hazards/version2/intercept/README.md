Compare some methods to remove intercept from models
========================================================


These are possible methods to remove the intercept from models fitting bole volume growth to prior bole volume:
* Recenter data around zero
* Fit linear/polynomial/B-spline models and explicitly specify no intercept

Data was centered by subtracting the means.





Test data from install `inst`, plot `plt`, and year `tm`.



```r
ind <- tst[, "bv"]
dep <- tst[, "bvgrowth"]

## Centered data
indCent <- ind - mean(ind)
depCent <- dep - mean(dep)

## Poly fits
fit1 <- lm(depCent ~ poly(indCent, degree = 3))
fit2 <- lm(depCent ~ poly(indCent, degree = 3) - 1)
fit3 <- lm(dep ~ poly(ind, degree = 3))
fit4 <- lm(dep ~ poly(ind, degree = 3) - 1)

## bs fits
fit5 <- lm(dep ~ bs(ind, degree = 3))
fit6 <- lm(dep ~ bs(ind, degree = 3) - 1)
fit7 <- lm(depCent ~ bs(indCent, degree = 3))
fit8 <- lm(depCent ~ bs(indCent, degree = 3) - 1)

## linear
fit9 <- lm(dep ~ ind)
fit10 <- lm(depCent ~ indCent)
fit11 <- lm(dep ~ ind - 1)
fit12 <- lm(depCent ~ indCent - 1)
```



## Poly
### centered

```r
summary(fit1)
```

```
## 
## Call:
## lm(formula = depCent ~ poly(indCent, degree = 3))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00207 -0.00001  0.00305  0.02995 
## 
## Coefficients:
##                             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                 3.90e-19   2.51e-03    0.00    1.000    
## poly(indCent, degree = 3)1  3.75e-01   1.55e-02   24.20   <2e-16 ***
## poly(indCent, degree = 3)2  3.41e-02   1.55e-02    2.20    0.035 *  
## poly(indCent, degree = 3)3 -9.91e-03   1.55e-02   -0.64    0.527    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0155 on 34 degrees of freedom
## Multiple R-squared:  0.946,	Adjusted R-squared:  0.941 
## F-statistic:  197 on 3 and 34 DF,  p-value: <2e-16
```

```r
summary(fit2)
```

```
## 
## Call:
## lm(formula = depCent ~ poly(indCent, degree = 3) - 1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00207 -0.00001  0.00305  0.02995 
## 
## Coefficients:
##                            Estimate Std. Error t value Pr(>|t|)    
## poly(indCent, degree = 3)1  0.37506    0.01528   24.55   <2e-16 ***
## poly(indCent, degree = 3)2  0.03407    0.01528    2.23    0.032 *  
## poly(indCent, degree = 3)3 -0.00991    0.01528   -0.65    0.521    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0153 on 35 degrees of freedom
## Multiple R-squared:  0.946,	Adjusted R-squared:  0.941 
## F-statistic:  203 on 3 and 35 DF,  p-value: <2e-16
```


### untransformed

```r
summary(fit3)
```

```
## 
## Call:
## lm(formula = dep ~ poly(ind, degree = 3))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00207 -0.00001  0.00305  0.02995 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)             0.04582    0.00251   18.22   <2e-16 ***
## poly(ind, degree = 3)1  0.37506    0.01550   24.20   <2e-16 ***
## poly(ind, degree = 3)2  0.03407    0.01550    2.20    0.035 *  
## poly(ind, degree = 3)3 -0.00991    0.01550   -0.64    0.527    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0155 on 34 degrees of freedom
## Multiple R-squared:  0.946,	Adjusted R-squared:  0.941 
## F-statistic:  197 on 3 and 34 DF,  p-value: <2e-16
```

```r
summary(fit4)
```

```
## 
## Call:
## lm(formula = dep ~ poly(ind, degree = 3) - 1)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.0183  0.0437  0.0458  0.0489  0.0758 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## poly(ind, degree = 3)1  0.37506    0.05013    7.48  9.2e-09 ***
## poly(ind, degree = 3)2  0.03407    0.05013    0.68     0.50    
## poly(ind, degree = 3)3 -0.00991    0.05013   -0.20     0.84    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0501 on 35 degrees of freedom
## Multiple R-squared:  0.617,	Adjusted R-squared:  0.585 
## F-statistic: 18.8 on 3 and 35 DF,  p-value: 1.92e-07
```


## B-spline
### centered

```r
summary(fit7)
```

```
## 
## Call:
## lm(formula = depCent ~ bs(indCent, degree = 3))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00207 -0.00001  0.00305  0.02995 
## 
## Coefficients:
##                          Estimate Std. Error t value Pr(>|t|)    
## (Intercept)              -0.04581    0.00452  -10.13  8.3e-12 ***
## bs(indCent, degree = 3)1  0.04705    0.02184    2.15    0.038 *  
## bs(indCent, degree = 3)2  0.15903    0.02315    6.87  6.5e-08 ***
## bs(indCent, degree = 3)3  0.25204    0.01594   15.81  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0155 on 34 degrees of freedom
## Multiple R-squared:  0.946,	Adjusted R-squared:  0.941 
## F-statistic:  197 on 3 and 34 DF,  p-value: <2e-16
```

```r
summary(fit8)
```

```
## 
## Call:
## lm(formula = depCent ~ bs(indCent, degree = 3) - 1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.05597 -0.03829 -0.01694  0.00839  0.03397 
## 
## Coefficients:
##                          Estimate Std. Error t value Pr(>|t|)    
## bs(indCent, degree = 3)1  -0.1125     0.0299   -3.76  0.00061 ***
## bs(indCent, degree = 3)2   0.2062     0.0448    4.60  5.3e-05 ***
## bs(indCent, degree = 3)3   0.1928     0.0293    6.58  1.3e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0306 on 35 degrees of freedom
## Multiple R-squared:  0.781,	Adjusted R-squared:  0.763 
## F-statistic: 41.7 on 3 and 35 DF,  p-value: 1.2e-11
```


### untransformed

```r
summary(fit5)
```

```
## 
## Call:
## lm(formula = dep ~ bs(ind, degree = 3))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00207 -0.00001  0.00305  0.02995 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          1.06e-05   4.52e-03    0.00    0.998    
## bs(ind, degree = 3)1 4.70e-02   2.18e-02    2.15    0.038 *  
## bs(ind, degree = 3)2 1.59e-01   2.32e-02    6.87  6.5e-08 ***
## bs(ind, degree = 3)3 2.52e-01   1.59e-02   15.81  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0155 on 34 degrees of freedom
## Multiple R-squared:  0.946,	Adjusted R-squared:  0.941 
## F-statistic:  197 on 3 and 34 DF,  p-value: <2e-16
```

```r
summary(fit6)
```

```
## 
## Call:
## lm(formula = dep ~ bs(ind, degree = 3) - 1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.06416 -0.00206 -0.00001  0.00306  0.02995 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## bs(ind, degree = 3)1   0.0471     0.0149    3.16   0.0033 ** 
## bs(ind, degree = 3)2   0.1590     0.0224    7.11  2.7e-08 ***
## bs(ind, degree = 3)3   0.2521     0.0146   17.24  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0153 on 35 degrees of freedom
## Multiple R-squared:  0.964,	Adjusted R-squared:  0.961 
## F-statistic:  317 on 3 and 35 DF,  p-value: <2e-16
```


## Linear
### centered

```r
summary(fit10)
```

```
## 
## Call:
## lm(formula = depCent ~ indCent)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07074 -0.00347  0.00287  0.00520  0.03220 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -4.41e-18   2.63e-03     0.0        1    
## indCent      4.56e-02   1.97e-03    23.2   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0162 on 36 degrees of freedom
## Multiple R-squared:  0.937,	Adjusted R-squared:  0.935 
## F-statistic:  537 on 1 and 36 DF,  p-value: <2e-16
```

```r
summary(fit12)
```

```
## 
## Call:
## lm(formula = depCent ~ indCent - 1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07074 -0.00347  0.00287  0.00520  0.03220 
## 
## Coefficients:
##         Estimate Std. Error t value Pr(>|t|)    
## indCent  0.04556    0.00194    23.5   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.016 on 37 degrees of freedom
## Multiple R-squared:  0.937,	Adjusted R-squared:  0.935 
## F-statistic:  552 on 1 and 37 DF,  p-value: <2e-16
```


### untransformed

```r
summary(fit9)
```

```
## 
## Call:
## lm(formula = dep ~ ind)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07074 -0.00347  0.00287  0.00520  0.03220 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -0.00660    0.00347    -1.9    0.065 .  
## ind          0.04556    0.00197    23.2   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0162 on 36 degrees of freedom
## Multiple R-squared:  0.937,	Adjusted R-squared:  0.935 
## F-statistic:  537 on 1 and 36 DF,  p-value: <2e-16
```

```r
summary(fit11)
```

```
## 
## Call:
## lm(formula = dep ~ ind - 1)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.07141 -0.00673 -0.00302 -0.00124  0.03438 
## 
## Coefficients:
##     Estimate Std. Error t value Pr(>|t|)    
## ind  0.04311    0.00154      28   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.0167 on 37 degrees of freedom
## Multiple R-squared:  0.955,	Adjusted R-squared:  0.954 
## F-statistic:  783 on 1 and 37 DF,  p-value: <2e-16
```


