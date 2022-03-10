Treatment Effect Heterogeneity Methods
================

The purpose of this repository is to make recent tools for understanding
treatment effect heterogeneity more accessible. This first commit
contains R code for average treatment effect (ATE) and variance of
treatment effect estimation (VTE), using both one-step bias correction,
and targeted maximum likelihood estimation (TMLE). Example code is given
below.

Source required files

``` r
require(tidyverse)
source("R/aipw.R")
source("R/tmle.R")
source("R/vte.R")
source("R/example_helpers.R") #Used for the current examples
```

## Obtain intitial fits

Construct some data I use the data in the example helpers file

``` r
set.seed(1234)
N <- 500 #size of generated data
df <- generate_data_simple(N)
```

Next we will fit the working models We use the function `fitMods` in the
example helpers file. At the moment the best way of modifying this code
(e.g. to use other machine learning methods) is to write a new `fitMods`
function. The existing `fitMods` function uses the `mgcv` package to fit
GAM models.

``` r
df_fit <- fitMods(df)
print(df_fit)
#> # A tibble: 500 × 5
#>        Y     A pi_hat mu1_hat mu0_hat
#>    <dbl> <int>  <dbl>   <dbl>   <dbl>
#>  1 0.641     0  0.603    1.45  0.624 
#>  2 0.731     0  0.420    1.96  0.825 
#>  3 2.16      0  0.416    2.63  1.16  
#>  4 1.24      0  0.489    3.34  1.10  
#>  5 0.465     1  0.468    1.31 -0.0306
#>  6 1.27      0  0.415    2.53  1.09  
#>  7 0.970     1  0.649    2.16  0.849 
#>  8 1.05      1  0.548    2.36  1.15  
#>  9 2.14      1  0.469    2.29  0.559 
#> 10 0.542     1  0.445    1.87  0.830 
#> # … with 490 more rows
```

The `fitmods` function returns a data frame that has 5 necessary
columns:

| Column    | Value                                                             |
|-----------|-------------------------------------------------------------------|
| `Y`       | outcome of interest                                               |
| `A`       | treatment indicator                                               |
| `pi_hat`  | predicted propensity score *E*(*A*\|*X*)                          |
| `mu1_hat` | predicted outcome regression in the treated *E*(*Y*\|*A*=1,*X*)   |
| `mu0_hat` | predicted outcome regression in the untreated *E*(*Y*\|*A*=0,*X*) |

Often it is better to use cross fitted models (see e.g. CV-TMLE notes in
Levy’s paper) This is implemented in the `crossFit` function which
essentially calls `fitmods` multiple times

``` r
foldIDs <- getFoldIDs(N,Nfolds=5) #5 fold cross validation in this example
df_xfit <- crossFit(df,foldIDs=foldIDs)
print(df_xfit)
#> # A tibble: 500 × 6
#>        Y     A pi_hat mu1_hat mu0_hat FoldID
#>    <dbl> <int>  <dbl>   <dbl>   <dbl>  <int>
#>  1 0.641     0  0.587    1.43  0.634       5
#>  2 0.731     0  0.469    2.00  0.883       3
#>  3 2.16      0  0.407    2.63  1.19        4
#>  4 1.24      0  0.570    3.20  1.13        3
#>  5 0.465     1  0.465    1.25 -0.0644      1
#>  6 1.27      0  0.417    2.63  1.17        5
#>  7 0.970     1  0.669    2.51  0.793       1
#>  8 1.05      1  0.542    2.48  1.02        1
#>  9 2.14      1  0.442    2.27  0.541       1
#> 10 0.542     1  0.425    2.00  0.881       5
#> # … with 490 more rows
```

Note there is an additional `FoldID` column. This is currently not used
for anything but may be useful in future.

## Pass initial fits to estimators

These fitted data frames can be passed used to infer ATE and VTEs We
also infer the square-root of the VTE which is on the same scale as the
ATE Two methods are provided: “AIPW” uses one step bias correction
estimators “TMLE” uses TMLE

``` r
aipw_nocv <- VTE(df_fit,method="AIPW")
tmle_nocv <- VTE(df_fit,method="TMLE")
aipw_cv   <- VTE(df_xfit,method="AIPW")
tmle_cv   <- VTE(df_xfit,method="TMLE")

print(tmle_cv) #example of displaying the output
#> 
#> Call:
#> VTE(data = df_xfit, method = "TMLE")
#> 
#> Estimator: TMLE
#> Estimand Values:
#>         Estimate Std.Error Wald.value  Wald.pval    
#> ATE     1.371375  0.100013   188.0186 < 2.22e-16 ***
#> rootVTE 0.917479  0.209974    19.0924 1.2454e-05 ***
#> VTE     0.841768  0.192647    19.0924 1.2454e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The method also returns the CATE model. e.g. For the AIPW method the
CATE is equivalent to the T-learner (by default)

``` r
tibble(
  TMLE_CATE = tmle_cv$CATE,
  AIPW_CATE = aipw_cv$CATE, 
  T_learner = df_xfit$mu1_hat-df_xfit$mu0_hat) 
#> # A tibble: 500 × 3
#>    TMLE_CATE AIPW_CATE T_learner
#>        <dbl>     <dbl>     <dbl>
#>  1     0.833     0.796     0.796
#>  2     1.14      1.12      1.12 
#>  3     1.45      1.44      1.44 
#>  4     2.06      2.07      2.07 
#>  5     1.33      1.32      1.32 
#>  6     1.47      1.46      1.46 
#>  7     1.72      1.72      1.72 
#>  8     1.47      1.45      1.45 
#>  9     1.73      1.73      1.73 
#> 10     1.14      1.11      1.11 
#> # … with 490 more rows
```

For the TMLE the CATE is a targeted T-learner so that the following are
equal

``` r
mean(tmle_cv$CATE)
#> [1] 1.371375
tmle_cv$coef["ATE"]
#>      ATE 
#> 1.371375

mean(tmle_cv$CATE^2) - mean(tmle_cv$CATE)^2
#> [1] 0.8417677
tmle_cv$coef["VTE"]
#>       VTE 
#> 0.8417677
```
