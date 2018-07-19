# PanelMatch: Matching Methods for Causal Inference with Time-Series Cross-Section Data [![Build Status](https://travis-ci.org/insongkim/PanelMatch.svg?branch=master)](https://travis-ci.org/insongkim/PanelMatch)

This R package provides a set of methodological tools that enable
researchers to apply matching methods to time-series cross-section
data.  Imai, Kim, and Wang (2018) proposes a nonparametric
generalization of difference-in-differences estimator, which does not
rely on the linearity assumption as often done in
practice. Researchers first select a method of matching each treated
observation from a given unit in a particular time period with control
observations from other units in the same time period that have a
similar treatment and covariate history.  These methods include
standard matching and weighting methods based on propensity score and Mahalanobis
distance.
Once matching is done, both short-term and long-term average treatment
effects for the treated can be estimated with standard errors.  The
package also offers a visualization technique that allows researchers
to assess the quality of matches by examining the resulting covariate
balance.

Installation Instructions
-------------------------

<!-- `panelMatch` is available on CRAN and can be installed using: -->

<!-- ``` r -->
<!-- install.packages("panelMatch") -->
<!-- ``` -->

You can install the most recent development version of `PanelMatch` using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `PanelMatch`:

``` r
library(devtools)
install_github("insongkim/PanelMatch", dependencies=TRUE)
```
If you encounter problems during installation, please consult [the wiki page](https://github.com/insongkim/PanelMatch/wiki/Installation-Troubleshooting) that has some ideas for handling common issues. 


Usage Examples
-------------------------

### Treatment Variation Plot

Users can visualize the variation of treatment across space and
time. This will help users build an intuition about how comparison of
treated and control observation can be made.

```r
library(PanelMatch)
DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem)
```
![](http://web.mit.edu/insong/www/pdf/varPlot.png)

### PanelMatch

`PanelMatch` identifies a matched set for each treated
 observation. Specifically, for a given treated unit, the matched set
 consists of control observations that have the identical treatment
 history up to a certain number of `lag` years. Researchers must
 specify `lag`. Researchers may also specify `covars.lagged`, a vector
 of covariates to be automatically lagged and included in the model formula.
 A further refinement of the matched set will be
 possible by setting the size of the matched set `M`, and adjusting
 for other confounders such as past outcomes and covariates via
 `formula`. Various matching and weighting methods such as `Mahalanobis distance`
 matching and `CBPS` and `Propensity score` matching and weighting can 
 be used.

``` r
matches.cbps <- PanelMatch(lag = 4, max.lead = 4, time.id = "year",
                           unit.id = "wbcode2", treatment = "dem",
                           covars.lagged = c("logpop", "tradewb", "nfagdp", 
                               "Populationages1564oftota", "unrest"),
                           formula =  y ~ dem, method = "CBPS",
                           weighting = FALSE,  qoi = "ate",  M = 5, data = dem)
```							

Users should closely examine the matched sets, and check the balance
between treated and control units in terms of their observable
pre-treatment characteristics.

### PanelEstimate

Once proper matched sets are attained by `PanelMatch`, users can
estimate the causal quantity of interest such as the average
treatment effect using `PanelEstimate`. Either bootstrap or weighted
fixed effects methods can be used for standard error
calculation. Users can estimate the contemporaneous effect as well as
long-term effects. In this example, we illustrate the use of
`PanelEstimate` to estimate the ATT of treatment at time `t` on the
outcomes on `t` through `t+4`.

```r

mod.bootSE <- PanelEstimate(lead = 0:4, inference = "bootstrap",
                            matched_sets = matches.cbps,
                            qoi = "att", CI = .95,
                            ITER = 500)

summary(mod.bootSE)

Weighted Difference-in-Differences with Covariate Balancing Propensity Score
Matches created with 4 lags

Standard errors computed with 500 Weighted bootstrap samples

Estimate of Average Treatment Effect on the Treated (ATT) by Period:

|                                                       |        t+0|        t+1|        t+2|        t+3|        t+4|
|:------------------------------------------------------|----------:|----------:|----------:|----------:|----------:|
|Point Estimate(s)                                      |  0.3697332|  1.9657985|  3.0537125|  4.1564464|  4.1361189|
|Standard Error(s)                                      |  0.8580289|  1.3487225|  1.7604525|  2.1510657|  2.5481439|
|Lower Limit of 95 % Regular Confidence Interval        | -1.3428403| -0.6863441| -0.4301161|  0.0415114| -0.5597593|
|Upper Limit of 95 % Regular Confidence Interval        |  1.9864734|  4.5182535|  6.8135019|  8.6022466|  9.6319008|
|Bias-corrected Estimate(s)                             |  0.3890818|  1.9662642|  3.0584442|  4.1538121|  4.1586525|
|Lower Limit of 95 % Bias-corrected Confidence Interval | -1.2470070| -0.5866564| -0.7060769| -0.2893537| -1.3596630|
|Upper Limit of 95 % Bias-corrected Confidence Interval |  2.0823067|  4.6179412|  6.5375412|  8.2713815|  8.8319971|
```
