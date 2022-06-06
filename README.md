# PanelMatch: Matching Methods for Causal Inference with Time-Series Cross-Section Data 
[![R build
status](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions) [![Build Status](https://travis-ci.org/insongkim/PanelMatch.svg?branch=master)](https://travis-ci.org/insongkim/PanelMatch) ![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/PanelMatch)  [![CRAN status](https://www.r-pkg.org/badges/version/PanelMatch)](https://CRAN.R-project.org/package=PanelMatch) 

Authors: In Song Kim (insong@mit.edu), Adam Rauh (amrauh@umich.edu), Erik Wang (haixiaow@Princeton.edu), Kosuke Imai (imai@harvard.edu)

PanelMatch is an R package implementing a set of methodological tools proposed by Imai, Kim, and Wang (2021) that enables researchers to apply matching methods for causal inference on time-series cross-sectional data with binary treatments. The package includes implementations of matching methods based on propensity scores and Mahalanobis distance, as well as weighting methods. PanelMatch enables users to easily calculate a variety of possible quantities of interest, along with standard errors. The software is flexible, allowing users to tune the matching, refinement, and estimation procedures with a large number of parameters. The package also offers a variety of visualization and diagnostic tools for researchers to better understand their data and assess their results.

Installation Instructions
-------------------------

`PanelMatch` is available on CRAN and can be installed using:

``` r
install.packages("PanelMatch")
```

You can install the most recent development version of `PanelMatch` using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `PanelMatch`:

``` r
library(devtools)
install_github("insongkim/PanelMatch", dependencies=TRUE, ref = "se_comparison")
```
If you encounter problems during installation, please consult [the wiki page](https://github.com/insongkim/PanelMatch/wiki/Installation-Troubleshooting) that has some ideas for handling common issues. 


Usage Examples
-------------------------

### Treatment Variation Plot

Users can visualize the variation of treatment across space and
time. This will help users build an intuition about how comparison of
treated and control observations can be made.

```r
library(PanelMatch)
DisplayTreatment(unit.id = "wbcode2",
                 time.id = "year", legend.position = "none",
                 xlab = "year", ylab = "Country Code",
                 treatment = "dem", data = dem)
```
![](https://github.com/insongkim/repo-data/blob/master/panelmatch/DT_plot.png)

### PanelMatch

`PanelMatch` identifies a matched set for each treated
 observation. Specifically, for a given treated unit, the matched set
 consists of control observations that have an identical treatment
 history up to a chosen number (`lag`) of years. This number corresponds with the `lag` parameter, which must
 be chosen by the user. Users must also consider various parameters regarding the refinement of created matched sets. Please consult the function documentation for a full set of descriptions, but some important arguments are described below:
 1) `refinement.method` -- Users may choose between standard propensity score weighting or matching (`ps.weight`, `ps.match`), covariate balanced propensity score weighting or matching (`CBPS.weight`, `CBPS.match`),  and mahalanobis distance matching (`mahalanobis`). Users may also opt to apply the idea of marginal structural models with the `CBPS.msm.weight` and `ps.msm.weight` methods. Alternatively users can do no refinement by setting this parameter to `none`.
 2) `size.match` -- This sets the maximum number of control units that can be included in a matched set.
 3) `covs.formula` -- This parameter defines which variables are considered in measuring the similarities/distances between units. These will then affect which control units are included/excluded during refinement. This can be set to include lagged versions of any variable as well. See the `PanelMatch` documentation for more information about this parameter.
 4) `match.missing` -- Should matches between treatment and control units with identical patterns of missingness in the treatment variable be considered? If set to FALSE, missing data is not permitted in the lag window of the treatment variable in either treated or control units. 
``` r
PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "mahalanobis", 
                         data = dem, match.missing = TRUE, 
                         covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                         size.match = 5, qoi = "att" ,outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE)

```							
The `PanelMatch` function will return an object of class "PanelMatch". This is a list that contains a few specific elements: First, a matched.set object(s) that has the same name as the provided qoi -- if the qoi is "att", "atc". If qoi = "ate" then two matched.set objects will be attached, named "att" and "atc." Users can extract information about individual matched sets as well as statistics about all created matched sets from this object. Consult the [Wiki page on Matched Set Objects](https://github.com/insongkim/PanelMatch/wiki/Matched-Set-Objects) for a more detailed walk through and description of these objects. Put simply, `matched.set` objects are merely lists with some assumed structure and special attributes.

The `PanelMatch` object also has some additional attributes: "qoi", "lead", "forbid.treatment.reversal" (a logical value that is the same as what was specified in the function call), and "outcome.var" (character value that is the same as what was specified in the function call)

You can check covariate balance using the `get_covariate_balance` function:

```{r}
get_covariate_balance(PM.results$att, dem, covariates = c("tradewb"), plot = FALSE, ylim = c(-2,2))
        tradewb
t_4  0.05459705
t_3 -0.03101839
t_2 -0.01828529
t_1  0.07784846
```
See the documentation for more information about this function.

### PanelEstimate

Once proper matched sets are attained by `PanelMatch`, users can
estimate the causal quantity of interest such as the average
treatment effect using `PanelEstimate`. Either bootstrap or weighted
fixed effects methods can be used for standard error
calculation. Users can estimate the contemporaneous effect as well as
long-term effects. In this example, we illustrate the use of
`PanelEstimate` to estimate the average treatment effect on treated units (att) at time `t` on the outcomes from time `t+0` to `t+4`.

```r
PE.results <- PanelEstimate(sets = PM.results, data = dem)
```

The `PanelEstimate` function returns a `PanelEstimate` object, which is a named list. This object will contain the point estimates, standard errors and other information about the calculations. See the wiki page about PanelEstimate objects for more information. 

Users can easily obtain and visualize important information about esimtates and standard errors using the `summary` and `plot` methods for PanelEstimate objects

```r
summary(PE.results)
Weighted Difference-in-Differences with Mahalanobis Distance
Matches created with 4 lags

Standard errors computed with 1000 Weighted bootstrap samples

Estimate of Average Treatment Effect on the Treated (ATT) by Period:
$summary
      estimate std.error      2.5%    97.5%
t+0 -0.5349572 0.9231828 -2.434598 1.314614
t+1 -0.2396204 1.4565119 -2.970407 2.506002
t+2  0.5532550 1.8562746 -2.949847 4.059749
t+3  1.8425824 2.1679704 -2.421909 5.953401
t+4  1.9920680 2.3756352 -2.648517 6.569985

$lag
[1] 4

$iterations
[1] 1000

$qoi
[1] "att"
```

```r
plot(PE.results)
```
![](https://github.com/insongkim/repo-data/blob/master/panelmatch/pe_plot_6_22.png)
