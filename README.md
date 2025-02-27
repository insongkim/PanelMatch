# PanelMatch: Matching Methods for Causal Inference with Time-Series Cross-Section Data 
[![R build
status](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions) [![Build Status](https://travis-ci.org/insongkim/PanelMatch.svg?branch=master)](https://travis-ci.org/insongkim/PanelMatch) ![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/PanelMatch)  [![CRAN status](https://www.r-pkg.org/badges/version/PanelMatch)](https://CRAN.R-project.org/package=PanelMatch) 

Authors: In Song Kim (insong@mit.edu), Adam Rauh (amrauh@umich.edu), Erik Wang (haixiaow@Princeton.edu), Kosuke Imai (imai@harvard.edu)

PanelMatch is an R package implementing a set of methodological tools proposed by Imai, Kim, and Wang (2023) that enables researchers to apply matching methods for causal inference on time-series cross-sectional data with binary treatments. The package includes implementations of matching methods based on propensity scores and Mahalanobis distance, as well as weighting methods. PanelMatch enables users to easily calculate a variety of possible quantities of interest, along with standard errors. The software is flexible, allowing users to tune the matching, refinement, and estimation procedures with a large number of parameters. The package also offers a variety of visualization and diagnostic tools for researchers to better understand their data and assess their results.

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
install_github("insongkim/PanelMatch", dependencies=TRUE, ref = "version3")
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
dem.panel <- PanelData(panel.data = dem,
               unit.id = "wbcode2",
               time.id = "year",
               treatment = "dem",
               outcome = "y")
DisplayTreatment(panel.data = dem.panel, legend.position = "none",
                 xlab = "year", ylab = "Country Code")
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
PM.results <- PanelMatch(panel.data = dem.panel, lag = 4, 
                          refinement.method = "ps.match", 
                          match.missing = TRUE, 
                          covs.formula = ~ tradewb + I(lag(tradewb, 1:4) + I(lag(y, 1:4))),
                          size.match = 5, qoi = "att",
                          lead = 0:4, 
                          forbid.treatment.reversal = FALSE)

```							
The `PanelMatch` function will return an object of class "PanelMatch". This is a list that contains a few specific elements: First, a matched.set object(s) that has the same name as the provided qoi -- if the qoi is "att", "atc". If qoi = "ate" then two matched.set objects will be attached, named "att" and "atc." Users can extract information about individual matched sets as well as statistics about all created matched sets from this object. 


You can check covariate balance using the `get_covariate_balance` function:

```{r}
PM.results <- PanelMatch(panel.data = dem.panel, lag = 4, 
                          refinement.method = "ps.match", 
                          match.missing = TRUE, 
                          covs.formula = ~ I(lag(tradewb, 1:4)),
                          size.match = 5, qoi = "att",
                          lead = 0:4, 
                          forbid.treatment.reversal = FALSE)

get_covariate_balance(PM.results, panel.data = dem.panel, covariates = c("tradewb"))

$att
       tradewb
t_4 0.27372313
t_3 0.17210534
t_2 0.12176551
t_1 0.11029284
t_0 0.07768758
```

### PanelEstimate

Once proper matched sets are attained by `PanelMatch`, users can
estimate the causal quantity of interest such as the average
treatment effect using `PanelEstimate`. Users can estimate the contemporaneous effect as well as
long-term effects. In this example, we illustrate the use of
`PanelEstimate` to estimate the average treatment effect on treated units (att) at time `t` on the outcomes from time `t+0` to `t+4`.

```r
PE.results <- PanelEstimate(sets = PM.results, panel.data = dem.panel, se.method = "bootstrap")
```

The `PanelEstimate` function returns a `PanelEstimate` object, which is a named list. This object will contain the point estimates, standard errors and other information about the calculations.

Users can easily obtain and visualize important information about estimates and standard errors using the `summary` and `plot` methods for PanelEstimate objects

```r
summary(PE.results)
      estimate std.error      2.5%     97.5%
t+0 -0.8760766 0.9135571 -2.676931 0.9031642
t+1 -1.5782041 1.6476270 -4.925601 1.6090211
t+2 -1.1811249 2.2342724 -5.745626 3.0826224
t+3 -0.6618564 2.6058876 -6.125138 4.2956798
t+4 -1.2512349 2.8367039 -7.107129 4.1198900
```

```r
plot(PE.results)
```
![](https://github.com/insongkim/repo-data/blob/master/panelmatch/pe_plot_6_22.png)
