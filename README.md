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
standard matching and weighting methods based on propensity score and Mahalanobis distance.
Once matching is done, both short-term and long-term average treatment
effects for the treated can be estimated with standard errors. The
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
treated and control observations can be made.

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
 consists of control observations that have an identical treatment
 history up to a chosen number (`lag`) of years. This number corresponds with the `lag` parameter, which must
 be chosen by the user. Users must also consider various parameters regarding the refinement of created matched sets:
 1) `refinement.method` -- Users may choose between standard propensity score weighting or matching (`ps.weight`, `ps.match`), covariate balanced propensity score weighting or matching (`CBPS.weight`, `CBPS.match`),  and mahalanobis distance matching (`mahalanobis`).
 2) `size.match` -- This sets the maximum number of control units that can be included in a matched set.
 3) `covs.formula` -- This sets which variables are considered in creating the refined matched sets. This can be set to include lagged versions of any variable as well. See the `PanelMatch` documentation for more information about this parameter.

``` r
PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
                         treatment = "dem", refinement.method = "mahalanobis", 
                         data = dem, match.missing = T, 
                         covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), size.match = 5, qoi = "att"
                         ,outcome.var = "y",
                         lead = 0:4, forbid.treatment.reversal = FALSE)

```							
The `PanelMatch` function will return an object of class "PanelMatch". This is a list that contains a few specific elements: First, a matched.set object(s) that has the same name as the provided qoi if the qoi is "att", "atc". If qoi = "ate" then two matched.set objects will be attached, named "att" and "atc." Users can extract information about individual matched sets as well as statistics about all created matched sets from this object. Consult the [Wiki page on Matched Set Objects](https://github.com/insongkim/PanelMatch/wiki/Matched-Set-Objects) for a much more detailed walkthrough and description of these objects.

The `PanelMatch` object also has some additional attributes: "qoi", "lead", "forbid.treatment.reversal" (a logical value that is the same as what was specified in the function call), and "outcome.var" (character value that is the same as what was specified in the function call)

You can check covariate balance using the `get_covariate_balance` function:

```{r}
get_covariate_balance(PM.results$att, dem, covariates = c("tradewb"), plot = F, ylim = c(-2,2))
      tradewb
t_4 0.2095181
t_3 0.1711747
t_2 0.1234327
t_1 0.1172117
t_0 0.1189729
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
PE.results <- PanelEstimate(inference = "bootstrap", sets = PM.results, 
                            data = dem)
```

The `PanelEstimate` function returns a `PanelEstimate` object, which is a named list. This object will contain the point estimates, standard errors and other information about the calculations. See the wiki page about PanelEstimate objects for more information. 

Users can easily obtain and visualize important information about esimtates and standard errors using the `summary` and `plot` methods for PanelEstimate objects

```r
summary(PE.results)
Weighted Difference-in-Differences with Mahalanobis Distance
Weighted Difference-in-Differences with Covariate Balancing Propensity Score
Matches created with 4 lags

Standard errors computed with 1000 Weighted bootstrap samples

Estimate of Average Treatment Effect on the Treated (ATT) by Period:
$summary
      estimate std.error      2.5%     97.5%
t+0 -0.7932684 0.6233449 -2.029409 0.4221848
t+1 -0.3090152 1.0489291 -2.310423 1.8106643
t+2  0.8810066 1.3762887 -1.847383 3.5910296
t+3  1.8564952 1.6584296 -1.692733 5.1152676
t+4  1.8219732 1.8112170 -1.906964 5.4912805

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
![](https://github.com/insongkim/repo-data/blob/master/panelmatch/peplot5-2.png)
