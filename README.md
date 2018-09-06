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
                           covars.lagged = c("tradewb"),
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

|                                                       |        t+0|        t+1|       t+2|       t+3|       t+4|
|:------------------------------------------------------|----------:|----------:|---------:|---------:|---------:|
|Point Estimate(s)                                      |  0.6882587|  1.1526035|  1.529355|  1.911578|  2.060706|
|Standard Error(s)                                      |  0.6874394|  1.0909965|  1.441706|  1.812096|  2.229936|
|Lower Limit of 95 % Regular Confidence Interval        | -0.5811716| -0.9481110| -1.305701| -1.642694| -2.544959|
|Upper Limit of 95 % Regular Confidence Interval        |  1.9722803|  3.0828392|  4.275508|  5.302505|  6.379517|
|Bias-corrected Estimate(s)                             |  0.6889523|  1.1386178|  1.509317|  1.921500|  2.091011|
|Lower Limit of 95 % Bias-corrected Confidence Interval | -0.5957630| -0.7776321| -1.216798| -1.479350| -2.258106|
|Upper Limit of 95 % Bias-corrected Confidence Interval |  1.9576889|  3.2533181|  4.364411|  5.465850|  6.666370|
```

## 9/6 Update Placeholder Title

Let's create a small, easy to work with form of our example data set. We'll do this by getting data for only the first 10 countries. 

```{r}
uid <-unique(dem$wbcode2)[1:10]
subdem <- dem[dem$wbcode2 %in% uid, ]
DisplayTreatment(unit.id = "wbcode2", time.id = "year", treatment = 'dem', data = subdem)
```
We can use the `findAllTreated` function to figure out which units received a treatment over a time period and the time at which the treatment was "given." 
We provide a matrix of data we want to search for control units, and the column names that correspond to the treatment, time, and unit identifier variables.

These units might have corresponding sets of matched control units. We can use the `get.matchedsets` function to find this set of matched control units. Each matched set will therefore correspond with a unique identifier of the treated unit, and time (ie. the time at which the treated unit received treatment). You must also specify a lag window, which will be used to match treated units to control units. Treated units receive treatment at time T, whereas control units do not. However, besides this, the control and treated units have identical treatment histories over the time specified by the lag variable, from time T-L, to T-1.

The `get.matchedsets` function will return a "matched.set" object, which is just a named list with some other attributes saved. By default, it prints out like a data frame, but its subsetting behaviors work like those of a list. The names of each entry in the matched.set object will specify relevant information about the treated unit in the format `[unit id variable].[time at which treatment was received]`


```{r}
treateds <- findAllTreated(dmat = subdem, treatedvar = "dem", time.var = "year", unit.var = "wbcode2")


msets <- get.matchedsets(t = treateds$year, id = treateds$wbcode2, data = subdem, L = 4, t.column = "year", id.column = "wbcode2", treatedvar = "dem")

names(msets)

[1] "4.1992"  "4.1997"  "6.1973"  "6.1983"  "7.1991"  "7.1998"  "12.1992" "13.2003"
```

```{r}
#data frame printing view: better as a summary view
print(msets)

  wbcode2 year matched.set.size
1       4 1992                3
2       4 1997                1
3       6 1973                5
4       6 1983                6
5       7 1991                5
6       7 1998                0
7      12 1992                3
8      13 2003                3
```
```{r}
#prints as a list, shows all data at once
print(msets, verbose = T)

$`4.1992`
[1]  2  3 13

$`4.1997`
[1] 7

$`6.1973`
[1]  2  4  7 12 13

$`6.1983`
[1]  2  3  4  7 12 13

$`7.1991`
[1]  2  3  4 12 13

$`7.1998`
numeric(0)

$`12.1992`
[1]  2  3 13

$`13.2003`
[1]  2  3 12

attr(,"lag")
[1] 4
attr(,"t.var")
[1] "year"
attr(,"id.var")
[1] "wbcode2"
attr(,"treated.var")
[1] "dem"
```

```{r}
#returns a "matched.set" object (list) of length 1
print(msets[1])

  wbcode2 year matched.set.size
1       4 1992                3
```

```{r}
#prints the control units in this matched set
print(msets[[1]])
[1]  2  3 13
```

Calling `plot` on a `matched.set` object will display a histogram of the sizes of the matched sets. the `summary` function provides a variety of information about the sizes of matched sets, number of empty sets, lag size, and a compact "overview" data frame showing matched set sizes. The `summary` function also has an option to print only the overview data frame. Toggle this by setting the `verbose` argument to `FALSE`.
```{r}
plot(msets)
```

```{r}
print(summary(msets))
$overview
  wbcode2 year matched.set.size
1       4 1992                3
2       4 1997                1
3       6 1973                5
4       6 1983                6
5       7 1991                5
6       7 1998                0
7      12 1992                3
8      13 2003                3

$set.size.summary
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00    2.50    3.00    3.25    5.00    6.00 

$number.of.treated.units
[1] 8

$num.units.empty.set
[1] 1

$lag
[1] 4

print(summary(msets, verbose = FALSE))

$overview
  wbcode2 year matched.set.size
1       4 1992                3
2       4 1997                1
3       6 1973                5
4       6 1983                6
5       7 1991                5
6       7 1998                0
7      12 1992                3
8      13 2003                3

```
Passing a matched set (one treated unit at a particular time and its corresponding set of controls) to the `DisplayTreatment` function will  highlight the treatment histories used to create the matched set. If you set the `show.set.only` argument to `TRUE`, then only units from the matched set (and the treated unit) will be shown on the plot. This is useful when working with larger data sets.

```{r}
DisplayTreatment(unit.id = "wbcode2", time.id = "year", treatment = 'dem', data = subdem, matched.set = msets[1])
DisplayTreatment(unit.id = "wbcode2", time.id = "year", treatment = 'dem', data = subdem, matched.set = msets[1], show.set.only = T)
```
