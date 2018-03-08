# panelMatch: Matching Methods for Causal Inference with Time-Series Cross-Section Data [![Build Status](https://travis-ci.org/insongkim/PanelMatch.svg?branch=master)](https://travis-ci.org/insongkim/panelMatch)

This R package provides a set of methodological tools that enable
researchers to apply matching methods to time-series cross-section
data.  Imai, Kim, and Wang (2018) proposes a nonparametric
generalization of difference-in-differences estimator, which does not
rely on the linearity assumption as often done in
practice. Researchers first select a method of matching each treated
observation from a given unit in a particular time period with control
observations from other units in the same time period that have a
similar treatment and covariate history.  These methods include
standard matching methods based on propensity score and Mahalanobis
distance as well as weighting methods such as synthetic controls.
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

You can also install the most recent development version of `panelMatch` using the `devtools` package. First you have to install `devtools` using the following code. Note that you only have to do this once:

``` r
if(!require(devtools)) install.packages("devtools")
```

Then, load `devtools` and use the function `install_github()` to install `panelMatch`:

``` r
library(devtools)
install_github("insongkim/panelMatch", dependencies=TRUE)
```
