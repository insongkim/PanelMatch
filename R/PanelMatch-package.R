#' Matching Methods for Causal Inference with Time-Series Cross-Section Data
#'
#' \code{PanelMatch} provides a set of methodological tools that
#' enable researchers to apply matching methods to time-series
#' cross-section data. Imai, Kim, and Wang (2018) proposes a
#' nonparametric generalization of difference-in-differences
#' estimator, which does not rely on the linearity assumption as often
#' done in practice. Researchers first select a method of matching
#' each treated observation from a given unit in a particular time
#' period with control observations from other units in the same time
#' period that have a similar treatment and covariate history. These
#' methods include standard matching methods based on propensity score
#' and Mahalanobis distance as well as weighting methods such as
#' synthetic controls. Once matching is done, both short-term and
#' long-term average treatment effects for the treated can be
#' estimated with standard errors. The package also offers a
#' visualization technique that allows researchers to assess the
#' quality of matches by examining the resulting covariate balance.
#'
#' \tabular{ll}{ Package: \tab PanelMatch\cr Type: \tab Package\cr Version: \tab 0.0.1-\cr
#' Date: \tab 2018-03-01\cr License: \tab GPL (>= 3)\cr }
#'
#' @name PanelMatch-package
#' @useDynLib PanelMatch
#' @aliases PanelMatch-package 
#' @docType package
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#' 
#' Maintainer: In Song Kim \email{insong@mit.edu}
#' @references Imai, Kosuke, In Song Kim and Erik Wang. (2018)
#' "Matching Methods for Causal Inference with Time-Series
#' Cross-Section Data." Working paper.
#' @keywords package
#' @import ggplot2 Synth stats MASS
#' @importFrom Rcpp sourceCpp
#' @importFrom DataCombine slide MoveFront
#' @importFrom CBPS CBPS
#' @importFrom lasso2 merge.formula
#' @importFrom data.table rbindlist
#' @importFrom plyr rbind.fill
#' @importFrom reshape2 dcast
#' @importFrom utils capture.output
#' @importFrom MCMCpack rdirichlet
#' @importFrom nloptr nloptr
#' @importFrom Rsolnp solnp
#' @importFrom matrixStats colSds colQuantiles
#' @importFrom knitr kable
#' @importFrom Matrix cBind rBind tcrossprod crossprod Diagonal Matrix drop0
#' @importFrom methods as
#' @importFrom utils flush.console tail object.size
NULL


