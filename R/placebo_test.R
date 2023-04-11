#' placebo_test
#'
#' Calculates results for a placebo test
#'
#'
#' Calculate the results of a placebo test, looking at the change in outcome at time = t-1, compared to other pre-treatment periods in the lag window.
#' @param pm.obj an object of class \code{PanelMatch}
#' @param data data.frame with the original data
#' @param lag.in integer indicating earliest the time period(s) in the future for which the placebo test change in outcome will be calculated. Calculations will be made over the period t - max(lag) to t-2, where t is the time of treatment. The results are similar to those returned by PanelEstimate(), except t-1 is used as the period of comparison, rather than the lead window.
#' @param number.iterations integer specifying the number of bootstrap iterations
#' @param confidence.level confidence level for the calculated standard error intervals
#' @param plot logical indicating whether or not a plot should be generated, or just return the raw data from the calculations
#' @param se.method character string describing the type of standard error to be used. Valid inputs include "bootstrap", "conditional" and "unconditional". When the QOI is ATE, only bootstrap can be used.
#' @param ... extra arguments to be passed to plot()
#'
#' @return list with 2 or 3 elements: "estimates", which contains the point estimates for the test, "standard.errors" which has the standard errors for each period and optionally "bootstrapped.estimates", containing the bootstrapped point estimates for the test for each specified lag window period.
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
#'                          placebo.test = TRUE)
#' placebo_test(PM.results, data = dem, number.iterations = 100, plot = FALSE)
#' 
#' 
#' @export
#'
#'
placebo_test <- function(pm.obj,
                         data,
                         lag.in = NULL,
                         number.iterations = 1000,
                         confidence.level = .95,
                         plot = FALSE,
                         se.method = "bootstrap",
                         ...)
{
  df.adjustment <- FALSE
  qoi.in <- attr(pm.obj, "qoi")
  if (is.null(lag.in))
  {
    if (identical(qoi.in, "ate"))
    {
      #just pick one
      matchedsets <- pm.obj[["att"]]
      lag.in <- attr(matchedsets, "lag")
    } else {
      matchedsets <- pm.obj[[qoi.in]]
      lag.in <- attr(matchedsets, "lag")
    }
  }
  
  if (identical(qoi.in, "ate"))
  {
    #just pick one for boundary check
    matchedsets <- pm.obj[["att"]]
  }
  
  if (lag.in == 1) stop("placebo test cannot be conducted for lag = 1")
  if (length(lag.in) >1) stop("lag.in should be a single integer")
  if (lag.in > attr(matchedsets, "lag")) stop("provided lag.in value exceeds lag parameter from matching stage. Please specify a valid lag.in value, such that lag.in < lag")
  
  lag.in <- lag.in:2
  
  
  placebo.results.raw <- panel_estimate(sets = pm.obj,
                                        data = data,
                                        number.iterations = number.iterations,
                                        df.adjustment = df.adjustment,
                                        placebo.test = TRUE,
                                        placebo.lead = lag.in,
                                        confidence.level = confidence.level,
                                        se.method = se.method)
  
  if (plot)
  {
    plot(placebo.results.raw, ...)
  } else {
    if (identical(se.method,"bootstrap"))
    {
      colnames(placebo.results.raw$bootstrapped.estimates) <- 
        names(placebo.results.raw$estimates)
      ses <- apply(placebo.results.raw$bootstrapped.estimates, 
                   2, 
                   sd, 
                   na.rm = TRUE)
      ret.results <- list(estimates = placebo.results.raw$estimates,
                          bootstrapped.estimates = placebo.results.raw$bootstrapped.estimates,
                          standard.errors = ses)
    } else
    {
      ret.results <- list(estimates = placebo.results.raw$estimates,
                          standard.errors = placebo.results.raw$standard.error)
    }
    
    return(ret.results)
  }
}