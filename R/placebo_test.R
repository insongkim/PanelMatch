#' Conduct a placebo test
#'
#'
#' Calculate the results of a placebo test, looking at the change in outcome at time = t-1, compared to other pre-treatment periods in the lag window.
#' @param pm.obj an object of class \code{PanelMatch}
#' @param panel.data \code{PanelData} object
#' @param lag.in integer indicating earliest the time period(s) in the future for which the placebo test change in outcome will be calculated. Calculations will be made over the period t - max(lag) to t-2, where t is the time of treatment. The results are similar to those returned by \code{PanelEstimate()}, except t-1 is used as the period of comparison, rather than the lead window. If not specified, the placebo test is conducted for periods from t - max(lag) to t-2.
#' @param number.iterations integer specifying the number of bootstrap iterations. This argument only has an effect if standard errors are calculated with the bootstrap.
#' @param confidence.level confidence level for the calculated standard error intervals. Should be specified as a numeric between 0 and 1.
#' @param plot logical indicating whether or not a plot should be generated, or just return the raw data from the calculations
#' @param se.method character string describing the type of standard error to be used. Valid inputs include "bootstrap", "conditional" and "unconditional". When the QOI is ATE, only bootstrap can be used. See the documentation of this argument in \code{PanelEstimate()} for more.
#' @param parallel Logical. If TRUE and \code{se.method = "bootstrap"}, bootstrap procedure will be parallelized. Default is FALSE. If \code{se.method} is not set to \code{bootstrap}, this option does nothing."
#' @param num.cores Integer. Specifies the number of cores to use for parallelization. If \code{se.method = "bootstrap"} and \code{parallel = TRUE}, then this option will take effect. Otherwise, it will do nothing. 
#' @param ... extra arguments to be passed to \code{plot()}
#'
#' @return list with 2 or 3 elements: "estimate", which contains the point estimates for the test, "standard.errors" which has the standard errors for each period and optionally "bootstrapped.estimates", containing the bootstrapped point estimates for the test for each specified lag window period.
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
#' # create subset of data for simplicity
#' PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4, 
#'                          refinement.method = "ps.match", 
#'                          match.missing = TRUE, 
#'                          covs.formula = ~ tradewb,
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, 
#'                          forbid.treatment.reversal = FALSE, placebo.test = TRUE)
#' placebo_test(PM.results, panel.data = dem.sub.panel, se.method = "unconditional", plot = FALSE)
#' 
#' 
#' @export
placebo_test <- function(pm.obj,
                         panel.data,
                         lag.in = NULL,
                         number.iterations = 1000,
                         confidence.level = .95,
                         plot = FALSE,
                         se.method = "bootstrap",
                         parallel = FALSE,
                         num.cores = 1,
                         ...)
{
  if (!inherits(panel.data, "PanelData")) stop("Please provide a PanelData object.")
  if (attr(pm.obj, "placebo.test") == FALSE)
  {
    stop("Placebo test cannot be executed. Please ensure placebo.test = TRUE in PanelMatch()")
  }
  
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
  } else
  {
    if (identical(qoi.in, "ate"))
    {
      #just pick one
      matchedsets <- pm.obj[["att"]]
      
    } else {
      matchedsets <- pm.obj[[qoi.in]]
    }
  }
  
  
  if (lag.in == 1) stop("placebo test cannot be conducted for lag = 1")
  if (length(lag.in) >1) stop("lag.in should be a single integer")
  if (lag.in > attr(matchedsets, "lag")) stop("provided lag.in value exceeds lag parameter from matching stage. Please specify a valid lag.in value, such that lag.in < lag")
  
  lag.in <- lag.in:2
  
  
  placebo.results.raw <- panel_estimate(sets = pm.obj,
                                        data = panel.data,
                                        number.iterations = number.iterations,
                                        df.adjustment = df.adjustment,
                                        placebo.test = TRUE,
                                        placebo.lead = lag.in,
                                        confidence.level = confidence.level,
                                        se.method = se.method,
                                        parallel = parallel,
                                        num.cores = num.cores)
  
  if (plot)
  {
    plot(placebo.results.raw, ...)
  } else {
    if (identical(se.method,"bootstrap"))
    {
      colnames(placebo.results.raw$bootstrapped.estimates) <- 
        names(placebo.results.raw$estimate)
      ses <- apply(placebo.results.raw$bootstrapped.estimates, 
                   2, 
                   sd, 
                   na.rm = TRUE)
      ret.results <- list(estimates = placebo.results.raw$estimate,
                          bootstrapped.estimates = placebo.results.raw$bootstrapped.estimates,
                          standard.errors = ses)
    } else
    {
      ret.results <- list(estimates = placebo.results.raw$estimate,
                          standard.errors = placebo.results.raw$standard.error)
    }
    
    return(ret.results)
  }
}
