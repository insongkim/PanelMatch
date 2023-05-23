#' Get summaries of PanelEstimate objects/calculations
#' 
#' 
#' \code{summary.PanelEstimate} takes an object returned by
#' \code{PanelEstimate}, and returns a summary table of point
#' estimates and confidence intervals
#'
#' @param object A PanelEstimate object
#' @param verbose logical indicating whether or not output should be printed in an expanded form. Default is TRUE
#' @param bias.corrected logical indicating whether or not bias corrected estimates should be provided. Default is FALSE. This argument only applies for standard errors calculated with the bootstrap. 
#' @param ... optional additional arguments. Currently, no additional arguments are supported. 
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "none", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem, number.iterations = 500)
#' summary(PE.results)
#' 
#'
#' 
#' @method summary PanelEstimate
#' @export
summary.PanelEstimate <- function(object, verbose = TRUE, 
                                  bias.corrected = FALSE, ...) {
  
  if (verbose)
  {
    if (object$qoi == "ate")
    {
      refinement.method <- attr(object$matched.set[[1]], "refinement.method")
      lag <- attr(object$matched.set[[1]], "lag")
    }
    else
    {
      refinement.method <- attr(object$matched.set, "refinement.method")
      lag <- attr(object$matched.set, "lag")
    }
    if (refinement.method == "mahalanobis")
    {
      cat("Weighted Difference-in-Differences with Mahalanobis Distance\n")
    } 
    if (refinement.method == "ps.weight" || 
        refinement.method == "ps.match") 
    {
      cat("Weighted Difference-in-Differences with Propensity Score\n")
    } 
    if (refinement.method == "CBPS.weight" || 
        refinement.method == "CBPS.match") 
    {
      cat("Weighted Difference-in-Differences with Covariate Balancing Propensity Score\n")
    }
    cat("Matches created with", attr(object$matched.sets, "lag"), "lags\n")
    if (!is.null(object$bootstrap.iterations))
    {
      cat("\nStandard errors computed with", 
          object$bootstrap.iterations, "Weighted bootstrap samples\n")  
    }
    
    if (identical(object$se.method, "conditional") || 
        identical(object$se.method, "unconditional"))
    {
      cat("\nStandard errors computed with", 
          object$se.method, "method\n")
    }
    
    if (object$qoi == "att")
    {
      qoi <- "Average Treatment Effect on the Treated (ATT)"
    }
    if (object$qoi == "atc")
    {
      qoi <- "Average Treatment Effect on the Control (ATC)"
    }
    if (object$qoi == "ate")
    {
      qoi <- "Average Treatment Effect (ATE)"
    }
    
    if (object$qoi == "art")
    {
      qoi <- "Average Effect of Treatment Reversal on the Reversed (ART)"
    }
    
    cat("\nEstimate of", qoi, "by Period:\n")
  }
  if (bias.corrected && identical(object$se.method, "bootstrap"))
  {
    if (!identical(object$se.method, "bootstrap")) stop("bias corrected estimates only available for bootstrap method currently")
    a1 <- t(as.data.frame(object$estimates))
    a2 <- apply(object$bootstrapped.estimates, 2, sd, na.rm = TRUE)
    a3 <- apply(object$bootstrapped.estimates, 2, 
          quantile, 
          probs = c( (1 - object$confidence.level)/2, 
                     object$confidence.level + (1 - object$confidence.level)/2 ), 
          na.rm = TRUE)
    a4 <- 2*object$estimates - colMeans(object$bootstrapped.estimates, na.rm = TRUE)
    a5 <- apply( (2*matrix(nrow = object$bootstrap.iterations, 
                     ncol = length(object$estimates), 
                     object$estimates, byrow = TRUE) - object$bootstrapped.estimates), 
           2, 
           quantile, 
           probs = c((1 - object$confidence.level)/2, 
                     object$confidence.level + (1 - object$confidence.level)/2), 
           na.rm = TRUE)
    df <- rbind(a1, # point estimate
                a2, # bootstrap se
                # Efron & Tibshirani 1993 p170 - 171
                a3, # percentile confidence.level
                # Efron & Tibshirani 1993 p138
                a4, # bc point estimate
                a5) # bc percentile confidence.level)
    
    rownames(df) <- c("estimate", "std.error", 
                      paste0((1 - object$confidence.level)/2 * 100, "%"),
                      paste0( (object$confidence.level + (1 - object$confidence.level)/2) * 100, "%"),
                      "estimate(bias corrected)", 
                      paste0((1 - object$confidence.level)/2 * 100, "%", "(bias corrected)"),
                      paste0((object$confidence.level + (1 - object$confidence.level)/2) * 100, "%", "(bias corrected)"))
  }
  else
  {
    if ( identical(object$se.method, "bootstrap") )
    {
      a1 <- t(as.data.frame(object$estimates))
      a2 <- apply(object$bootstrapped.estimates, 2, sd, na.rm = TRUE)
      a3 <- apply(object$bootstrapped.estimates, 2, quantile, 
            probs = c( (1 - object$confidence.level)/2, 
                       object$confidence.level + (1 - object$confidence.level)/2 ), 
            na.rm = TRUE)
      df <- rbind(a1, # point estimate
                  a2, # bootstrap se
                  a3) # Efron & Tibshirani 1993 p170 - 171
      
      rownames(df) <- c("estimate", "std.error", 
                        paste0((1 - object$confidence.level)/2 * 100, "%"),
                        paste0( (object$confidence.level + (1 - object$confidence.level)/2) * 100, "%"))
      tdf <- t(df)
      if (!verbose) return(t(df))
      return(list("summary" = tdf, 
                  "lag" = lag, 
                  "iterations" = object$bootstrap.iterations, 
                  "qoi" = object$qoi) )   
    }
    else if (identical(object$se.method, "conditional") || 
             identical(object$se.method, "unconditional")) 
    {
      
      critical.vals <- rep(qnorm( (1 - object$confidence.level) / 2), 2) * c(1, -1)
      
      quants <- mapply(FUN = function(x, y) x + critical.vals * y,
                       x = object$estimates, 
                       y = object$standard.error,
                       SIMPLIFY = FALSE)
      qts <- do.call(rbind, quants)
      df <- data.frame(estimate = object$estimates,
                       std.error = object$standard.error)
      
      tdf <- cbind(df, qts)
      colnames(tdf) <- c("estimate", "std.error", 
                         paste0((1 - object$confidence.level)/2 * 100, "%"),
                         paste0( (object$confidence.level + (1 - object$confidence.level)/2) * 100, "%"))
      rownames(tdf) <- names(object$estimates)
      if (!verbose) return(tdf)
      return(list("summary" = tdf, 
                  "lag" = lag, 
                  "qoi" = object$qoi) )
    } else {
      stop("se.method not specified correctly")
    }
    
  }
  
}


#' Plot point estimates and standard errors from a PanelEstimate calculation.
#' 
#' 
#' The \code{plot.PanelEstimate} method takes an object returned by the \code{PanelEstimate} function and plots the calculated 
#' point estimates and standard errors over the specified \code{lead} time period. 
#' The only mandatory argument is an object of the \code{PanelEstimate} class.
#'
#' @param x a \code{PanelEstimate} object
#' @param ylab default is "Estimated Effect of Treatment." This is the same argument as the standard argument for \code{plot()}
#' @param xlab default is "Time". This is the same argument as the standard argument for \code{plot()}
#' @param main default is "Estimated Effects of Treatment Over Time". This is the same argument as the standard argument for \code{plot}
#' @param ylim default is NULL. This is the same argument as the standard argument for \code{plot()}
#' @param pch default is NULL. This is the same argument as the standard argument for \code{plot()}
#' @param cex default is NULL. This is the same argument as the standard argument for \code{plot()}
#' @param include.placebo Logical. If TRUE, then the results of the placebo test are included on the plot. Note that, by definition, the results for t-1 will be 0. This requires that \code{PanelEstimate()} was run with the placebo test option set to \code{TRUE}. If FALSE, the results are not included. Default is FALSE. Please see \code{placebo_test()} and \code{PanelEstimate()} for more information. 
#' @param ... Additional optional arguments to be passed to \code{plot()}.
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem)
#' plot(PE.results)
#' 
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = TRUE, 
#'                          placebo.test = TRUE)
#' 
#' PE.results <- PanelEstimate(sets = PM.results, 
#'                             data = dem, 
#'                             number.iterations = 100,
#'                             include.placebo.test = TRUE)
#' 
#' plot(PE.results, include.placebo = TRUE)
#'
#' @method plot PanelEstimate
#' @export
plot.PanelEstimate <- function(x, 
                               ylab = "Estimated Effect of Treatment", 
                               xlab = "Time", 
                               main = "Estimated Effects of Treatment Over Time",
                               ylim = NULL, 
                               pch = NULL,
                               cex = NULL,
                               include.placebo = FALSE,
                               ...)
{
  
  
  if (include.placebo)
  {
    if (is.null(x[["placebo.test"]]))
    {
      stop("placebo test data is missing! Please re-run PanelEstimate() with include.placebo.test = TRUE.")
    }
    pt.data <- summarize_placebo_data(x[["placebo.test"]], 
                                      x[["confidence.level"]], 
                                      x[["se.method"]])
    pt.data <- rbind(pt.data, 0)
    rownames(pt.data)[nrow(pt.data)] <- "t-1"
  }
  pe.object <- x
  plot.data <- summary(pe.object, 
                       verbose = FALSE, bias.corrected = FALSE)
  
  if (include.placebo)
  {
    plot.data <- rbind(pt.data, plot.data)
  }
  
  if (is.null(ylim))
  {
    ylim <- c(min(plot.data[, 3]) - abs(mean(plot.data[, 3])),
              max(plot.data[, 4]) + abs(mean(max(plot.data[, 4]))))
  }
  if (is.null(pch))
  {
    pch <- 16
    
  }
  if (is.null(cex))
  {
    cex <- 1.5
  }
  graphics::plot(x = 1:(nrow(plot.data)),y = plot.data[, 1],
                 xaxt = "n", ylab = ylab, xlab = xlab, 
                 main = main, ylim = ylim, pch = pch, cex = cex, ...)
  graphics::axis(side = 1, 
                 at = 1:nrow(plot.data), 
                 labels = rownames(plot.data))
  graphics::segments(1:(nrow(plot.data)), 
                     plot.data[,3], 
                     1:(nrow(plot.data)), 
                     plot.data[,4])
  graphics::abline(h = 0, lty = "dashed")
}
