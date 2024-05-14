#' calculate_estimates
#'
#' Mid-level function that helps with estimation process. Calls lower level helper functions
#' @param qoi.in String specifying qoi
#' @param data.in data.frame object with the data
#' @param lead integer specifying lead window
#' @param number.iterations integer. specifies number of bootstrap iterations
#' @param att.treated.unit.ids Integer vector specifying the treated units for the att or art
#' @param atc.treated.unit.ids Integer vector specifying the "treated" units under the atc definition
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable
#' @param confidence.level double. specifies confidence level for confidence interval
#' @param att.sets matched.set object specifying the att or art sets
#' @param atc.sets matched.set object specifying the atc sets
#' @param lag integer vector specifying size of the lag.
#' @param se.method string specifying which method should be used for standard error calculation
#' @param pooled bool. specifies whether or not estimates should be calculated for each lead period, or pooled across all lead periods
#' @param parallel bool. Specifies whether or not parallelization should be used
#' @param num.cores Integer. specifies how many cores to use for parallelization
#' @return Returns PanelEstimate object.
#' @keywords internal
calculate_estimates <- function(qoi.in, data.in, lead,
                               number.iterations,
                               att.treated.unit.ids,
                               atc.treated.unit.ids,
                               outcome.variable,
                               unit.id.variable,
                               confidence.level,
                               att.sets,
                               atc.sets,
                               placebo.test = FALSE,
                               lag,
                               se.method,
                               pooled = FALSE,
                               parallel = FALSE,
                               num.cores = 1)
{
  
  pt.estimates  <- calculate_point_estimates(qoi.in, data.in, 
                                             lead, outcome.variable, pooled)
  if (identical(se.method, "bootstrap"))
  {
    
    if (parallel) {
      coefs <- handle_bootstrap_parallel(qoi.in = qoi.in, 
                                         data.in = data.in, 
                                         lead = lead,
                                         number.iterations = number.iterations,
                                         att.treated.unit.ids = att.treated.unit.ids,
                                         atc.treated.unit.ids = atc.treated.unit.ids,
                                         outcome.variable = outcome.variable,
                                         unit.id.variable = unit.id.variable,
                                         confidence.level = confidence.level,
                                         lag = lag,
                                         pooled = pooled,
                                         num.cores = num.cores)
    } else {
      coefs <- handle_bootstrap(qoi.in = qoi.in, 
                                data.in = data.in, 
                                lead = lead,
                                number.iterations = number.iterations,
                                att.treated.unit.ids = att.treated.unit.ids,
                                atc.treated.unit.ids = atc.treated.unit.ids,
                                outcome.variable = outcome.variable,
                                unit.id.variable = unit.id.variable,
                                confidence.level = confidence.level,
                                lag = lag,
                                pooled = pooled) 
    }
    
    
    if (identical(qoi.in, "att") ||
        identical(qoi.in, "art"))
    {
      sets <- att.sets
    } else if (identical(qoi.in, "atc")) {
      sets <- atc.sets
    }
    ses <- apply(coefs, 2, sd, na.rm = TRUE)
    
    if (!pooled)
    {
      names(ses) <- paste0("t+", lead)
    }
    if (identical(qoi.in, "ate"))
    {
      z <- list("estimates" = pt.estimates,
                "bootstrapped.estimates" = coefs, 
                "bootstrap.iterations" = number.iterations, 
                "standard.error" = ses,
                "lead" = lead, "confidence.level" = confidence.level, 
                "qoi" = "ate", "matched.sets" = list(att = att.sets, 
                                                     atc = atc.sets),
                "se.method" = se.method,
                "pooled" = pooled)
    } else {
      z <- list(
        "estimates" = pt.estimates,
        "bootstrapped.estimates" = coefs,
        "bootstrap.iterations" = number.iterations,
        "standard.error" = ses,
        "lag" = lag,
        "lead" = lead,
        "confidence.level" = confidence.level,
        "qoi" = qoi.in,
        "matched.sets" = sets,
        "se.method" = se.method,
        "pooled" = pooled
      )
    }
    
    class(z) <- "PanelEstimate"
    return(z)
    
  } else if (identical(se.method, "conditional"))
  {
    if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
    if (pooled) stop("Analytical SE's not available for pooled estimates")
    
    se.s <- handle_conditional_se(qoi.in, data.in, 
                                  lead, outcome.variable, 
                                  unit.id.variable)
    
    if (is.null(atc.sets))
    {
      sets <- att.sets
    } else if (is.null(att.sets)) {
      sets <- atc.sets
    } else {
      stop("missing sets")
    }
    z <- list("estimates" = pt.estimates,
              "standard.error" = se.s,
              "lag" = lag,
              "lead" = lead,
              "confidence.level" = confidence.level,
              "qoi" = qoi.in,
              "matched.sets" = sets,
              "se.method" = se.method,
              "pooled" = pooled)
    class(z) <- "PanelEstimate"
    return(z)
    
  } else if (identical(se.method, "unconditional")) 
  {
    
    se.s <- handle_unconditional_se(qoi.in, data.in, 
                                    lead, outcome.variable, 
                                    unit.id.variable)
    
    if (is.null(atc.sets))
    {
      sets <- att.sets
    } else if (is.null(att.sets)) {
      sets <- atc.sets
    } else {
      stop("missing sets")
    }
    z <- list("estimates" = pt.estimates,
              "standard.error" = se.s,
              "lag" = lag,
              "lead" = lead,
              "confidence.level" = confidence.level,
              "qoi" = qoi.in,
              "matched.sets" = sets,
              "se.method" = se.method,
              "pooled" = pooled)
    class(z) <- "PanelEstimate"
    return(z)
    
  } else {
    stop("invalid standard error method")
  }
}