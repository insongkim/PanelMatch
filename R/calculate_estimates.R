# calculates estimates and standard errors within PanelEstimate()
# calls the appropriate helper functions for each step
# creates and returns a PanelEstimate object
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
      coefs <- handle_bootstrap_parallel(qoi.in, 
                                         data.in, 
                                         lead,
                                         number.iterations,
                                         att.treated.unit.ids,
                                         atc.treated.unit.ids,
                                         outcome.variable,
                                         unit.id.variable,
                                         confidence.level,
                                         lag,
                                         se.method,
                                         pooled,
                                         num.cores = num.cores)
    } else {
      coefs <- handle_bootstrap(qoi.in, 
                                data.in, 
                                lead,
                                number.iterations,
                                att.treated.unit.ids,
                                atc.treated.unit.ids,
                                outcome.variable,
                                unit.id.variable,
                                confidence.level,
                                lag,
                                se.method,
                                pooled) 
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