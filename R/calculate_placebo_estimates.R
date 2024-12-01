#' calculate_placebo_estimates
#'
#' @description Handles the procedures for calculating point estimates and standard errors for the placebo test. Code is structured very similarly to the calculate_estimates() code, but with appropriate modifications for the placebo test. See that function for description of arguments. Bootstrap SEs are available for any specification. Conditional, unconditional standard errors only available for att, art, atc. 
#' @return Returns a PanelEstimate object
#' @keywords internal
calculate_placebo_estimates <- function(qoi.in, data.in, lead,
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
                                      placebo.lead,
                                      se.method = "bootstrap", 
                                      parallel = FALSE, 
                                      num.cores = 1)
{
  
  if (se.method == "bootstrap")
  {
    if ( identical(qoi.in, "att") ||
         identical(qoi.in, "atc") ||
         identical(qoi.in, "art"))
    {
      
      col.idx <- sapply(placebo.lead - 2, function(x) paste0("Wit_", qoi.in, x))
      x.in <- data.in[, col.idx, drop = FALSE]
      
      create.lagged.dfs <- function(d, dv, idx, k)
      {
        d[, paste0(dv, "l", idx)] <- 
          lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
        return(d)
      }
      y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                 FUN = create.lagged.dfs,
                 dv = outcome.variable,
                 idx = placebo.lead - 1,
                 k = placebo.lead - 1,
                 simplify = FALSE)
      
      data.in <- do.call(rbind, y.in)
      y.in <- data.in[, paste0(outcome.variable, "l", 
                               (placebo.lead - 1)), drop = FALSE]
      z.in <- data.in[, paste0("dits_", qoi.in)]
      
      o.coefs <- equality_four_placebo(x.in, y.in, z.in)
      
      #do coefficient flip for atc
      if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
      
      if (isTRUE(parallel)) {
        coefs <- handle_bootstrap_placebo_parallel(qoi.in = qoi.in,
                                          data.in = data.in,
                                          placebo.lead = placebo.lead,
                                          number.iterations = number.iterations,
                                          att.treated.unit.ids = att.treated.unit.ids,
                                          atc.treated.unit.ids = atc.treated.unit.ids,
                                          outcome.variable = outcome.variable,
                                          unit.id.variable = unit.id.variable,
                                          confidence.level = confidence.level,
                                          lag = lag, 
                                          num.cores = num.cores)
      } else {
        coefs <- handle_bootstrap_placebo(qoi.in = qoi.in,
                                          data.in = data.in,
                                          placebo.lead = placebo.lead,
                                          number.iterations = number.iterations,
                                          att.treated.unit.ids = att.treated.unit.ids,
                                          atc.treated.unit.ids = atc.treated.unit.ids,
                                          outcome.variable = outcome.variable,
                                          unit.id.variable = unit.id.variable,
                                          confidence.level = confidence.level,
                                          lag = lag)
      }
      
      
      if (identical(qoi.in, "att") || identical(qoi.in, "art"))
      {
        sets <- att.sets
      } else {
        sets <- atc.sets
      }
      
      names(o.coefs) <- paste0("t-", placebo.lead)
      
      z <- list("estimate" = o.coefs,
                "bootstrapped.estimates" = coefs,
                "bootstrap.iterations" = number.iterations,
                "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lag" = lag,
                "lead" = lead,
                "confidence.level" = confidence.level,
                "qoi" = qoi.in,
                "matched.sets" = sets,
                "se.method" = se.method)
      class(z) <- "PanelEstimate"
      return(z)
      
    } else if (identical(qoi.in, "ate")) 
    {
      col.idx <- sapply(placebo.lead - 2, 
                        function(x) paste0("Wit_", "att", x))
      x.in <- data.in[, col.idx, drop = FALSE]
      
      create.lagged.dfs <- function(d, dv, idx, k)
      {
        d[, paste0(dv, "l", idx)] <- 
          lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
        return(d)
      }
      y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                 FUN = create.lagged.dfs,
                 dv = outcome.variable,
                 idx = placebo.lead - 1,
                 k = placebo.lead - 1,
                 simplify = FALSE)
      
      data.in <- do.call(rbind, y.in)
      y.in <- data.in[, paste0(outcome.variable, "l", 
                               (placebo.lead - 1)), drop = FALSE]
      z.in <- data.in[, paste0("dits_", "att")]
      
      att.coefs <- equality_four_placebo(x.in, y.in, z.in)
      
      col.idx <- sapply(placebo.lead - 2, 
                        function(x) paste0("Wit_", "atc", x))
      x.in <- data.in[, col.idx, drop = FALSE]
      
      
      
      create.lagged.dfs <- function(d, dv, idx, k)
      {
        d[, paste0(dv, "l", idx)] <- 
          lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
        return(d)
      }
      y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                 FUN = create.lagged.dfs,
                 dv = outcome.variable,
                 idx = placebo.lead - 1,
                 k = placebo.lead - 1,
                 simplify = FALSE)
      
      data.in <- do.call(rbind, y.in)
      y.in <- data.in[, paste0(outcome.variable, "l", 
                               (placebo.lead - 1)), drop = FALSE]
      z.in <- data.in[, paste0("dits_", "atc")]
      
      atc.coefs <- equality_four_placebo(x.in, y.in, z.in)
      atc.coefs <- -atc.coefs
      
      
      o.coefs_ate <- (att.coefs*sum(data.in$dits_att) + 
                        atc.coefs*sum(data.in$dits_atc))/
        (sum(data.in$dits_att) + sum(data.in$dits_atc))
      
      
      if (isTRUE(parallel)) {
        coefs <- handle_bootstrap_placebo_parallel(qoi.in = qoi.in,
                                          data.in = data.in,
                                          placebo.lead = placebo.lead,
                                          number.iterations = number.iterations,
                                          att.treated.unit.ids = att.treated.unit.ids,
                                          atc.treated.unit.ids = atc.treated.unit.ids,
                                          outcome.variable = outcome.variable,
                                          unit.id.variable = unit.id.variable,
                                          confidence.level = confidence.level,
                                          lag = lag,
                                          num.cores = num.cores)
      } else {
        coefs <- handle_bootstrap_placebo(qoi.in = qoi.in,
                                          data.in = data.in,
                                          placebo.lead = placebo.lead,
                                          number.iterations = number.iterations,
                                          att.treated.unit.ids = att.treated.unit.ids,
                                          atc.treated.unit.ids = atc.treated.unit.ids,
                                          outcome.variable = outcome.variable,
                                          unit.id.variable = unit.id.variable,
                                          confidence.level = confidence.level,
                                          lag = lag)
      }
      
      
      
    
      names(o.coefs_ate) <- paste0("t-", placebo.lead)
      
      colnames(coefs) <- names(o.coefs_ate)
      z <- list("estimate" = o.coefs_ate,
                "bootstrapped.estimates" = coefs,
                "bootstrap.iterations" = number.iterations,
                "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lag" = lag,
                "lead" = lead,
                "confidence.level" = confidence.level,
                "qoi" = qoi.in,
                "matched.sets" = list(att = att.sets,
                                      atc = atc.sets),
                "se.method" = se.method)
      class(z) <- "PanelEstimate"
      return(z)
    } else {
      stop("invalid qoi")
      
    }
  } else if (se.method %in% c("conditional", "unconditional") && 
             qoi.in %in% c("att", "art", "atc")) {
    
    col.idx <- sapply(placebo.lead - 2, function(x) paste0("Wit_", qoi.in, x))
    x.in <- data.in[, col.idx, drop = FALSE]
    
    create.lagged.dfs <- function(d, dv, idx, k)
    {
      d[, paste0(dv, "l", idx)] <- 
        lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
      return(d)
    }
    y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
               FUN = create.lagged.dfs,
               dv = outcome.variable,
               idx = placebo.lead - 1,
               k = placebo.lead - 1,
               simplify = FALSE)
    
    data.in <- do.call(rbind, y.in)
    y.in <- data.in[, paste0(outcome.variable, "l", 
                             (placebo.lead - 1)), drop = FALSE]
    z.in <- data.in[, paste0("dits_", qoi.in)]
    
    o.coefs <- equality_four_placebo(x.in, y.in, z.in)
    
    #do coefficient flip for atc
    if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
    
    if (identical(se.method, "conditional"))
    {
      ll <- list()
      for (f in placebo.lead) {
        
        l.outcome <- paste0(outcome.variable, "l", f -1)
        data.in[, l.outcome][is.na(data.in[, l.outcome])] <- 0
        ll[[paste0("X", f)]] <- 
          handle_conditional_se(qoi.in = qoi.in,
                                data.in = data.in, lead = f - 2, 
                                outcome.variable = paste0(outcome.variable, "l", f -1),
                                unit.id.variable = unit.id.variable)
      }
      ses <- unlist(ll)
      names(ses) <- paste0("t-", placebo.lead)
    }
    
    if (identical(se.method, "unconditional"))
    {
      ll <- list()
      for (f in placebo.lead) {
        
        l.outcome <- paste0(outcome.variable, "l", f -1)
        data.in[, l.outcome][is.na(data.in[, l.outcome])] <- 0
        ll[[paste0("X", f)]] <- 
          handle_unconditional_se(qoi.in = qoi.in, data.in = data.in,
                                  lead = f - 2, 
                                  outcome.variable = paste0(outcome.variable, "l", f -1),
                                  unit.id.variable = unit.id.variable)
      }
      
      ses <- unlist(ll)
      names(ses) <- paste0("t-", placebo.lead)
      
    }
    if (identical(qoi.in, "att") || 
        identical(qoi.in, "art"))
    {
      sets <- att.sets
    } else {
      sets <- atc.sets
    }
    
    names(o.coefs) <- paste0("t-", placebo.lead)
    
    z <- list("estimate" = o.coefs,
              "standard.error" = ses,
              "lag" = lag,
              "lead" = lead,
              "confidence.level" = confidence.level,
              "qoi" = qoi.in,
              "matched.sets" = sets,
              "se.method" = se.method)
    class(z) <- "PanelEstimate"
    return(z)
  } else {
    stop("qoi and/or se.method not well specified")
  }
  
}

#' equality_four_placebo
#' @description Small helper function implementing estimation function from Imai, Kim, and Wang (2023)
#' @return Returns numeric vector of results.
#' @keywords internal
equality_four_placebo <- function(x, y, z){
  
  y[is.na(y)] <- 0
  res <- colSums(x * y) / sum(z)
  return(res)
}
