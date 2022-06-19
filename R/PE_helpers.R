prepareData <- function(data.in, lead, sets.att = NULL,
                        sets.atc = NULL, continuous.treatment,
                        qoi.in, dependent.variable)
{
  if ( identical(qoi.in, "att") || identical(qoi.in, "art") || identical(qoi.in, "ate"))
  {
    if (!identical(qoi.in, "ate")) qoi.t <- qoi.in
    if (identical(qoi.in, "ate")) qoi.t <- "att"
    for (j in lead)
    {
      
      dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.att,
                            continuous.treatment = continuous.treatment)
      data.in = merge(x = data.in, y = dense.wits, all.x = TRUE,
                   by.x = colnames(data.in)[1:2], by.y = c("id", "t"))
      colnames(data.in)[length(data.in)] <- paste0("Wit_", qoi.t, j)
      data.in[is.na(data.in[, length(data.in)]), length(data.in)] <- 0 #replace NAs with zeroes
    }

    data.in[, paste0("dit_", qoi.t)] <- getDits(matched_sets = sets.att,
                                                data = data.in)
    colnames(data.in)[length(data.in)] <- paste0("dits_", qoi.t)
    data.in[, paste0("Wit_", qoi.t, "-1")] <- 0

  }
  if (identical(qoi.in, "atc") || identical(qoi.in, "ate"))
  {

    for (j in lead)
    {
      
      dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.atc,
                            continuous.treatment = continuous.treatment)
      data.in = merge(x = data.in, y = dense.wits, all.x = TRUE,
                   by.x = colnames(data.in)[1:2], by.y = c("id", "t"))
      colnames(data.in)[length(data.in)] <- paste0("Wit_atc", j)
      data.in[is.na(data.in[, length(data.in)]), length(data.in)] <- 0 #replace NAs with zeroes
    }

    data.in$dit_atc <- getDits(matched_sets = sets.atc, data = data.in)
    colnames(data.in)[length(data.in)] <- "dits_atc"
    data.in$`Wit_atc-1` <- 0

  }
  #NOTE THE COMMENT/ASSUMPTION

  data.in[, dependent.variable][is.na(data.in[, dependent.variable])] <- 0 #replace the NAs with zeroes.
  #I think this is ok because the dits should always be zero for these, so the value is irrelevant.
  return(data.in)

}


calculateEstimates <- function(qoi.in, data.in, lead,
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
                                pooled = FALSE)
{
  
  pt.estimates  <- calculate_point_estimates(qoi.in, data.in, 
                                             lead, outcome.variable, pooled)
  if (identical(se.method, "bootstrap"))
  {
    
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
    
    if (identical(qoi.in, "att") ||
        identical(qoi.in, "art"))
    {
      sets <- att.sets
    } else if (identical(qoi.in, "atc")) {
      sets <- atc.sets
    }
    ses <- apply(coefs, 2, sd, na.rm = T)
    
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
                "qoi" = "ate", "matched.sets" = list(att = att.sets, atc = atc.sets),
                "se.method" = se.method)
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
        "se.method" = se.method
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
              "se.method" = se.method)
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
              "se.method" = se.method)
    class(z) <- "PanelEstimate"
    return(z)
    
  } else {
    stop("invalid standard error method")
  }
}
  
calculatePlaceboEstimates <- function(qoi.in, data.in, lead,
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
                                        se.method = "bootstrap")
{
    
    if (se.method == "bootstrap")
    {
      if ( identical(qoi.in, "att") ||
           identical(qoi.in, "atc") ||
           identical(qoi.in, "art"))
      {
        
        col.idx <- sapply(placebo.lead - 2, function(x) paste0("Wit_", qoi.in, x))
        x.in <- data.in[, col.idx, drop = FALSE]
        
        #y.in <- data.in[c(outcome.variable)][,1]
        
        create.lagged.dfs <- function(d, dv, idx, k)
        {
          d[, paste0(dv, "l", idx)] <- lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
          return(d)
        }
        y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                   FUN = create.lagged.dfs,
                   dv = outcome.variable,
                   idx = placebo.lead - 1,
                   k = placebo.lead - 1,
                   simplify = FALSE)
        
        data.in <- do.call(rbind, y.in)
        y.in <- data.in[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
        z.in <- data.in[, paste0("dits_", qoi.in)]
        
        o.coefs <- equality_four_placebo(x.in, y.in, z.in)
        
        #do coefficient flip for atc
        if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
        
        
        coefs <- matrix(NA, nrow = number.iterations, ncol = length(placebo.lead))
        
        for (k in 1:number.iterations)
        {
          # make new data
          clusters <- unique(data.in[, unit.id.variable])
          units <- sample(clusters, size = length(clusters), replace = TRUE)
          if (identical(qoi.in, "att") || identical(qoi.in, "art"))
          {
            treated.unit.ids <- att.treated.unit.ids
          } else {
            treated.unit.ids <- atc.treated.unit.ids
          }
          while (all(!units %in% treated.unit.ids)) #while none of the units are treated units, resample
          {
            units <- sample(clusters, size = length(clusters), replace = TRUE)
          }
          
          df.bs <- lapply(units, function(x) which(data.in[, unit.id.variable] == x))
          d.sub1 <- data.in[unlist(df.bs),]
          
          y.in <- d.sub1[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
          z.in <- d.sub1[, paste0("dits_", qoi.in)]
          
          
          at__new <- equality_four_placebo(d.sub1[, col.idx, drop = FALSE],
                                           y = y.in,
                                           z = z.in)
          if (identical(qoi.in, "atc")) at__new <- -at__new
          coefs[k,] <- at__new
        }
        
        if (identical(qoi.in, "att") || identical(qoi.in, "art"))
        {
          sets <- att.sets
        } else {
          sets <- atc.sets
        }
        
        names(o.coefs) <- paste0("t-", placebo.lead)
        #o.coefs <- rev(o.coefs) # for ease
        z <- list("estimates" = o.coefs,
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
        col.idx <- sapply(placebo.lead - 2, function(x) paste0("Wit_", "att", x))
        x.in <- data.in[, col.idx, drop = FALSE]
        
        #y.in <- data.in[c(outcome.variable)][,1]
        
        create.lagged.dfs <- function(d, dv, idx, k)
        {
          d[, paste0(dv, "l", idx)] <- lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
          return(d)
        }
        y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                   FUN = create.lagged.dfs,
                   dv = outcome.variable,
                   idx = placebo.lead - 1,
                   k = placebo.lead - 1,
                   simplify = FALSE)
        
        data.in <- do.call(rbind, y.in)
        y.in <- data.in[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
        z.in <- data.in[, paste0("dits_", "att")]
        
        att.coefs <- equality_four_placebo(x.in, y.in, z.in)
        
        col.idx <- sapply(placebo.lead - 2, function(x) paste0("Wit_", "atc", x))
        x.in <- data.in[, col.idx, drop = FALSE]
        
        #y.in <- data.in[c(outcome.variable)][,1]
        
        create.lagged.dfs <- function(d, dv, idx, k)
        {
          d[, paste0(dv, "l", idx)] <- lapply(k, function(x) data.table::shift(d[, dv], n = x, type = "lag"))
          return(d)
        }
        y.in <- by(data.in, as.factor(data.in[, unit.id.variable]),
                   FUN = create.lagged.dfs,
                   dv = outcome.variable,
                   idx = placebo.lead - 1,
                   k = placebo.lead - 1,
                   simplify = FALSE)
        
        data.in <- do.call(rbind, y.in)
        y.in <- data.in[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
        z.in <- data.in[, paste0("dits_", "atc")]
        
        atc.coefs <- equality_four_placebo(x.in, y.in, z.in)
        atc.coefs <- -atc.coefs
        
        
        o.coefs_ate <- (att.coefs*sum(data.in$dits_att) + atc.coefs*sum(data.in$dits_atc))/
          (sum(data.in$dits_att) + sum(data.in$dits_atc))
       
        
        
        coefs <- matrix(NA, nrow = number.iterations, ncol = length(placebo.lead))
        
        for (k in 1:number.iterations)
        {
          # make new data
          clusters <- unique(data.in[, unit.id.variable])
          units <- sample(clusters, size = length(clusters), replace = TRUE)
          
          while(all(!units %in% att.treated.unit.ids) || 
                all(!units %in% atc.treated.unit.ids)) #while none of the units are treated units (att and atc), resample
          {
            units <- sample(clusters, size = length(clusters), replace=TRUE)
          }
          
          
          
          df.bs <- lapply(units, function(x) which(data.in[, unit.id.variable] == x))
          d.sub1 <- data.in[unlist(df.bs),]
          
          y.in <- d.sub1[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
          z.in <- d.sub1[, paste0("dits_", "att")]
          
          
          att_new <- equality_four_placebo(d.sub1[, sapply(placebo.lead - 2, function(x) paste0("Wit_", "att", x)), drop = FALSE],
                                           y = y.in,
                                           z = z.in)
          
          
          y.in <- d.sub1[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
          z.in <- d.sub1[, paste0("dits_", "atc")]
          
          
          atc_new <- equality_four_placebo(d.sub1[, sapply(placebo.lead - 2, function(x) paste0("Wit_", "atc", x)), drop = FALSE],
                                           y = y.in,
                                           z = z.in)
          
          atc_new <- -atc_new
          coefs[k,] <- (att_new*sum(d.sub1$dits_att) + atc_new*sum(d.sub1$dits_atc))/
            (sum(d.sub1$dits_att) + sum(d.sub1$dits_atc))
        }
        
        
        
        names(o.coefs_ate) <- paste0("t-", placebo.lead)
        #o.coefs <- rev(o.coefs) # for ease
        colnames(coefs) <- names(o.coefs_ate)
        z <- list("estimates" = o.coefs_ate,
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
    } else if (identical(se.method, "conditional") || 
               identical(se.method, "unconditional")) 
    {
      stop("not currently implemented")
    } else {
    stop("placebo tests only currently available with bootstrap!")
  }
  
  
}
