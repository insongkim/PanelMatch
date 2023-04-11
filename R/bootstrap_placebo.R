# calculates standard errors for each period in the placebo lead window using
# the bootstrap

# returns a matrix of bootstrapped coefficients
handle_bootstrap_placebo <- function(qoi.in, 
                             data.in, 
                             placebo.lead,
                             number.iterations,
                             att.treated.unit.ids,
                             atc.treated.unit.ids,
                             outcome.variable,
                             unit.id.variable,
                             confidence.level,
                             lag) 
{

    if ( identical(qoi.in, "att") ||
         identical(qoi.in, "atc") ||
         identical(qoi.in, "art"))
    {
      
      coefs <- matrix(NA, 
                      nrow = number.iterations, 
                      ncol = length(placebo.lead))
      
      for (k in 1:number.iterations)
      {
        # make new data
        clusters <- unique(data.in[, unit.id.variable])
        units <- sample(clusters, size = length(clusters), replace = TRUE)
        if (identical(qoi.in, "att") || 
            identical(qoi.in, "art"))
        {
          treated.unit.ids <- att.treated.unit.ids
        } else {
          treated.unit.ids <- atc.treated.unit.ids
        }
        while (all(!units %in% treated.unit.ids)) #while none of the units are treated units, resample
        {
          units <- sample(clusters, size = length(clusters), replace = TRUE)
        }
        # build the bootstrap samples
        df.bs <- lapply(units,
                        function(x) which(data.in[, unit.id.variable] == x))
        d.sub1 <- data.in[unlist(df.bs),]
        y.in <- d.sub1[, paste0(outcome.variable, "l", 
                                (placebo.lead - 1)), drop = FALSE]
        z.in <- d.sub1[, paste0("dits_", qoi.in)]
        
        col.idx <- sapply(placebo.lead - 2, 
                          function(x) paste0("Wit_", qoi.in, x))
        # calculate the point estimate(s) using the bootstrap sample
        at__new <- equality_four_placebo(d.sub1[, col.idx, drop = FALSE],
                                         y = y.in,
                                         z = z.in)
        if (identical(qoi.in, "atc")) at__new <- -at__new
        coefs[k,] <- at__new
      }
      
    } else if (identical(qoi.in, "ate")) 
    { #ate special because it requires both att and atc calculation
      coefs <- matrix(NA, 
                      nrow = number.iterations, 
                      ncol = length(placebo.lead))
      
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
        ## build data for doing att and atc
        df.bs <- lapply(units, 
                        function(x) which(data.in[, unit.id.variable] == x))
        d.sub1 <- data.in[unlist(df.bs),]
        y.in <- d.sub1[, paste0(outcome.variable, "l", 
                                (placebo.lead - 1)), drop = FALSE]
        z.in <- d.sub1[, paste0("dits_", "att")]
        tidx <- sapply(placebo.lead - 2, 
                       function(x) paste0("Wit_", "att", x))
        att_new <- equality_four_placebo(d.sub1[,tidx, drop = FALSE],
                                         y = y.in, z = z.in)
        y.in <- d.sub1[, paste0(outcome.variable, "l", 
                                (placebo.lead - 1)), drop = FALSE]
        z.in <- d.sub1[, paste0("dits_", "atc")]
        tidx <- sapply(placebo.lead - 2, 
                       function(x) paste0("Wit_", "atc", x))
        atc_new <- equality_four_placebo(d.sub1[, tidx, drop = FALSE], 
                                         y = y.in, z = z.in)
        atc_new <- -atc_new
        # construct the ate estimates from the att and atc estimates
        coefs[k,] <- (att_new*sum(d.sub1$dits_att) + 
                        atc_new*sum(d.sub1$dits_atc))/
          (sum(d.sub1$dits_att) + sum(d.sub1$dits_atc))
      }
      
    } else {
      stop("invalid qoi")
      
    }
  return(coefs)
   
}

