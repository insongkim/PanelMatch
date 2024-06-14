#' handle_bootstrap_placebo
#'
#' Helper function for calculating bootstrapped estimates for the placebo test. This version is not parallelized.
#' @param qoi.in String specifying qoi
#' @param data.in data.frame object with the data
#' @param number.iterations integer. specifies number of bootstrap iterations
#' @param att.treated.unit.ids Integer vector specifying the treated units for the att or art
#' @param atc.treated.unit.ids Integer vector specifying the "treated" units under the atc definition
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable
#' @param confidence.level double. specifies confidence level for confidence interval
#' @param lag integer vector specifying size of the lag.
#' @return Returns a matrix of bootstrapped QOI estimate values.
#'
#' @keywords internal
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


#' handle_bootstrap_placebo_parallel
#'
#' Helper function for calculating bootstrapped estimates for the placebo test. This version is parallelized.
#' @param qoi.in String specifying qoi
#' @param data.in data.frame object with the data
#' @param number.iterations integer. Specifies number of bootstrap iterations
#' @param att.treated.unit.ids Integer vector specifying the treated units for the att or art
#' @param atc.treated.unit.ids Integer vector specifying the "treated" units under the atc definition
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable
#' @param confidence.level double. specifies confidence level for confidence interval
#' @param lag integer vector specifying size of the lag.
#' @param num.cores number of cores to be used for parallelization
#' @return Returns a matrix of bootstrapped QOI estimate values.
#'
#' @keywords internal
handle_bootstrap_placebo_parallel <- function(qoi.in, 
                                     data.in, 
                                     placebo.lead,
                                     number.iterations,
                                     att.treated.unit.ids,
                                     atc.treated.unit.ids,
                                     outcome.variable,
                                     unit.id.variable,
                                     confidence.level,
                                     lag, 
                                     num.cores = 1) 
{
  if (identical(num.cores, 1))
  {
    warning("Only 1 core specified for paralleization. Did you mean to specify more?")
  }
  
  if ( identical(qoi.in, "att") ||
       identical(qoi.in, "atc") ||
       identical(qoi.in, "art"))
  {
    
    coefs <- matrix(NA, 
                    nrow = number.iterations, 
                    ncol = length(placebo.lead))
    
    doParallel::registerDoParallel(num.cores)
    #initialize k
    k <- NA
    coefs <- foreach::foreach(k = 1:number.iterations, .combine = rbind) %dopar% {
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
    
    doParallel::registerDoParallel(num.cores)
    k <- NA    
    coefs <- foreach::foreach(k = 1:number.iterations, .combine = rbind) %dopar% {
      # make new data
      clusters <- unique(data.in[, unit.id.variable])
      units <- sample(clusters, size = length(clusters), replace = TRUE)
      
      while(all(!units %in% att.treated.unit.ids) || 
            all(!units %in% atc.treated.unit.ids)) 
        #while none of the units are treated units (att and atc), resample
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
  rownames(coefs) <- NULL
  colnames(coefs) <- NULL
  return(coefs)
  
}