#' perunitSum
#' This is a low level function that is used to calculate a value associated with each unit. This value is a weighted summation of the dependent variable, based on the Wit values discussed in Imai et al. (2023)
#' @param udf data.frame
#' @param lead.in integer. A particular lead value
#' @param dependent.in string specifying the dependent variable name
#' @param qoi_in string specifying the QOI
#'
#' @return Named vector containing the per-unit sums. 
#' @keywords internal
perunitSum <- function(udf,
                       lead.in,
                       dependent.in,
                       qoi_in) {
  w.it.stars <- udf[, sapply(lead.in, 
                             function(x) paste0("Wit_", qoi_in, x)), 
                    drop = FALSE]
  
  w.it.stars[is.na(w.it.stars)] <- 0
  return(colSums(apply(w.it.stars, 
                       MARGIN = 2, 
                       FUN = function(j) return(j * udf[, dependent.in]))))
}


#' perunitSum_Dit
#' Similar to perunitSum, this is a low level helper function for calculating specific values defined in Imai et al. (2023). This focuses on Dit rather than Wit
#' @param udf data.frame
#' @param qoi_in string specifying the QOI
#'
#' @return Named vector containing the per-unit sums.
#' @keywords internal
perunitSum_Dit <- function(udf, qoi_in) {
  
  d.it <- udf[, paste0("dits_",qoi_in)] #should always return a vector
  d.it[is.na(d.it)] <- 0
  return(sum(d.it))
}

#' handle_bootstrap_parallel
#'
#' Helper function for calculating bootstrapped estimates for the QOI. This version is parallelized.
#' @param qoi.in String specifying qoi
#' @param data.in data.frame object with the data
#' @param number.iterations integer. Specifies number of bootstrap iterations
#' @param att.treated.unit.ids Integer vector specifying the treated units for the att or art
#' @param atc.treated.unit.ids Integer vector specifying the "treated" units under the atc definition
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable
#' @param confidence.level double. specifies confidence level for confidence interval
#' @param lag integer vector specifying size of the lag.
#' @param pooled logical. Specifies whether or not to calculate point estimates for each specified lead value, or a single pooled estimate.
#' @param num.cores number of cores to be used for parallelization
#' @return Returns a matrix of bootstrapped QOI estimate values.
#'
#' @keywords internal
handle_bootstrap_parallel <- function(qoi.in, 
                             data.in, 
                             lead,
                             number.iterations,
                             att.treated.unit.ids,
                             atc.treated.unit.ids,
                             outcome.variable,
                             unit.id.variable,
                             confidence.level,
                             lag,
                             pooled,
                             num.cores = 1) 
{
  coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
  
  
  if (qoi.in %in% c("att", "art", "atc"))
  {
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable,
                        qoi_in = qoi.in)
    
    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    units.id <- names(per.unit.sums)
    tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf) <- sapply(lead, function(x) paste0("Wit_",qoi.in, x))
    tdf$unit.id <- NA
    tdf$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    ######### Do the bootstrapping
    if (identical(num.cores, 1))
    {
      warning("Only 1 core specified for paralleization. Did you mean to specify more?")
    }
    doParallel::registerDoParallel(num.cores)
    k <- NA
    coefs <- foreach::foreach(k = 1:number.iterations, .combine = rbind) %dopar% {
      # make new data
      clusters <- unique(data.in[, unit.id.variable])
      units <- sample(clusters, size = length(clusters), replace = TRUE)
      if (identical(qoi.in, "att") || identical(qoi.in, "art")) {
        treated.unit.ids <- att.treated.unit.ids
      } else {
        treated.unit.ids <- atc.treated.unit.ids
      }
      # while none of the units are treated units, resample:
      while (all(!units %in% treated.unit.ids)) {
        units <- sample(clusters, size = length(clusters), replace = TRUE)
      }
      
      fdf <- as.data.frame(table(units))
      colnames(fdf) <- c('unit.id', 'frequency')
      bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
      bdf$frequency[is.na(bdf$frequency)] <- 0
      if (identical(qoi.in, "atc")) {
        multiplier <- -1
      } else {
        multiplier <- 1
      }
      multiplier * colSums(bdf[, 
                               sapply(lead, 
                                      function(x) paste0("Wit_", 
                                                         qoi.in, x)), 
                               drop = FALSE] * bdf[, 'frequency'], 
                           na.rm = FALSE) / (sum(bdf[, 'frequency'] * bdf[, 'Dit']))
    }
    
    rownames(coefs) <- NULL
    colnames(coefs) <- NULL
    
    ########
    if (pooled)
    {
      coefs <- as.matrix(rowMeans(coefs), ncol = 1)
    }
    return(coefs)
  } else if (identical(qoi.in, "ate"))
  {
    qoi.in <- "att"
    
    
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable, 
                        qoi_in = qoi.in)
    
    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    units.id <- names(per.unit.sums)
    tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf) <- sapply(lead, function(x) paste0("Wit_","att", x))
    tdf$unit.id <- NA
    tdf$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    qoi.in <- "atc"
    
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable,
                        qoi_in = qoi.in)
    
    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    
    units.id <- names(per.unit.sums)
    tdf.atc <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf.atc) <- sapply(lead, function(x) paste0("Wit_","atc", x))
    tdf.atc$unit.id <- NA
    tdf.atc$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf.atc$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    
    doParallel::registerDoParallel(num.cores)
    k <- NA
    coefs <- foreach::foreach(k = 1:number.iterations, .combine = rbind) %dopar% {
      # make new data
      clusters <- unique(data.in[, unit.id.variable])
      units <- sample(clusters, size = length(clusters), replace=TRUE)
      while(all(!units %in% att.treated.unit.ids) || 
            all(!units %in% atc.treated.unit.ids)) 
      { #while none of the units are treated units (att and atc), resample
        units <- sample(clusters, size = length(clusters), replace=TRUE)
      }
      
      fdf <- as.data.frame(table(units))
      colnames(fdf) <- c('unit.id', 'frequency')
      bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
      bdf$frequency[is.na(bdf$frequency)] <- 0
      
      att_new <- colSums(bdf[, sapply(lead, function(x) paste0("Wit_att", x)),
                             drop = FALSE] * bdf[, 'frequency'],
                         na.rm = FALSE) / (sum(bdf[, 
                                                   'frequency'] * bdf[, 'Dit']))
      
      
      bdf.atc <- merge(x = tdf.atc, y = fdf, all.x = TRUE)
      bdf.atc$frequency[is.na(bdf.atc$frequency)] <- 0
      
      atc_new <- -colSums(bdf.atc[, sapply(lead, 
                                           function(x) paste0("Wit_atc", x)),
                                  drop = FALSE] * bdf.atc[, 'frequency'], 
                          na.rm = FALSE) / (sum(bdf.atc[, 'frequency'] * 
                                                  bdf.atc[, 'Dit']))
      
      coefs[k,] <- (att_new*sum(bdf$Dit * bdf$frequency) +
                      atc_new*sum(bdf.atc$Dit * bdf.atc$frequency)) /
        (sum(bdf$Dit * bdf$frequency) + sum(bdf.atc$Dit * bdf.atc$frequency))
      
      
      
    }
    rownames(coefs) <- NULL
    colnames(coefs) <- NULL
    
    if (pooled)
    {
      coefs <- as.matrix(rowMeans(coefs), ncol = 1)
    }
    
    return(coefs)
    
    
  } else {
    stop("qoi misspecified")
  }
  
}


#' handle_bootstrap
#'
#' Helper function for calculating bootstrapped estimates for the QOI. This version is not parallelized.
#' @param qoi.in String specifying qoi
#' @param data.in data.frame object with the data
#' @param number.iterations integer. Specifies number of bootstrap iterations
#' @param att.treated.unit.ids Integer vector specifying the treated units for the att or art
#' @param atc.treated.unit.ids Integer vector specifying the "treated" units under the atc definition
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable
#' @param confidence.level double. specifies confidence level for confidence interval
#' @param lag integer vector specifying size of the lag.
#' @param pooled logical. Specifies whether or not to calculate point estimates for each specified lead value, or a single pooled estimate.
#' @return Returns a matrix of bootstrapped QOI estimate values.
#'
#' @keywords internal
handle_bootstrap <- function(qoi.in, 
                            data.in, 
                            lead,
                            number.iterations,
                            att.treated.unit.ids,
                            atc.treated.unit.ids,
                            outcome.variable,
                            unit.id.variable,
                            confidence.level,
                            lag,
                            pooled) 
{
  coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
  
  
  if (qoi.in %in% c("att", "art", "atc"))
  {
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable,
                        qoi_in = qoi.in)
    
    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    units.id <- names(per.unit.sums)
    tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf) <- sapply(lead, function(x) paste0("Wit_",qoi.in, x))
    tdf$unit.id <- NA
    tdf$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    ######### Do the bootstrapping
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
      } #while none of the units are treated units, resample:
      while (all(!units %in% treated.unit.ids)) 
      {
        units <- sample(clusters, size = length(clusters), replace = TRUE)
      }
      
      fdf <- as.data.frame(table(units))
      colnames(fdf) <- c('unit.id', 'frequency')
      bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
      bdf$frequency[is.na(bdf$frequency)] <- 0
      if (identical(qoi.in, "atc"))
      {
        multiplier <- -1
      } else {
        multiplier <- 1
      }
      coefs[k,] <- multiplier * 
        colSums(bdf[, 
                    sapply(lead, 
                           function(x) paste0("Wit_", qoi.in, x)),
                                            drop = FALSE] * bdf[, 'frequency'], 
                na.rm = FALSE) / (sum(bdf[, 'frequency'] * bdf[, 'Dit']))
      
      
    }
    
    if (pooled)
    {
      coefs <- as.matrix(rowMeans(coefs), ncol = 1)
    }
    return(coefs)
  } else if (identical(qoi.in, "ate"))
  {
    qoi.in <- "att"
    
    
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable, 
                        qoi_in = qoi.in)
    
    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    units.id <- names(per.unit.sums)
    tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf) <- sapply(lead, function(x) paste0("Wit_","att", x))
    tdf$unit.id <- NA
    tdf$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    qoi.in <- "atc"
    
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable,
                        qoi_in = qoi.in)

    per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                            FUN = perunitSum_Dit,
                            qoi_in = qoi.in)
    
    
    units.id <- names(per.unit.sums)
    tdf.atc <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
    colnames(tdf.atc) <- sapply(lead, function(x) paste0("Wit_","atc", x))
    tdf.atc$unit.id <- NA
    tdf.atc$unit.id <- as.character(unlist(units.id))
    # should be in the same order...
    tdf.atc$Dit <- as.numeric(unlist(per.unit.dit.sums))
    
    
    for (k in 1:number.iterations) {
      # make new data
      clusters <- unique(data.in[, unit.id.variable])
      units <- sample(clusters, size = length(clusters), replace=TRUE)
      while(all(!units %in% att.treated.unit.ids) || 
            all(!units %in% atc.treated.unit.ids)) 
      { #while none of the units are treated units (att and atc), resample
        units <- sample(clusters, size = length(clusters), replace=TRUE)
      }
      
      fdf <- as.data.frame(table(units))
      colnames(fdf) <- c('unit.id', 'frequency')
      bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
      bdf$frequency[is.na(bdf$frequency)] <- 0
      
      att_new <- colSums(bdf[, sapply(lead, function(x) paste0("Wit_att", x)),
                             drop = FALSE] * bdf[, 'frequency'],
                         na.rm = FALSE) / (sum(bdf[, 
                                                   'frequency'] * bdf[, 'Dit']))
      
      
      bdf.atc <- merge(x = tdf.atc, y = fdf, all.x = TRUE)
      bdf.atc$frequency[is.na(bdf.atc$frequency)] <- 0
      
      atc_new <- -colSums(bdf.atc[, sapply(lead, 
                                           function(x) paste0("Wit_atc", x)),
                                  drop = FALSE] * bdf.atc[, 'frequency'], 
                          na.rm = FALSE) / (sum(bdf.atc[, 'frequency'] * 
                                                  bdf.atc[, 'Dit']))
      
      coefs[k,] <- (att_new*sum(bdf$Dit * bdf$frequency) +
                      atc_new*sum(bdf.atc$Dit * bdf.atc$frequency)) /
        (sum(bdf$Dit * bdf$frequency) + sum(bdf.atc$Dit * bdf.atc$frequency))
    }
    
    if (pooled)
    {
      coefs <- as.matrix(rowMeans(coefs), ncol = 1)
    }
    
    return(coefs)
    
    
  } else {
    stop("qoi misspecified")
  }
  
}

#' handle_conditional_se
#' Calculates conditional standard errors analytically, as defined in Imai et al. (2023). See PanelEstimate() for a more complete description of the standard error types. 
#' @param qoi.in string specifying the QOI
#' @param data.in data.frame specifying the data
#' @param lead See PanelMatch() documentation
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable.
#'
#' @return Named vector with standard error estimates
#' @keywords internal
handle_conditional_se <- function(qoi.in, data.in, lead,
                                  outcome.variable,
                                  unit.id.variable)
{
  
  if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
  per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                      FUN = perunitSum,
                      lead.in = lead,
                      dependent.in = outcome.variable,
                      qoi_in = qoi.in)
  tdf <- do.call(rbind, as.list(per.unit.sums))
  
  vdf <- apply(tdf, 2, var, na.rm = TRUE) #should return a number or vector
  D.it <- sum(data.in[, paste0("dits_", qoi.in)])
  D.it.denom <- D.it^2
  
  checkWits <- function(udf,
                        lead.in,
                        qoi_in) {
    w.it.stars <- udf[, sapply(lead.in, 
                               function(x) paste0("Wit_", qoi_in, x)), 
                      drop = FALSE]
    w.it.stars[is.na(w.it.stars)] <- 0
    apply(w.it.stars, 2, FUN = function(x) all(x == 0))
  }
  
  check.vecs <- by(data.in, as.factor(data.in[, unit.id.variable]),
                   FUN = checkWits,
                   lead.in = lead,
                   qoi_in = qoi.in)
  
  ndf <- do.call(rbind, as.list(check.vecs))
  N.nums <- apply(ndf, 2, function(x) sum(!x))
  estimator.var <- (N.nums * vdf) / D.it.denom
  names(estimator.var) <- paste0("t+",lead)
  return(sqrt(estimator.var))
  
}

#' handle_conditional_se
#' Calculates conditional standard errors analytically, as defined in Imai et al. (2023). See PanelEstimate() for a more complete description of the standard error types. 
#' @param qoi.in string specifying the QOI
#' @param data.in data.frame specifying the data
#' @param lead See PanelMatch() documentation
#' @param outcome.variable string specifying the name of the outcome variable
#' @param unit.id.variable string specifying the name of the unit id variable.
#'
#' @return Named vector with standard error estimates
#' @keywords internal
handle_unconditional_se <- function(qoi.in, data.in, lead,
                                 outcome.variable,
                                 unit.id.variable) 
{
  if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
  Ais <- by(data.in, as.factor(data.in[, unit.id.variable]),
            FUN = perunitSum,
            lead.in = lead,
            dependent.in = outcome.variable,
            qoi_in = qoi.in)
  
  tdf <- do.call(rbind, as.list(Ais))
  As <- colSums(tdf, na.rm = TRUE)
  
  
  perunitDits <- function(udf, qoi_in) {
    dits <- udf[, paste0("dits_", qoi_in), drop = FALSE]
    dits[is.na(dits)] <- 0
    return(sum(dits, na.rm = TRUE))
  }
  
  Bi <- as.numeric(by(data.in, as.factor(data.in[, unit.id.variable]),
                      FUN = perunitDits,
                      qoi_in = qoi.in))
  
  N <- length(unique(data.in[, unit.id.variable]))
  
  EB <- mean(Bi) * N
  VB <- var(Bi, na.rm = TRUE) * N
  vdf <- apply(tdf, 2, var, na.rm = TRUE) #should return a number or vector
  VA <- N * vdf
  EA <- N * colMeans(tdf, na.rm = TRUE)
  covAB <- apply(tdf, 2, 
                 FUN = function(x) return(cov(x, Bi)))
  
  estimator.var <- (1 / (EB^2)) * 
    (VA - (2 * (EA / EB) * covAB) + ( (EA^2 / EB^2) * VB) )
  
  names(estimator.var) <- paste0("t+",lead)
  return(sqrt(estimator.var))
}