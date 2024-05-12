################################################################################################
# Functions for calculating standard errors and associated helper functions########################
################################################################################################
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


perunitSum_Dit <- function(udf, qoi_in) {
  d.it <- udf[, paste0("dits_",qoi_in)] #should always return a vector
  d.it[is.na(d.it)] <- 0
  return(sum(d.it))
}


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
                             se.method,
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
                            se.method,
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

handle_unconditional_se <- function(qoi.in, data.in, lead,
                                 outcome.variable,
                                 unit.id.variable) 
{
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

