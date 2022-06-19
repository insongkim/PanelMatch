calculate_point_estimates <- function(qoi.in, data.in, lead,
                                      outcome.variable,
                                      pooled = FALSE)
{
  if ( identical(qoi.in, "att") ||
       identical(qoi.in, "atc") ||
       identical(qoi.in, "art"))
  {
    
    col.idx <- sapply(lead, function(x) paste0("Wit_", qoi.in, x))
    x.in <- data.in[, col.idx, drop = FALSE]
    y.in <- data.in[c(outcome.variable)][,1]
    z.in <- data.in[, paste0("dits_", qoi.in)]
    
    o.coefs <- sapply(x.in, equality_four,
                      y = y.in,
                      z = z.in)
    
    #do coefficient flip for atc
    if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
    
    
    
    if (all(lead >= 0)) 
    {
      if (length(lead[lead < 0]) > 1)
      {
        names(o.coefs)[(length(o.coefs) - max(lead[lead >= 0])):length(o.coefs)] <-
          sapply(lead[lead >= 0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs) - length(lead) + 1):length(lead[lead < 0])] <-
          sapply(lead[lead < 0], function(x) paste0("t", x))
        
      } else
      {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
    }
    
    if (pooled)
    {
      o.coefs <- mean(o.coefs, na.rm = TRUE)
      names(o.coefs) <- NULL
    }
  } else if (identical(qoi.in, "ate")) {
    o.coefs_att <-  sapply(data.in[, sapply(lead, function(x) paste0("Wit_att", x)),
                                   drop = FALSE],
                           equality_four,
                           y = data.in[c(outcome.variable)][,1],
                           z = data.in$dits_att)
    
    o.coefs_atc <-  -sapply(data.in[, sapply(lead, function(x) paste0("Wit_atc", x)),
                                    drop = FALSE],
                            equality_four,
                            y = data.in[c(outcome.variable)][,1],
                            z = data.in$dits_atc)
    
    o.coefs <- (o.coefs_att*sum(data.in$dits_att) + o.coefs_atc*sum(data.in$dits_atc))/
      (sum(data.in$dits_att) + sum(data.in$dits_atc))
    
    
    
    if (length(lead[lead<0]) > 1)
    {
      names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                           length(o.coefs)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
      names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                           length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
      
    } else
    {
      names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
    }
    
    if (pooled)
    {
      o.coefs <- mean(o.coefs, na.rm = TRUE)
      names(o.coefs) <- NULL
    }
  }
  return(o.coefs)
}

perunitSum <- function(udf,
                       lead.in,
                       dependent.in,
                       qoi_in) {
  w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_", qoi_in, x)), drop = FALSE]
  
  w.it.stars[is.na(w.it.stars)] <- 0
  return(colSums(apply(w.it.stars, MARGIN = 2, FUN = function(j) return(j * udf[, dependent.in]))))
}


perunitSum_Dit <- function(udf, qoi_in) {
  d.it <- udf[, paste0("dits_",qoi_in)] #should always return a vector
  d.it[is.na(d.it)] <- 0
  return(sum(d.it))
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
  ##### ** precompute some values for the bootstrap iterations
  
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
      }
      while (all(!units %in% treated.unit.ids)) #while none of the units are treated units, resample
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
      coefs[k,] <- multiplier * colSums(bdf[,
                                            sapply(lead, function(x) paste0("Wit_", qoi.in, x)),
                                            drop = FALSE] * bdf[, 'frequency'], na.rm = FALSE) / (sum(bdf[, 'frequency'] * bdf[, 'Dit']))
      
      
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
            all(!units %in% atc.treated.unit.ids)) #while none of the units are treated units (att and atc), resample
      {
        units <- sample(clusters, size = length(clusters), replace=TRUE)
      }
      
      fdf <- as.data.frame(table(units))
      colnames(fdf) <- c('unit.id', 'frequency')
      bdf <- merge(x = tdf, y = fdf, all.x = TRUE)
      bdf$frequency[is.na(bdf$frequency)] <- 0
      
      att_new <- colSums(bdf[, sapply(lead, function(x) paste0("Wit_att", x)),
                             drop = FALSE] * bdf[, 'frequency'], na.rm = FALSE) / (sum(bdf[, 'frequency'] * bdf[, 'Dit']))
      
      
      bdf.atc <- merge(x = tdf.atc, y = fdf, all.x = TRUE)
      bdf.atc$frequency[is.na(bdf.atc$frequency)] <- 0
      
      atc_new <- -colSums(bdf.atc[, sapply(lead, function(x) paste0("Wit_atc", x)),
                                  drop = FALSE] * bdf.atc[, 'frequency'], na.rm = FALSE) / (sum(bdf.atc[, 'frequency'] * bdf.atc[, 'Dit']))
      
      coefs[k,] <- (att_new*sum(bdf$Dit * bdf$frequency) + atc_new*sum(bdf.atc$Dit * bdf.atc$frequency)) /
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
    w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_", qoi_in, x)), drop = FALSE]
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
  covAB <- apply(tdf, 2, FUN = function(x) return(cov(x, Bi)))
  
  estimator.var <- (1 / (EB^2)) * (VA - (2 * (EA / EB) * covAB) + ( (EA^2 / EB^2) * VB) )
  
  names(estimator.var) <- paste0("t+",lead)
  return(sqrt(estimator.var))
}

handle_moderating_variable <- function(ordered.data, att.sets, atc.sets, PM.object,
                                       moderator, unit.id, time.id, qoi.in)
{
  .reconstruct_pm_objects <- function(att.set = NULL, 
                                      atc.set = NULL, 
                                      PM.object_)
  {
    #TODO: maybe add somthing about the qoi in? 
    # needs to handle both att and art here..
    t.pm.object <- list()
    if(!is.null(att.set))
    {
      #t.pm.object[["att"]] <- att.set
      t.pm.object[[qoi.in]] <- att.set
    }
    if(!is.null(atc.set))
    {
      t.pm.object[["atc"]] <- atc.set
    }
    attrib <- names(attributes(PM.object_))[names(attributes(PM.object_)) != "names"]
    for(tatt in attrib)
    {
      attr(t.pm.object, tatt) <- attr(PM.object_, tatt)
    }
    return(t.pm.object)
  }
  
  ref.names <- paste0(ordered.data[, unit.id], ".", ordered.data[, time.id])
  moderator.vector <- ordered.data[, moderator]
  names(moderator.vector) <- ref.names
  moderated.sets.att <- list()
  moderated.sets.atc <- list()
  subset.list <- list()
  moderating.values <- unique(ordered.data[, moderator])
  for(val in as.vector(na.omit(moderating.values)))
  {
    #make sure we handle empty set situation
    if(!is.null(att.sets))
    {
      indx.set <- moderator.vector[names(att.sets)] == val
      t.set <- att.sets[indx.set]
      if(length(t.set) > 0)
      {
        moderated.sets.att[[make.names(val)]] <- t.set
      }
      else
      {
        moderated.sets.att[[make.names(val)]] <- NULL
      }
    }
    if(!is.null(atc.sets))
    {
      indx.set <- moderator.vector[names(atc.sets)] == val
      t.set <- atc.sets[indx.set]
      if(length(t.set) > 0)
      {
        moderated.sets.atc[[make.names(val)]] <- t.set
      }
      else
      {
        moderated.sets.atc[[make.names(val)]] <- NULL
      }
      
    }
  }
  if(!is.null(att.sets) & !is.null(atc.sets))
  {
    if(!identical(names(moderated.sets.att), names(moderated.sets.atc)))
    {
      stop("Data insufficient for calculating ATE with moderating variables")
    }
  }
  
  if(is.null(atc.sets))
  {
    #do att
    ret.obj <-lapply(moderated.sets.att,
                     FUN = .reconstruct_pm_objects,
                     PM.object_ = PM.object, atc.set = NULL)
  }
  else if(is.null(att.sets))
  {
    #do atc
    ret.obj <-lapply(moderated.sets.atc,
                     FUN = .reconstruct_pm_objects,
                     PM.object_ = PM.object, att.set = NULL)
  }
  else
  {
    
    ret.obj <- mapply(FUN = .reconstruct_pm_objects, SIMPLIFY = FALSE,
                      att.set = moderated.sets.att, atc.set = moderated.sets.atc,
                      MoreArgs = list(PM.object_ = PM.object))
  }
  names(ret.obj) <- as.character(as.vector(na.omit(moderating.values)))
  return(ret.obj)
}




#pcs and pts create data frames with the time/id combinations--that need to be found so that they can be easily looked up in the data frame via a hash table. The data frame also contains information
# about the weight of that unit at particular times, so we use the hashtable to look up where to put this data so that we can easily assign the appropriate weights in the original data frame containing the problem data.
# pcs does this for all control units in a matched set
# ("Prepare Control unitS)
pcs <- function(sets, lead.in,
                continuous.treatment = FALSE)
{
  
  L <- attr(sets, "lag")
  ts <- as.numeric(sub(".*\\.", "", names(sets)))
  
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  lsets <- sapply(sets, length)
  ts <- unlist(mapply(FUN = make.years,
                      t = ts, lead = lead.in,
                      repnum = lsets, SIMPLIFY = FALSE))
  ids <- rep(unlist(sets), rep((2), length(unlist(sets)) ) )
  names(ids) <- NULL
  if (continuous.treatment)
  {
    wts <- unlist(sapply(sets,
                         function(s){return(attr(s, "weights") / attr(s, 'treatment.change'))}))
  } else {
    wts <- unlist(sapply(sets,
                         function(s){return(attr(s, "weights"))}))
  }
  
  names(wts) <- NULL
  
  wts <- rep(wts, rep((2), length(unlist(sets)))) * c(1,-1)
  
  
  num.empty <- sum(!sapply(sets, length) > 0)
  set.nums <- rep(0:(length(sets) - num.empty - 1), (lsets[lsets != 0] * (2)))
  # if (lead.in < 0)
  # {
  #   wts <- wts * -1
  # }
  dtf <- data.frame(t = ts, id = ids, weight = wts, set.number = set.nums)
  return(dtf)
}
# refer to the description above -- pts works on treated units ("Prepare Treated unitS)
pts <- function(sets, lead.in,
                continuous.treatment = FALSE)
{
  include <- sapply(sets, length) > 0
  num.empty <- sum(!include)
  
  ts <- as.numeric(sub(".*\\.", "", names(sets)))[include]
  tids <- as.numeric(sub("\\..*", "", names(sets)))[include]
  tids <- rep(tids, rep(2, length(tids)))
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  
  ts <- unlist(mapply(FUN = make.years, t = ts,
                      lead = lead.in, repnum = 1, SIMPLIFY = FALSE))
  
  
  wts <- rep(c(-1, 1), length(sets) - num.empty)
  
  set.nums <- rep(0:(length(sets) - num.empty - 1 ), rep(2, length(sets) - num.empty ))
  # if (lead.in < 0)
  # {
  #   wts <- wts * -1
  # }
  ldf = data.frame(t = ts, id = tids, weight = wts, set.number = set.nums)
  if (continuous.treatment)
  {
    existing.sets <- sets[sapply(sets, length) > 0]
    treatment.changes <- unlist(lapply(existing.sets,
                                       function(x) return(attr(x, "treatment.change"))))
    x.std <- unlist(sapply(treatment.changes,
                           function(x) rep(x, 2), simplify = FALSE))
    names(x.std) <- NULL
    ldf$weight <- ldf$weight / x.std
  }
  return(ldf)
  
  
}
#returns a vector of Wits, as defined in the paper (equation 25 or equation 23). They should be in the same order as the data frame containing the original problem data. The pts, pcs, and getWits functions act for a specific
# lead. So, for instance if our lead window is 0,1,2,3,4, these function must be called for each of those -- so for 0, then for 1, etc.
getWits <- function(matched_sets, lead, data, 
                    estimation.method = "bootstrap",
                    continuous.treatment = FALSE)
{
  
  #sort the data
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  #include <- sapply(matched_sets, length) > 0
  #num.empty <- sum(!include)
  
  
  #prep control sets, prep treatment sets for search/summation vector
  p.df <- pcs(matched_sets, lead,
              continuous.treatment = continuous.treatment)
  t.df <- pts(matched_sets, lead,
              continuous.treatment = continuous.treatment)
  
  ##to solve cran check note about unbound variables
  . <- weight <- id <- NULL
  w.it.df <- rbind(p.df, t.df)
  w.it.df <- data.table::as.data.table(w.it.df[w.it.df$weight != 0,])
  summarized.Wits <- w.it.df[,.(Wit = sum(weight)), by = .(t,id)]
  return(summarized.Wits)
  
}


# returns a vector of dit values, as defined in the paper. They should be in the same order as the data frame containing the original problem data.
getDits <- function(matched_sets, data)
{
  
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  include <- sapply(matched_sets, length) > 0
  msets <- matched_sets[include]
  nms <- names(msets)
  refnames <- paste0(data[, id.var], ".", data[, t.var])
  dit.vect <- get_dits(refnames, nms)
  
  return(dit.vect)
}

#function that ultimately calculates the point estimate values
equality_four <- function(x, y, z){
  
  return(sum(x*y)/sum(z))
}

equality_four_placebo <- function(x, y, z){
  
  y[is.na(y)] <- 0
  res <- colSums(x * y) / sum(z)
  return(res)
}



