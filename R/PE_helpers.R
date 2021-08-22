handle_moderating_variable <- function(ordered.data, att.sets, atc.sets, PM.object,
                                       moderator, unit.id, time.id)
{
  .reconstruct_pm_objects <- function(att.set = NULL, atc.set = NULL, PM.object_)
  {

    t.pm.object <- list()
    if(!is.null(att.set))
    {
      t.pm.object[["att"]] <- att.set
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
  w.it.df <- as.data.table(w.it.df[w.it.df$weight != 0,])
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
                               se.method)
{
  if (identical(se.method, "bootstrap"))
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
      
      coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
      ##### ** precompute some values for the bootstrap iterations
      perunitSum <- function(udf,
                             lead.in,
                             dependent.in) {
        w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_", qoi.in, x)), drop = FALSE]
        
        w.it.stars[is.na(w.it.stars)] <- 0
        return(colSums(apply(w.it.stars, MARGIN = 2, FUN = function(j) return(j * udf[, dependent.in]))))
      }
      
      per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                          FUN = perunitSum,
                          lead.in = lead,
                          dependent.in = outcome.variable)
      
      
      perunitSum_Dit <- function(udf) {
        d.it <- udf[, paste0("dits_",qoi.in)] #should always return a vector
        d.it[is.na(d.it)] <- 0
        return(sum(d.it))
      }
      
      per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                              FUN = perunitSum_Dit)
      
      units.id <- names(per.unit.sums)
      tdf <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
      colnames(tdf) <- sapply(lead, function(x) paste0("Wit_",qoi.in, x))
      tdf$unit.id <- NA
      tdf$unit.id <- as.character(unlist(units.id))
      # should be in the same order...
      tdf$Dit <- as.numeric(unlist(per.unit.dit.sums))
      #####***************
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
      
      if (identical(qoi.in, "att") || identical(qoi.in, "art"))
      {
        sets <- att.sets
      } else if (identical(qoi.in, "atc")) {
        sets <- atc.sets
      }
      ses <- apply(coefs, 2, sd, na.rm = T)
      names(ses) <- paste0("t+",lead)
      z <- list("estimates" = o.coefs,
                "bootstrapped.estimates" = coefs,
                "bootstrap.iterations" = number.iterations,
                "standard.error" = ses,
                "lag" = lag,
                "lead" = lead,
                "confidence.level" = confidence.level,
                "qoi" = qoi.in,
                "matched.sets" = sets,
                "se.method" = se.method)
      class(z) <- "PanelEstimate"
      return(z)
      
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
      
      o.coefs_ate <- (o.coefs_att*sum(data.in$dits_att) + o.coefs_atc*sum(data.in$dits_atc))/
        (sum(data.in$dits_att) + sum(data.in$dits_atc))
      
      if (length(lead[lead<0]) > 1)
      {
        names(o.coefs_ate)[(length(o.coefs_ate)-max(lead[lead>=0])):
                             length(o.coefs_ate)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs_ate)[(length(o.coefs_ate)-length(lead) + 1):
                             length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else
      {
        names(o.coefs_ate) <- sapply(lead, function(x) paste0("t+", x))
      }
      
      coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
      
      ##### ** precompute some values for the bootstrap iterations
      qoi.in <- "att"
      perunitSum <- function(udf,
                             lead.in,
                             dependent.in) {
        w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_", qoi.in, x)), drop = FALSE]
        
        w.it.stars[is.na(w.it.stars)] <- 0
        return(colSums(apply(w.it.stars, MARGIN = 2, FUN = function(j) return(j * udf[, dependent.in]))))
      }
      
      perunitSum_Dit <- function(udf) {
        d.it <- udf[, paste0("dits_",qoi.in)] #should always return a vector
        d.it[is.na(d.it)] <- 0
        return(sum(d.it))
      }
      
      per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                          FUN = perunitSum,
                          lead.in = lead,
                          dependent.in = outcome.variable)
      
      

      
      per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                              FUN = perunitSum_Dit)
      
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
                          dependent.in = outcome.variable)
      
      
      
      
      per.unit.dit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                              FUN = perunitSum_Dit)
      
      
      units.id <- names(per.unit.sums)
      tdf.atc <- as.data.frame(do.call(rbind, as.list(per.unit.sums)))
      colnames(tdf.atc) <- sapply(lead, function(x) paste0("Wit_","atc", x))
      tdf.atc$unit.id <- NA
      tdf.atc$unit.id <- as.character(unlist(units.id))
      # should be in the same order...
      tdf.atc$Dit <- as.numeric(unlist(per.unit.dit.sums))
      #####***************
      
      for (k in 1:number.iterations) {
        # make new data
        clusters <- unique(data.in[, unit.id.variable])
        units <- sample(clusters, size = length(clusters), replace=TRUE)
        while(all(!units %in% att.treated.unit.ids) || all(!units %in% atc.treated.unit.ids)) #while none of the units are treated units (att and atc), resample
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
      ses <- apply(coefs, 2, sd, na.rm = T)
      names(ses) <- paste0("t+",lead)
      
      z <- list("estimates" = o.coefs_ate,
                "bootstrapped.estimates" = coefs, 
                "bootstrap.iterations" = number.iterations, 
                "standard.error" = ses,
                "lead" = lead, "confidence.level" = confidence.level, 
                "qoi" = "ate", "matched.sets" = list(att = att.sets, atc = atc.sets),
                "se.method" = se.method)
      class(z) <- "PanelEstimate"
      return(z)
      
    } else {
      stop("invalid qoi")
    }
    
    
  } else if (identical(se.method, "conditional"))
  {
    if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
    o.coefs <- sapply(data.in[, sapply(lead, function(x) paste0("Wit_",qoi.in, x)), drop = FALSE],
                      equality_four,
                      y = data.in[c(outcome.variable)][,1],
                      z = data.in[, paste0("dits_", qoi.in)])
    if (qoi.in == "atc")
    {
      o.coefs <- o.coefs * -1
    }
    
    ## analytical v1
    perunitSum <- function(udf,
                           lead.in,
                           dependent.in) {
      w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_",qoi.in, x)), drop = FALSE]
      w.it.stars[is.na(w.it.stars)] <- 0
      return(colSums(apply(w.it.stars, MARGIN = 2, FUN = function(j) return(j * udf[, dependent.in]))))
    }
    
    per.unit.sums <- by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitSum,
                        lead.in = lead,
                        dependent.in = outcome.variable)
    
    
    tdf <- do.call(rbind, as.list(per.unit.sums))
    
    vdf <- apply(tdf, 2, var, na.rm = TRUE) #should return a number or vector
    D.it <- sum(data.in[, paste0("dits_", qoi.in)])
    D.it.denom <- D.it^2
    
    checkWits <- function(udf,
                          lead.in) {
      w.it.stars <- udf[, sapply(lead.in, function(x) paste0("Wit_", qoi.in, x)), drop = FALSE]
      w.it.stars[is.na(w.it.stars)] <- 0
      apply(w.it.stars, 2, FUN = function(x) all(x == 0))
    }
    
    check.vecs <- by(data.in, as.factor(data.in[, unit.id.variable]),
                     FUN = checkWits,
                     lead.in = lead)
    
    ndf <- do.call(rbind, as.list(check.vecs))
    N.nums <- apply(ndf, 2, function(x) sum(!x))
    
    #N.units <- length(unique(data[, unit.id]))
    #browser()
    
    estimator.var <- (N.nums * vdf) / D.it.denom
    names(estimator.var) <- paste0("t+",lead)
    names(o.coefs) <- paste0("t+", lead)
    
    if (is.null(atc.sets))
    {
      sets <- att.sets
    } else if (is.null(att.sets)) {
      sets <- atc.sets
    } else {
      stop("missing sets")
    }
    z <- list("estimates" = o.coefs,
              "standard.error" = sqrt(estimator.var),
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
    if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
    
    if (identical(qoi.in, "ate")) stop("analytical standard errors not available for ATE")
    o.coefs <- sapply(data.in[, sapply(lead, function(x) paste0("Wit_",qoi.in, x)), drop = FALSE],
                      equality_four,
                      y = data.in[c(outcome.variable)][,1],
                      z = data.in[, paste0("dits_", qoi.in)])
    if (qoi.in == "atc")
    {
      o.coefs <- o.coefs * -1
    }
    
    perunitSum <- function(udf,
                           lead.in,
                           dependent.in) {
      w.it.stars <- udf[, sapply(lead.in,
                                 function(x) paste0("Wit_", qoi.in ,x)),
                        drop = FALSE]
      w.it.stars[is.na(w.it.stars)] <- 0
      return(colSums(apply(w.it.stars, MARGIN = 2,
                           FUN = function(j) return(j * udf[, dependent.in]))))
    }

    Ais <- by(data.in, as.factor(data.in[, unit.id.variable]),
              FUN = perunitSum,
              lead.in = lead,
              dependent.in = outcome.variable)

    tdf <- do.call(rbind, as.list(Ais))
    As <- colSums(tdf, na.rm = TRUE)

    ###

    perunitDits <- function(udf) {
      dits <- udf[, paste0("dits_", qoi.in), drop = FALSE]
      dits[is.na(dits)] <- 0
      return(sum(dits, na.rm = TRUE))
    }

    Bi <- as.numeric(by(data.in, as.factor(data.in[, unit.id.variable]),
                        FUN = perunitDits))

    N <- length(unique(data.in[, unit.id.variable]))

    EB <- mean(Bi) * N
    VB <- var(Bi, na.rm = TRUE) * N
    vdf <- apply(tdf, 2, var, na.rm = TRUE) #should return a number or vector
    VA <- N * vdf
    EA <- N * colMeans(tdf, na.rm = TRUE)
    covAB <- apply(tdf, 2, FUN = function(x) return(cov(x, Bi)))

    estimator.var <- (1 / (EB^2)) * (VA - (2 * (EA / EB) * covAB) + ( (EA^2 / EB^2) * VB) )

    names(estimator.var) <- paste0("t+",lead)
    names(o.coefs) <- paste0("t+", lead)
    
    if (is.null(atc.sets))
    {
      sets <- att.sets
    } else if (is.null(att.sets)) {
      sets <- atc.sets
    } else {
      stop("missing sets")
    }
    z <- list("estimates" = o.coefs,
              "standard.error" = sqrt(estimator.var),
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
          
          while(all(!units %in% att.treated.unit.ids) || all(!units %in% atc.treated.unit.ids)) #while none of the units are treated units (att and atc), resample
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
    } else if (identical(se.method, "conditional") || identical(se.method, "unconditional")) 
    {
      stop("not currently implemented")
    } else {
    stop("placebo tests only currently available with bootstrap!")
  }
  
  
}
