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
getWits <- function(matched_sets, lead, 
                    data, continuous.treatment = FALSE)
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
  if ( identical(qoi.in, "att") || identical(qoi.in, "art")) 
  {
    if (!identical(qoi.in, "ate")) qoi.t <- qoi.in
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
  if (identical(qoi.in, "atc")) 
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


preparePlaceboTestData <- function(data.in, lead, sets.att = NULL,
                                   sets.atc = NULL, continuous.treatment,
                                   qoi.in, dependent.variable,
                                   time.var,
                                   id.var,
                                   placebo.test = TRUE,
                                   placebo.lag = NULL)
{
  
  # lag DV by one
  browser()
  data.in$lag.dv <- c(NA, data.in[, dependent.variable][-nrow(data.in)])
  data.in$lag.dv[which(!duplicated(data.in[,id.var]))] <- NA
  
  data.in[, dependent.variable] <- data.in$lag.dv
  data.in <- data.in[, -ncol(data.in)] #remove lag.dv column
  
  if ( identical(qoi.in, "att") || identical(qoi.in, "art")) 
  {
    if (!identical(qoi.in, "ate")) qoi.t <- qoi.in
    
    for (lag.val in 2:placebo.lag) 
    {
      j <- lag.val - 2 #set lead
      dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.att, 
                            continuous.treatment = continuous.treatment)
      data.in = merge(x = data.in, y = dense.wits, all.x = TRUE, 
                      by.x = colnames(data.in)[1:2], by.y = c("id", "t"))
      colnames(data.in)[length(data.in)] <- paste0("Wit_", qoi.t, j)
      data.in[is.na(data.in[, length(data.in)]), length(data.in)] <- 0 #replace NAs with zeroes
      
      
      data.in$lag.dv <- c(NA, data.in[, dependent.variable][-nrow(data.in)])
      data.in$lag.dv[which(!duplicated(data.in[,id.var]))] <- NA
      
      data.in[, dependent.variable] <- data.in$lag.dv
      data.in <- data.in[, -ncol(data.in)] #remove lag.dv column
      
    }
    
    data.in[, paste0("dit_", qoi.t)] <- getDits(matched_sets = sets.att, 
                                                data = data.in)
    colnames(data.in)[length(data.in)] <- paste0("dits_", qoi.t)
    data.in[, paste0("Wit_", qoi.t, "-1")] <- 0
    
  } 
  if (identical(qoi.in, "atc")) 
  {
    stop("atc not supported for placebo tests")
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
                               lag)
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
    
    if (all(lead >= 0)) #not placebo/regular case
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
      
      y.in <- d.sub1[,outcome.variable]
      z.in <- d.sub1[, paste0("dits_", qoi.in)]
      
      at__new <-  sapply(d.sub1[, col.idx, 
                                drop = FALSE],
                         equality_four,
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
    
    z <- list("estimates" = o.coefs,
              "bootstrapped.estimates" = coefs, 
              "bootstrap.iterations" = number.iterations, 
              "standard.error" = apply(coefs, 2, sd, na.rm = T),
              "lag" = lag,
              "lead" = lead, 
              "confidence.level" = confidence.level, 
              "qoi" = qoi.in, 
              "matched.sets" = sets)
    class(z) <- "PanelEstimate"
    return(z)
    
  } else {
    stop("invalid qoi")
   
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
                                      placebo.lead)
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
      d[, paste0(dv, "l", idx)] <- lapply(k, function(x) data.table::shift(d[, dv], n = k, type = "lag"))
      return(d)
    }
    y.in <- by(data.in, as.factor(data.in[, unit.id.variable]), 
               FUN = create.lagged.dfs, 
               dv = outcome.variable,
               idx = placebo.lead - 1,
               k = placebo.lead - 1, 
               simplify = FALSE)
    browser()
    data.in <- do.call(rbind, y.in)
    y.in <- data.in[, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
    z.in <- data.in[, paste0("dits_", qoi.in)]
    
    o.coefs <- equality_four_placebo(x.in, y.in, z.in)
    
    #do coefficient flip for atc
    if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
    
    
    coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
    
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
      
      y.in <- d.sub1[,, paste0(outcome.variable, "l", (placebo.lead - 1)), drop = FALSE]
      z.in <- d.sub1[, paste0("dits_", qoi.in)]
      
      # at__new <- sapply(d.sub1[, col.idx, 
      #                           drop = FALSE],
      #                    equality_four,
      #                    y = y.in,
      #                    z = z.in)
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
    
    z <- list("estimates" = o.coefs,
              "bootstrapped.estimates" = coefs, 
              "bootstrap.iterations" = number.iterations, 
              "standard.error" = apply(coefs, 2, sd, na.rm = T),
              "lag" = lag,
              "lead" = lead, 
              "confidence.level" = confidence.level, 
              "qoi" = qoi.in, 
              "matched.sets" = sets)
    class(z) <- "PanelEstimate"
    return(z)
    
  } else {
    stop("invalid qoi")
    
  }
}

