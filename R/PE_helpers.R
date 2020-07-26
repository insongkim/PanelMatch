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
    ret.obj <-lapply(moderated.sets.att, FUN = .reconstruct_pm_objects, PM.object_ = PM.object, atc.set = NULL)
  }
  else if(is.null(att.sets))
  {
    #do atc
    ret.obj <-lapply(moderated.sets.atc, FUN = .reconstruct_pm_objects, PM.object_ = PM.object, att.set = NULL)
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
pcs <- function(sets, lead.in, method)
{
  L <- attr(sets, "lag")
  ts <- as.numeric(sub(".*\\.", "", names(sets)))
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  lsets <- sapply(sets, length)
  ts <- unlist(mapply(FUN = make.years, t = ts, lead = lead.in, repnum = lsets, SIMPLIFY = F))
  ids <- rep(unlist(sets), rep((2), length(unlist(sets)) ) )
  names(ids) <- NULL
  wts <- unlist(sapply(sets, function(s){return(attr(s, "weights"))}))
  names(wts) <- NULL
  if(method == "bootstrap")
  {
    wts <- rep(wts, rep((2), length(unlist(sets)))) * c(1,-1)
  }
  if(method == "wfe")
  {
    wts <- rep(wts, rep((2), length(unlist(sets)))) * c(-1,1)
  }
  num.empty <- sum(!sapply(sets, length) > 0)
  set.nums <- rep(0:(length(sets) - num.empty - 1), (lsets[lsets != 0] * (2)))
  data.frame(t = ts, id = ids, weight = wts, set.number = set.nums)
}
# refer to the description above -- pts works on treated units ("Prepare Treated unitS)
pts <- function(sets, lead.in, method)
{
  include <- sapply(sets, length) > 0
  num.empty <- sum(!include)
  tids <- as.numeric(sub("\\..*", "", names(sets)))[include]
  tids <- rep(tids, rep(2, length(tids)))
  ts <- as.numeric(sub(".*\\.", "", names(sets)))[include]
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  ts <- unlist(mapply(FUN = make.years, t = ts, lead = lead.in, repnum = 1, SIMPLIFY = F))
  if(method == "bootstrap")
  {
    wts <- rep(c(-1, 1), length(sets) - num.empty)
  }
  if(method == "wfe")
  {
    wts <- rep(c(1, 1), length(sets) - num.empty)
  }
  set.nums <- rep(0:(length(sets) - num.empty -1 ), rep(2, length(sets) - num.empty ))
  data.frame(t = ts, id = tids, weight = wts, set.number = set.nums)
}
#returns a vector of Wits, as defined in the paper (equation 25 or equation 23). They should be in the same order as the data frame containing the original problem data. The pts, pcs, and getWits functions act for a specific 
# lead. So, for instance if our lead window is 0,1,2,3,4, these function must be called for each of those -- so for 0, then for 1, etc.
getWits <- function(matched_sets, lead, data, estimation.method)
{
  
  #sort the data
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  include <- sapply(matched_sets, length) > 0
  num.empty <- sum(!include)
  
  
  #prep control sets, prep treatment sets for search/summation vector
  p.df <- pcs(matched_sets, lead, estimation.method)
  t.df <- pts(matched_sets, lead, estimation.method)
  
  t.idvector <- paste0(c(p.df$id, t.df$id), ".", c(p.df$t, t.df$t))
  setnums <- c(p.df$set.num, t.df$set.num)
  
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
}

#function that ultimately calculates the point estimate values
equality_four <- function(x, y, z){
    return(sum(x*y)/sum(z))
} 

  