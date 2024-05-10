################################################################################
# Point estimates and standard errors first require the computation of a number of
# weights, parameters, indicator variables, as specified in Imai et al. (2021)
# This file has a variety of helper functions for calculating these
# specifically, v,w,d in the paper
################################################################################
#pcs and pts create data frames with the time/id combinations--that need to be found so that they can be easily looked up in the data frame via a hash table. The data frame also contains information
# about the weight of that unit at particular times, so we use the hashtable to look up where to put this data so that we can easily assign the appropriate weights in the original data frame containing the problem data.
# pcs does this for all control units in a matched set
# ("Prepare Control unitS")
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
  dtf <- data.frame(t = ts, id = ids, weight = wts, set.number = set.nums)
  return(dtf)
}
# refer to the description above -- pts works on treated units ("Prepare Treated unitS")
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
  set.nums <- rep(0:(length(sets) - num.empty - 1 ), 
                  rep(2, length(sets) - num.empty ))
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
# returns a data.table object
getWits <- function(matched_sets, lead, data, 
                    continuous.treatment = FALSE)
{
  
  #sort the data
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  
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