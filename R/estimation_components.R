################################################################################
# Point estimates and standard errors first require the computation of a number of
# weights, parameters, indicator variables, as specified in Imai et al. (2023)
# This file has a variety of helper functions for calculating these
# Specifically, v,w,d in the paper
################################################################################

#' Prepare Control Units
#' pcs and pts create data frames with the time/id combinations--that need to be found so that they can be easily looked up in the data frame via a hash table. The data frame also contains information about the weight of that unit at particular times, so we use the hash table to look up where to put this data so that we can easily assign the appropriate weights in the original data frame containing the problem data. pcs does this for all control units in a matched set. pts does this for all treated units.
#' @param sets object describing the matched sets
#' @param lead.in integer describing a particular lead value.
#'
#' @return data.frame object with time-id combinations
#' @keywords internal
pcs <- function(sets, lead.in)
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
  wts <- unlist(sapply(sets, 
                       function(s){return(attr(s, "weights"))}))
  names(wts) <- NULL
  wts <- rep(wts, rep((2), length(unlist(sets)))) * c(1,-1)
  num.empty <- sum(!sapply(sets, length) > 0)
  set.nums <- rep(0:(length(sets) - num.empty - 1), (lsets[lsets != 0] * (2)))
  dtf <- data.frame(t = ts, id = ids, weight = wts, set.number = set.nums)
  return(dtf)
}

# refer to the description above 
pts <- function(sets, lead.in)
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
                      lead = lead.in, repnum = 1, 
                      SIMPLIFY = FALSE))
  wts <- rep(c(-1, 1), length(sets) - num.empty)
  set.nums <- rep(0:(length(sets) - num.empty - 1 ), 
                  rep(2, length(sets) - num.empty ))
  ldf = data.frame(t = ts, id = tids, weight = wts, 
                   set.number = set.nums)
  return(ldf)
  
  
}

#' getWits
#' returns a vector of Wits, as defined in the paper (equation 25 or equation 23). They should be in the same order as the data frame containing the original problem data. The pts, pcs, and getWits functions act for a specific lead. So, for instance if our lead window is 0,1,2,3,4, these function must be called for each of those -- so for 0, then for 1, etc.
#'
#' @param matched_sets matched.set object
#' @param lead integer providing a specific lead value
#' @param data data.frame object
#' @param estimation.method method of estimation for calculating standard errors.
#'
#' @return data.table of Wits, as described above
#' @keywords internal
getWits <- function(matched_sets, lead, data, 
                    estimation.method = "bootstrap")
{
  
  #sort the data
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  
  #prep control sets, prep treatment sets for search/summation vector
  p.df <- pcs(matched_sets, lead)
  t.df <- pts(matched_sets, lead)
  
  ##to solve cran check note about unbound variables
  . <- weight <- id <- NULL
  w.it.df <- rbind(p.df, t.df)
  w.it.df <- data.table::as.data.table(w.it.df[w.it.df$weight != 0,])
  summarized.Wits <- w.it.df[,.(Wit = sum(weight)), by = .(t,id)]
  return(summarized.Wits)
  
}


#' getDits
#' returns a vector of Dit values, as defined in the paper. They should be in the same order as the data frame containing the original problem data.
#'
#' @param matched_sets matched.set object
#' @param data data.frame object
#'
#' @return vector of Dits, as described in Imai et al. (2023)
#' @keywords internal
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