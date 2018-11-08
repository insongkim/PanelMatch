#' @export
pcs <- function(sets, lead.in)
{
  L <- attr(sets, "lag")
  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])
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
  wts <- rep(wts, rep((2), length(unlist(sets)))) * c(1,-1)
  num.empty <- sum(!sapply(sets, length) > 0)
  set.nums <- rep(0:(length(sets) - num.empty - 1), (lsets[lsets != 0] * (2)))
  data.frame(t = ts, id = ids, weight = wts, set.number = set.nums)
}

#' @export
pts <- function(sets, lead.in)
{
  include <- sapply(sets, length) > 0
  num.empty <- sum(!include)
  tids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])[include]
  tids <- rep(tids, rep(2, length(tids)))
  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])[include]
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  ts <- unlist(mapply(FUN = make.years, t = ts, lead = lead.in, repnum = 1, SIMPLIFY = F))
  wts <- rep(c(-1, 1), length(sets) - num.empty)
  set.nums <- rep(0:(length(sets) - num.empty -1 ), rep(2, length(sets) - num.empty ))
  data.frame(t = ts, id = tids, weight = wts, set.number = set.nums)
}

#' @export
getWits <- function(matched_sets, lead, data)
{
  #sort the data
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  include <- sapply(matched_sets, length) > 0
  num.empty <- sum(!include)
  
  vit.vect <- numeric(nrow(data) * (length(matched_sets) - num.empty))
  #prep control sets, prep treatment sets for search/summation vector
  p.df <- pcs(matched_sets, lead)
  t.df <- pts(matched_sets, lead)
  
  t.idvector <- paste0(c(p.df$id, t.df$id), ".", c(p.df$t, t.df$t))
  setnums <- c(p.df$set.num, t.df$set.num)
  #function to figure out where to put each weight
  idxes <- get_vit_index(paste0(data[, id.var], ".", data[, t.var]), t.idvector, setnums)
  
  vit.vect[idxes] <- c(p.df$weight, t.df$weight)
  #summing the weights
  Wits <- sumwits(nrow(data), vit_vect = vit.vect) 
  #Wits <- sapply(1:nrow(data), function(i){sum(vit.vect[seq(i, length(vit.vect), by = nrow(data))])}) #also to be converted to c++ later more than likely
  
  return(Wits)
  
}

#' @export
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


#' @export
  prep_for_leads <- function(matched_sets, ordered.data, max.lead, t.var, id.var, outcome.var) #add other columns to check here
{
  #CHECK TO MAKE SURE COLUMNS ARE IN ORDER
  #browser()
  ordered.data <- ordered.data[order(ordered.data[,id.var], ordered.data[,t.var]), ]
  compmat <- data.table::dcast(data.table::as.data.table(ordered.data), formula = paste0(id.var, "~", t.var), value.var = outcome.var)
  ts <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(F,T)])
  tids <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(T,F)])
  class(matched_sets) <- "list" #so that Rcpp::List is accurate when we pass it into cpp functions
  compmat <- data.matrix(compmat)
  idx <- check_treated_units(compmat = compmat, compmat_row_units = as.numeric(compmat[, 1]), compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]), lead = max.lead, treated_ids = tids, treated_ts = ts)
  if(any(!idx))
  {
    class(matched_sets) <- c("matched.set", "list") #to get the matched.set subsetting with attributes
    matched_sets <- matched_sets[idx]
    ts <- ts[idx]
    
  }
  #colnames must be numeric in some form because we must be able to sort them into an ascending column order
  class(matched_sets) <- "list" # for rcpp reasons again
  ll <- re_norm_index(compmat = compmat, compmat_row_units = as.numeric(compmat[, 1]), compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]), lead = max.lead, sets = matched_sets, control_start_years = ts)
  idx <- needs_renormalization(ll)
  class(matched_sets) <- c("matched.set", "list")
  if(any(idx))
  {
    sub.index <- ll[idx]
    sub.set <- matched_sets[idx]
    matched_sets[idx] <- renormalize(sub.index, sub.set) #utilize the [.matched.set operator
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}
