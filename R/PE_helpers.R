#pcs and pts create data frames with the time/id combinations--that need to be found so that they can be easily looked up in the data frame via a hash table. The data frame also contains information
# about the weight of that unit at particular times, so we use the hashtable to look up where to put this data so that we can easily assign the appropriate weights in the original data frame containing the problem data.
# pcs does this for all control units in a matched set
# ("Prepare Control unitS)
pcs <- function(sets, lead.in, method)
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
  tids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])[include]
  tids <- rep(tids, rep(2, length(tids)))
  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])[include]
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

# this will return a new matched.set object containing only treated/controls that can be used to calculate a point estimate. In particular, this function looks from time = t - 1 to t + lead
# and verifies that the necessary data is present. If not, those units are removed. The matched.set object that comes out of this function will need to be reweighted
prep_for_leads <- function(matched_sets, ordered.data, max.lead, t.var, id.var, outcome.var) 
{
  #CHECK TO MAKE SURE COLUMNS ARE IN ORDER
  
  ordered.data <- ordered.data[order(ordered.data[,id.var], ordered.data[,t.var]), ]
  compmat <- data.table::dcast(data.table::as.data.table(ordered.data), formula = paste0(id.var, "~", t.var), value.var = outcome.var)
  ts <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(F,T)])
  tids <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(T,F)])
  class(matched_sets) <- "list" #so that Rcpp::List is accurate when we pass it into cpp functions
  compmat <- data.matrix(compmat)
  
  idx <- check_treated_units(compmat = compmat, compmat_row_units = as.numeric(compmat[, 1]), compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]), lead = max.lead, treated_ids = tids, treated_ts = ts)
  if(all(!idx)) stop("estimation not possible: All treated units are missing data necessary for the calculations to proceed")
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
    create_new_sets <- function(set, index)
    {
      return(set[index])
    }
    sub.set.new <- mapply(FUN = create_new_sets, sub.set, sub.index, SIMPLIFY = FALSE)
    attributes(sub.set.new) <- attributes(sub.set)
    all.gone.counter <- sapply(sub.set.new, function(x){sum(x)})
    if(sum(all.gone.counter == 0) > 0) #case in which all the controls in a particular group were dropped
    {
      #warning("all controls in a particular matched set were removed due to missing data")
      
      idx[all.gone.counter == 0] <- FALSE
      sub.index <- ll[idx]
      sub.set <- matched_sets[idx]
      create_new_sets <- function(set, index)
      {
        return(set[index])
      }
      sub.set.new <- mapply(FUN = create_new_sets, sub.set, sub.index, SIMPLIFY = FALSE)
      attributes(sub.set.new) <- attributes(sub.set)
    }
    if(all(sapply(sub.set.new, length) == 0)) stop('estimation not possible: none of the matched sets have viable control units due to a lack of necessary data')
    pm2 <- perform_refinement(ordered.data = ordered.data, mset.object = sub.set.new)
    
    
    matched_sets[idx] <- pm2
    
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}
#function that ultimately calculates the point estimate values
equality_four <- function(x, y, z){
    return(sum(x*y)/sum(z))
} 

  