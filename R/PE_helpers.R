#pcs and pts create data frames with the time/id combinations--that need to be found so that they can be easily looked up in the data frame via a hash table. The data frame also contains information
# about the weight of that unit at particular times, so we use the hashtable to look up where to put this data so that we can easily assign the appropriate weights in the original data frame containing the problem data.
# pcs does this for all control units in a matched set
# ("Prepare Control unitS)
#' @export
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
#' @export
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
#' @export
getWits <- function(matched_sets, lead, data, estimation.method)
{
  #sort the data
  
  t.var <- attr(matched_sets, "t.var")
  id.var <- attr(matched_sets, "id.var")
  data <- data[order(data[,id.var], data[,t.var]), ]
  include <- sapply(matched_sets, length) > 0
  num.empty <- sum(!include)
  
  vit.vect <- numeric(nrow(data) * (length(matched_sets) - num.empty))
  #prep control sets, prep treatment sets for search/summation vector
  p.df <- pcs(matched_sets, lead, estimation.method)
  t.df <- pts(matched_sets, lead, estimation.method)
  
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
# returns a vector of dit values, as defined in the paper. They should be in the same order as the data frame containing the original problem data.
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

# this will return a new matched.set object containing only treated/controls that can be used to calculate a point estimate. In particular, this function looks from time = t - 1 to t + lead
# and verifies that the necessary data is present. If not, those units are removed. The matched.set object that comes out of this function will need to be reweighted
#' @export
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
      stop("all controls in a particular matched set were removed due to missing data")
      
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
    
    pm2 <- reweight(sub.set.new, ordered.data, outcome.var)
    
    matched_sets[idx] <- pm2
    #matched_sets[idx] <- renormalize(sub.index, sub.set) #utilize the [.matched.set operator
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}
#function that ultimately calculates the point estimate values
equality_four <- function(x, y, z){
    return(sum(x*y)/sum(z))
} 



#' functions nearly identically to PM
#' This is currently only used inside PE2, where it is used to update the weights of a matched set object after control units that are missing data in the necessary future or past periods have been identified and taken out.
#' Function assumes that the units that are included are the only units it needs to consider. Pulls argument information from matched.set object 
#' Major difference between this function and the actual PanelMatch function is that this one does not identify matched sets and/or treated units. It merely updates and refines the matched set object provided to it
#' @export
reweight <- function(mset.object, data, outcome.var)
{
  lag = attr(mset.object, "lag")
  time.id = attr(mset.object, "t.var")
  unit.id = attr(mset.object, "id.var")
  treatment = attr(mset.object, "treated.var")
  outcome = outcome.var
  refinement.method = attr(mset.object, "refinement.method")
  size.match = attr(mset.object, "max.match.size")
  covs.formula = attr(mset.object, "covs.formula")
  match.missing <- attr(mset.object, "match.missing")
  verbose = FALSE
  msets <- mset.object
  
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    stop("panel data is not balanced")
  }
  
  #order and ensure order does not get violated
  #convert to data.table?
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment, outcome)]
  data <- data[, c(unit.id, time.id, treatment, outcome, othercols)] #reorder columns 
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  
  
  msets <- msets[sapply(msets, length) > 0 ]
  treated.ts <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(T,F)])
  
  ordered.data <- as.matrix(parse_and_prep(formula = covs.formula, data = ordered.data, unit.id = unit.id)) #every column > 4 at this point should be used in distance/refinement calculation
  ordered.data <- as.matrix(handle.missing.data(ordered.data, 5:ncol(ordered.data)))
  #RE IMPLEMENT RESTRICTED OR NAIVE?
  if(refinement.method == "mahalanobis")
  {
    
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, idx =  idxlist)
    weighted.mset <- handle_mahalanobis_calculations(mahalmats, msets, size.match, verbose)
    attr(weighted.mset, "covs.formula") <- covs.formula
    attr(weighted.mset, "match.missing") <- match.missing
    attr(weighted.mset, "max.match.size") <- size.match
    return(weighted.mset)
  }
  
  else
  {
    
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- rbindlist(expanded.sets.t0)
    pooled <- pre.pooled[complete.cases(pre.pooled), ]
    
    # Because we need certain data to be present (most notably, outcome variable data over t-1 to t + f), we will likely need to remove
    # some added columns about missing data because they will all be the same value (0 for instance because all of the data is guaranteed to be present)
    # Theoretically, we could have a situation where every value in the column is identical and thus cant be used for finding propensity scores
    # In these situations we end up with NA coefficients which causes a bunch of other code to break. This is a bit of a hack, but works for now.
    
    cols.to.remove <- which(unlist(lapply(pooled, function(x){all(x[1] == x)})))
    cols.to.remove <- unique(c(cols.to.remove, which(!colnames(pooled) %in% colnames(t(unique(t(pooled))))))) #removing columns that are identical to another column 
    if(length(cols.to.remove) > 0)
    {
      class(pooled) <- c("data.frame")
      pooled <- pooled[, -cols.to.remove]
      rmv <- function(x, cols.to.remove_)
      {
        return(x[, -cols.to.remove_])
      }
      expanded.sets.t0 <- lapply(expanded.sets.t0, rmv, cols.to.remove_ = cols.to.remove)
    }
    
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match")
    {
      fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:4)]), 
                                          family = binomial(link = "logit"), data = pooled))
    }
    if(refinement.method == "ps.weight" | refinement.method == "ps.match")
    {
      fit0 <- glm(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:4)]), 
                  family = binomial(link = "logit"), data = pooled)
    }
    
    just.ps.sets <- find_ps(expanded.sets.t0, fit0)
    
    if(refinement.method == "CBPS.weight" | refinement.method == "ps.weight")
    {
      msets <- handle_ps_weighted(just.ps.sets, msets, refinement.method)
    }
    if(refinement.method == "CBPS.match" | refinement.method == "ps.match")
    {
      msets <- handle_ps_match(just.ps.sets, msets, refinement.method, verbose, size.match)
    }
    
  }
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "match.missing") <- match.missing
  attr(msets, "max.match.size") <- size.match
  return(msets)    
}

  