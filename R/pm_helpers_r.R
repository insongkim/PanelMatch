# File contains helper functions written in R for PanelMatch functionality
perform_refinement <- function(lag, time.id, unit.id, treatment, refinement.method, size.match, 
                               ordered.data, match.missing, covs.formula, verbose, mset.object = NULL, lead = NULL, outcome.var = NULL, restricted = FALSE)
{
  if(!is.null(mset.object))
  {
    lag = attr(mset.object, "lag")
    time.id = attr(mset.object, "t.var")
    unit.id = attr(mset.object, "id.var")
    treatment = attr(mset.object, "treated.var")
    refinement.method = attr(mset.object, "refinement.method")
    size.match = attr(mset.object, "max.match.size")
    covs.formula = attr(mset.object, "covs.formula")
    match.missing <- attr(mset.object, "match.missing")
    verbose = FALSE
    msets <- mset.object
  }
  else
  {
    temp.treateds <- findAllTreated(ordered.data, treatedvar = treatment, time.var = time.id, unit.var = unit.id, hasbeensorted = TRUE)
    if(nrow(temp.treateds) == 0) stop("no treated units")
    msets <- get.matchedsets(temp.treateds[, time.id], temp.treateds[, unit.id], ordered.data, lag, time.id, unit.id, treatment, hasbeensorted = TRUE)
    msets <- msets[sapply(msets, length) > 0 ]
    if(max(lead) > 0)
    {
      msets <- clean_leads(msets, ordered.data, max(lead), time.id, unit.id, outcome.var)  
    }
    if(restricted)
    {
      msets <- enforce_lead_restrictions(msets, ordered.data, max(lead), time.id, unit.id, treatment.var = treatment)
    }
  }
  treated.ts <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(T,F)])
  ordered.data <- as.matrix(parse_and_prep(formula = covs.formula, data = ordered.data, unit.id = unit.id)) #every column > 3 at this point should be used in distance/refinement calculation
  ordered.data <- as.matrix(handle.missing.data(ordered.data, 4:ncol(ordered.data)))
  
  #RE IMPLEMENT RESTRICTED OR NAIVE?
  if(refinement.method == "mahalanobis")
  {
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, idx =  idxlist)
    msets <- handle_mahalanobis_calculations(mahalmats, msets, size.match, verbose)
  }
  if(refinement.method == "ps.msm.match" | refinement.method == "CBPS.msm.match" | 
     refinement.method == "ps.msm.weight" | refinement.method == "CBPS.msm.weight")
  {
    store.msm.data <- list()
    for(i in 1:length(lead))
    {
        #tf = t + f
        f <- lead[i]
        tf <- expand.treated.ts(lag, treated.ts = treated.ts + f)
        tf.index <- get_yearly_dmats(ordered.data, treated.ids, tf, paste0(ordered.data[,unit.id], ".", 
                                                                           ordered.data[, time.id]), matched_sets = msets, lag)
        expanded.sets.tf <- build_ps_data(tf.index, ordered.data, lag)
        pre.pooled <- ordered.data[ordered.data[, time.id] %in% (treated.ts + f), ]
        pooled <- pre.pooled[complete.cases(pre.pooled), ]
        pooled <- as.data.frame(pooled)
        #do the column removal thing
        cols.to.remove <- which(unlist(lapply(pooled, function(x){all(x[1] == x)}))) #checking for columns that only have one value
        cols.to.remove <- unique(c(cols.to.remove, which(!colnames(pooled) %in% colnames(t(unique(t(pooled))))))) #removing columns that are identical to another column 
        if(length(cols.to.remove) > 0)
        {
          class(pooled) <- c("data.frame")
          pooled <- pooled[, -cols.to.remove]
          rmv <- function(x, cols.to.remove_)
          {
            return(x[, -cols.to.remove_])
          }
          expanded.sets.tf <- lapply(expanded.sets.tf, rmv, cols.to.remove_ = cols.to.remove)
        }
        if(qr(pooled)$rank != ncol(pooled)) stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
        
        if(refinement.method == "CBPS.msm.match" | refinement.method == "CBPS.msm.weight") #obviously update these conditionals
        {
          fit.tf <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                                                family = binomial(link = "logit"), data = pooled))
        }
        else
        {
          fit.tf <- glm(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                        family = binomial(link = "logit"), data = pooled)
        }
        store.msm.data[[i]] <- find_ps(expanded.sets.tf, fit.tf)
        
      }
    
    msm.sets <- gather_msm_sets(store.msm.data)
    if(refinement.method == "CBPS.msm.match" | refinement.method == "ps.msm.match")
    {
      msets <- handle_ps_match(msm.sets, msets = msets, refinement.method, verbose = verbose, max.set.size = size.match)
      attr(msets, "max.match.size") <- size.match
    }
    if(refinement.method == "CBPS.msm.weight" | refinement.method == "ps.msm.weight")
    {
      msets <- handle_ps_weighted(msm.sets, msets, refinement.method)
    }
  } else #not msm
  {
    if(!all(refinement.method %in% c("CBPS.weight", "CBPS.match", "ps.weight", "ps.match"))) stop("please choose valid refinement method")
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- rbindlist(expanded.sets.t0)
    pooled <- pre.pooled[complete.cases(pre.pooled), ]
  
    cols.to.remove <- which(unlist(lapply(pooled, function(x){all(x[1] == x)}))) #checking for columns that only have one value
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
    if(qr(pooled)$rank != ncol(pooled)) stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match")
    {
      fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                                          family = binomial(link = "logit"), data = pooled))
    }
    if(refinement.method == "ps.weight" | refinement.method == "ps.match")
    {
      fit0 <- glm(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
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
      attr(msets, "max.match.size") <- size.match
    }
  }
  
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "match.missing") <- match.missing
  return(msets)
}
# builds a list that contains all times in a lag window that correspond to a particular treated unit. This is structured as a list of vectors. Each vector is lag + 1 units long. The overall list will 
# be the same length as the number of matched sets
expand.treated.ts <- function(lag, treated.ts)
{
  helper <- function(treated.t)
  {
    return(seq(from =  (treated.t - lag), to = treated.t, by = 1))
  }
  lapply(treated.ts, helper)
}

#builds the matrices that we will then use to calculate the mahalanobis distances for each matched set
build_maha_mats <- function(idx, ordered_expanded_data)
{
  subset.per.matchedset <- function(sub.idx)
  {
    ordered_expanded_data[sub.idx,]
  }
  unnest <- function(mset.idx)
  {
    lapply(mset.idx, subset.per.matchedset)
  }
  result <- lapply(idx, unnest)
  return(result)
}

# Will need to be updated if the syntax/implementation of the covs.formula argument is changed
# Applies necessary transformation to the data based on the specified parameters, including using the covs.formula argument to apply the necessary lags to particular columns.
# It also will structure the data frame such that every column after the 3rd column is used in the following calculations
parse_and_prep <- function(formula, data, unit.id)
{
  #check if formula empty or no covariates provided -- error checks
  terms <- attr(terms(formula),"term.labels")
  terms <- gsub(" ", "", terms) #remove whitespace
  lag.calls <- terms[grepl("lag(*)", terms)] #regex to get calls to lag function
  other.terms <- terms[!grepl("lag(*)", terms)]
  sub.data <- data[, c(1:3, which(colnames(data) %in% other.terms) )] #including only what is specified in the formula
  if(any(grepl("=", lag.calls))) stop("fix lag calls to use only unnamed arguments in the correct positions")
  data <- data.table::as.data.table(data) #check sorting
  if(length(lag.calls) > 0)
  {
    results.unmerged <- mapply(FUN = handle.calls, call.as.string = lag.calls, MoreArgs =  list(.data = data, .unitid = unit.id), SIMPLIFY = FALSE)
    names(results.unmerged) <- NULL
    full.data <- cbind(sub.data, do.call("cbind", results.unmerged))
    return(full.data)
  }
  else
  {
    return(sub.data)
  }
  
}
#helper function that deals with the covs.formula argument. This is what specifically will need to be modified if the format of that argument is changed. It handles the lagging of the various covariates
# in the process of preparing the data frame for the remaining calculations.
handle.calls <- function(call.as.string, .data, .unitid)
{
  lag <- function(var.name, lag.window, data = .data, unitid = .unitid) #want to make sure its being called appropriately
  {
    tr <- data[, data.table::shift(.SD, n = lag.window, fill = NA, type = "lag", give.names = TRUE), .SDcols=var.name, by = unitid]
    tr <- tr[, -1]
    colnames(tr) <- paste0(var.name, "_l", lag.window)
    return(tr)
  }
  return(eval(parse(text = call.as.string)))
}

#use col.index to determine which columns we want to "scan" for missing data
# Note that in earlier points in the code, we rearrange the columns and prepare the data frame such that cols 1-4 are bookkeeping (unit id, time id, treated variable, unlagged outcome variable)
# and all remaining columns are used in the calculations after going through parse_and_prep function, so col.index should usually be 5:ncol(data)
# In practice, this function just looks over the data in the specified columns in the "data" data frame for missing data. Then it creates columns with indicator variables about the missingness of those variables
# 1 for missing data, 0 for present
handle.missing.data <- function(data, col.index)
{
  
  new.names <- paste0(colnames(data)[col.index], "_NA")
  if(length(col.index) == 1)
  {
    new.col <- is.na(data[, col.index]) * 1
    data <- cbind(data, new.col)
    colnames(data)[ncol(data)] <- new.names
  }
  else
  {
    missing.mat <- as.data.frame(apply(data[, col.index], 2, is.na)) * 1
    colnames(missing.mat) <- new.names
    data <- cbind(data, missing.mat)  
  }
  
  data[is.na(data)] <- 0
  return(data)
}

# prepares the data for calculating propensity scores. Will return a list of length equal to the number of matched sets. Each item is a data frame and each data frame contains information at time = t + 0
# for each treated unit and their corresponding controls. 
build_ps_data <- function(idxlist, data, lag)
{
  obtain.t.rows <- function(idx)
  {
    return(idx[length(idx)])
  }
  unnest <- function(subidxlist,  lag)
  {
    temp <- sapply(subidxlist[[lag + 1]], obtain.t.rows)
    return(data.frame(data[temp, ]))
  }
  results <- lapply(idxlist, unnest, lag = lag)
  #results <- rbindlist(results)
  #results <- results[complete.cases(results), ]
  return(results)
}

#returns a list of data frames with propensity scores for each unit in a matched set. Each element in the list is a data frame which corresponds to a matched set of 1 treatment and all
# matched control units #NOTE: NOT THE PROPENSITY SCORE? ACTUALLY THE WEIGHTS...RENAME THIS AT SOME POINT
find_ps <- function(sets, fitted.model)
{
  
  apply_formula <- function (x, B) 
  {
    xx <- cbind(1, as.matrix(x[, 4:ncol(x)]))
    x[, (ncol(x) + 1)] <- 1 - 1/(1+exp(xx %*% B))
    names(x)[ncol(x)] <- "ps"
    return(x[, c(1:3, ncol(x))])
  }
  sets_with_ps <- lapply(sets, apply_formula, B = fitted.model$coefficients)
  return(sets_with_ps)
}

# Each of the following similarly named functions carry out the refinement procedures -- using either propensity scores or mahalanobis distances according to the paper. Each of these functions
# will return a matched.set object with the appropriate weights for control units assigned and new, additional attributes (such as refinement method). 
# These functions could use some work to be better optimized. For instance, when the "verbose" argument is set to true, they will essentially do all of the refinement calculations twice.
handle_mahalanobis_calculations <- function(mahal.nested.list, msets, max.size, verbose)
{
  
  do.calcs <- function(year.df)
  { 
    if(nrow(year.df) == 2)
    {
      return(1)
    }
    
    cov.data <- year.df[1:(nrow(year.df) - 1), 4:ncol(year.df)]
    cov.matrix <- cov(cov.data)
    center.data <- year.df[nrow(year.df), 4:ncol(year.df)]
    
    if( isTRUE(all.equal(det(cov.matrix), 0, tolerance = .00001)) ) #we might want to make this smaller, but had some errors here about computationally infeasible problems because of values very close to zero
    {
      cov.matrix <- ginv(cov.matrix)
      return(mahalanobis(x = cov.data, center = center.data, cov = cov.matrix, inverted = TRUE))
    }
    else
    {
      
      return(mahalanobis(x = cov.data, center = center.data, cov = cov.matrix))
    }
    
  }
  handle_set <- function(sub.list, max.set.size, idx)
  {
    results.temp <- lapply(sub.list, do.calcs)
    tmat <- do.call(rbind, results.temp)
    colnames(tmat) <- NULL
    dists <- colMeans(tmat)
    n.dists <- dists[dists > 0]
    if(length(n.dists) == 0 ) stop('matched set contains only identical units')
    if(length(n.dists) < max.set.size)
    {
      w <- 1 / length(n.dists)
      newdists <- dists
      newdists[newdists > 0 ] <- w
    }
    else
    {
      ordered.dists <- sort(n.dists)
      scoretobeat <- max(head(ordered.dists, n = max.set.size + 1))
      # might have situation where the Mth largest distance is the same as the Mth - 1 distance. This means that we either choose to leave out both and have a matched set smaller than the max, 
      # or include both of them and relax the size of our maximum set size
      if(sum(dists < scoretobeat & dists > 0) < max.set.size) #change this if we want to be more strict about max.set.size enforcements
      {
        new.denom <- sum(dists <= scoretobeat & dists > 0)
        newdists <- ifelse(dists <= scoretobeat & dists > 0, 1 / new.denom, 0)
      }
      else
      {
        newdists <- ifelse(dists < scoretobeat & dists > 0, 1 / max.set.size, 0)  
      }
      
    }
    names(newdists) <- NULL
    return(newdists)
    
  }
  scores <- mapply(FUN = handle_set, sub.list = mahal.nested.list, idx = 1:length(msets), MoreArgs = list(max.set.size = max.size), SIMPLIFY = FALSE)
  for(i in 1:length(msets))
  {
    names(scores[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- scores[[i]]
  }
  if(verbose) #in future versions, avoid doing the same calculations twice
  {
    handle_set_verbose <- function(sub.list)
    {
      results.temp <- lapply(sub.list, do.calcs)
      dists <- colMeans(do.call(rbind, results.temp))
      names(dists) <- NULL
      return(dists)
    }
    
    full.scores <- mapply(FUN = handle_set_verbose, sub.list = mahal.nested.list, SIMPLIFY = FALSE)
    
    for(i in 1:length(msets))
    {
      names(full.scores[[i]]) <- msets[[i]]
      attr(msets[[i]], "distances") <- full.scores[[i]]
    }
  }
  
  attr(msets, "refinement.method") <- "mahalanobis"
  return(msets)
}


handle_ps_weighted <- function(just.ps.sets, msets, refinement.method)
{
  
  handle_set <- function(set)
  {
    control.ps.set <- set[1:(nrow(set) - 1), ncol(set)]
    if(length(control.ps.set) == 1)
    {
      return(1)
    }
    vec.ratio <- control.ps.set / (1 - control.ps.set) #just for clarity
    wts <- ( vec.ratio ) / sum( vec.ratio )
    return(as.vector(wts))
  }
  wts <- lapply(just.ps.sets, handle_set)
  for(i in 1:length(msets))
  {
    names(wts[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- wts[[i]]
  }
  attr(msets, "refinement.method") <- refinement.method
  return(msets)
}

handle_ps_match <- function(just.ps.sets, msets, refinement.method, verbose, max.set.size)
{
  handle_set <- function(set, max.size)
  {
    treated.ps <- as.numeric(set[nrow(set), "ps"])
    control.ps.set <- as.numeric(set[1:(nrow(set) - 1), "ps"])
    if(length(control.ps.set) == 1)
    {
      return(1)
    }
    dists <- abs(treated.ps - control.ps.set)
    dists.to.consider <- dists[dists > 0]
    if(length(dists.to.consider) < max.size)
    {
      dists[ dists > 0 ] <- 1 / length(dists.to.consider)
      wts <- dists
    }
    else
    {
      dist.to.beat <- max(head(sort(dists.to.consider), max.size + 1))
      if(sum(dists < dist.to.beat & dists > 0) < max.set.size)
      {
        new.denom <- sum(dists <= dist.to.beat & dists > 0)
        wts <- ifelse(dists <= dist.to.beat & dists > 0, 1 / new.denom, 0)
        
      }
      else
      {
        wts <- ifelse(dists < dist.to.beat & dists > 0, (1 / max.size), 0)    
      }
      
      #wts <- ifelse(dists < dist.to.beat & dists > 0, (1 / max.size), 0)  
    }
    return(wts)
  }
  wts <- lapply(just.ps.sets, handle_set, max.size = max.set.size)
  for(i in 1:length(msets))
  {
    names(wts[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- wts[[i]]
  }
  if(verbose) #again this is not well designed, would want to avoid having to do everything again, but works for now.
  {
    handle_set <- function(set, max.size)
    {
      treated.ps <- as.numeric(set[nrow(set), "ps"])
      control.ps.set <- as.numeric(set[1:(nrow(set) - 1), "ps"])
      if(length(control.ps.set) == 1)
      {
        return(1)
      }
      dists <- abs(treated.ps - control.ps.set)
      return(dists)
    }
    dts <- lapply(just.ps.sets, handle_set, max.size = max.set.size)
    for(i in 1:length(msets))
    {
      names(dts[[i]]) <- msets[[i]]
      attr(msets[[i]], "distances") <- dts[[i]]
    }
  }
  
  attr(msets, "refinement.method") <- refinement.method
  return(msets)
}
#right now this function just checks outcome data and cleans up based on that, but when msm is implemented, we will also need to check reversion of treatment
clean_leads <- function(matched_sets, ordered.data, max.lead, t.var, id.var, outcome.var) 
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
      warning("all controls in a particular matched set were removed due to missing data")
      
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
    #pm2 <- perform_refinement(ordered.data = ordered.data, mset.object = sub.set.new)
    #matched_sets[idx] <- pm2
    matched_sets[idx] <- sub.set.new
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}


enforce_lead_restrictions <- function(matched_sets, ordered.data, max.lead, t.var, id.var, treatment.var) 
{
  #CHECK TO MAKE SURE COLUMNS ARE IN ORDER
  
  ordered.data <- ordered.data[order(ordered.data[,id.var], ordered.data[,t.var]), ]
  compmat <- data.table::dcast(data.table::as.data.table(ordered.data), formula = paste0(id.var, "~", t.var), value.var = treatment.var)
  ts <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(F,T)])
  tids <- as.numeric(unlist(strsplit(names(matched_sets), split = "[.]"))[c(T,F)])
  class(matched_sets) <- "list" #so that Rcpp::List is accurate when we pass it into cpp functions
  compmat <- data.matrix(compmat)
  
  idx <- check_treated_units_for_treatment_reversion(compmat = compmat, compmat_row_units = as.numeric(compmat[, 1]), 
                                                     compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]),
                                                     lead = max.lead, treated_ids = tids, treated_ts = ts)
  if(all(!idx)) stop("estimation not possible: All treated units are missing data necessary for the calculations to proceed")
  if(any(!idx))
  {
    class(matched_sets) <- c("matched.set", "list") #to get the matched.set subsetting with attributes
    matched_sets <- matched_sets[idx]
    ts <- ts[idx]
    
  }
  #colnames must be numeric in some form because we must be able to sort them into an ascending column order
  class(matched_sets) <- "list" # for rcpp reasons again
  ll <- check_control_units_for_treatment_restriction(compmat = compmat, 
                                                      compmat_row_units = as.numeric(compmat[, 1]), 
                                                      compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]), 
                                                      lead = max.lead, sets = matched_sets, control_start_years = ts)
  #probably should rename this function, but working in a similar context here so seeing if it works
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
      warning("all controls in a particular matched set were removed due to missing data")
      
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
    #pm2 <- perform_refinement(ordered.data = ordered.data, mset.object = sub.set.new)
    #matched_sets[idx] <- pm2
    matched_sets[idx] <- sub.set.new
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}

gather_msm_sets <- function(lead.data.list)
{
  number.of.sets <- sapply(lead.data.list, length)
  if(length(unique(number.of.sets)) != 1) stop("error with matched sets in msm calculations")
  number.of.sets <- unique(number.of.sets)
  long.data.lead.list <- unlist(lead.data.list, recursive = F)
  
  long.weights.list <- lapply(long.data.lead.list, function(x){return(as.vector(x[, 4]))})
  multiplied.weights <-  multiply_weights_msm(long.weights.list, number.of.sets)
  reassembled.sets <- long.data.lead.list[1:number.of.sets]
  reassemble.weights <- function(set, weights)
  {
    set[, "ps"] <- weights #again this ps is misleading but for consistency with the other functions lets go with it
    return(set)
  }
  reassembled.sets <- mapply(FUN = reassemble.weights, set = reassembled.sets, weights = multiplied.weights, SIMPLIFY = F)
  
  return(reassembled.sets)
}
