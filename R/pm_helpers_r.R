#' @export
subset_expanded.data <- function(expanded_data, matched_set_indices)
{
  create_subset <- function(index, expanded.data)
  {
    return(data.frame(expanded.data[index, ]))
  }
  res <- lapply(matched_set_indices, FUN = create_subset, expanded.data = expanded_data)
  return(res)
}

expand.treated.ts <- function(lag, treated.ts)
{
  helper <- function(treated.t)
  {
    return(seq(from =  (treated.t - lag), to = treated.t, by = 1))
  }
  lapply(treated.ts, helper)
}

##OPTIMIZE THIS
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

#needs more error checking, details will also depend on the syntax/mechanics of the covs.formula argument to be determined later.
parse_and_prep <- function(formula, data, unit.id)
{
  #check if formula empty or no covariates provided -- error checks
  terms <- attr(terms(formula),"term.labels")
  terms <- gsub(" ", "", terms) #remove whitespace
  lag.calls <- terms[grepl("lag(*)", terms)] #regex to get calls to lag function
  other.terms <- terms[!grepl("lag(*)", terms)]
  sub.data <- data[, c(1:4, which(other.terms == colnames(data)))] #including only what is specified in the formula
  
  if(any(grepl("=", lag.calls))) stop("fix lag calls to use only unnamed arguments in the correct positions")
  data <- data.table::as.data.table(data) #check sorting
  results.unmerged <- mapply(FUN = handle.calls, call.as.string = lag.calls, MoreArgs =  list(.data = data, .unitid = unit.id), SIMPLIFY = FALSE)
  names(results.unmerged) <- NULL
  full.data <- cbind(sub.data, do.call("cbind", results.unmerged))
  return(full.data)
}

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
handle.missing.data <- function(data, col.index)
{
  new.names <- paste0(colnames(data)[col.index], "_NA")  
  missing.mat <- as.data.frame(apply(data[, col.index], 2, is.na)) * 1
  colnames(missing.mat) <- new.names
  data <- cbind(data, missing.mat)
  data[is.na(data)] <- 0
  return(data)
}

#OPTIMIZATION OPPORTUNITIES HERE
handle_mahalanobis_calculations <- function(mahal.nested.list, msets, max.size, verbose)
{
  #browser()
  do.calcs <- function(year.df)
  { 
    if(nrow(year.df) == 2)
    {
      return(1)
    }
    
    cov.data <- year.df[1:(nrow(year.df) - 1), 5:ncol(year.df)]
    cov.matrix <- cov(cov.data)
    center.data <- year.df[nrow(year.df), 5:ncol(year.df)]
    
    if( isTRUE(all.equal(det(cov.matrix), 0, tolerance = .00001)) ) #using the default tolerance value
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
    #if(idx == 55) browser()
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
      newdists <- ifelse(dists < scoretobeat & dists > 0, 1 / max.set.size, 0)
    }
    names(newdists) <- NULL
    return(newdists)
    
  }
  scores <- mapply(FUN = handle_set, sub.list = mahal.nested.list, idx = 1:length(msets), MoreArgs = list(max.set.size = max.size))
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
    
    full.scores <- mapply(FUN = handle_set_verbose, sub.list = mahal.nested.list)
    
    for(i in 1:length(msets))
    {
      names(full.scores[[i]]) <- msets[[i]]
      attr(msets[[i]], "distances") <- full.scores[[i]]
    }
  }
  
  attr(msets, "refinement.method") <- "mahalanobis"
  return(msets)
}

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

build_expanded_sets_for_coef_mult <- function(idxlist, data)
{
  unnest <- function(subidxlist)
  {
    temp.index <- unlist(subidxlist) #again, need to be sure that the first four columns can be thrown out
    x <- as.data.frame(cbind(1, data[temp.index, 5:ncol(data)]))
    return(x)
  }
  results <- lapply(idxlist, unnest)
}

find_ps <- function(sets, fitted.model)
{
  apply_formula <- function (x, B) 
  {
    xx <- cbind(1, as.matrix(x[, 5:ncol(x)]))
    x[, (ncol(x) + 1)] <- 1 - 1/(1+exp(xx %*% B))
    names(x)[ncol(x)] <- "ps"
    return(x[, c(1:4, ncol(x))])
  }
  sets_with_ps <- lapply(sets, apply_formula, B = fitted.model$coefficients)
  return(sets_with_ps)
}

handle_ps_weighted <- function(just.ps.sets, msets, refinement.method)
{
  #browser()
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
        wts <- ifelse(dists < dist.to.beat & dists > 0, (1 / max.size), 0)  
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

