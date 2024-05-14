#' build_maha_mats
#' Builds the matrices that we will then use to calculate the mahalanobis distances for each matched set
#' @param idx List of vectors specifying which observations should be extracted
#' @param ordered_expanded_data data.frame of prepared/parsed input data
#'
#' @returns List of parsed distance matrices, with elements corresponding to each matched set
#' @keywords internal
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

#' handle_mahalanobis_calculations
#' Returns a matched.set object with weights for control units, along with some other metadata
#'
#' @param mahal.nested.list Output from build_maha_mats function 
#' @param msets matched.set object -- list containing the treated observations and matched controls
#' @param max.size maximum number of control units that will receive non-zero weights within a matched set
#' @param verbose Logical. See PanelMatch() documentation
#' @param use.diagonal.covmat Logical. See PanelMatch() documentation
#'
#' @return matched.set object with weights for control units, along with some other metadata
#' @keywords internal
handle_mahalanobis_calculations <- function(mahal.nested.list, 
                                            msets,
                                            max.size, 
                                            verbose, 
                                            use.diagonal.covmat)
{
  do.calcs <- function(year.df)
  {
    if(nrow(year.df) == 2)
    {
      return(1)
    }
    cov.data <- year.df[1:(nrow(year.df) - 1), 
                        4:ncol(year.df), drop = FALSE]
    if(use.diagonal.covmat)
    {
      cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
    } else
    {
      cov.matrix <- cov(cov.data)
    }
    #### try to preemptively address some simple computational problems
    center.data <- year.df[nrow(year.df), 4:ncol(year.df), drop = FALSE]
    if(isTRUE(all.equal(det(cov.matrix), 0, tolerance = .00001))) #setting tolerance here seemed to reduce some computational/instability issues.
    {
      cols.to.remove <- which(apply(cov.data, 2, 
                                    function(x) isTRUE(length(unique(x)) == 1))) #checking for columns that only have one value
      cols.to.remove <- unique(c(cols.to.remove, 
                                 which(!colnames(cov.data) %in% colnames(t(unique(t(cov.data))))))) #removing columns that are identical to another column
      if(length(cols.to.remove) > 0 & length(cols.to.remove) < ncol(cov.data))
      {
        cov.data <- cov.data[, -cols.to.remove, drop = FALSE]
        center.data <- center.data[-cols.to.remove, drop = FALSE]
        if(use.diagonal.covmat)
        {
          cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
        } else {
          cov.matrix <- cov(cov.data)
        }
      }
      
    }
    
    result = tryCatch({
      mahalanobis(x = cov.data, center = center.data, cov = cov.matrix)
    }, warning = function(w) {
      
    }, error = function(e) {
      cov.matrix <- cov(cov.data)
      cov.matrix <- ginv(cov.matrix)
      mahalanobis(x = cov.data, center = center.data, cov = cov.matrix, inverted = TRUE)
    }, finally = {
    })
    
    return(result)
    
    
  }
  handle_set <- function(sub.list, max.set.size, idx)
  {
    
    results.temp <- lapply(sub.list, do.calcs)
    tmat <- do.call(rbind, results.temp)
    colnames(tmat) <- NULL
    dists <- colMeans(tmat)
    
    n.dists <- dists
    
    
    if(length(n.dists) < max.set.size) #case where total number of units in matched set < max.set size
    {
      w <- 1 / length(n.dists)
      newdists <- dists
      newdists <- rep(w, length(newdists))
      
    }
    else
    {
      ordered.dists <- sort(n.dists)
      scoretobeat <- max(utils::head(ordered.dists, 
                                     n = max.set.size + 1))
      # might have situation where the Mth largest distance is the same as the Mth - 1 distance. This means that we either choose to leave out both and have a matched set smaller than the max,
      # or include both of them and relax the size of our maximum set size
      
      # rounding distances to avoid computational bugs/inconsistencies that popped up
      dists <- round(dists, 5)
      scoretobeat <- round(scoretobeat, 5)
      if(sum(dists < scoretobeat) < max.set.size) 
      {
        new.denom <- sum(dists <= scoretobeat)
        newdists <- ifelse(dists <= scoretobeat, 1 / new.denom, 0)
      }
      else
      {
        newdists <- ifelse(dists < scoretobeat, 1 / max.set.size, 0)
      }
      
    }
    names(newdists) <- NULL
    return(newdists)
    
  }
  
  scores <- mapply(FUN = handle_set, sub.list = mahal.nested.list,
                   idx = 1:length(msets),
                   MoreArgs = list(max.set.size = max.size), SIMPLIFY = FALSE)
  for(i in 1:length(msets))
  {
    names(scores[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- scores[[i]]
  }
  if(verbose) 
  {
    handle_set_verbose <- function(sub.list)
    {
      results.temp <- lapply(sub.list, do.calcs)
      dists <- colMeans(do.call(rbind, results.temp))
      names(dists) <- NULL
      return(dists)
    }
    
    full.scores <- mapply(FUN = handle_set_verbose, 
                          sub.list = mahal.nested.list,
                          SIMPLIFY = FALSE)
    
    for(i in 1:length(msets))
    {
      names(full.scores[[i]]) <- msets[[i]]
      attr(msets[[i]], "distances") <- full.scores[[i]]
    }
  }
  
  attr(msets, "refinement.method") <- "mahalanobis"
  return(msets)
}


#' handle_ps_match
#' Returns a matched.set object with weights for control units, along with some other metadata
#'
#' @param just.ps.sets Output from find_ps() function 
#' @param msets matched.set object -- list containing the treated observations and matched controls
#' @param max.set.size maximum number of control units that will receive non-zero weights within a matched set
#' @param verbose Logical. See PanelMatch() documentation
#'
#' @return matched.set object with weights for control units, along with some other metadata
#' @keywords internal
handle_ps_match <- function(just.ps.sets, msets,
                            refinement.method,
                            verbose, max.set.size)
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
      dist.to.beat <- max(utils::head(sort(dists.to.consider), max.size + 1))
      if(sum(dists < dist.to.beat & dists > 0) < max.set.size)
      {
        new.denom <- sum(dists <= dist.to.beat & dists > 0)
        wts <- ifelse(dists <= dist.to.beat & dists > 0, 1 / new.denom, 0)
        
      }
      else
      {
        wts <- ifelse(dists < dist.to.beat & dists > 0, (1 / max.size), 0)
      }
    }
    return(wts)
  }
  wts <- lapply(just.ps.sets, handle_set, max.size = max.set.size)
  for(i in 1:length(msets))
  {
    names(wts[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- wts[[i]]
  }
  if(verbose) 
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