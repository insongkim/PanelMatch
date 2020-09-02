
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

# # Each of the following similarly named functions carry out the refinement procedures -- using either propensity scores or mahalanobis distances according to the paper. Each of these functions
# # will return a matched.set object with the appropriate weights for control units assigned and new, additional attributes (such as refinement method).
# # These functions could use some work to be better optimized. For instance, when the "verbose" argument is set to true, they will essentially do all of the refinement calculations twice.
# handle_mahalanobis_calculations <- function(mahal.nested.list, msets, 
#                                             max.size, verbose, use.diagonal.covmat)
# {
#   do.calcs <- function(year.df)
#   {
#     if(nrow(year.df) == 2)
#     {
#       return(1)
#     }
#     cov.data <- year.df[1:(nrow(year.df) - 1), 4:ncol(year.df), drop = FALSE]
#     if(use.diagonal.covmat)
#     {
#       cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
#     } else
#     {
#       cov.matrix <- cov(cov.data)
#     }
#     
#     center.data <- year.df[nrow(year.df), 4:ncol(year.df), drop = FALSE]
#     if(isTRUE(all.equal(det(cov.matrix), 0, tolerance = .00001))) #might not be the conditions we want precisely
#     {
#       cols.to.remove <- which(apply(cov.data, 2, function(x) isTRUE(length(unique(x)) == 1))) #checking for columns that only have one value
#       cols.to.remove <- unique(c(cols.to.remove, 
#                                  which(!colnames(cov.data) %in% colnames(t(unique(t(cov.data))))))) #removing columns that are identical to another column
#       if(length(cols.to.remove) > 0 & length(cols.to.remove) < ncol(cov.data))
#       {
#         cov.data <- cov.data[, -cols.to.remove, drop = FALSE]
#         center.data <- center.data[-cols.to.remove, drop = FALSE]
#         if(use.diagonal.covmat)
#         {
#           cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
#         } else {
#           cov.matrix <- cov(cov.data)
#         }
#       }
#       
#     }
#     
#     result = tryCatch({
#       mahalanobis(x = cov.data, center = center.data, cov = cov.matrix)
#     }, warning = function(w) {
#       
#     }, error = function(e) {
#       cov.matrix <- cov(cov.data)
#       cov.matrix <- ginv(cov.matrix)
#       mahalanobis(x = cov.data, center = center.data, cov = cov.matrix, inverted = TRUE)
#     }, finally = {
#       #(mahalanobis(x = cov.data, center = center.data, cov = cov.matrix))
#     })
#     
#     return(result)
#     
#     
#   }
#   handle_set <- function(sub.list, max.set.size)
#   {
#     
#     #results.temp <- lapply(sub.list, do.calcs)
#     #tmat <- do.call(rbind, results.temp)
#     
#     tmat <- do.calcs(sub.list)
#     colnames(tmat) <- NULL
#     #dists <- colMeans(tmat)
#     #n.dists <- dists[dists > 0]
#     dists <- tmat
#     #n.dists <- dists
#     if(length(dists) < max.set.size) #case where total number of units in matched set < max.set size
#     {
#       w <- 1 / length(dists)
#       newdists <- dists
#       newdists <- rep(w, length(newdists))
#       #newdists[newdists > 0 ] <- w
#     }
#     else
#     {
#       ordered.dists <- sort(dists)
#       scoretobeat <- max(utils::head(ordered.dists, n = max.set.size + 1))
# 
#       if(sum(dists < scoretobeat) < max.set.size) #change this if we want to be more strict about max.set.size enforcements
#       {
#         new.denom <- sum(dists <= scoretobeat)
#         newdists <- ifelse(dists <= scoretobeat, 1 / new.denom, 0)
#       }
#       else
#       {
#         newdists <- ifelse(dists < scoretobeat, 1 / max.set.size, 0)
#       }
#       
#     }
#     names(newdists) <- NULL
#     return(newdists)
#     
#   }
#   
#   # scores <- mapply(FUN = handle_set, 
#   #                  sub.list = mahal.nested.list, 
#   #                  idx = 1:length(msets), 
#   #                  MoreArgs = list(max.set.size = max.size),
#   #                  SIMPLIFY = FALSE)
#   
#   #scores <- handle_set(sub.list = mahal.nested.list[[1]], max.set.size = max.size)
#   scores <- handle_set(sub.list = mahal.nested.list, max.set.size = max.size)
#   # for(i in 1:length(msets))
#   # {
#   #   names(scores[[i]]) <- msets[[i]]
#   #   attr(msets[[i]], "weights") <- scores[[i]]
#   # }
#   attr(msets, "weights") <- scores
#   if(verbose) #in future versions, avoid doing the same calculations twice
#   {
#     warning("verbose mode temporarily deprecated")
#     # handle_set_verbose <- function(sub.list)
#     # {
#     #   results.temp <- lapply(sub.list, do.calcs)
#     #   dists <- colMeans(do.call(rbind, results.temp))
#     #   names(dists) <- NULL
#     #   return(dists)
#     # }
#     # 
#     # full.scores <- mapply(FUN = handle_set_verbose, sub.list = mahal.nested.list, SIMPLIFY = FALSE)
#     # 
#     # for(i in 1:length(msets))
#     # {
#     #   names(full.scores[[i]]) <- msets[[i]]
#     #   attr(msets[[i]], "distances") <- full.scores[[i]]
#     # }
#   }
#   
#   return(msets)
# }

# Each of the following similarly named functions carry out the refinement procedures -- using either propensity scores or mahalanobis distances according to the paper. Each of these functions
# will return a matched.set object with the appropriate weights for control units assigned and new, additional attributes (such as refinement method).
# These functions could use some work to be better optimized. For instance, when the "verbose" argument is set to true, they will essentially do all of the refinement calculations twice.
handle_mahalanobis_calculations <- function(mahal.nested.list, msets, max.size, verbose, use.diagonal.covmat)
{
  do.calcs <- function(year.df)
  {
    if(nrow(year.df) == 2)
    {
      return(1)
    }
    cov.data <- year.df[1:(nrow(year.df) - 1), 4:ncol(year.df), drop = FALSE]
    if(use.diagonal.covmat)
    {
      cov.matrix <- diag(apply(cov.data, 2, var), ncol(cov.data), ncol(cov.data))
    } else
    {
      cov.matrix <- cov(cov.data)
    }
    
    center.data <- year.df[nrow(year.df), 4:ncol(year.df), drop = FALSE]
    if(isTRUE(all.equal(det(cov.matrix), 0, tolerance = .00001))) #might not be the conditions we want precisely
    {
      cols.to.remove <- which(apply(cov.data, 2, function(x) isTRUE(length(unique(x)) == 1))) #checking for columns that only have one value
      cols.to.remove <- unique(c(cols.to.remove, which(!colnames(cov.data) %in% colnames(t(unique(t(cov.data))))))) #removing columns that are identical to another column
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
      #(mahalanobis(x = cov.data, center = center.data, cov = cov.matrix))
    })
    
    return(result)
    
    
  }
  handle_set <- function(sub.list, max.set.size, idx)
  {
    #browser()
    results.temp <- lapply(sub.list, do.calcs)
    tmat <- do.call(rbind, results.temp)
    colnames(tmat) <- NULL
    dists <- colMeans(tmat)
    #n.dists <- dists[dists > 0]
    n.dists <- dists
    #if(length(n.dists) == 0)
    #{
    #  w <- 1 / length(dists)
    #  newdists <- dists
    #  newdists <- rep(w, length(newdists))
    #}
    #if(length(n.dists) == 0 ) browser()#stop("a matched set contain only identical units. Please examine the data and remove this set.")
    #else 
    
    if(length(n.dists) < max.set.size) #case where total number of units in matched set < max.set size
    {
      w <- 1 / length(n.dists)
      newdists <- dists
      newdists <- rep(w, length(newdists))
      #newdists[newdists > 0 ] <- w
    }
    else
    {
      ordered.dists <- sort(n.dists)
      scoretobeat <- max(utils::head(ordered.dists, n = max.set.size + 1))
      # might have situation where the Mth largest distance is the same as the Mth - 1 distance. This means that we either choose to leave out both and have a matched set smaller than the max,
      # or include both of them and relax the size of our maximum set size
      # if(sum(dists < scoretobeat & dists > 0) < max.set.size) #change this if we want to be more strict about max.set.size enforcements
      # {
      #   new.denom <- sum(dists <= scoretobeat & dists > 0)
      #   newdists <- ifelse(dists <= scoretobeat & dists > 0, 1 / new.denom, 0)
      # }
      # else
      # {
      #   newdists <- ifelse(dists < scoretobeat & dists > 0, 1 / max.set.size, 0)
      # }
      
      if(sum(dists < scoretobeat) < max.set.size) #change this if we want to be more strict about max.set.size enforcements
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
    if(sum(vec.ratio) == 0)
    {
      wts <- rep(1 / length(control.ps.set), length(control.ps.set))
    }
    else
    {
      wts <- ( vec.ratio ) / sum( vec.ratio )
    }
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
