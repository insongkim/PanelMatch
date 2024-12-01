####################################################################################
####################################################################################
####### Performing listwise deletion requires specialized refinement procedures.
####### This file contains helper functions for carrying out those procedures.
####### They are largely modified/adapted versions of other refinement code
####################################################################################
####################################################################################

#' lwd_refinement
#' master function that performs refinement with listwise deletion = TRUE
#'
#' @param msets 
#' @param global.data data.frame. needs to be fully prepped/parsed data set that is internally balanced, full of NAs likely
#' @param treated.ts vector of the times of treatment for treated observations
#' @param treated.ids vector of unit identifiers of treated observations
#' @param lag 
#' @param time.id string specifying
#' @param unit.id 
#' @param lead vector of lead values
#' @param refinement.method string specifying refinement method
#' @param treatment string specifying treatment variable
#' @param size.match maximum number of units to give non-zero weight to when using matching refinement method
#' @param match.missing logical. indicates whether or not to allow the package to match units on missingness in treatment history
#' @param covs.formula see PanelMatch documentation for descriptions
#' @param verbose see PanelMatch documentation for descriptions
#' @param outcome.var string specifying outcome variable
#' @param e.sets empty sets (treated observations with no matched controls)
#' @param use.diag.covmat see PanelMatch documentation for descriptions
#'
#' @return matched.set object with refined matched sets.
#' @keywords internal
lwd_refinement <- function(msets, global.data, treated.ts, 
                           treated.ids, lag, time.id, unit.id, 
                           lead, refinement.method, treatment, size.match,
                           match.missing, covs.formula, verbose, 
                           outcome.var, e.sets,use.diag.covmat)
{
  #extract other attributes about the msets object, attach to individual mset for consistency
  
  terms <- attr(terms(covs.formula),"term.labels")
  terms <- gsub(" ", "", terms) #remove whitespace
  lag.calls <- terms[grepl("lag(*)", terms)] #regex to get calls to lag function
  if(any(grepl("=", lag.calls))) stop("fix lag calls to use only unnamed arguments in the correct positions")
  lag.calls <- gsub(pattern = "\"", replacement = "", 
                    lag.calls, fixed = TRUE)
  lag.calls <- gsub(pattern = "\'", replacement = "", 
                    lag.calls, fixed = TRUE)
  lag.calls <- gsub(pattern = "lag(", replacement = "", 
                    lag.calls, fixed = TRUE)
  lag.vars <- unlist(strsplit(lag.calls, 
                              split = ","))[c(TRUE,FALSE)]
  if(length(lag.calls) > 0)
  {
    lag.nums <- unlist(strsplit(lag.calls, 
                                split = ","))[c(FALSE,TRUE)]
    lag.nums <- as.numeric(gsub(")", "", 
                                unlist(strsplit(lag.nums, 
                                                split = ":"))[c(FALSE,TRUE)]))
    max.lag <- max(lag.nums)
    stm <- max.lag + 1
    global.data <- do.call(rbind, 
                           by(global.data, as.factor(global.data[, unit.id]), 
                              function(x) x[stm:nrow(x), ] ))
    rownames(global.data) <- NULL
    global.data <- as.matrix(global.data[order(global.data[,unit.id], 
                                               global.data[, time.id]), ] )
  }  
  if(length(msets) == 0) stop("There are no matched sets!")
  new.msets <- list()
  for(i in 1:length(msets))
  {
    
    if(length(msets[[i]]) > 0)
    {
      time <- treated.ts[i]
      uid <- treated.ids[i]
      localdata <- global.data[ global.data[, time.id]  %in% ( (time - lag):(time + max(lead) )), ]
      localdata <- lwd_units(localdata, unit.id)
      viable.units <- unique(localdata[, unit.id])
      if(uid %in% viable.units)
      {
        mset <- msets[i]
        controls <- mset[[1]]
        controls <- controls[controls %in% viable.units]
        mset[[1]] <- controls
        
        if(length(mset[[1]]) > 0 ) # do something else if we delete everything because of listwise deletion
        {
          if(refinement.method == "mahalanobis")
          {
            tset <- set_lwd_refinement(mset, localdata, time, uid, lag, 
                                       refinement.method, lead, 
                                       verbose, size.match, unit.id, time.id, 
                                       covs.formula, match.missing, treatment,
                                       use.diag.covmat = use.diag.covmat)  
            new.msets[[i]] <- tset
          } else
          {
            new.msets[[i]] <- mset
          }
           
        } else
        {
          new.msets[[i]] <- NA
        }
      } else #case where treated unit does not have complete data
      {
        new.msets[[i]] <- NA
      }
    }
    
  }
  t.newsets <- unlist(new.msets, recursive = FALSE)
  idx <- sapply(t.newsets, function(x) !any(is.na(x)))
  t.newsets <- t.newsets[idx]
  if(length(t.newsets) == 0) stop("There are no matched sets!")
  treated.ts <- as.numeric(sub(".*\\.", "", names(t.newsets)))
  treated.ids <- as.numeric(sub("\\..*", "", names(t.newsets)))
  if(refinement.method != "mahalanobis")
  {
    t.newsets <- set_lwd_refinement(t.newsets, global.data, 
                                    treated.ts, treated.ids, 
                                    lag, refinement.method, lead, 
                                    verbose, size.match, unit.id, 
                                    time.id, covs.formula, 
                                    match.missing, treatment,
                                    use.diag.covmat = use.diag.covmat)
  }
  
  class(e.sets) <- 'list'
  treated.ts <- as.numeric(sub(".*\\.", "", names(e.sets)))
  treated.ids <- as.numeric(sub("\\..*", "", names(e.sets)))
  
  esetlist <- logical(length(e.sets))
  if(length(e.sets) > 0)
  {
    for(i in 1:length(e.sets))
    {
      
      time <- treated.ts[i]
      uid <- treated.ids[i]
      
      localdata <- global.data[ global.data[, time.id]  %in% ((time - lag):time), ]
      localdata <- lwd_units(localdata, unit.id)
      viable.units <- unique(localdata[, unit.id])
      
      if(uid %in% viable.units)
      {
        esetlist[i] <- TRUE
      }
    }
    e.sets <- e.sets[esetlist]
    t.newsets <- c(t.newsets, e.sets)
  }
  
  attrib <- names(attributes(msets))[names(attributes(msets)) != "names"]
  for(tatt in attrib)
  {
    attr(t.newsets, tatt) <- attr(msets, tatt)
  }
  attr(t.newsets, "refinement.method") <- refinement.method
  attr(t.newsets, 'max.match.size') <- size.match
  attr(t.newsets, "covs.formula") <- covs.formula
  attr(t.newsets, "match.missing") <- match.missing
  return(t.newsets)
}

#' set_lwd_refinement
#' Performs the set-level operations for refinement with listwise deletion. See documentation for lwd_refinement for descriptions of most parameters.
#' @param mset individual matched set
#' @param local.data data.frame containing the data relevant for set level refinement
#' @param time time of treated observation
#' @param id id of treated observation
#' @return an individual matched set
#' @keywords internal
set_lwd_refinement <- function(mset, local.data, time, id, 
                               lag, refinement.method, lead, verbose, 
                               size.match, unit.id, time.id, 
                               covs.formula, match.missing,
                               treatment, use.diag.covmat)
{
  
  treated.ts <- time
  treated.ids <- id
  ordered.data <- as.matrix(local.data)
  msets <- mset
  if(refinement.method == "mahalanobis")
  {
    old.lag <- lag
    lag <- 0
    tlist <- expand_treated_ts(lag, treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, 
                                treated.ids, 
                                tlist, 
                                msets, 
                                lag)
    mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, 
                                 idx = idxlist)
    msets <- handle_mahalanobis_calculations(mahalmats, 
                                             msets, 
                                             size.match, 
                                             verbose, 
                                             use.diag.covmat)
    lag <- old.lag
  }
  if(refinement.method == "ps.msm.weight" | refinement.method == "CBPS.msm.weight")
  {
    
    store.msm.data <- list()
    for(i in 1:length(lead))
    {
      f <- lead[i]
      tf <- expand_treated_ts(lag, treated.ts = treated.ts + f)
      tf.index <- get_yearly_dmats(ordered.data, treated.ids, tf, matched_sets = msets, lag)
      expanded.sets.tf <- build_ps_data(tf.index, ordered.data, lag)
  
      pre.pooled <- rbindlist(expanded.sets.tf)
      pooled <- pre.pooled[complete.cases(pre.pooled), ]
      pooled <- unique(as.data.frame(pooled))
      cols.to.remove <- which(unlist(lapply(pooled, function(x){all(x[1] == x)}))) #checking for columns that only have one value
      cols.to.remove <- unique(c(cols.to.remove, which(!colnames(pooled) %in% colnames(t(unique(t(pooled))))))) #removing columns that are identical to another column 
      cols.to.remove <- cols.to.remove[cols.to.remove > 3] #leave the first three columns alone
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
      if(qr(pooled)$rank != ncol(pooled)) 
      {
        stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
        #}
        
      }
      if(refinement.method == "CBPS.msm.weight") #obviously update these conditionals
      {
        dummy <- capture.output(fit.tf <- (CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                                                      family = binomial(link = "logit"), data = pooled)))
      }
      if(refinement.method == "ps.msm.weight")
      {
        fit.tf <- glm(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                      family = binomial(link = "logit"), data = pooled)
      }
      store.msm.data[[i]] <- find_ps(expanded.sets.tf, fit.tf)
      
    }
    
    msm.sets <- gather_msm_sets(store.msm.data)
    #can only have weighting in these situations
    msets <- handle_ps_weighted(msm.sets, msets, refinement.method)
  } 
  if(all(refinement.method %in% c("CBPS.weight", "CBPS.match", 
                                  "ps.weight", "ps.match")))
  {
    if(!all(refinement.method %in% c("CBPS.weight", "CBPS.match", 
                                     "ps.weight", "ps.match"))) stop("please choose valid refinement method")
    tlist <- expand_treated_ts(lag, treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, 
                                treated.ids, 
                                tlist, 
                                msets, 
                                lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- data.table::rbindlist(expanded.sets.t0)
    pooled <- unique(pre.pooled[complete.cases(pre.pooled), ])
    
    cols.to.remove <- which(unlist(lapply(pooled, 
                                          function(x){all(x[1] == x)}))) #checking for columns that only have one value
    cols.to.remove <- unique(c(cols.to.remove, 
                               which(!colnames(pooled) %in% colnames(t(unique(t(pooled))))))) #removing columns that are identical to another column 
    cols.to.remove <- cols.to.remove[cols.to.remove > 3] #leave the first three columns alone
    if(length(cols.to.remove) > 0)
    {
      class(pooled) <- c("data.frame")
      pooled <- pooled[, -cols.to.remove]
      rmv <- function(x, cols.to.remove_)
      {
        return(x[, -cols.to.remove_])
      }
      expanded.sets.t0 <- lapply(expanded.sets.t0, rmv, 
                                 cols.to.remove_ = cols.to.remove)
    }
    if(qr(pooled)$rank != ncol(pooled)) 
    {
      stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
      
    }
    if(refinement.method == "CBPS.weight" || refinement.method == "CBPS.match")
    {
      dummy <-capture.output(fit0 <- (CBPS::CBPS(reformulate(response = treatment, 
                                                             termlabels = colnames(pooled)[-c(1:3)]), 
                                          family = binomial(link = "logit"), data = pooled)))
    }
    if(refinement.method == "ps.weight" || refinement.method == "ps.match")
    {
      fit0 <- glm(reformulate(response = treatment, 
                              termlabels = colnames(pooled)[-c(1:3)]), 
                  family = binomial(link = "logit"), data = pooled)
    }
    
    just.ps.sets <- find_ps(expanded.sets.t0, fit0)
    
    if(refinement.method == "CBPS.weight" || refinement.method == "ps.weight")
    {
      msets <- handle_ps_weighted(just.ps.sets, msets, refinement.method)
    }
    if(refinement.method == "CBPS.match" || refinement.method == "ps.match")
    {
      msets <- handle_ps_match(just.ps.sets, msets, 
                               refinement.method, verbose, size.match)
      attr(msets, "max.match.size") <- size.match
    }
  }
  t.attributes <- attributes(msets)[names(attributes(msets)) != "names"]
  
  for(idx in names(t.attributes))
  {
    attr(msets, idx) <- t.attributes[[idx]]
  }
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "match.missing") <- match.missing
  return(msets)
}

#' lwd_units
#' helper function that actually subsets sets down to contain units with complete data
#' @param full.local.data data.frame containing the data to be used in set-level refinement, but containing missing data
#' @param unit.id 
#' @return data.frame with the missing data removed to be used for set-level refinement.
#' @keywords internal
lwd_units <- function(full.local.data, unit.id)
{
  # can assume the structure of the data such that columns 4 and higher are relevant covariate data
  ld <- full.local.data[, c(1, 4:ncol(full.local.data)), drop = FALSE]
  idx <- unlist(by(ld, 
                   as.factor(ld[, unit.id]), 
                   FUN = function(x) any(!complete.cases(x))))
  units.to.remove <- as.numeric(names(idx)[idx])
  sub.data <- subset(full.local.data, 
                     !full.local.data[, unit.id] %in% units.to.remove)
  return(as.matrix(sub.data))
}