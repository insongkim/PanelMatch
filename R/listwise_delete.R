lwd_units <- function(full.local.data, unit.id)
{
  ld <- full.local.data[, c(1, 4:ncol(full.local.data)), drop = F]
  idx <- unlist(by(ld, as.factor(ld[, unit.id]), FUN = function(x) any(!complete.cases(x))))
  units.to.remove <- as.numeric(names(idx)[idx])
  sub.data <- subset(full.local.data, !full.local.data[, unit.id] %in% units.to.remove)
  return(as.matrix(sub.data))
}

#global.data needs to be fully prepped/parsed data set that is internally balanced, full of NAs likely
lwd_refinement <- function(msets, global.data, treated.ts, 
                           treated.ids, lag, time.id, unit.id, lead, refinement.method, treatment, size.match,
                           match.missing, covs.formula, verbose, outcome.var, e.sets,use.diag.covmat)
{
  #print(refinement.method)
  #extract other attributes about the msets object, attach to individual mset for consistency
  
  terms <- attr(terms(covs.formula),"term.labels")
  terms <- gsub(" ", "", terms) #remove whitespace
  lag.calls <- terms[grepl("lag(*)", terms)] #regex to get calls to lag function
  if(any(grepl("=", lag.calls))) stop("fix lag calls to use only unnamed arguments in the correct positions")
  lag.calls <- gsub(pattern = "\"", replacement = "", lag.calls, fixed = T)
  lag.calls <- gsub(pattern = "\'", replacement = "", lag.calls, fixed = T)
  lag.calls <- gsub(pattern = "lag(", replacement = "", lag.calls, fixed = T)
  lag.vars <- unlist(strsplit(lag.calls, split = ","))[c(T,F)]
  if(length(lag.calls) > 0)
  {
    lag.nums <- unlist(strsplit(lag.calls, split = ","))[c(F,T)]
    lag.nums <- as.numeric(gsub(")", "", unlist(strsplit(lag.nums, split = ":"))[c(F,T)]))
    max.lag <- max(lag.nums)
    stm <- max.lag + 1
    global.data <- do.call(rbind, by(global.data, as.factor(global.data[, unit.id]), function(x) x[stm:nrow(x), ] ))
    rownames(global.data) <- NULL
    global.data <- as.matrix(global.data[order(global.data[,unit.id], global.data[, time.id]), ] )
  }  

  
  
  if(length(msets) == 0) stop("There are no matched sets!")
  new.msets <- list()
  for(i in 1:length(msets))
  {
    #probably some check on non-empty sets or something
    if(length(msets[[i]]) > 0)
    {
      time <- treated.ts[i]
      uid <- treated.ids[i]
      # if(refinement.method == "ps.msm.weight" | refinement.method == "CBPS.msm.weight")
      # {
      #   localdata <- global.data[ global.data[, time.id]  %in% ( (time - lag):(time + max(lead) )), ]
      #   localdata <- lwd_units(localdata, unit.id)
      # } else
      # {
      #   localdata <- global.data[ global.data[, time.id]  %in% ((time - lag):time), ]
      #   localdata <- lwd_units(localdata, unit.id)  
      # }
      localdata <- global.data[ global.data[, time.id]  %in% ( (time - lag):(time + max(lead) )), ]
      localdata <- lwd_units(localdata, unit.id)
      viable.units <- unique(localdata[, unit.id])
      if(uid %in% viable.units)
      {
        mset <- msets[i]
        #print(mset)
        controls <- mset[[1]]
        controls <- controls[controls %in% viable.units]
        mset[[1]] <- controls
        #mset <- mset[mset %in% viable.units]
        if(length(mset[[1]]) > 0 ) # do something else if we delete everything because of listwise deletion
        {
          #print(mset)
          if(refinement.method == "mahalanobis")
          {
            tset <- set_lwd_refinement(mset, localdata, time, uid, lag, refinement.method, lead, 
                                       verbose, size.match, unit.id, time.id, covs.formula, match.missing, treatment,
                                       use.diag.covmat = use.diag.covmat)  
            new.msets[[i]] <- tset
          } else
          {
            new.msets[[i]] <- mset
          }
          
          #maybe add things/attributes? idk
           
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
  
  t.newsets <- unlist(new.msets, recursive = F)
  idx <- sapply(t.newsets, function(x) !any(is.na(x)))
  t.newsets <- t.newsets[idx]
  if(length(t.newsets) == 0) stop("There are no matched sets!")
  treated.ts <- as.numeric(unlist(strsplit(names(t.newsets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(t.newsets), split = "[.]"))[c(T,F)])
  if(refinement.method != "mahalanobis")
  {
    t.newsets <- set_lwd_refinement(t.newsets, global.data, treated.ts, treated.ids, lag, refinement.method, lead, 
                                              verbose, size.match, unit.id, time.id, covs.formula, match.missing, treatment,
                                    use.diag.covmat = use.diag.covmat)
  }
  
  class(e.sets) <- 'list'
  
  treated.ts <- as.numeric(unlist(strsplit(names(e.sets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(e.sets), split = "[.]"))[c(T,F)])
  
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

set_lwd_refinement <- function(mset, local.data, time, id, 
                               lag, refinement.method, lead, verbose, size.match, unit.id, time.id, covs.formula, match.missing,
                               treatment, use.diag.covmat)
{
  treated.ts <- time
  treated.ids <- id
  ordered.data <- local.data 
  msets <- mset
  if(refinement.method == "mahalanobis")
  {
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, idx =  idxlist)
    msets <- handle_mahalanobis_calculations(mahalmats, msets, size.match, verbose, use.diagonal.covmat = use.diag.covmat)
  }
  if(refinement.method == "ps.msm.weight" | refinement.method == "CBPS.msm.weight")
  {
    
    store.msm.data <- list()
    for(i in 1:length(lead))
    {
      f <- lead[i]
      tf <- expand.treated.ts(lag, treated.ts = treated.ts + f)
      tf.index <- get_yearly_dmats(ordered.data, treated.ids, tf, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
      expanded.sets.tf <- build_ps_data(tf.index, ordered.data, lag)
      #pre.pooled <- ordered.data[ordered.data[, time.id] %in% (treated.ts + f), ]
      pre.pooled <- rbindlist(expanded.sets.tf)
      pooled <- pre.pooled[complete.cases(pre.pooled), ]
      pooled <- unique(as.data.frame(pooled))
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
      if(qr(pooled)$rank != ncol(pooled)) 
      {
        # print("Data used to generate propensity scores is not linearly independent. Calculations cannot be completed.
        #       Would you like to save the problematic matrix to file for manual inspection? File and variable will be saved as 'problematic_matrix.rda'. ")
        # inkey <- readline("Press 'y' to save and any other key to do nothing: ")
        # if(inkey == "y")
        # {
        #   problematic_matrix <- pooled
        #   save(problematic_matrix, file = "problematic_matrix.rda")
        #   stop("PanelMatch terminated")
        # }
        # else
        # {
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
  }  #not msm
  if(all(refinement.method %in% c("CBPS.weight", "CBPS.match", "ps.weight", "ps.match")))
  {
    if(!all(refinement.method %in% c("CBPS.weight", "CBPS.match", "ps.weight", "ps.match"))) stop("please choose valid refinement method")
    tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                         ordered.data[, time.id]), matched_sets = msets, lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- rbindlist(expanded.sets.t0)
    pooled <- unique(pre.pooled[complete.cases(pre.pooled), ])
    
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
    if(qr(pooled)$rank != ncol(pooled)) 
    {
      # print("Data used to generate propensity scores is not linearly independent. Calculations cannot be completed.
      #       Would you like to save the problematic matrix to file for manual inspection? File and variable will be saved as 'problematic_matrix.rda'. ")
      # inkey <- readline("Press 'y' to save and any other key to do nothing: ")
      # if(inkey == "y")
      # {
      #   problematic_matrix <- pooled
      #   save(problematic_matrix, file = "problematic_matrix.rda")
      #   stop("PanelMatch terminated")
      # }
      # else
      # {
      stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
      #}
      
    }
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match")
    {
      dummy <-capture.output(fit0 <- (CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]), 
                                          family = binomial(link = "logit"), data = pooled)))
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
  t.attributes <- attributes(msets)[names(attributes(msets)) != "names"]
  # msets <- c(msets, e.sets)
  for(idx in names(t.attributes))
  {
    attr(msets, idx) <- t.attributes[[idx]]
  }
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "match.missing") <- match.missing
  return(msets)
}
