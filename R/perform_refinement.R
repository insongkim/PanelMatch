# Performs refinement on matched sets
# The function mostly calls lower level functions that apply the specified 
# refinement processes
# Returns refined matched sets
perform_refinement <- function(lag, time.id, unit.id, treatment, 
                               refinement.method, size.match,
                               ordered.data, match.missing, covs.formula, 
                               verbose,
                               mset.object = NULL, lead, outcome.var = NULL,
                               forbid.treatment.reversal = FALSE, qoi = "",
                               matching = TRUE, exact.matching.variables = NULL,
                               listwise.deletion,
                               use.diag.covmat = FALSE, 
                               placebo.test = FALSE,
                               restrict.control.period = NULL,
                               caliper.formula = NULL,
                               continuous.treatment.info = NULL)
{

  if (inherits(ordered.data[, unit.id], "numeric"))
  {
    warning("converting unit id variable data to integer")
    class(ordered.data[, unit.id]) <- "integer"
  }
  if (is.null(continuous.treatment.info))
  {
    temp.treateds <- findBinaryTreated(ordered.data, qoi.in = qoi,
                                       treatedvar = treatment, 
                                       time.var = time.id,
                                       unit.var = unit.id, 
                                       hasbeensorted = TRUE)
  } else #continuous
  {
    temp.treateds <- findContinuousTreated(dmat = ordered.data, 
                                           treatedvar = treatment, 
                                           time.var = time.id,
                                           unit.var = unit.id, qoi = qoi,
                                           continuous.treatment.info = continuous.treatment.info)
    
    
  }
  

  idx <- !((temp.treateds[, time.id] - lag) < min(ordered.data[, time.id]))
  temp.treateds <- temp.treateds[idx, ]
  if (nrow(temp.treateds) == 0)
  {
    warn.str <- paste0("no viable treated units for ", qoi, " specification")
    stop(warn.str)
  }
  msets <- get.matchedsets(temp.treateds[, time.id],
                             temp.treateds[, unit.id],
                             data = ordered.data,
                             L = lag,
                             t.column = time.id,
                             id.column = unit.id,
                             treatedvar = treatment,
                             hasbeensorted = TRUE,
                             match.on.missingness = match.missing,
                             matching = TRUE,
                             qoi.in = qoi,
                             restrict.control.period = restrict.control.period,
                             continuous.treatment.info = continuous.treatment.info)
    

  
  e.sets <- msets[sapply(msets, length) == 0]
  msets <- msets[sapply(msets, length) > 0 ]
  if (length(msets) == 0)
  {
    t.attributes <- attributes(e.sets)[names(attributes(e.sets)) != "names"]
    msets <- e.sets
    for (idx in names(t.attributes))
    {
      attr(msets, idx) <- t.attributes[[idx]]
    }
    attr(msets, "covs.formula") <- covs.formula
    attr(msets, "match.missing") <- match.missing
    return(msets)
  }
  
  if (length(msets) > 0)
  {
    msets <- identifyDirectionalChanges(msets, ordered.data,
                                        unit.id, time.id, treatment, qoi)
    # add attribute to track treatment levels in treated observations at time t-1.
    msets <- extract.baseline.treatment(matched.sets = msets,
                                        data.in = ordered.data,
                                        id.variable = unit.id,
                                        time.variable = time.id,
                                        treatment.variable = treatment)
  }
  if (length(e.sets) > 0)
  {
    e.sets <- identifyDirectionalChanges(e.sets, ordered.data,
                                         unit.id, time.id, treatment, qoi)
    e.sets <- extract.baseline.treatment(matched.sets = e.sets,
                                         data.in = ordered.data,
                                         id.variable = unit.id,
                                         time.variable = time.id,
                                         treatment.variable = treatment)
  }
  
  if (!is.null(continuous.treatment.info))
  {
    setlist <- filterContinuousTreated(msets, e.sets, qoi, 
                                       continuous.treatment.info[["control.threshold"]])
    msets <- setlist[["sets"]]
    e.sets <- setlist[["empty.sets"]]
  }
  

  msets <- clean_leads(msets, 
                       ordered.data, 
                       max(lead), 
                       time.id, 
                       unit.id, 
                       outcome.var)


  
  if (placebo.test)
  {
    
    treated.ts <- as.integer(sub(".*\\.", "", names(msets)))
    treated.ids <- as.integer(sub("\\..*", "", names(msets)))

    rownames(ordered.data) <- paste0(ordered.data[, unit.id],
                                     ".",
                                     ordered.data[, time.id])
    msets2 <- filter_placebo_results(expanded_data = as.matrix(ordered.data[, c(unit.id, time.id)]),
                                     ordered_outcome_data = ordered.data[, outcome.var],
                                     treated_ids = treated.ids,
                                     treated_ts = treated.ts,
                                     sets = msets,
                                     lag = lag)
    
    
    msets2 <- msets2[sapply(msets2, length) > 0 ]
    
    attr(msets2, "lag") <- lag
    attr(msets2, "t.var") <- time.id
    attr(msets2, "id.var" ) <- unit.id
    attr(msets2, "treatment.var") <- treatment
    class(msets2) <- "matched.set"
    
  }

  if (forbid.treatment.reversal)
  {
    msets <- enforce_lead_restrictions(msets, 
                                       ordered.data, 
                                       max(lead), 
                                       time.id, 
                                       unit.id, 
                                       treatment.var = treatment)
  }
  if (length(msets) == 0)
  {
    warn.str <- paste0("no matched sets for ", qoi, " specification")
    stop(warn.str)
  }

  if (!is.null(exact.matching.variables))
  {
    msets <- do_exact_matching(msets, 
                               ordered.data, 
                               exact.matching.variables)

    e.sets <- c(e.sets, msets[sapply(msets, length) == 0])
    msets <- msets[sapply(msets, length) > 0 ]
  }


  

  
  
  ####apply calipers here
  if(!is.null(caliper.formula))
  {
    msets <- handle_calipers(plain.ordered.data = ordered.data, caliper.formula,
                             matched.sets = msets, lag.window = 0:lag)
    
  }
  

  
  if(refinement.method == "none")
  {
    for(i in 1:length(msets))
    {
      attr(msets[[i]], "weights") <- rep(1/length(msets[[i]]),
                                         length(msets[[i]]))
      names(attr(msets[[i]], "weights")) <- msets[[i]]
    }
    attr(msets, "refinement.method") <- refinement.method
    t.attributes <- attributes(msets)[names(attributes(msets)) != "names"]
    msets <- c(msets, e.sets)
    for(idx in names(t.attributes))
    {
      attr(msets, idx) <- t.attributes[[idx]]
    }
    attr(msets, "covs.formula") <- covs.formula
    attr(msets, "match.missing") <- match.missing
    attr(msets, "caliper.formula") <- caliper.formula
    return(msets)
  }

  treated.ts <- as.integer(sub(".*\\.", "", names(msets)))
  treated.ids <- as.integer(sub("\\..*", "", names(msets)))


  ordered.data <- parse_and_prep(formula = covs.formula, 
                                 data = ordered.data)
  
  
  if (!is.null(continuous.treatment.info))
  { #make the treatment variable binary in case PS based method is used
    
    idx <- paste0(ordered.data[, unit.id], ".", 
                  ordered.data[, time.id]) %in% paste0(treated.ids, ".", treated.ts)
    ordered.data[, treatment] <- ifelse(idx, 1, 0)
  }
  
  if (any(c("character", "factor") %in% sapply(ordered.data, class)))
  {
    stop("please convert covs.formula variables to numerical data")
  }
  if(any(apply(ordered.data, 2, FUN = function(x) any(is.infinite(x)))))
  {
    stop("Data needed for refinement contains infinite values. Code cannot proceed!")
  }
  ################################################################################################
  if(listwise.deletion) #code will just return from here when listwise.deletion = TRUE
  {
    msets <- lwd_refinement(msets, ordered.data, 
                            treated.ts, treated.ids, lag,
                            time.id, unit.id, 
                            lead, refinement.method,
                            treatment, size.match,
                            match.missing, covs.formula,
                            verbose, outcome.var, e.sets,
                            use.diag.covmat = use.diag.covmat)
    attr(msets, "covs.formula") <- covs.formula
    attr(msets, "match.missing") <- match.missing
    attr(msets, "caliper.formula") <- caliper.formula
    return(msets)
  }
  ################################################################################################
  if(!listwise.deletion)
  {
    ordered.data <- as.matrix(handle.missing.data(ordered.data, 
                                                  4:ncol(ordered.data)))
  }

  if (refinement.method == "mahalanobis")
  {

    old.lag <- lag
    lag <- 0
    tlist <- expand.treated.ts(lag, treated.ts)

    idxlist <- get_yearly_dmats(ordered.data, 
                                treated.ids,
                                tlist, 
                                msets, 
                                lag)

    mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, 
                                 idx =  idxlist)

    msets <- handle_mahalanobis_calculations(mahalmats, msets, 
                                             size.match, verbose, 
                                             use.diag.covmat)

    lag <- old.lag
  }
  if(all(refinement.method %in% c("CBPS.weight", "CBPS.match", 
                                  "ps.weight", "ps.match")))
  {

    if(!all(refinement.method %in% c("CBPS.weight", "CBPS.match", 
                                     "ps.weight", "ps.match")))
    {
      stop("please choose valid refinement method")
    }
      
      
    tlist <- expand.treated.ts(lag, treated.ts)
    idxlist <- get_yearly_dmats(ordered.data, 
                                treated.ids, 
                                tlist, 
                                msets, lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- data.table::rbindlist(expanded.sets.t0)
    pooled <- unique(pre.pooled[complete.cases(pre.pooled), ])
    # Do what we can to minimize computational errors
    cols.to.remove <- which(unlist(lapply(pooled, 
                                          function(x){all(x[1] == x)}))) 
    #checking for columns that only have one value
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
      expanded.sets.t0 <- lapply(expanded.sets.t0, 
                                 rmv, 
                                 cols.to.remove_ = cols.to.remove)
    }
    if(qr(pooled)$rank != ncol(pooled))
    {
      stop("Error: Provided data is not linearly independent so calculations cannot be completed. Please check the data set for any redundant, unnecessary, or problematic information.")
    }
    if(refinement.method == "CBPS.weight" || 
       refinement.method == "CBPS.match")
    {
      dummy <- capture.output(fit0 <- (CBPS::CBPS(reformulate(response = treatment, 
                                                              termlabels = colnames(pooled)[-c(1:3)]),
                                                  family = binomial(link = "logit"), 
                                                  data = pooled)))
    }
    if(refinement.method == "ps.weight" || 
       refinement.method == "ps.match")
    {
      fit0 <- glm(reformulate(response = treatment, 
                              termlabels = colnames(pooled)[-c(1:3)]),
                  family = binomial(link = "logit"), data = pooled)
    }

    just.ps.sets <- find_ps(expanded.sets.t0, fit0)

    if(refinement.method == "CBPS.weight" || 
       refinement.method == "ps.weight")
    {
      msets <- handle_ps_weighted(just.ps.sets, 
                                  msets, 
                                  refinement.method)
    }
    if(refinement.method == "CBPS.match" || 
       refinement.method == "ps.match")
    {
      msets <- handle_ps_match(just.ps.sets, msets, 
                               refinement.method, verbose, size.match)
      attr(msets, "max.match.size") <- size.match
    }
  }
  if(refinement.method == "mahalanobis")
  {
    attr(msets, "max.match.size") <- size.match
  }
  t.attributes <- attributes(msets)[names(attributes(msets)) != "names"]
  msets <- c(msets, e.sets)
  for(idx in names(t.attributes))
  {
    attr(msets, idx) <- t.attributes[[idx]]
  }
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "caliper.formula") <- caliper.formula
  attr(msets, "match.missing") <- match.missing
  return(msets)
}