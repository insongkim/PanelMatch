# File contains helper functions written in R for PanelMatch functionality
perform_refinement <- function(lag, time.id, unit.id, treatment, refinement.method, size.match,
                               ordered.data, match.missing, covs.formula, verbose,
                               mset.object = NULL, lead, outcome.var = NULL, forbid.treatment.reversal = FALSE, qoi = "",
                               matching = TRUE, exact.matching.variables = NULL, listwise.deletion,
                               use.diag.covmat = FALSE, caliper.formula = NULL, calipers.in.refinement = FALSE,
                               continuous.treatment = FALSE, continuous.treatment.formula = NULL, restrict.control.period = NULL)
{


  
  if (continuous.treatment)
  {
    temp.treateds <- findContinuousTreated(dmat = ordered.data, treatedvar = treatment, time.var = time.id,
                          unit.var = unit.id)
  } else
  {

    temp.treateds <- findBinaryTreated(ordered.data, treatedvar = treatment, time.var = time.id,
                                       unit.var = unit.id, hasbeensorted = TRUE)

  }

  idx <- !((temp.treateds[, time.id] - lag) < min(ordered.data[, time.id]))
  temp.treateds <- temp.treateds[idx, ]
  if (nrow(temp.treateds) == 0)
  {
    warn.str <- paste0("no viable treated units for ", qoi, " specification")
    stop(warn.str)
  }
  
  msets <- get.matchedsets(temp.treateds[, time.id], temp.treateds[, unit.id], data = ordered.data,
                             L = lag, t.column = time.id, id.column = unit.id,
                             treatedvar = treatment, hasbeensorted = TRUE,
                             match.on.missingness = match.missing, matching = TRUE,
                             continuous = continuous.treatment, continuous.treatment.formula = NULL,
                           restrict.control.period = restrict.control.period)
  e.sets <- msets[sapply(msets, length) == 0]
  msets <- msets[sapply(msets, length) > 0 ]


  msets <- clean_leads(msets, ordered.data, max(lead), time.id, unit.id, outcome.var)

  if(forbid.treatment.reversal)
  {
    msets <- enforce_lead_restrictions(msets, ordered.data, max(lead), time.id, unit.id, treatment.var = treatment)
  }
  if(length(msets) == 0)
  {
    warn.str <- paste0("no matched sets for ", qoi, " specification")
    stop(warn.str)
  }

  print("refinement process can now begin")
  if(!is.null(exact.matching.variables))
  {
    msets <- do_exact_matching(msets, ordered.data, exact.matching.variables)

    e.sets <- c(e.sets, msets[sapply(msets, length) == 0])
    msets <- msets[sapply(msets, length) > 0 ]
  }

  ####apply calipers here
  if(!is.null(caliper.formula))
  {
    msets <- handle_calipers(plain.ordered.data = ordered.data, caliper.formula,
                             matched.sets = msets, lag.window = 0:lag)
    if(calipers.in.refinement)
    {
      covs.formula <- merge_formula(covs.formula, caliper.formula)
    }
  }

  if(refinement.method == "none")
  {
    for(i in 1:length(msets))
    {
      attr(msets[[i]], "weights") <- rep(1/length(msets[[i]]), length(msets[[i]]))
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
    return(msets)
  }

  treated.ts <- as.integer(sub(".*\\.", "", names(msets)))
  treated.ids <- as.integer(sub("\\..*", "", names(msets)))



  print("data reformatting beginning")
  # print("covs.formula")
  # print(covs.formula)
  # print("preview of data:")
  # print(head(ordered.data))
  # print(class(ordered.data))
  
  ordered.data <- parse_and_prep(formula = covs.formula, 
                                 data = ordered.data, 
                                 unit.id = unit.id,
                                 time.id = time.id,
                                 treatment.variable = treatment)
  print("data reformatting complete")
  
  if (any(apply(ordered.data, 2, FUN = function(x) any(is.infinite(x)))))
  {
    stop("Data needed for refinement contains infinite values. Code cannot proceed!")
  }
  ################################################################################################
  if (listwise.deletion) #code will just return from here when listwise.deletion = T
  {
    msets <- lwd_refinement(msets, ordered.data, treated.ts, treated.ids, lag,
                            time.id, unit.id, lead, refinement.method, treatment, size.match,
                            match.missing, covs.formula, verbose, outcome.var, e.sets,
                            use.diag.covmat = use.diag.covmat)
    return(msets)
  }
  ################################################################################################
  if (!listwise.deletion)
  {
    ordered.data <- as.matrix(handle.missing.data(ordered.data, 4:ncol(ordered.data)))
  }


  if (refinement.method == "mahalanobis")
  {

    old.lag <- lag
    lag <- 0
    ## new version
    print("mahalanobis refinement starting")
    
    rownames(ordered.data) <- paste0(ordered.data[, unit.id], ".", ordered.data[, time.id])
    ##now we assume that rownames have been given
    msets <- handle_distance_matrices_maha(ordered_expanded_data = ordered.data,
                                           matched.sets = msets, id.var = unit.id,
                                           time.var = time.id, lag.in = lag, maxSize = size.match,
                                           useDiagonalCovmat = use.diag.covmat, 
                                           verbose.in = verbose, treat.var = treatment)
    
    ###old code
    # tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
    # 
    # idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, matched_sets = msets, lag)
    # 
    # mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, idx =  idxlist)
    # 
    # msets <- handle_mahalanobis_calculations(mahalmats, msets, size.match, 
    #                                          verbose, use.diagonal.covmat = use.diag.covmat)

    lag <- old.lag
    print("Refinement complete!")
  }
  if(refinement.method == "ps.msm.weight" | refinement.method == "CBPS.msm.weight")
  {
    store.msm.data <- list()
    for(i in 1:length(lead))
    {
        f <- lead[i]
        tf <- expand.treated.ts(lag, treated.ts = treated.ts + f)
        tf.index <- get_yearly_dmats(ordered.data, treated.ids, tf, matched_sets = msets, lag)
        expanded.sets.tf <- build_ps_data(tf.index, ordered.data, lag)
        #pre.pooled <- ordered.data[ordered.data[, time.id] %in% (treated.ts + f), ]
        pre.pooled <- rbindlist(expanded.sets.tf)
        pooled <- pre.pooled[complete.cases(pre.pooled), ]
        pooled <- unique(as.data.frame(pooled))
        #do the column removal thing
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
    idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, matched_sets = msets, lag)
    expanded.sets.t0 <- build_ps_data(idxlist, ordered.data, lag)
    pre.pooled <- rbindlist(expanded.sets.t0)
    pooled <- unique(pre.pooled[complete.cases(pre.pooled), ])

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
      dummy <- capture.output(fit0 <- (CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:3)]),
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
  msets <- c(msets, e.sets)
  for(idx in names(t.attributes))
  {
    attr(msets, idx) <- t.attributes[[idx]]
  }
  attr(msets, "covs.formula") <- covs.formula
  attr(msets, "match.missing") <- match.missing
  return(msets)
}

#data has unit, time, treatment, everything else column order at this point
parse_and_prep <- function(formula, data, unit.id, 
                           time.id, treatment.variable)
{
  internal.lag <- function (x, n = 1L, default = NA)
  {
    if (n == 0) return(x)
    xlen <- length(x)
    n <- pmin(n, xlen)
    out <- c(rep(default, n), x[seq_len(xlen - n)])
    attributes(out) <- attributes(x)
    out
  }

  lag <- function(y, lwindow)
  {
    
    sapply(lwindow, internal.lag, x = y)
  }

  apply_formula <- function(x, form)
  {
    attr(form, ".Environment") <- environment()
    # print("class of x:")
    # print(class(x))
    # print("preview of x:")
    # print(head(x))
    # print("dimensions of x:")
    # print(dim(x))
    tdf <- model.frame(form, x, na.action = NULL)
    #tdf2 <- model.frame(form, x, na.action = NULL)
    
    #return(cbind(x[, c(1, 2, 3)], model.matrix(form, tdf)[, -1]))
    
    return(as.matrix(tdf))
    #cbind(x[, c(1, 2)], model.matrix(form, tdf)[, -1])
  }
  ###assuming everything is ordered when it comes in!!!!!!
  #by(data, as.factor(data[, unit.id]), FUN = tfunc, form = formula)
  #original.id <- sub("\\..*", "", unit.id)
  #DT <- data.table(data)
  # DT[,(original.id):=NULL]
  t.data <- do.call(rbind, by(data, as.factor(data[, 1]), FUN = apply_formula, form = formula))
  #idk.dt <- DT[, := apply_formula(x = .SD, form = formula), by = unit.id]
  #may not be necessary?
  #idk.dt[, !original.id]
  #t.data <- as.data.frame(idk.dt)
  
  t.data <- cbind(data[, c(unit.id, time.id, treatment.variable)], t.data)
  #othercols <- colnames(t.data)[!colnames(t.data) %in% c(time.id, unit.id, treatment.variable)]
  #t.data <- t.data[, c(unit.id, time.id, treatment.variable, othercols)]
  #t.data <- t.data[order(t.data[,1], t.data[,2]), ]
  rownames(t.data) <- NULL
  
  return(t.data)
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

#right now this function just checks outcome data and cleans up based on that, but when msm is implemented, we will also need to check reversion of treatment
clean_leads <- function(matched_sets, ordered.data, max.lead, t.var, id.var, outcome.var)
{

  old.attributes <- attributes(matched_sets)[names(attributes(matched_sets)) != "names"]
  print('using regex')
  ts <- as.numeric(sub(".*\\.", "", names(matched_sets)))
  tids <- as.numeric(sub("\\..*", "", names(matched_sets)))
  class(matched_sets) <- "list" #so that Rcpp::List is accurate when we pass it into cpp functions
  print('starting cpp cleaning process')

  sub_sets <- clean_leads_cpp(subset_data = as.matrix(ordered.data[, c(id.var,t.var,outcome.var)]),
                                       sets = matched_sets, treated_tid_pairs = names(matched_sets),
                                       treated_ids = tids, lead =  as.integer(max.lead), times = ts)
  print('cleaning complete')
 
  print('all units checked')
  if(all(sapply(sub_sets, length) == 0)) stop('estimation not possible: none of the matched sets have viable control units due to a lack of necessary data')
  matched_sets <- sub_sets[sapply(sub_sets, length) > 0]
  print("subsetting process complete")
  for(idx in names(old.attributes))
  {
    attr(matched_sets, idx) <- old.attributes[[idx]]
  }
  class(matched_sets) <- c("matched.set")
  print("bookkeeping completed")
  return(matched_sets)


}


enforce_lead_restrictions <- function(matched_sets, ordered.data, max.lead, t.var, id.var, treatment.var)
{
  #CHECK TO MAKE SURE COLUMNS ARE IN ORDER
  ordered.data <- ordered.data[order(ordered.data[,id.var], ordered.data[,t.var]), ]
  compmat <- data.table::dcast(data.table::as.data.table(ordered.data), formula = paste0(id.var, "~", t.var), value.var = treatment.var)
  ts <- as.numeric(sub(".*\\.", "", names(matched_sets)))
  tids <- as.numeric(sub("\\..*", "", names(matched_sets)))
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
      #warning("all controls in a particular matched set were removed due to missing data")
      #browser()
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

merge_formula <- function(form1, form2)
{
  # assuming one sided formulae
  # form1 should be the covs.formula argument
  # get character strings of the right hand sides
  #rhs1 <- strsplit(deparse(form1[[2]]), " \\+ ")[[1]]
  rhs1 <- trimws(unlist(strsplit(as.character(form1)[2], "\\+")))
  rhs2 <- strsplit(deparse(form2[[2]]), " \\+ ")[[1]]

  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)

  # put the two sides together with the amazing
  # reformulate function
  out <- reformulate(rhs)

  # set the environment of the formula (i.e. where should
  # R look for variables when data aren't specified?)
  environment(out) <- environment(form1)

  return(out)
}
