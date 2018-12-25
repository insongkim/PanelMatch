#' PanelMatch2
#' 
#' \code{PanelMatch2} identifies a matched set for each treated
#' observation. Specifically, for a given treated unit, the matched
#' set consists of control observations that have the identical
#' treatment history up to a certain number of \code{lag}
#' years. Researchers must specify \code{lag}. A further refinement of
#' the matched set will be possible by setting the size of the matched
#' set \code{M}, \code{weighting}, and adjusting for other confounders
#' such as past outcomes and covariates via \code{formula}.
#'
#' @param lag An integer value indicating the length of history to be matched
#' @param time.id A character string indicating the name of time identifier
#' variable in the \code{data}. currently must be integer.
#' @param unit.id A character string indicating the name of unit identifier in the data. Currently must be integer.
#' @param treatment A character string indicating the name of
#' treatment variable in the \code{data}. The treatment should be a
#' binary indicator (integer with 0 for the control group and 1 for
#' the treatment group).
#' @param outcome Character string of outcome variable. 
#' @param refinement.method character string of matching method used for refining the matched sets. Currently only mahalanobis is supported
#' @param match.missing Logical variable indicating whether or not units can be matched on the patterns of missingness in their treatment histories
#' @param data Balanced panel data set
#' @param size.match Maximum size of the matched sets after refinement
#' @param covs.formula One sided formula indicating which variables should be used for matching and refinement.
#' Can specify lags using a function "lag" which takes two positional arguments. The first is the name of the variable which you wish to lag, specified as a string second is the lag window, specified as an integer sequence
#' See the example below. 
#' @param verbose option to include more information about the matched.set object calculations
#' @return \code{PanelMatch} returns a list of class `matched.set'
#' containing the following components:
#' \item{treatment}{treatment variable name}
#'
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples \dontrun{
#' results <- PanelMatch2(4, "year", "wbcode2", "dem", "y", refinement.method = "mahalanobis", 
#'                       data = dem, match.missing = T, covs.formula = ~ lag("tradewb", 1:4) + lag("y", 1:4), size.match = 5, verbose = T)
#' }
#' @export
PanelMatch2 <- function(lag, time.id, unit.id, treatment, outcome,
                       refinement.method = c(NULL, "ps.weight", "ps.match", "mahalanobis", "CBPS.weight", "CBPS.match"),
                       size.match = 10,
                       data,
                       match.missing = TRUE,
                       covs.formula,
                       verbose = FALSE
                       ) 
{
  #browser()
  #refinement.method <- "mahalanobis"
  # set covariates and dependent
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    stop("panel data is not balanced")
  }
  
  #order and ensure order does not get violated
  #convert to data.table?
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment, outcome)]
  data <- data[, c(unit.id, time.id, treatment, outcome, othercols)] #reorder columns 
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  temp.treateds <- findAllTreated(ordered.data, treatedvar = treatment, time.var = time.id, unit.var = unit.id, hasbeensorted = TRUE)
  msets <- get.matchedsets(temp.treateds[, time.id], temp.treateds[, unit.id], ordered.data, lag, time.id, unit.id, treatment, hasbeensorted = TRUE)
  msets <- msets[sapply(msets, length) > 0 ]
  treated.ts <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(T,F)])
  
  ordered.data <- as.matrix(parse_and_prep(formula = covs.formula, data = ordered.data, unit.id = unit.id)) #every column > 4 at this point should be used in distance/refinement calculation
  #browser()
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
    
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match")
    {
      fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:4)]), 
                                          family = binomial(link = "logit"), data = pooled))
    }
    if(refinement.method == "ps.weight" | refinement.method == "ps.score")
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

#' @export
PanelMatch2.matched.set <- function(mset.object, data, outcome.var)
{
  
  #browser()
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
    
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match")
    {
      fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = colnames(pooled)[-c(1:4)]), 
                                          family = binomial(link = "logit"), data = pooled))
    }
    if(refinement.method == "ps.weight" | refinement.method == "ps.score")
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





























