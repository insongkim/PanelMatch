#' PanelMatch
#' 
#' \code{PanelMatch} identifies a matched set for each treated
#' observation. Specifically, for a given treated unit, the matched
#' set consists of control observations that have the identical
#' treatment history up to a certain number of \code{lag}
#' years. Researchers must specify \code{lag}. A further refinement of
#' the matched set will be possible by setting a maximum size of the matched
#' set \code{size.match} (the maximum number of control units that can be matched to a treated unit) and adjusting for other confounders
#' such as past outcomes and covariates via \code{covs.formula} and specifying a method of refinement via \code{refinement.method}.
#'
#' @param lag An integer value indicating the length of history to be matched
#' @param time.id A character string indicating the name of time identifier
#' variable in the \code{data}. This data currently must be integer.
#' @param unit.id A character string indicating the name of unit identifier in the data. This data currently must be integer.
#' @param treatment A character string indicating the name of
#' treatment variable in the \code{data}. The treatment should be a
#' binary indicator (integer with 0 for the control group and 1 for
#' the treatment group).
#' @param outcome Character string of outcome variable. 
#' @param refinement.method character string of matching or weighting method used for refining the matched sets. The user can choose "mahalanobis", "ps.match", "CBPS.match", "ps.weight", "CBPS.weight". The first three methods will use the \code{size.match} argument to create sets of at most \code{size.match} closest control units.
#' @param match.missing Logical variable indicating whether or not units should be matched on the patterns of missingness in their treatment histories
#' @param data A data.frame object containing time series cross sectional data
#' @param size.match Maximum size of the matched sets after refinement
#' @param covs.formula One sided formula indicating which variables should be used for matching and refinement.
#' The user can specify lags using a function "lag" which takes two, unnamed, positional arguments. The first is the name of the variable which you wish to lag, specified as a string. The second is the lag window, specified as an integer sequence
#' See the example below. 
#' @param verbose option to include more information about the matched.set object calculations, like the distances used to create the refined sets and weights
#' @return \code{PanelMatch} returns a list of class `matched.set'. Each element in the list is a vector of integers corresponding to the control unit ids in a matched set. 
#' Additionally, these vectors might have additional attributes -- "weights" or "distances". These correspond to the weights or distances corresponding to each control unit, as determined by the specified refinement method.
#' Each element also has a name, which corresponds to the unit id and time variable of the treated unit and time of treatment, concatenated together and separated by a period.
#'  
#' matched.set objects also have a number of other potential attributes:
#' \item{lag}{same as lag parameter -- an integer value indicating the length of treatment history to be}
#' \item{t.var}{time variable name}
#' \item{id.var}{unit id variable name}
#' \item{treated.var}{treated variable name}
#' \item{class}{class of the object: should always be "matched.set"} 
#' \item{refinement.method}{method used to refine and/or weight the control units in each set.}
#' \item{covs.formula}{see covs.formula argument}
#' \item{match.missing}{see match.missing argument}
#' \item{max.match.size}{same as size.match argument}
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples \dontrun{
#' results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", treatment = "dem", outcome = "y", refinement.method = "mahalanobis", 
#'                       data = dem, match.missing = T, covs.formula = ~ lag("tradewb", 1:4) + lag("y", 1:4), size.match = 5)
#' results[[1]] #to see the control units matched to the first treated units.                     
#' }
#' @export
PanelMatch <- function(lag, time.id, unit.id, treatment, outcome,
                       refinement.method = c(NULL, "ps.weight", "ps.match", "mahalanobis", "CBPS.weight", "CBPS.match"),
                       size.match = 10,
                       data,
                       match.missing = TRUE,
                       covs.formula,
                       verbose = FALSE
                       ) 
{
  if(!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment, outcome)]
  data <- data[, c(unit.id, time.id, treatment, outcome, othercols)] #reorder columns 
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  temp.treateds <- findAllTreated(ordered.data, treatedvar = treatment, time.var = time.id, unit.var = unit.id, hasbeensorted = TRUE)
  if(nrow(temp.treateds) == 0) stop("no treated units")
  msets <- get.matchedsets(temp.treateds[, time.id], temp.treateds[, unit.id], ordered.data, lag, time.id, unit.id, treatment, hasbeensorted = TRUE)
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