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
PanelMatch <- function(lag, time.id, unit.id, treatment,
                       refinement.method = c(NULL, "ps.weight", "ps.match", "mahalanobis", "CBPS.weight", "CBPS.match"),
                       size.match = 10,
                       data,
                       match.missing = TRUE,
                       covs.formula,
                       verbose = FALSE,
                       qoi
                       ) 
{
  if(!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  check_time_data(data, time.id)
  # can probably add some checks to avoid doing all this stuff when already integer??
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  ordered.data[, paste0(unit.id, ".int")] <- as.integer(as.factor(data[, unit.id]))
  #ordered.data[, paste0(time.id,".int")] <- as.integer(factor(x = as.character(ordered.data[, time.id]), levels = as.character(sort(unique(ordered.data[, time.id]))), 
  #                  labels = as.character(1:length(unique(ordered.data[, time.id]))), ordered = T))
  unit.index.map <- data.frame(original.id = make.names(as.character(unique(ordered.data[, unit.id]))), new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  #time.index.map <- data.frame(original.time.id = make.names(as.character(unique(ordered.data[, time.id]))), new.time.id = unique(ordered.data[, paste0(time.id, ".int")]), stringsAsFactors = F)
  og.unit.id <- unit.id
  #og.time.id <- time.id
  unit.id <- paste0(unit.id, ".int")
  #time.id <- paste0(time.id, ".int")
  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  
  if(qoi == "atc")
  {
    ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables 
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, match.missing, covs.formula, verbose)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    pm.obj <- list("atc" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    return(msets)
  }
  else if(qoi == "att")
  {
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, match.missing, covs.formula, verbose)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    pm.obj <- list("att" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    return(pm.obj)
  }
  else if(qoi == "ate")
  {
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, match.missing, covs.formula, verbose)
    ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables 
    msets2 <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, match.missing, covs.formula, verbose)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    msets2 <- decode_index(msets2, unit.index.map, og.unit.id)
    pm.obj <- list("att" = msets, "atc" = msets2)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    return(pm.obj)
  }
  
} 