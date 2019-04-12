#' PanelMatch
#' 
#' Create refined/weighted sets of treated and control units
#' 
#' \code{PanelMatch} identifies a matched set for each treated
#' observation. Specifically, for a given treated unit, the matched
#' set consists of control observations that have the identical
#' treatment history up to a certain number of \code{lag}
#' years. Researchers must specify \code{lag}. A further refinement of
#' the matched set will be possible by setting a maximum size of the matched
#' set \code{size.match} (the maximum number of control units that can be matched to a treated unit) and adjusting for other confounders
#' such as past outcomes and covariates via \code{covs.formula} and specifying a method of refinement via \code{refinement.method}.
#' @param lag An integer value indicating the length of history to be matched
#' @param time.id A character string indicating the name of time identifier
#' variable in the \code{data}. This data currently must be integers that increase by one. 
#' @param unit.id A character string indicating the name of unit identifier in the data. This data must be character, integer, or numeric. However, it is recommended to use integers as ids.
#' @param treatment A character string indicating the name of treatment variable in the \code{data}. The treatment must be a binary indicator (integer with 0 for the control group and 1 for the treatment group).
#' @param outcome.var Character string of the outcome variable.
#' @param refinement.method character string of matching or weighting method used for refining the matched sets. The user can choose "mahalanobis", "ps.match", "CBPS.match", "ps.weight", "CBPS.weight", "ps.msm.weight", "CBPS.msm.weight", or "none". The first three methods will use the \code{size.match} argument to create sets of at most \code{size.match} closest control units. Choosing "none" will assign equal weights to all control units in each matched sets.
#' @param match.missing Logical variable indicating whether or not units should be matched on the patterns of missingness in their treatment histories
#' @param data A data.frame object containing time series cross sectional data. Time data must be integers that increase by 1.
#' @param size.match Maximum size of the matched sets after refinement. This argument only affects results when using a matching method (any of the refinement methods that end in .match). This argument is not needed and will have no impact if included on a weighting method.
#' @param covs.formula One sided formula indicating which variables should be used for matching and refinement. Argument is optional if \code{refinement.method} is set to "none"
#' If the user wants to include lagged variables, this can be done using a function, "lag()", which takes two, unnamed, positional arguments. The first is the name of the variable which you wish to lag, specified as a string. The second is the lag window, specified as an integer sequence
#' For instance, lag("x", 1:4) will then add new columns to the data for variable "x" for time t-1, t-2, t-3, and t-4.
#' @param verbose option to include more information about the matched.set object calculations, like the distances used to create the refined sets and weights.
#' @param qoi quantity of interest: att (average treatment effect on treated units), atc (average treatment effect on control units), ate (average treatment effect), or ade (average controlled direct effect). ade should only be used with the ps.msm.weight or CBPS.msm.weight methods.
#' @param lead integer sequence specifying the lead window for which qoi estimates will ultimately be produced. Default is 0.
#' @param matching logical indicating whether or not any matching on treatment history should be performed. This is used for diagnostic purposes. Default is TRUE.
#' @param restricted Logical indicating whether or not it is permissible for treatment to reverse. This must be set to TRUE for msm methods. When set to TRUE, only matched sets where treatment is applied continuously are included.
#' @return \code{PanelMatch} returns an object of class "PanelMatch". This is a list that contains a few specific elements: First, a matched.set object(s) that has the same name as the provided qoi if the qoi is "att", "atc", or "ade". 
#' If qoi = "ate" then two matched.set objects will be attached, named "att" and "atc." This object also has some additional attributes:
#' \item{qoi}{The qoi specified in the original function call}
#' \item{lead}{the lead window specified in the original function call}
#' \item{restricted}{logial value matching the restricted parameter provided in the function call.}
#' \item{outcome.var}{character string matching the outcome variable provided in the original function call.}
#' For more information about matched.set objects, see documentation for the "matched_set" function
#' @author Adam Rauh <adamrauh@mit.edu>, In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples \dontrun{
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = T, 
#'                          covs.formula = ~ lag("tradewb", 1:4) + lag("y", 1:4), 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, restricted = FALSE)
#' }
#' @export
PanelMatch <- function(lag, time.id, unit.id, treatment,
                       refinement.method,
                       size.match = 10,
                       data,
                       match.missing = TRUE,
                       covs.formula = NULL,
                       verbose = FALSE,
                       qoi,
                       lead = 0,
                       outcome.var,
                       restricted = FALSE,
                       matching = TRUE
                       ) 
{
  if(!matching & match.missing)
  {
    lag <- 1
  }
  if(lag < 1) stop("please specify a lag value >= 1")
  if(!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  if(!all(refinement.method %in% c("mahalanobis", "ps.weight", "ps.match", "CBPS.weight", "CBPS.match", "ps.msm.weight", "CBPS.msm.weight", "none"))) stop("please choose a valid refinement method")
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  check_time_data(data, time.id)
  if(!all(qoi %in% c("att", "atc", "ate", "ade"))) stop("please choose a valid qoi")
  
  if(qoi == "ade" & !all(refinement.method %in% c("CBPS.msm.weight", "ps.msm.weight")))
  {
    stop("ade must have one of the following refinement methods: CBPS.msm.weight, ps.msm.weight")
  }
  if(!restricted & (qoi == 'ade' | all(refinement.method %in% c("CBPS.msm.weight", "ps.msm.weight"))))
  {
    stop("please set restricted to TRUE for msm methods")
  }
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  
  ordered.data[, paste0(unit.id, ".int")] <- as.integer(as.factor(data[, unit.id]))
  if(class(data[, unit.id]) == "character") {
    unit.index.map <- data.frame(original.id = make.names(as.character(unique(ordered.data[, unit.id]))), 
                                 new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(data[, unit.id]) == "integer") {
    unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))),
                                 new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(data[, unit.id]) == "numeric") {
    if(all(unique(ordered.data[, unit.id]) == as.integer(unique(ordered.data[, unit.id])))) #actually integers
    {
      unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))),
                                   new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
    }
  }
  else {
  stop("Unit ID Data is not integer, numeric, or character.")
  }
  og.unit.id <- unit.id
  unit.id <- paste0(unit.id, ".int")
  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 

  if(qoi == "atc")
  {
    ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables 
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, 
                                match.missing, covs.formula, verbose, lead= lead, outcome.var = outcome.var, 
                                restricted = restricted, qoi = qoi, matching = matching)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    pm.obj <- list("atc" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "restricted") <- restricted
    return(pm.obj)
  } else if(qoi == "att")
  {
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data,
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                restricted = restricted, qoi = qoi, matching = matching)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    pm.obj <- list("att" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "restricted") <- restricted
    return(pm.obj)
  } else if(qoi == "ate")
  {
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, 
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                restricted = restricted, qoi = qoi, matching = matching)
    ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables 
    msets2 <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data,
                                 match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                 restricted = restricted, qoi = qoi, matching = matching)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    msets2 <- decode_index(msets2, unit.index.map, og.unit.id)
    pm.obj <- list("att" = msets, "atc" = msets2)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "restricted") <- restricted
    return(pm.obj)
  } else #ade
  {
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data,
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                restricted = restricted, qoi = qoi, matching = matching)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    pm.obj <- list("ade" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "restricted") <- restricted
    return(pm.obj)
  }
  
} 