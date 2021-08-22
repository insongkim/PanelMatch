#' PanelMatch
#' 
#' Create refined/weighted sets of treated and control units
#' 
#' \code{PanelMatch} identifies a matched set for each treated
#' observation. Specifically, for a given treated unit, the matched
#' set consists of control observations that have an identical
#' treatment history up to a number of \code{lag}
#' time periods. Researchers must specify \code{lag}. A further refinement of
#' the matched set may be performed by setting a maximum size of each matched
#' set, \code{size.match} (the maximum number of control units that can be matched to a treated unit). Users can 
#' also specify covariates that should be used to identify
#' similar control units and a method for defining similarity/distance between units. This is done 
#' via the \code{covs.formula} and \code{refinement.method} arguments, respectively, which are explained in more detail below.
#' @param lag An integer value indicating the length of treatment history periods to be matched on
#' @param time.id A character string indicating the name of the time 
#' variable in the \code{data}. This data currently must be formatted as sequential integers. 
#' @param unit.id A character string indicating the name of unit identifier in the data. This data must be integer.
#' @param treatment A character string indicating the name of the treatment variable in the \code{data}. 
#' The treatment must be a binary indicator variable (integer with 0 for the control group and 1 for the treatment group).
#' @param outcome.var A character string identifying the outcome variable.
#' @param refinement.method A character string specifying the matching or weighting method to be used for refining the matched sets. The user can choose "mahalanobis", "ps.match", "CBPS.match", "ps.weight", "CBPS.weight", "ps.msm.weight", "CBPS.msm.weight", or "none". The first three methods will use the \code{size.match} argument to create sets of at most \code{size.match} closest control units. Choosing "none" will assign equal weights to all control units in each matched set.
#' @param match.missing Logical variable indicating whether or not units should be matched on the patterns of missingness in their treatment histories. Default is TRUE. When FALSE, neither treated nor control units are allowed to have missing treatment data in the lag window.
#' @param data A \code{data.frame} object containing time series cross sectional data. 
#' Time data must be sequential integers that increase by 1. Unit identifiers must be integers. Treatment data must be binary.
#' @param size.match An integer dictating the number of permitted closest control units in a matched set after refinement. 
#' This argument only affects results when using a matching method ("mahalanobis" or any of the refinement methods that end in ".match").
#' This argument is not needed and will have no impact if included when a weighting method is specified (any \code{refinement.method} that includes "weight" in the name).
#' @param covs.formula One sided formula object indicating which variables should be used for matching and refinement. 
#' Argument is not needed if \code{refinement.method} is set to "none"
#' If the user wants to include lagged variables, this can be done using a function, "lag()", which takes two, unnamed, 
#' positional arguments. The first is the name of the variable which you wish to lag. The second is the lag window, 
#' specified as an integer sequence in increasing order.
#' For instance, I(lag(x, 1:4)) will then add new columns to the data for variable "x" for time t-1, t-2, t-3, and t-4 internally
#' and use them for defining/measuring similarity between units. 
#' Other transformations using the I() function, such as I(x^2) are also permitted.
#' The variables specified in this formula are used to define the similarity/distances between units.
#' @param verbose option to include more information about the \code{matched.set} object calculations, 
#' like the distances used to create the refined sets and weights.
#' @param qoi quantity of interest, provided as a string: \code{att} (average treatment effect on treated units), \code{atc} (average treatment effect of treatment on the control units) \code{art} (average effect of treatment reversal for units that experience treatment reversal), or \code{ate} (average treatment effect). 
#' Note that the qoi for MSM methods will give the estimated average treatment effect of being treated for a chosen \code{lead} 
#' time periods. This differs slightly from the non-MSM methods, where treatment reversal is permitted.
#' @param lead integer sequence specifying the lead window, for which qoi point estimates (and standard errors) will 
#' ultimately be produced. Default is 0 (which corresponds to contemporaneous treatment effect).
#' @param matching logical indicating whether or not any matching on treatment history should be performed. 
#' This is primarily used for diagnostic purposes, and most users will never need to set this to FALSE. Default is TRUE.
#' @param forbid.treatment.reversal Logical indicating whether or not it is permissible for treatment to reverse in the specified lead window. 
#' This must be set to TRUE for MSM methods. When set to TRUE, only matched sets for treated units where treatment is 
#' applied continuously in the lead window are included in the results. Default is FALSE.
#' @param exact.match.variables character vector giving the names of variables to be exactly matched on. These should be time invariant variables. 
#' Exact matching for time varying covariates is not currently supported. 
#' @param listwise.delete TRUE/FALSE indicating whether or not missing data should be handled using listwise deletion or the package's default missing data handling procedures. Default is FALSE.
#' @param use.diagonal.variance.matrix TRUE/FALSE indicating whether or not a regular covariance matrix should be used in 
#' mahalanobis distance calculations during refinement, 
#' or if a diagonal matrix with only covariate variances should be used instead. 
#' In many cases, setting this to TRUE can lead to better covariate balance, especially when there is 
#' high correlation between variables. Default is FALSE. This argument is only necessary when 
#' \code{refinement.method = mahalanobis} and will have no impact otherwise.
#' @param placebo.test logical TRUE/FALSE. indicates whether or not you want to be able to run a placebo test. This will add additional requirements on the data -- specifically, it requires that no unit included in the matching/refinement process can having missing outcome data over the lag window. Additionally, you should not use the outcome variable in refinement when \code{placebo.test = TRUE}.
#' @return \code{PanelMatch} returns an object of class "PanelMatch". This is a list that contains a few specific elements: 
#' First, a \code{matched.set} object(s) that has the same name as the provided qoi if the qoi is "att", "art", or "atc". 
#' If qoi = "ate" then two \code{matched.set} objects will be attached, named "att" and "atc." Please consult the documentation for
#' \code{matched_set} to read more about the structure and usage of \code{matched.set} objects. Also, see the vignette page about matched.set objects for 
#' more information about these objects: \code{vignette("matched_set_objects", package = "PanelMatch")}.
#' The \code{PanelMatch} object also has some additional attributes:
#' \item{qoi}{The qoi specified in the original function call}
#' \item{lead}{the lead window specified in the original function call}
#' \item{forbid.treatment.reversal}{logial value matching the forbid.treatment.reversal parameter provided in the function call.}
#' \item{outcome.var}{character string matching the outcome variable provided in the original function call.}
#' 
#' @references Imai, Kosuke, In Song Kim, and Erik Wang (2021)
#' @author Adam Rauh <amrauh@umich.edu>, In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' #not including any lagged variables
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ tradewb, 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#'
#'
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
                       exact.match.variables = NULL,
                       forbid.treatment.reversal = FALSE,
                       matching = TRUE,
                       listwise.delete = FALSE,
                       use.diagonal.variance.matrix = FALSE,
                       placebo.test = FALSE
                       ) 
{
  #since these features are not done yet, force them to be set to defaults
  adjacency.matrix = NULL
  neighborhood.degree = NULL
  caliper.formula = NULL
  network.caliper.info = NULL
  network.refinement.info = NULL 
  continuous.treatment.info = NULL
  
  if (placebo.test) warning("when placebo.test = TRUE, using the dependent variable in refinment is invalid")
  if(class(lag) == "list" & class(time.id) == "list" & class(unit.id) == "list" & class(treatment) == "list" & 
     class(refinement.method) == "list" & class(size.match) == "list" & class(match.missing) == "list" & 
     class(covs.formula) == "list" & class(verbose) == "list" & class(qoi) == "list" & class(lead) == "list" & 
     class(outcome.var) == "list" & class(forbid.treatment.reversal) == "list" &
     class(matching) == "list" & class(listwise.delete) == "list" & class(use.diagonal.variance.matrix) == "list") #everything but data must be provided explicitly
  {
    #stop("looped version of PanelMatch currently unavailable")
    if(length(unique(length(lag) , length(time.id)  , length(unit.id)  , length(treatment)  , 
                     length(refinement.method)  , length(size.match)  , length(match.missing)  , 
                     length(covs.formula)  , length(verbose)  , length(qoi)  , length(lead)  , 
                     length(outcome.var)  , length(exact.match.variables)  , length(forbid.treatment.reversal)  ,
                     length(matching)  , length(listwise.delete)  , length(use.diagonal.variance.matrix) )) == 1)
    {
      
     list.res <- mapply(FUN = panel_match,
             lag = lag, time.id = time.id, unit.id = unit.id, treatment = treatment,
             refinement.method = refinement.method,
             size.match = size.match,
             match.missing = match.missing,
             covs.formula = covs.formula,
             verbose = verbose,
             qoi = qoi,
             lead = lead,
             outcome.var  = outcome.var,
             exact.match.variables = exact.match.variables,
             forbid.treatment.reversal = forbid.treatment.reversal,
             matching = matching,
             listwise.delete = listwise.delete,
             use.diagonal.variance.matrix = use.diagonal.variance.matrix,
             #adjacency.matrix = adjacency.matrix,
             #neighborhood.degree = neighborhood.degree,
             #caliper.formula = caliper.formula,
             #network.caliper.info = network.caliper.info,
             #network.refinement.info = network.refinement.info,
             #continuous.treatment.info = continuous.treatment.info,
             placebo.test = placebo.test,
             MoreArgs = list(data = data)
             , SIMPLIFY = FALSE)
     return(list.res)
    }
    else {
      stop("arguments are not provided in equal length lists")
    }
  }
  else
  {
    panel_match(lag, time.id, unit.id, treatment,
                refinement.method,
                size.match,
                data,
                match.missing,
                covs.formula,
                verbose,
                qoi,
                lead,
                outcome.var,
                exact.match.variables,
                forbid.treatment.reversal,
                matching,
                listwise.delete,
                use.diagonal.variance.matrix,
                adjacency.matrix,
                neighborhood.degree,
                caliper.formula,
                network.caliper.info,
                network.refinement.info,
                continuous.treatment.info,
                placebo.test)
  }
  
}

panel_match <- function(lag, time.id, unit.id, treatment,
                        refinement.method,
                        size.match,
                        data,
                        match.missing,
                        covs.formula,
                        verbose,
                        qoi,
                        lead,
                        outcome.var,
                        exact.match.variables,
                        forbid.treatment.reversal,
                        matching,
                        listwise.delete,
                        use.diagonal.variance.matrix,
                        adjacency.matrix = NULL,
                        neighborhood.degree = NULL,
                        caliper.formula = NULL,
                        network.caliper.info = NULL,
                        network.refinement.info = NULL, 
                        continuous.treatment.info = NULL,
                        placebo.test)
{
  if (!matching & match.missing)
  {
    old.lag <- lag
    lag <- 1
  }
  if (refinement.method == "CBPS.msm.weight" | refinement.method == "ps.msm.weight")
  {
    warning("Note that for msm methods, PanelMatch will attempt to find the estimated average treatment effect of being treated for the entire specified 'lead' time periods.")
  }
  if (listwise.delete & match.missing) stop("set match.missing = FALSE when listwise.delete = TRUE")
  if (lag < 1) stop("please specify a lag value >= 1")
  if (class(data) != "data.frame") stop("please convert data to data.frame class")
  if (!all(refinement.method %in% c("mahalanobis", "ps.weight", "ps.match", "CBPS.weight", "CBPS.match", "ps.msm.weight", "CBPS.msm.weight", "none"))) stop("please choose a valid refinement method")
  if (any(duplicated(data[, c(unit.id, time.id)]))) stop("Time, unit combinations should uniquely identify rows. Please remove duplicates")
  if (!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if (class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  if ( !all(c(time.id, unit.id, treatment, outcome.var)  %in% colnames(data)) ) stop("time id, unit id, outcome, or treatment column name invalid")
  
  if (!is.null(continuous.treatment.info))
  {
    if (!(all(c("treatment.threshold", "units", "matching.threshold", "control.threshold") %in% names(continuous.treatment.info))))
    {
      
      stop("Missing parameter in continuous matching specification.
           Please include all of the treatment.threshold, 
           units, matching.threshold, control.threshold parameters")
    }
    if (continuous.treatment.info[["treatment.threshold"]] == 0)
    {
      stop("treatment.threshold must be > 0")
    }
  }
  
  if (any(c("character", "factor") %in% sapply(data, class)))
  {
    warning("non-numeric data exists. Only numeric (including binary) data can be used for refinement and calculations")
  }
  
  #######take this out when negative lead is implemented:
  if (any(lead < 0)) stop("Please provide positive lead values. Please see the placeboTest function for more.")
  
  if (any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data), formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table(testmat), id = unit.id, variable = time.id, value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table(d), data.table(data), all.x = TRUE, by = c(unit.id, time.id))
    data <- as.data.frame(data)
    # data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  check_time_data(data, time.id)
  if (!all(qoi %in% c("att", "atc", "ate", "art"))) stop("please choose a valid qoi")
  if (!is.null(continuous.treatment.info) && !identical(qoi, "att")) stop('only att is supported for continuous treatment')
  
  if (!forbid.treatment.reversal & all(refinement.method %in% c("CBPS.msm.weight", "ps.msm.weight")))
  {
    stop("please set forbid.treatment.reversal to TRUE for msm methods")
  }
  
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  
  
  ########################################################################### NETWORK CODE ###########################################################################
  ######################################################################################################################################################
  if(!is.null(network.caliper.info) || 
     !is.null(network.refinement.info)) 
  {
    
    ordered.data <- calculate_neighbor_treatment(ordered.data, adjacency.matrix,
                                                 neighborhood.degree, unit.id, 
                                                 time.id, treatment)
    # ordered.data <- handle_network_covariates(ordered.data, adjacency.matrix,
    #                                              neighborhood.degree, unit.id, 
    #                                              time.id, c("logpop", "y"))
    
    ll <- handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, 
                                                ordered.data,
                                                adjacency.matrix, neighborhood.degree, 
                                                unit.id, time.id, treatment,
                                                covs.formula, caliper.formula)
    covs.formula <- ll[[1]]
    caliper.formula <- ll[[2]]
  }

  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  if(!is.null(exact.match.variables))
  {
    for(variable in exact.match.variables)
    {
      ordered.data[, variable] <- as.numeric(as.factor(ordered.data[, variable]))
    }  
  }
  
  if (identical(qoi,"atc") || identical(qoi,"art"))
  {
    if (is.null(continuous.treatment.info))
    {
      ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables   
    }
    if (identical(qoi, "art") && !is.null(continuous.treatment.info)) {
      continuous.treatment.info[["direction"]] <- "negative"
      # in the art case, we focus on negative changes in treatment
      # this parameter of the list does not get modified by the user externally
    } 
    
    msets <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                treatment = treatment, 
                                refinement.method = refinement.method,
                                size.match = size.match, ordered.data = ordered.data, 
                                match.missing = match.missing, covs.formula = covs.formula,
                                verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, 
                                qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, 
                                listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                caliper.formula = caliper.formula, 
                                continuous.treatment.info = continuous.treatment.info,
                                placebo.test = placebo.test)
    #msets <- decode_index(msets, unit.index.map, og.unit.id)
    if (!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }

    pm.obj <- list()
    pm.obj[[qoi]] <- msets
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    if (!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    return(pm.obj)
  } else if (qoi == "att")
  { #note that ordered.data at this point is in column order: unit, time, treatment, everything else
    
    if (identical(qoi, "att") && !is.null(continuous.treatment.info)) {
      continuous.treatment.info[["direction"]] <- "positive"
      # in the art case, we focus on negative changes in treatment
      # this parameter of the list does not get modified by the user externally
    } 
    
    msets <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                treatment = treatment, 
                                refinement.method = refinement.method,
                                size.match = size.match, ordered.data = ordered.data, 
                                match.missing = match.missing, covs.formula = covs.formula,
                                verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, 
                                qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, 
                                listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                caliper.formula = caliper.formula, 
                                continuous.treatment.info = continuous.treatment.info,
                                placebo.test = placebo.test)
    
    
    if (!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    
    pm.obj <- list("att" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    if (!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    attr(pm.obj, "placebo.test") <- placebo.test
    return(pm.obj)
  } else if (identical(qoi, "ate") && is.null(continuous.treatment.info))
  {
    msets <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                treatment = treatment, 
                                refinement.method = refinement.method,
                                size.match = size.match, ordered.data = ordered.data, 
                                match.missing = match.missing, covs.formula = covs.formula,
                                verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, 
                                qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, 
                                listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                caliper.formula = caliper.formula, 
                                continuous.treatment.info = continuous.treatment.info,
                                placebo.test = placebo.test)
    ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables 
    
    msets2 <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                 treatment = treatment, 
                                 refinement.method = refinement.method,
                                 size.match = size.match, ordered.data = ordered.data, 
                                 match.missing = match.missing, covs.formula = covs.formula,
                                 verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                 forbid.treatment.reversal = forbid.treatment.reversal, 
                                 qoi = qoi, matching = matching,
                                 exact.matching.variables = exact.match.variables, 
                                 listwise.deletion = listwise.delete,
                                 use.diag.covmat = use.diagonal.variance.matrix, 
                                 caliper.formula = caliper.formula, 
                                 continuous.treatment.info = continuous.treatment.info,
                                 placebo.test = placebo.test)
    
    if(!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    
    if(!matching & match.missing)
    {
      attr(msets2, "lag") <- old.lag
    }
    pm.obj <- list("att" = msets, "atc" = msets2)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    return(pm.obj)
    
    
  } else {
    stop("qoi not specified correctly")
  }
  
  
}
  
  
