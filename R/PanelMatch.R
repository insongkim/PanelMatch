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
#' @param refinement.method A character string specifying the matching or weighting method to be used for refining the matched sets. The user can choose "mahalanobis", "ps.match", "CBPS.match", "ps.weight", "CBPS.weight", or "none". The first three methods will use the \code{size.match} argument to create sets of at most \code{size.match} closest control units. Choosing "none" will assign equal weights to all control units in each matched set.
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
#' @param lead integer sequence specifying the lead window, for which qoi point estimates (and standard errors) will 
#' ultimately be produced. Default is 0 (which corresponds to contemporaneous treatment effect).
#' @param matching logical indicating whether or not any matching on treatment history should be performed. 
#' This is primarily used for diagnostic purposes, and most users will never need to set this to FALSE. Default is TRUE.
#' @param forbid.treatment.reversal Logical indicating whether or not it is permissible for treatment to reverse in the specified lead window. 
#' When set to TRUE, only matched sets for treated units where treatment is 
#' applied continuously in the lead window are included in the results. Default is FALSE.
#' @param exact.match.variables character vector giving the names of variables to be exactly matched on. These should be time invariant variables. 
#' Exact matching for time varying covariates is not currently supported. 
#' @param listwise.delete TRUE/FALSE indicating whether or not missing data should be handled using listwise deletion or the package's default missing data handling procedures. Default is FALSE.
#' @param use.diagonal.variance.matrix TRUE/FALSE indicating whether or not a regular covariance matrix should be used in mahalanobis distance calculations during refinement, 
#' or if a diagonal matrix with only covariate variances should be used instead. 
#' In many cases, setting this to TRUE can lead to better covariate balance, especially when there is 
#' high correlation between variables. Default is FALSE. This argument is only necessary when 
#' \code{refinement.method = mahalanobis} and will have no impact otherwise.
#' @param restrict.control.period (optional) integer specifying the number of pre-treatment periods that treated units and potentially matched control units should be non-NULL and in the control state. For instance, specifying 4 would mean that the treatment history cannot contain any missing data or treatment from t-4 to t. 
#' @param placebo.test logical TRUE/FALSE. indicates whether or not you want to be able to run a placebo test. This will add additional requirements on the data -- specifically, it requires that no unit included in the matching/refinement process can having missing outcome data over the lag window. Additionally, you should not use the outcome variable in refinement when \code{placebo.test = TRUE}.
#' @param caliper.formula a formula object that specifies the caliper to be applied to the data. The caliper is applied after units are matched on treatment history, but before the refinement process. Each caliper is specified in a format
#' similar to the format of the \code{lag} function in the \code{covs.formula} argument. Specifically, it takes the form of a custom function with five positional arguments. 
#' One must specify the name of the variable to use, the caliper method ("average", or "max"), the threshold value, whether the data is categorical or numeric ("categorical" or "numeric"), and whether the threshold should be applied 
#' on a standardized scale (ie. use standard deviations as units) or the original units ("sd" or "raw"). The caliper method parameter specifies whether the threshold should be applied to the distance between
#' a treated unit and each possible control unit over each period in the lag window or to the average distance over the entire lag window (as specified by L). For instance, 
#' \code{caliper.formula = ~ I(caliper(cal.data,"average", .5, "categorical", "raw")) would filter matched sets such that they only contain control units that are <= .5 units of the cal.data variable
#' of the treated unit, on average. In contrast, ~ I(caliper(cal.data,"max", .5, "categorical", "raw"))} would include only control units that are <= .5 units of the cal.data at every period in the lag window. Please see the vignette for more.
#' @param continuous.treatment.info a named list with elements \code{treatment.threshold}, \code{type}, \code{units},, \code{matching.threshold}, \code{control.treshold}, and optionally \code{minimum.treatment.value} and/or \code{maximum.treatment.value}.
#' The treatment threshold corresponds to the minimum of the magnitude of the change in the treatment variable from time \code{t-1} to time \code{t}. It must be a positive number. When the qoi is set to be "att", then treatment will be defined as a positive change (ie. an increase) in the treatment variable. When the qoi is set to be "art", treatment is defined as a negative change
#' The type is either "raw" or "sd" (specified as a character string) and corresponds to the units of the treatment and matching threshold. If "raw", then thresholds will be applied to the raw treatment data in the original units.
#' If "sd", then the thresholds are interpreted to be provided in standard deviations, and the treatment variable data will be standardized. 
#' The \code{matching.threshold} parameter corresponds to the maximum permissible absolute value of the difference in the treatment variable 
#' between treated units and control units. This is similar to the threshold provided for numerical data in the caliper formula. \code{control.threshold} is a number that specifies the maximum permissable deviation from t-1 to t to still be considered a treated unit.
#' When specified, the minimum (maximum) treatment value specifies the minimum (maximum) acceptable value for a treated unit at time t to be considered a viable treatment unit.
#' @return \code{PanelMatch()} returns an object of class "PanelMatch". This is a list that contains a few specific elements: 
#' First, a \code{matched.set} object(s) that has the same name as the provided qoi if the qoi is "att", "art", or "atc". 
#' If qoi = "ate" then two \code{matched.set} objects will be attached, named "att" and "atc." Please consult the documentation for
#' \code{matched_set()} to read more about the structure and usage of \code{matched.set} objects. Also, see the vignette page about matched.set objects for 
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
#'                          treatment = "dem", refinement.method = "ps.match", 
#'                          data = dem_small, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#'
#'
#' @export
PanelMatch <- function(lag, time.id, unit.id, 
                       treatment, outcome.var,
                       refinement.method, data,qoi,
                       size.match = 10,
                       match.missing = TRUE,
                       covs.formula = NULL,
                       lead = 0,
                       verbose = FALSE,
                       exact.match.variables = NULL,
                       forbid.treatment.reversal = FALSE,
                       matching = TRUE,
                       listwise.delete = FALSE,
                       use.diagonal.variance.matrix = FALSE,
                       restrict.control.period = NULL,
                       placebo.test = FALSE,
                       caliper.formula = NULL,
                       continuous.treatment.info = NULL) 
{
 
  if (placebo.test) warning("when placebo.test = TRUE, using the dependent variable in refinment is invalid")
  if (!matching && match.missing)
  {
    old.lag <- lag
    lag <- 1
  }
  ##############################error checking##############################
  if (listwise.delete & match.missing) stop("set match.missing = FALSE when listwise.delete = TRUE")
  if (lag < 1) stop("please specify a lag value >= 1")
  if (any(class(data) != "data.frame")) stop("please convert data to data.frame class")
  
  if (!all(refinement.method %in% c("mahalanobis", "ps.weight", "ps.match", "CBPS.weight", "CBPS.match", "none"))) stop("please choose a valid refinement method")
  if (any(duplicated(data[, c(unit.id, time.id)]))) stop("Time, unit combinations should uniquely identify rows. Please remove duplicates")
  if (!inherits(data[, unit.id], "integer") && !inherits(data[, unit.id], "numeric")) stop("please convert unit id column to integer or numeric")
  if ( !all(c(time.id, unit.id, treatment, outcome.var)  %in% colnames(data)) ) stop("time id, unit id, outcome, or treatment column name invalid")
  if (forbid.treatment.reversal && !identical(qoi, "att"))
  {
    stop("forbid.treatment.reversal = TRUE only valid for qoi = att")
  }
  
  if(!is.null(restrict.control.period))
  {
    if(restrict.control.period < 1) stop("restricted control period specification must be >=1")
    if(restrict.control.period > lag) stop("restricted control period specification cannot be greater than lag")
  }
  if (any(lead < 0)) stop("Please provide positive lead values. Please see the placebo_test function for more.")
  if (!all(qoi %in% c("att", "atc", "ate", "art"))) stop("please choose a valid qoi")
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  
  
  if (!is.null(continuous.treatment.info))
  {
    if (!(all(c("treatment.threshold",
                "units", 
                "matching.threshold",
                "control.threshold") %in% names(continuous.treatment.info))))
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
  if (!is.null(continuous.treatment.info) && 
      !(qoi %in% c("att", "art")))
  {
    stop("Only ATT and ART are valid for continuous treatment.")
  }
  ##############################error checking##############################
  
  ##############################balance the panel############################## 
  if (any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data), 
                                 formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table(testmat), 
                          id = unit.id, 
                          variable = time.id, 
                          value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table::data.table(d), 
                  data.table::data.table(data), 
                  all.x = TRUE, 
                  by = c(unit.id, time.id))
    data <- as.data.frame(data)
    
  }
  ##############################balance the panel##############################
  
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  ordered.data <- check_time_data(ordered.data, time.id)
  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  
  # remove this when caliper functionality is implemented
  if(!is.null(exact.match.variables))
  {
    for(variable in exact.match.variables)
    {
      ordered.data[, variable] <- as.numeric(as.factor(ordered.data[, variable]))
    }  
  }
  
  
  if (identical(qoi,"art"))
  {
    # flip the treatment variables only in the case of the binary ART
    if (is.null(continuous.treatment.info))
    {
      ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables     
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
                                restrict.control.period = restrict.control.period,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                placebo.test = placebo.test,
                                caliper.formula = caliper.formula,
                                continuous.treatment.info = continuous.treatment.info)
    
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
    attr(pm.obj, "placebo.test") <- placebo.test
    if (!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    
    return(pm.obj)
  } else if (identical(qoi,"att") || identical(qoi,"atc"))
  { #note that ordered.data at this point is in column order: unit, time, treatment, everything else
    
    
    msets <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                treatment = treatment, 
                                refinement.method = refinement.method,
                                size.match = size.match, 
                                ordered.data = ordered.data, 
                                match.missing = match.missing, 
                                covs.formula = covs.formula,
                                verbose = verbose,
                                lead = lead, 
                                outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, 
                                qoi = qoi, 
                                matching = matching,
                                exact.matching.variables = exact.match.variables, 
                                listwise.deletion = listwise.delete,
                                restrict.control.period = restrict.control.period,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                placebo.test = placebo.test,
                                caliper.formula = caliper.formula,
                                continuous.treatment.info = continuous.treatment.info)
    
    
    if (!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    
    pm.obj <- list( msets)
    names(pm.obj) <- qoi
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    attr(pm.obj, "placebo.test") <- placebo.test
    if (!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    
    return(pm.obj)
  } else if (identical(qoi, "ate"))
  { # for ate, we have to calculate both att and atc
    msets <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                treatment = treatment, 
                                refinement.method = refinement.method,
                                size.match = size.match, ordered.data = ordered.data, 
                                match.missing = match.missing, covs.formula = covs.formula,
                                verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, 
                                qoi = "att", matching = matching,
                                exact.matching.variables = exact.match.variables, 
                                listwise.deletion = listwise.delete,
                                restrict.control.period = restrict.control.period,
                                use.diag.covmat = use.diagonal.variance.matrix, 
                                placebo.test = placebo.test,
                                caliper.formula = caliper.formula,
                                continuous.treatment.info = continuous.treatment.info)
    
    msets2 <- perform_refinement(lag = lag, time.id = time.id, unit.id = unit.id, 
                                 treatment = treatment, 
                                 refinement.method = refinement.method,
                                 size.match = size.match, ordered.data = ordered.data, 
                                 match.missing = match.missing, covs.formula = covs.formula,
                                 verbose = verbose, lead = lead, outcome.var = outcome.var, 
                                 forbid.treatment.reversal = forbid.treatment.reversal, 
                                 qoi = "atc", matching = matching,
                                 exact.matching.variables = exact.match.variables, 
                                 listwise.deletion = listwise.delete,
                                 restrict.control.period = restrict.control.period,
                                 use.diag.covmat = use.diagonal.variance.matrix, 
                                 placebo.test = placebo.test,
                                 caliper.formula = caliper.formula,
                                 continuous.treatment.info = continuous.treatment.info)
    
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
    attr(pm.obj, "placebo.test") <- placebo.test
    if (!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    
    return(pm.obj)
    
  } else {
    stop("qoi not specified correctly")
  }
  
}