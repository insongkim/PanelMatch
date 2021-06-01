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
#' @param qoi quantity of interest: \code{att} (average treatment effect on treated units), 
#' \code{atc} (average treatment effect on control units), \code{ate} (average treatment effect). 
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
#' @param adjacency.matrix a data.frame object describing an undirected network. The first two columns of the data frame should contain unit ids, and the third column should contain a 1 or 0 indicating adjacency. The unit ids in the first two columns
#' must match the unit ids in the column specified in the \code{unit.id} argument. The third column must be binary, as there are no edge weights currently incorporated into this feature.
#' @param neighborhood.degree an integer specifying the degree of the neighborhood to be considered for each node in the calculations. For instance, setting \code{neighborhood.degree = 1}means
#' that the code will consider the first degree neighbors of each node, but \code{neighborhood.degree = 2} will consider both first and second degree neighbors of each node. Please see the vignette for more information about what network information is calculated.
#' @param caliper.formula a formula object that specifies the caliper to be applied to the data. The caliper is applied after units are matched on treatment history, but before the refinement process. Each caliper is specified in a format
#' similar to the format of the \code{lag} function in the \code{covs.formula} argument. Specifically, it takes the form of a custom function with five positional arguments. 
#' One must specify the name of the variable to use, the caliper method ("average", or "max"), the threshold value, whether the data is categorical or numeric ("categorical" or "numeric"), and whether the threshold should be applied 
#' on a standardized scale (ie. use standard deviations as units) or the original units ("sd" or "raw"). The caliper method parameter specifies whether the threshold should be applied to the distance between
#' a treated unit and each possible control unit over each period in the lag window or to the average distance over the entire lag window (as specified by L). For instance, 
#' \code{caliper.formula = ~ I(caliper(cal.data,"average", .5, "categorical", "raw")) would filter matched sets such that they only contain control units that are <= .5 units of the cal.data variable
#' of the treated unit, on average. In contrast, ~ I(caliper(cal.data,"max", .5, "categorical", "raw"))} would include only control units that are <= .5 units of the cal.data at every period in the lag window. Please see the vignette for more.
#' @param network.caliper.info a named list specifying if/how network data should be used as a caliper. If not specified or set to \code{NULL}, no network information is used in the caliper process. 
#' The provided list must specify the following named elements:
#' \code{use.proportion.data}, \code{proportion.caliper.method}, \code{proportion.caliper.threshold}, \code{prop.unit.type}, \code{use.count.data}, \code{count.caliper.method},
#' \code{count.caliper.threshold}, \code{count.unit.type}. Please see the vignette for more about these parameters and a general discussion about how network information is calculated and used.
#' @param network.refinement.info a named list specifying if/how network data should be used in the refinment process. If not specified or set to \code{NULL}, no network information is used in the refinement process.
#' The provided list must specify the following named elements:
#' \code{use.proportion.data}, \code{proportion.lags}, \code{use.count.data}, \code{count.lags}. Please see the vignette for a discussion of these parameters and how network information is utilized in the package.
#' @param continuous.treatment.info a named list with elements \code{treatment.threshold}, \code{type}, \code{units}, \code{direction}, \code{matching.threshold}, \code{control.treshold}, and optionally \code{minimum.treatment.threshold} and/or \code{maximum.treatment.threshold}.
#' The treatment threshold corresponds to the minimum of the magnitude of the change in the treatment variable from time \code{t-1} to time \code{t}. It must be a positive number.
#' The type is either "raw" or "sd" (specified as a character string) and corresponds to the units of the treatment and mathing threshold. If "raw", then thresholds will be applied to the raw treatment data in the original units.
#' If "sd", then the thresholds are interpreted to be provided in standard deviations, and the treatment variable data will be standardized. 
#' The direction parameter specifies the direction of interest for the change in the treatment variable. It is specified as a character string: "positive", "negative", or "both".
#' So, for example, if the \code{qoi} is set to "att" and the direction is "positive", then only treated units with a positive change in the treatment variable larger than the specified \code{treatment.threshold}
#' will factor into the calculated estimates. The \code{matching.threshold} parameter corresponds to the maximum permissible absolute value of the difference in the treatment variable 
#' between treated units and control units. This is similar to the threshold provided for numerical data in the caliper formula. \code{control.threshold} is a number that specifies the maximum permissable deviation from t-1 to t to still be considered a treated unit.
#' When specified, the minimum (maximum) treatment threshold specifies the minimum (maximum) acceptable value for a treated unit at time t to be considered a viable treatment unit.
#' @return \code{PanelMatch} returns an object of class "PanelMatch". This is a list that contains a few specific elements: 
#' First, a \code{matched.set} object(s) that has the same name as the provided qoi if the qoi is "att" or "atc". See the vignette about continuous treatment for examples and more information.
#' If qoi = "ate" then two \code{matched.set} objects will be attached, named "att" and "atc." Please consult the documentation for
#' \code{matched_set} to read more about the structure and usage of \code{matched.set} objects. Also, see the wiki page for 
#' more information about these objects: \url{https://github.com/insongkim/PanelMatch/wiki/Matched-Set-Objects}.
#' The \code{PanelMatch} object also has some additional attributes:
#' \item{qoi}{The qoi specified in the original function call}
#' \item{lead}{the lead window specified in the original function call}
#' \item{forbid.treatment.reversal}{logial value matching the forbid.treatment.reversal parameter provided in the function call.}
#' \item{outcome.var}{character string matching the outcome variable provided in the original function call.}
#' 
#' @references Imai, Kosuke, In Song Kim, and Erik Wang (2018)
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
#' # Running multiple configurations at once REMOVED FOR NOW...
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
                       adjacency.matrix = NULL,
                       neighborhood.degree = NULL,
                       caliper.formula = NULL,
                       network.caliper.info = NULL,
                       network.refinement.info = NULL, 
                       continuous.treatment.info = NULL
                       ) 
{
  if(class(lag) == "list" & class(time.id) == "list" & class(unit.id) == "list" & class(treatment) == "list" & 
     class(refinement.method) == "list" & class(size.match) == "list" & class(match.missing) == "list" & 
     class(covs.formula) == "list" & class(verbose) == "list" & class(qoi) == "list" & class(lead) == "list" & 
     class(outcome.var) == "list" & class(forbid.treatment.reversal) == "list" &
     class(matching) == "list" & class(listwise.delete) == "list" & class(use.diagonal.variance.matrix) == "list") #everything but data must be provided explicitly
  {
    stop("looped version of PanelMatch currently unavailable")
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
             adjacency.matrix = adjacency.matrix,
             neighborhood.degree = neighborhood.degree,
             caliper.formula = caliper.formula,
             network.caliper.info = network.caliper.info,
             network.refinement.info = network.refinement.info,
             continuous.treatment.info = continuous.treatment.info,
             MoreArgs = list(data = data)
             , SIMPLIFY = F)
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
                continuous.treatment.info)
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
                        adjacency.matrix,
                        neighborhood.degree, 
                        caliper.formula,
                        network.caliper.info,
                        network.refinement.info, 
                        continuous.treatment.info)
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
  if (listwise.delete & match.missing) stop("set match.missing = F when listwise.delete = TRUE")
  if (lag < 1) stop("please specify a lag value >= 1")
  if (class(data) != "data.frame") stop("please convert data to data.frame class")
  if (!all(refinement.method %in% c("mahalanobis", "ps.weight", "ps.match", "CBPS.weight", "CBPS.match", "ps.msm.weight", "CBPS.msm.weight", "none"))) stop("please choose a valid refinement method")
  if (any(duplicated(data[, c(unit.id, time.id)]))) stop("Time, unit combinations should uniquely identify rows. Please remove duplicates")
  if (!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if (class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  if ( !all(c(time.id, unit.id, treatment, outcome.var)  %in% colnames(data)) ) stop("time id, unit id, outcome, or treatment column name invalid")
  
  if (!is.null(continuous.treatment.info))
  {
    if (!(all(c("treatment.threshold", "units", "matching.threshold", "direction", "control.threshold") %in% names(continuous.treatment.info))))
    {
      
      stop("Missing parameter in continuous matching specification.
           Please include all of the treatment.threshold, 
           units, matching.threshold, control.threshold, and direction parameters")
    }
    if (continuous.treatment.info[["treatment.threshold"]] == 0)
    {
      stop("treatment.threshold must be > 0")
    }
  }
  
  #######take this out when negative lead is implemented:
  if (any(lead < 0)) stop("Please provide positive lead values. Negative lead values will be supported in future versions")
  
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
  if(!all(qoi %in% c("att", "atc", "ate"))) stop("please choose a valid qoi")
  
  if(!forbid.treatment.reversal & all(refinement.method %in% c("CBPS.msm.weight", "ps.msm.weight")))
  {
    stop("please set forbid.treatment.reversal to TRUE for msm methods")
  }
  
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  
  #if(!is.null(edge.matrix) & !is.null(neighborhood.degree)) #do early to avoid the encoding/index change, should be safe? 
  ########################################################################### NETWORK CODE ###########################################################################
  ######################################################################################################################################################
  if(!is.null(network.caliper.info) || 
     !is.null(network.refinement.info)) 
  {
    browser()
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
  ###################################################### NETWORK CODE######################################################################################################################################################
  ###########################################################################
  ordered.data[, paste0(unit.id, ".int")] <- as.integer(as.factor(ordered.data[, unit.id]))
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
  if(!is.null(exact.match.variables))
  {
    for(variable in exact.match.variables)
    {
      ordered.data[, variable] <- as.numeric(as.factor(ordered.data[, variable]))
    }  
  }
  
  if(qoi == "atc")
  {
    if(is.null(continuous.treatment.info))
    {
      ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables   
    }
    
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, 
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, caliper.formula = caliper.formula, 
                                continuous.treatment.info = continuous.treatment.info)
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    if(!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    # ### THIS IS FOR ATC, NOT SURE IF IT REALLY MAKES SENSE
    # if (!is.null(continuous.treatment.info))
    # {
    #   if (continuous.treatment.info[["direction"]] == "positive")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) >= 0
    #     msets <- msets[idx]
    #   } else if (continuous.treatment.info[["direction"]] == "negative")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) <= 0
    #     msets <- msets[idx]
    #   } else if (continuous.treatment.info[["direction"]] != "both")
    #   {
    #     stop("direction not well specified")
    #   }  
    # }
    
    
    pm.obj <- list("atc" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    if(!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      attr(pm.obj, "continous.treatment.direction") <- continuous.treatment.info[["direction"]]
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    return(pm.obj)
  } else if(qoi == "att")
  { #note that ordered.data at this point is in column order: unit, time, treatment, everything else
    
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data,
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, caliper.formula = caliper.formula,
                                continuous.treatment.info = continuous.treatment.info)
    
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    if(!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    
    # if (!is.null(continuous.treatment.info))
    # {
    #   if (continuous.treatment.info[["direction"]] == "positive")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) >= 0
    #     msets <- msets[idx]
    #   } else if (continuous.treatment.info[["direction"]] == "negative")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) <= 0
    #     msets <- msets[idx]
    #   } else if (continuous.treatment.info[["direction"]] != "both")
    #   {
    #     stop("direction not well specified")
    #   }  
    # }
    # 
    
    pm.obj <- list("att" = msets)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    if(!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      attr(pm.obj, "continous.treatment.direction") <- continuous.treatment.info[["direction"]]
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    return(pm.obj)
  } else if(qoi == "ate")
  {
    if (!is.null(continuous.treatment.info))
    {
      qoi <- "att"
    }
    
    msets <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data, 
                                match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                forbid.treatment.reversal = forbid.treatment.reversal, qoi = qoi, matching = matching,
                                exact.matching.variables = exact.match.variables, listwise.deletion = listwise.delete,
                                use.diag.covmat = use.diagonal.variance.matrix, caliper.formula = caliper.formula,
                                continuous.treatment.info = continuous.treatment.info)
    
    if (is.null(continuous.treatment.info))
    {
      ordered.data[, treatment] <- ifelse(ordered.data[, treatment] == 1,0,1) #flip the treatment variables   
    }
    if (!is.null(continuous.treatment.info))
    {
      qoi <- "atc"
    }
    msets2 <- perform_refinement(lag, time.id, unit.id, treatment, refinement.method, size.match, ordered.data,
                                 match.missing, covs.formula, verbose, lead = lead, outcome.var = outcome.var, 
                                 forbid.treatment.reversal = forbid.treatment.reversal, qoi = qoi, matching = matching,
                                 exact.matching.variables = exact.match.variables, listwise.deletion = listwise.delete,
                                 use.diag.covmat = use.diagonal.variance.matrix, caliper.formula = caliper.formula,
                                 continuous.treatment.info = continuous.treatment.info)
    
    if (!is.null(continuous.treatment.info)) qoi <- "ate"
    
    msets <- decode_index(msets, unit.index.map, og.unit.id)
    
    if (!matching & match.missing)
    {
      attr(msets, "lag") <- old.lag
    }
    msets2 <- decode_index(msets2, unit.index.map, og.unit.id)
    if(!matching & match.missing)
    {
      attr(msets2, "lag") <- old.lag
    }
    
    # if (!is.null(continuous.treatment.info))
    # {
    #   if (continuous.treatment.info[["direction"]] == "positive")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) >= 0
    #     msets <- msets[idx]
    #     
    #     idx <- sapply(msets2, function(x) attr(x, "treatment.change")) >= 0
    #     msets2 <- msets2[idx]
    #     
    #     
    #   } else if (continuous.treatment.info[["direction"]] == "negative")
    #   {
    #     idx <- sapply(msets, function(x) attr(x, "treatment.change")) <= 0
    #     msets <- msets[idx]
    #     
    #     idx <- sapply(msets2, function(x) attr(x, "treatment.change")) <= 0
    #     msets2 <- msets2[idx]
    #     
    #   } else if (continuous.treatment.info[["direction"]] != "both")
    #   {
    #     stop("direction not well specified")
    #   }  
    # }
    
    pm.obj <- list("att" = msets, "atc" = msets2)
    class(pm.obj) <- "PanelMatch"
    attr(pm.obj, "qoi") <- qoi
    attr(pm.obj, "outcome.var") <- outcome.var
    attr(pm.obj, "lead") <- lead
    attr(pm.obj, "forbid.treatment.reversal") <- forbid.treatment.reversal
    
    if(!is.null(continuous.treatment.info))
    {
      attr(pm.obj, "continuous.treatment") <- TRUE
      attr(pm.obj, "continous.treatment.direction") <- continuous.treatment.info[["direction"]]
    } else {
      attr(pm.obj, "continuous.treatment") <- FALSE
    }
    return(pm.obj)
  } else
  {
    stop("qoi not specified correctly!")
  }
  
  
}
  
  
