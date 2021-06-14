#' PanelEstimate
#'
#' \code{PanelEstimate} estimates a causal quantity of interest, including the average treatment effect for 
#' treated or control units (att and atc, respectively), or average treatment effect (ate), as specified in \code{PanelMatch}.
#' This is done by estimating the counterfactual outcomes for each treated unit using
#' matched sets. Users will provide matched sets that were obtained by the
#' \code{PanelMatch} function and obtain point estimates via a
#' weighted average computation with weighted bootstrap standard errors. Point estimates and standard errors will be 
#' produced for each period in the lead window specified by the \code{lead} argument from \code{PanelMatch}. 
#' Users may run multiple estimations by providing lists of each argument to the function. 
#' However, in this format, every argument must be explicitly specified in each configuration 
#' and must adhere to the same data types/structures outlined below. See the included code examples for more about 
#' how this functionality works.
#' 
#' @param number.iterations An integer value indicating the number of bootstrap
#' iterations. The default is 1000.
#' @param sets A \code{PanelMatch} object attained via the
#' \code{PanelMatch} function.
#' @param df.adjustment A logical value indicating whether or not a
#' degree-of-freedom adjustment should be performed for the standard error
#' calculation. The default is \code{FALSE}.
#' @param confidence.level A numerical value specifying the confidence level and range of interval
#' estimates for statistical inference. The default is .95.
#' @param moderator The name of a moderating variable, provided as a character string. If a moderating variable is provided
#' the returned object will be a list of \code{PanelEstimate} objects. The names of the list will reflect the different values of the 
#' moderating variable. More specifically, the moderating variable values will be converted to syntactically proper names using 
#' \code{make.names}.
#' @param data The same time series cross sectional data set provided to the PanelMatch function used to produce 
#' the matched sets
#' @return \code{PanelEstimate} returns a list of class
#' `PanelEstimate' containing the following components:
#' \item{estimates}{the point estimates of the quantity of interest for the lead periods specified}
#' \item{bootstrapped.estimates}{the bootstrapped point estimate values}
#' \item{bootstrap.iterations}{the number of iterations used in bootstrapping}
#' \item{method}{refinement method used to create the matched sets from which the estimates were calculated}
#' \item{lag}{See PanelMatch argument \code{lag} for more information.}
#' \item{lead}{The lead window sequence for which \code{PanelEstimate} is producing point estimates and standard errors.}
#' \item{confidence.level}{the confidence level}
#' \item{qoi}{the quantity of interest}
#' \item{matched.sets}{the refined matched sets used to produce the estimations}
#' \item{standard.error}{the standard error(s) of the point estimates}
#' 
#' @references Imai, Kosuke, In Song Kim, and Erik Wang (2018)
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <amrauh@umich.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = TRUE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem, number.iterations = 500)
#'
#' 
#'
#' @export
PanelEstimate <- function(sets, data,
                          number.iterations = 1000,
                          df.adjustment = FALSE,
                          confidence.level = .95,
                          moderator = NULL) 
{
  
  
  if(class(number.iterations) == "list" & class(df.adjustment) == "list" & 
     class(confidence.level) == "list" & class(sets) == "list")
  {
    if(length(unique(length(number.iterations), length(df.adjustment), 
                     length(confidence.level), length(sets))) == 1)
    {
      if(!is.null(moderator))
      {
        
        handle.nesting <- function(data, sets.in, moderating.variable.in,
                                   number.iterations.in, df.adjustment.in, 
                                   confidence.level.in) 
        {
            if(attr(sets.in, "qoi") == "att")
            {
              s1 = sets.in[["att"]]
              unit.id <- attr(s1, "id.var")
              time.id <- attr(s1, "t.var")
            }
            if(attr(sets.in, "qoi") == "atc")
            {
              s1 = sets.in[["atc"]]
              unit.id <- attr(s1, "id.var")
              time.id <- attr(s1, "t.var")
            }
            if(attr(sets.in, "qoi") == "ate")
            { #can assume they are the same
              s1 <- sets.in[["att"]]
              # sets[["atc"]]
              unit.id <- attr(s1, "id.var")
              time.id <- attr(s1, "t.var")
            }
              
              ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
          
              set.list <- handle_moderating_variable(ordered.data = ordered.data, 
                                                     att.sets = sets.in[["att"]], 
                                                  atc.sets = sets.in[["atc"]],
                                                 moderator = moderating.variable.in, 
                                                 unit.id = unit.id, time.id = time.id,
                                                 PM.object = sets.in)
            
              res <- lapply(set.list, FUN = panel_estimate, 
                            number.iterations = number.iterations.in, 
                          df.adjustment = df.adjustment.in, 
                          confidence.level = confidence.level.in, data = data)
              return(res)
        }
        res <- mapply(FUN = handle.nesting, number.iterations.in = number.iterations, 
                     df.adjustment.in = df.adjustment, 
                     confidence.level.in= confidence.level, 
                     sets.in = sets, 
                     MoreArgs = list(data = data, 
                                     moderating.variable.in = moderator), 
                     SIMPLIFY = FALSE)
      }
      else
      {
        res = mapply(FUN = panel_estimate, number.iterations = number.iterations, 
                     df.adjustment = df.adjustment, 
                     confidence.level = confidence.level, sets = sets, 
                     MoreArgs = list(data = data),
                     SIMPLIFY = FALSE)  
      }
      
    }
    else {
      stop("arguments are not provided in equal length lists")
    }
  }
  else 
  {
    if (!is.null(moderator))
    {
      if (attr(sets, "qoi") == "att")
      {
        s1 = sets[["att"]]
        unit.id <- attr(s1, "id.var")
        time.id <- attr(s1, "t.var")
      }
      if (attr(sets, "qoi") == "atc")
      {
        s1 = sets[["atc"]]
        unit.id <- attr(s1, "id.var")
        time.id <- attr(s1, "t.var")
      }
      if (attr(sets, "qoi") == "ate")
      { #can assume they are the same
        s1 <- sets[["att"]]
        # sets[["atc"]]
        unit.id <- attr(s1, "id.var")
        time.id <- attr(s1, "t.var")
      }
      
      ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
      set.list <- handle_moderating_variable(ordered.data = ordered.data, 
                                             att.sets = sets[["att"]],
                                             atc.sets = sets[["atc"]],
                                             moderator = moderator, 
                                             unit.id = unit.id, time.id = time.id,
                                             PM.object = sets)
      
      res <- lapply(set.list, FUN = panel_estimate, 
                    number.iterations = number.iterations, 
                    df.adjustment = df.adjustment, 
                    confidence.level = confidence.level, data = data)
      
    }
    else
    {
      res = panel_estimate(number.iterations = number.iterations, 
                           df.adjustment = df.adjustment, 
                           confidence.level = confidence.level, 
                           sets = sets, data = data)  
    }
    
  }
  return(res)
}


panel_estimate <- function(sets,
                           data,
                           number.iterations = 1000,
                           df.adjustment = FALSE,
                           confidence.level = .95
                           )
{
  
  lead <- attr(sets, "lead")
  outcome.variable <- attr(sets, "outcome.var")
  continuous.treatment <- attr(sets,'continuous.treatment')
  direction.treatment <- attr(sets, "continuous.treatment.direction")
  if (class(sets) != "PanelMatch") stop("sets parameter is not a PanelMatch object")
  qoi <- attr(sets, "qoi")
  
  
  #this is just for meta data extraction
  if (qoi == "ate")
  {
    temp.sets <- sets
    sets <- sets[["att"]] #just picking one of the two because they should be the same
  }
  else
  {
    sets <- sets[[qoi]]  
  }
  sets <- sets[sapply(sets, length) > 0]
  

  dependent <- outcome.variable
  treatment <- attr(sets, "treatment.var")
  unit.id <- attr(sets, "id.var")
  time.id <- attr(sets, "t.var")
  
  
  if (!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  
  if (!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if (class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  
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
    
  }
  check_time_data(data, time.id)
  
  data <- data[order(data[,unit.id], data[,time.id]), ]
  if (any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment)]
  data <- data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  
  
  if (identical(qoi, "atc"))
  {
    sets.atc <- sets
    sets.att <- NULL
  }
  if (identical(qoi, "att"))
  {
    sets.att <- sets
    sets.atc <- NULL
  }
  if (qoi == "ate")
  {
    sets.att <- temp.sets$att 
    sets.atc <- temp.sets$atc
  }
  
  if (qoi == "att" | qoi == "ate") 
  {
    
    treated.unit.ids.att <- as.numeric(sub("\\..*", "", names(sets.att)))
    if (identical(qoi, "att")) treated.unit.ids.atc <- NULL
  } 
  if (qoi == "atc" | qoi == "ate") 
  {
    treated.unit.ids.atc <- as.numeric(sub("\\..*", "", names(sets.atc)))
    if (identical(qoi, "att")) treated.unit.ids.att <- NULL
  } 
  
  
  data <- prepareData(data.in = data, lead = lead,
                      sets.att = sets.att, sets.atc = sets.atc,
                      continuous.treatment = continuous.treatment,
                      qoi.in = qoi,
                      dependent.variable = dependent)
  
  pe.results <- calculateEstimates(qoi.in = qoi,
                                   data.in = data,
                                   lead = lead,
                                   number.iterations = number.iterations,
                                   att.treated.unit.ids = treated.unit.ids.att,
                                   atc.treated.unit.ids = treated.unit.ids.atc,
                                   outcome.variable = dependent,
                                   unit.id.variable = unit.id,
                                   confidence.level = confidence.level,
                                   att.sets = sets.att,
                                   atc.sets = sets.atc)
  return(pe.results)
  
  
}