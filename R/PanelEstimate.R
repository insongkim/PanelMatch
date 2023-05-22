#' PanelEstimate
#'
#' \code{PanelEstimate} estimates a causal quantity of interest, including the average treatment effect for
#' treated or control units (att and atc, respectively), the average effect of treatment reversal on reversed units (art), or average treatment effect (ate), as specified in \code{PanelMatch()}.
#' This is done by estimating the counterfactual outcomes for each treated unit using
#' matched sets. Users will provide matched sets that were obtained by the
#' \code{PanelMatch} function and obtain point estimates via a
#' weighted average computation with weighted bootstrap standard errors. Point estimates and standard errors will be
#' produced for each period in the lead window specified by the \code{lead} argument from \code{PanelMatch()}.
#' Users may run multiple estimations by providing lists of each argument to the function.
#' However, in this format, every argument must be explicitly specified in each configuration
#' and must adhere to the same data types/structures outlined below. 
#'
#' @param sets A \code{PanelMatch} object attained via the
#' \code{PanelMatch()} function.
#' @param data The same time series cross sectional data set provided to the \code{PanelMatch()} function used to produce
#' the matched sets.
#' @param se.method Method used for calculating standard errors, provided as a character string. Users must choose between "bootstrap", "conditional", and "unconditional" methods. Default is "bootstrap". "bootstrap" uses a block bootstrapping procedure to calculate standard errors. The conditional method calculates the variance of the estimator, assuming independence across units but not across time. The unconditional method also calculates the variance of the estimator analytically, but makes no such assumptions about independence across units. When the quantity of interest is "att", "atc", or "art", all methods are available. Only "bootstrap" is available for the ate. If \code{pooled} argument is TRUE, then only bootstrap is available. 
#' @param number.iterations If using bootstrapping for calculating standard errors, this is the number of bootstrap iterations. Provide as integer. If \code{se.method} is not equal to "bootstrap", this argument has no effect.
#' @param df.adjustment A logical value indicating whether or not a
#' degree-of-freedom adjustment should be performed for the standard error
#' calculation. The default is \code{FALSE}. This parameter is only available for the bootstrap method of standard error calculation.
#' @param confidence.level A numerical value specifying the confidence level and range of interval
#' estimates for statistical inference. The default is .95.
#' @param moderator The name of a moderating variable, provided as a character string. If a moderating variable is provided,the returned object will be a list of \code{PanelEstimate} objects. The names of the list will reflect the different values of the moderating variable. More specifically, the moderating variable values will be converted to syntactically proper names using \code{make.names()}.
#' @param pooled Logical. If TRUE, estimates and standard errors are returned for treatment effects pooled across the entire lead window. Only available for \code{se.method = ``bootstrap''}
#' @param include.placebo.test Logical. If TRUE, a placebo test is run and returned in the results. The placebo test uses the same specifications for calculating standard errors as the main results. That is, standard errors are calculated according to the user provided \code{se.method} and \code{confidence.level} arguments. If these are invalid for some reason, an error will be thrown. 
#' 
#' @return \code{PanelEstimate} returns a list of class
#' `PanelEstimate' containing the following components:
#' \item{estimates}{the point estimates of the quantity of interest for the lead periods specified}
#' \item{se.method}{The method used to calculate standard errors. This is the same as the argument provided to the function.}
#' \item{bootstrapped.estimates}{the bootstrapped point estimate values, when applicable}
#' \item{bootstrap.iterations}{the number of iterations used in bootstrapping, when applicable}
#' \item{method}{refinement method used to create the matched sets from which the estimates were calculated}
#' \item{lag}{See PanelMatch() argument \code{lag} for more information.}
#' \item{lead}{The lead window sequence for which \code{PanelEstimate()} is producing point estimates and standard errors.}
#' \item{confidence.level}{the confidence level}
#' \item{qoi}{the quantity of interest}
#' \item{matched.sets}{the refined matched sets used to produce the estimations}
#' \item{standard.error}{the standard error(s) of the point estimates}
#' \item{pooled}{Logical indicating whether or not estimates were calculated for individual lead periods or pooled.}
#' \item{placebo.test}{if \code{include.placebo.test = TRUE}, a placebo test is conducted using \code{placebo_test()} and returned as a list. See documentation for \code{placebo_test()} for more about each individual item.}
#' @references Imai, Kosuke, In Song Kim, and Erik Wang (2021)
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
#' PE.results <- PanelEstimate(sets = PM.results, data = dem, number.iterations = 100)
#'
#' @export
PanelEstimate <- function(sets, data,
                          number.iterations = 1000,
                          df.adjustment = FALSE,
                          confidence.level = .95,
                          moderator = NULL,
                          se.method = "bootstrap",
                          pooled = FALSE,
                          include.placebo.test = FALSE)
{
  
  if (pooled && !identical(se.method, "bootstrap"))
  {
    se.method = "bootstrap"
    warning("Pooled results only available with bootstrap SEs")
  }
  if (se.method == "wfe") stop("wfe is no longer supported. Please specify se.method = 'bootstrap', 'conditional', or 'unconditional'")
  if (inherits(number.iterations, "list") & 
      inherits(df.adjustment, "list") &
      inherits(confidence.level, "list")& 
      inherits(sets, "list") == "list")
  {
    if (length(unique(length(se.method), 
                      length(number.iterations), 
                      length(df.adjustment),
                      length(confidence.level), 
                      length(sets))) == 1)
    {
      if(!is.null(moderator))
      {
        
        handle.nesting <- function(data, sets.in, moderating.variable.in,
                                   se.method.in, number.iterations.in, 
                                   df.adjustment.in, confidence.level.in)
        {
          
          if (attr(sets.in, "qoi") == "att")
          {
            s1 = sets.in[["att"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
          }
          if (attr(sets.in, "qoi") == "atc")
          {
            s1 = sets.in[["atc"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
          }
          if (attr(sets.in, "qoi") == "ate")
          { #can assume they are the same
            s1 <- sets.in[["att"]]
            # sets[["atc"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
          } 
          if (attr(sets.in, "qoi") == "art")
          {
            s1 <- sets.in[["art"]]
            unit.id <- attr(s1, "id.var")
            time.id <- attr(s1, "t.var")
          }
          
          if ((attr(sets.in, "qoi") == "art") || 
              (attr(sets.in, "qoi") == "att"))
          {
            att.sets.in <- s1
          } else {
            att.sets.in <- NULL
          }
          
          ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
          set.list <- handle_moderating_variable(ordered.data = ordered.data,
                                                 att.sets = att.sets.in,
                                                 atc.sets = sets[["atc"]],
                                                 moderator = moderating.variable.in,
                                                 unit.id = unit.id, 
                                                 time.id = time.id,
                                                 PM.object = sets, 
                                                 qoi.in = attr(sets, "qoi"))
          res <- lapply(set.list, FUN = panel_estimate, 
                        se.method = se.method, 
                        number.iterations = number.iterations.in,
                        df.adjustment = df.adjustment.in, 
                        confidence.level = confidence.level.in, 
                        data = data,
                        pooled = pooled, 
                        include.placebo.test = include.placebo.test)
          return(res)
        }
        res <- mapply(FUN = handle.nesting, 
                      number.iterations.in = number.iterations,
                      df.adjustment.in = df.adjustment,
                      confidence.level.in = confidence.level, 
                      sets.in = sets,
                      MoreArgs = list(data = data, 
                                      se.method = se.method, 
                                      moderating.variable.in = moderator,
                                      pooled = pooled),
                      SIMPLIFY = FALSE)
      }
      else
      {
        res = mapply(FUN = panel_estimate, 
                     number.iterations = number.iterations,
                     df.adjustment = df.adjustment,
                     confidence.level = confidence.level, 
                     sets = sets,
                     MoreArgs = list(data = data, 
                                     se.method = se.method,
                                     pooled = pooled,
                                     include.placebo.test = include.placebo.test),
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
      if (attr(sets, "qoi") == "art")
      {
        s1 <- sets[["art"]]
        unit.id <- attr(s1, "id.var")
        time.id <- attr(s1, "t.var")
      }
      
      if ((attr(sets, "qoi") == "art") || 
          (attr(sets, "qoi") == "att"))
      {
        att.sets.in <- s1
      } else {
        att.sets.in <- NULL
      }
      # again, functionally treating art and att the same for moderating variable functionality. no need to create separate argument for art sets
      ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
      set.list <- handle_moderating_variable(ordered.data = ordered.data,
                                             att.sets = att.sets.in,
                                             atc.sets = sets[["atc"]],
                                             moderator = moderator,
                                             unit.id = unit.id, 
                                             time.id = time.id,
                                             PM.object = sets, 
                                             qoi.in = attr(sets, "qoi"))
      
      
      res <- lapply(set.list, 
                    FUN = panel_estimate, 
                    se.method = se.method, 
                    number.iterations = number.iterations,
                    df.adjustment = df.adjustment, 
                    confidence.level = confidence.level, 
                    data = data,
                    pooled = pooled, 
                    include.placebo.test = include.placebo.test)
      
    }
    else
    {
      res = panel_estimate(se.method = se.method, 
                           number.iterations = number.iterations,
                           df.adjustment = df.adjustment, 
                           confidence.level = confidence.level, 
                           sets = sets, 
                           data = data,
                           pooled = pooled, 
                           include.placebo.test = include.placebo.test)
    }
    
  }
  return(res)
}


panel_estimate <- function(sets,
                           data,
                           se.method = "bootstrap",
                           number.iterations = 1000,
                           df.adjustment = FALSE,
                           confidence.level = .95,
                           placebo.test = FALSE,
                           placebo.lead = NULL,
                           pooled = FALSE, 
                           include.placebo.test = FALSE)
{
  lead <- attr(sets, "lead")
  outcome.variable <- attr(sets, "outcome.var")
  if (include.placebo.test == TRUE)
  {
    old.data <- data
  }
  if (!inherits(sets, "PanelMatch")) stop("sets parameter is not a PanelMatch object")
  qoi <- attr(sets, "qoi")
  
  #this is just for metadata extraction
  if (qoi == "ate")
  {
    t.sets <- sets[["att"]]
    t.sets.2 <- sets[["atc"]]
    sets.x <- t.sets[sapply(t.sets, length) > 0]
    sets.y <- t.sets.2[sapply(t.sets.2, length) > 0]
    if (length(sets.x) == 0 && length(sets.y) == 0) stop("do not have adequate data to proceed")
  }
  else
  {
    t.sets <- sets[[qoi]]
  }
  
  if (length(t.sets) == 0)
  {
    return(NA)
  }
  
  lag.in <- attr(t.sets, "lag")
  dependent = outcome.variable
  treatment <- attr(t.sets, "treatment.var")
  unit.id <- attr(t.sets, "id.var")
  time.id <- attr(t.sets, "t.var")
  
  if (any(class(data) != "data.frame")){
    stop("please convert data to data.frame class")
  } 
  
  if (!inherits(data[, unit.id], "integer") && 
      !inherits(data[, unit.id], "numeric")){
    stop("please convert unit id column to integer or numeric")
  } 
  
  
  if (any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data),
                                 formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table(testmat), id = unit.id,
                          variable = time.id, value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table::data.table(d), 
                  data.table::data.table(data), 
                  all.x = TRUE, by = c(unit.id, time.id))
    data <- as.data.frame(data)
    
  }
  
  data <- data[order(data[,unit.id], data[,time.id]), ]
  data <- check_time_data(data, time.id)
  
  if (any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment)]
  data <- data[, c(unit.id, time.id, treatment, othercols)] #reorder columns
  
  
  if (identical(qoi, "ate"))
  {
    sets.att <- sets[["att"]]
    sets.att <- sets.att[sapply(sets.att, length) > 0]
    sets.atc <- sets[["atc"]]
    sets.atc <- sets.atc[sapply(sets.atc, length) > 0]
  }
  if (identical(qoi, "atc"))
  {
    sets.atc <- sets[["atc"]]
    sets.atc <- sets.atc[sapply(sets.atc, length) > 0]
    sets.att <- NULL
  }
  if (identical(qoi, "att") || identical(qoi, "art"))
  {
    sets.att <- sets[[qoi]] #art is functionally equivalent to att from here on out
    sets.atc <- NULL
    sets.att <- sets.att[sapply(sets.att, length) > 0]
  }
  
  if (qoi == "att" || 
      qoi == "art" || 
      qoi == "ate")
  {
    treated.unit.ids.att <- as.numeric(sub("\\..*", "", names(sets.att)))
    if (identical(qoi, "att") || identical(qoi, "art")) treated.unit.ids.atc <- NULL
  }
  if (qoi == "atc" || 
      qoi == "ate")
  {
    treated.unit.ids.atc <- as.numeric(sub("\\..*", "", names(sets.atc)))
    if (identical(qoi, "att")) treated.unit.ids.att <- NULL
    
    
  }
  
  data <- prepare_data(data.in = data, lead = lead,
                      sets.att = sets.att, sets.atc = sets.atc,
                      qoi.in = qoi,
                      dependent.variable = dependent)
  
  
  if (placebo.test)
  {
    pe.results <- calculate_placebo_estimates(qoi.in = qoi,
                                            data.in = data,
                                            lead = lead,
                                            number.iterations = number.iterations,
                                            att.treated.unit.ids = treated.unit.ids.att,
                                            atc.treated.unit.ids = treated.unit.ids.atc,
                                            outcome.variable = dependent,
                                            unit.id.variable = unit.id,
                                            confidence.level = confidence.level,
                                            att.sets = sets.att,
                                            atc.sets = sets.atc,
                                            lag = lag.in,
                                            placebo.lead = placebo.lead,
                                            se.method = se.method)
  } else {
    
    pe.results <- calculate_estimates(qoi.in = qoi,
                                     data.in = data,
                                     lead = lead,
                                     number.iterations = number.iterations,
                                     att.treated.unit.ids = treated.unit.ids.att,
                                     atc.treated.unit.ids = treated.unit.ids.atc,
                                     outcome.variable = dependent,
                                     unit.id.variable = unit.id,
                                     confidence.level = confidence.level,
                                     att.sets = sets.att,
                                     atc.sets = sets.atc,
                                     lag = lag.in,
                                     se.method = se.method,
                                     pooled = pooled)

   
    
  }
  
  if (include.placebo.test == TRUE)
  {
    pt.results <- placebo_test(sets,
                               old.data,
                               number.iterations = number.iterations,
                               confidence.level = confidence.level,
                               plot = FALSE,
                               se.method = se.method)
    pe.results[["placebo.test"]] <- pt.results
  }
  
  
  
  return(pe.results)
  
}
