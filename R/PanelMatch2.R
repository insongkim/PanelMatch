#' PanelMatch
#' 
#' \code{PanelMatch} identifies a matched set for each treated
#' observation. Specifically, for a given treated unit, the matched
#' set consists of control observations that have the identical
#' treatment history up to a certain number of \code{lag}
#' years. Researchers must specify \code{lag}. A further refinement of
#' the matched set will be possible by setting the size of the matched
#' set \code{M}, \code{weighting}, and adjusting for other confounders
#' such as past outcomes and covariates via \code{formula}.
#'
#' @param formula A symbolic description of the outcome variable and
#' only covariates with contemporaneous effects.
#' @param lag An integer value indicating the length of history to be matched
#' @param max.lead An optional integer value indicating the length of maximum leads
#' @param unit.id A character string indicating the name of unit identifier
#' variable in the \code{data}. The variable must be numeric. 
#' @param time.id A character string indicating the name of time identifier
#' variable in the \code{data}. This variable must be numeric. 
#' @param treatment A character string indicating the name of
#' treatment variable in the \code{data}. The treatment should be a
#' binary indicator (integer with 0 for the control group and 1 for
#' the treatment group).
#' @param covars.lagged A character vector with the names of the covariates
#' to be lagged.
#' @param qoi One of ``att'' (average treatment effect for the
#' treated) or ``ate'' (average treatment effect) or atc (average
#' treatment effect for the control). The default is ``att.''
#' @param data The data frame in data.frame class
#' @param weighting A logical value indicating whether refinement
#' through weighting should be performed. The default is \code{FALSE}.
#' @param M An integer value indicating the maximum number of the
#' control units in each matched set. The default is 3.
#' @param covariate.only A logical value indicating whether only covariates
#' are adjusted for when refining the matched sets, not lagged outcome.
#' The default is
#' \code{FALSE}.
#' @param method An optional character string indicating the matching
#' methods for refining the matched set based on covariates. One of
#' ``Maha'', ``Pscore'', ``CBPS'', ``SynthPscore'', and ``SynthCBPS''.
#' @param naive An logical value indicating whether
#' the user decides not to match on treatment history. By default 
#' it is FALSE. If TRUE, then \code{PanelMatch} will not match on
#' treatment history. 
#' @param restricted A logical value indicating whether the user decides 
#' to match with restricted treatment pattern
#' @param refinement A logical value indicating whether the user wants to 
#' refine matched sets based on lagged outcomes and/or covariates (lagged
#' or not lagged). If set to FALSE, then the function will give equal 
#' weight to all control units per matched set. By default it is TRUE.
#' 
#' @return \code{PanelMatch} returns a list of class `panelmatch'
#' containing the following components:
#' \item{treatment}{treatment variable name}
#' \item{qoi}{the quantity of interest}
#' \item{dependent}{the dependent variable}
#' \item{covariate}{the covariates}
#' \item{unit.id}{the unit id}
#' \item{time.id}{the time id}
#' \item{M}{the length of history to be matched}
#' \item{covariate.only}{the indicator of whether only covariates are adjusted for when refining the matched sets}
#' \item{lag}{the length of lags}
#' \item{max.lead}{the length of maximum lead}
#' \item{data}{the data frame}
#' \item{ATT_matches}{the matched set given the att as the quantity of interest}
#' \item{ATC_matches}{the matched set given the atc as the quantity of interest}
#' \item{covariate_names}{Names of covariates in the user-specified formula}
#'
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples \dontrun{
#' matches.cbps <- PanelMatch(lag = 4, max.lead = 4,
#'                            time.id = "year", unit.id = "wbcode2",
#'                            treatment = "dem", formula = y ~ dem,
#'                            method = "CBPS", weighting = FALSE,
#'                            qoi = "ate", M = 5, data = dem)
#' }
#' @export
PanelMatch2 <- function(lag, time.id, unit.id, treatment, outcome,
                       refinement.method = c(NULL, ps.weight, ps.match, mahalanobis, CBPS.weight, CBPS.match),
                       size.match = 10,
                       data,
                       match.missing = TRUE,
                       covs.formula) {
  browser()
  
  # set covariates and dependent
  
  #order and ensure order does not get violated
  #convert to data.table?
  #rename outcome variable to dependendent_l*?
  # EXPAND AND SORT WITH NAs
  ##ignoring the fill na's part for right now
  
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  temp.treateds <- findAllTreated(ordered.data, treatedvar = treatment, time.var = time.id, unit.var = unit.id, hasbeensorted = TRUE)
  msets <- get.matchedsets(temp.treateds[, time.id], temp.treateds[, unit.id], ordered.data, lag, time.id, unit.id, treatment, hasbeensorted = TRUE)
  msets <- msets[sapply(msets, length) > 0 ]
  treated.ts <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(msets), split = "[.]"))[c(T,F)])
  
  # compmat <- data.table::dcast(data.table::as.data.table(ordered.data), formula = paste0(unit.id, "~", time.id), value.var = outcome.var)
  
  ordered.data <- as.matrix(parse_and_prep(formula = covs.formula, data = ordered.data, unit.id = unit.id))
  ordered.data <- as.matrix(handle.missing.data(ordered.data, 5:ncol(ordered.data)))
  #result.index <- create_dmats_for_distance(expanded_data = ordered.data, treated_ids = treated.ids, treated_ts = treated.ts,
  #                          row_key = paste0(ordered.data[,unit.id], ".", ordered.data[, time.id]), matched_sets = msets)
  #might not need to use lapply + unlist -- maybe just unlist the indices and then store as one giant data frame? might require modification to subset_expanded_data function
  qoi <- "att"
  if (qoi == "att") {

    #RE IMPLEMENT RESTRICTED OR NAIVE?
    if(refinement.method == "mahalanobis")
    {
      tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
      idxlist <- get_yearly_dmats(ordered.data, treated.ids, tlist, paste0(ordered.data[,unit.id], ".", ordered.data[, time.id]), matched_sets = msets, 4)
      mahalmats <- build_maha_mats(ordered_expanded_data = ordered.data, idx =  idxlist)
      # idx <- create_dmats_for_distance(expanded_data = ordered.data, treated_ids = treated.ids, treated_ts = tlist,
      #                           row_key = paste0(ordered.data[,unit.id], ".", ordered.data[, time.id]), matched_sets = msets)
      #dmats <- subset_expanded.data(ordered.data, idxlist)
      results.maha <- handle_mahalanobis_calculations(mahalmats)
    }
    n.pooled <- subset_expanded.data(ordered.data, result.index)
      n.pooled <- rbindlist(n.pooled)
      n.pooled <- n.pooled[complete.cases(n.pooled), ]
      browser()
      colnames(pooled) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
      
      # get propensity scores
      
      if (covariate.only == TRUE|method == "SynthPscore"|method == "SynthCBPS") {
        if (method == "CBPS"|method == "SynthCBPS") {
          fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = covariate), 
                                              family = binomial(link = "logit"), data = pooled))
        } else {
          fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                      family = binomial(link = "logit"), data = pooled)
        }
        
      } else {
        if (method == "CBPS") {
          fit0 <- suppressMessages(CBPS::CBPS(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                                              family = binomial(link = "logit"), data = pooled))
        } else {
          fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                      family = binomial(link = "logit"), data = pooled)
        }
        
      }
      
      even_smaller1 <- lapply(even_smaller1, function (x) {
        if (covariate.only == TRUE|method == "SynthPscore"|method == "SynthCBPS") {
          x <- as.data.frame(x)
          colnames(x)[1:4] <- c("V1", "V2", "V3", "V4")
          x$ps <- 1 - 1/(1+exp(as.matrix(cbind(1, x[, covariate])) %*% fit0$coefficients))
          return(x)
        } else {
          x <- as.data.frame(x)
          colnames(x)[1:4] <- c("V1", "V2", "V3", "V4")
          x$ps <- 1 - 1/(1+exp(as.matrix(cbind(1, x[, 5:(length(x))])) %*% fit0$coefficients))
          return(x)
        }
        
      })
      
    }
    
    if (is.null(method)){
      return(list("treatment" = treatment, "qoi" = qoi,
                  "dependent" = dependent, "covariate" = covariate,
                  "unit.id" = unit.id, "time.id" = time.id, "M" = M,
                  "covariate.only" = covariate.only, "lag" = lag, 
                  "covariate_names" = covariate_names,
                  "restricted" = restricted,
                  "max.lead" = max.lead, "data" = d2, "ATT_matches" = even_smaller1))
    } else if (method == "Synth") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                  "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                  "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "restricted" = restricted,
                  "ATT_matches" = lapply(even_smaller1, 
                                         Panel_vit, lag = lag, 
                                         max.lead = max.lead, 
                                         M = M, method = method)))
    } else if (method == "SynthPscore"|method == "SynthCBPS") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                  "M" = M, "covariate.only" = FALSE, "lag" = lag, "max.lead" = max.lead, 
                  "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "restricted" = restricted,
                  "ATT_matches" = lapply(even_smaller1, 
                                         Panel_vit, lag = lag, 
                                         max.lead = max.lead, 
                                         M = M, method = method)))
    } else if (method == "Maha") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id,
                  "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                  "max.lead" = max.lead, "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "restricted" = restricted,
                  "ATT_matches" = delete.NULLs(lapply(even_smaller1, Panel_vit, lag = lag, 
                                                      max.lead = max.lead, M = M, method = method, 
                                                      covariate_names = covariate_names))))
    } else if (method == "Pscore"|method == "CBPS") {
      
      return(list("treatment" = treatment, "qoi" = qoi, 
                  "dependent" = dependent, "covariate" = covariate, 
                  "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                  "covariate.only" = covariate.only, "lag" = lag, 
                  "covariate_names" = covariate_names,
                  "restricted" = restricted,
                  "max.lead" = max.lead, "data" = d2, "method" = method,
                  "ATT_matches" = lapply(even_smaller1, Panel_vit, lag = lag,
                                         weighting = weighting,
                                         max.lead = max.lead, M = M, method = method)))
    } else {
      cat("Please select either NULL or one of the following
          estimation methods: Synth, Maha, Pscore, SynthPscore and SynthCBPS")
    }
  } 















