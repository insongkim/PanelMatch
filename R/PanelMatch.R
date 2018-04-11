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
#' covariates whose history should be used to refine the matched set.
#' @param lag An integer value indicating the length of history to be matched
#' @param max.lead An optional integer value indicating the length of maximum leads
#' @param unit.id A character string indicating the name of unit
#' variable in the \code{data}
#' @param time.id A character string indicating the name of time
#' variable in the \code{data}
#' @param treatment A character string indicating the name of
#' treatment variable in the \code{data}. The treatment should be a
#' binary indicator (integer with 0 for the control group and 1 for
#' the treatment group).
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
#' \item{NC_ATT}{...}
#' \item{NC_ATC}{...}
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
PanelMatch <- function(formula = y ~ treat, lag, max.lead,
                       unit.id, time.id, treatment, qoi = "att",
                       data, weighting = FALSE,
                       M = 3, covariate.only = FALSE,
                       naive = FALSE,
                       method = NULL) {
  
  ## Warning for missing unit & time index
  if (missing(unit.id))
    stop("'unit.id' should be provided")
  
  if (missing(time.id))
    stop("'time.id' should be provided")
  
  if (missing(data))
    stop("A data.frame should be provided")
  
  if (missing(treatment))
    stop("A treatment variable should be indicated")
  
  if (missing(formula))
    stop("A formula should be provided")
  
  
  # set covariates and dependent
  covariate <- attr(terms(formula),"term.labels")[!attr(terms(formula),"term.labels") == treatment]
  if(length(covariate) == 0) {
    covariate <- NULL # if there is no covariate then it's null
  }
  dependent <- all.vars(formula)[1]
  covariate_names <- covariate
  formula <- suppressWarnings(lasso2::merge.formula(reformulate(termlabels = c(time.id, unit.id), response = dependent),formula))
  
  
  d2 <- as.data.frame(model.matrix(formula, data = data))[,-1]
  d2[dependent] <- model.frame(formula, data=data)[,1]
  d2 <- DataCombine::MoveFront(d2, Var = c(time.id, unit.id, treatment, dependent))
  
  if(method == "Maha" & covariate.only == FALSE){
    
    d2 <- DataCombine::slide(data = d2, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -1,
                             NewVar = "dependent_l1")
    # to include ldvs in varnames
    # varnames <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
    
    d2 <- d2[is.na(c(time.id, unit.id, treatment, dependent, covariate)) == FALSE, ]
    
  } 
  
  
  if (method == "Pscore"|method == "CBPS") {
    
    dlist <- lapply(1:lag, 
                    function (i) DataCombine::slide(data = d2, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                                    NewVar = paste("dependent_l", i, sep="")))
    d2 <- Reduce(function(x, y) {merge(x, y)}, dlist)
    # to include ldvs in varnames
    # varnames <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
    
    d2 <- d2[is.na(c(time.id, unit.id, treatment, dependent, covariate)) == FALSE, ]
  }
  
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  if (method == "Pscore"|method == "CBPS"|method == "SynthPscore"|method == "SynthCBPS"|
      method == "Maha") {
    d2 <- d2[order(d2[,2], d2[,1]), ]
  }
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  if (qoi == "att") {
    # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
    
    ### finding matches using the cpp function ###
    # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
    ### finding matches using the cpp function ###
    
    ### cleaning the output from cpp ###
    # delete both higher level and lower level null entries
    if (naive == FALSE) {
      smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
    } else {
      smallerlist <- lapply(Filter(function (x) !is.null(x), findDDNaive(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
    }
    # further cleaning
    smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
    # use function dframelist.rb_dup to turn every list element into a data.frame
    even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
    
    if (length(even_smaller1) == 0)
      stop("There are no matches for ATT")
    
    if (method == "Pscore"|method == "CBPS"|method == "SynthPscore"|method == "SynthCBPS") {
      only.t0 <- lapply(even_smaller1, function (x) {
        x <- x[x$V1 == unique(x$V1)[lag+1], ]
        return(x)
      })
      
      # # take the forward periods from each subset:
      # # IMPORTANT
      # Fs <- lapply(even_smaller1, function(x) {
      #   x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
      #   return(x)
      # })
      # 
      # # to only include lag and the first treatment period
      # even_smaller1 <- lapply(even_smaller1, function(x) {
      #   x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
      #   return(x)
      # })
      
      # add varnames
      even_smaller1 <- lapply(even_smaller1, function (x) {
        colnames(x) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
        return(x)
      })
      
      # # add varnames to the forward-periods-subset
      # Fs <- lapply(Fs, function (x) {
      #   colnames(x) <- varnames
      #   return(x)
      # })
      
      pooled <- data.table::rbindlist(only.t0) # get a dataset for propensity score generation
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
                  "max.lead" = max.lead, "data" = d2, "ATT_matches" = even_smaller1, 
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Synth") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                  "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                  "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "ATT_matches" = lapply(even_smaller1, 
                                         Panel_vit, lag = lag, 
                                         max.lead = max.lead, 
                                         M = M, method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "SynthPscore"|method == "SynthCBPS") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                  "M" = M, "covariate.only" = FALSE, "lag" = lag, "max.lead" = max.lead, 
                  "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "ATT_matches" = lapply(even_smaller1, 
                                         Panel_vit, lag = lag, 
                                         max.lead = max.lead, 
                                         M = M, method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Maha") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id,
                  "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                  "max.lead" = max.lead, "data" = d2, "method" = method,
                  "covariate_names" = covariate_names,
                  "ATT_matches" = delete.NULLs(lapply(even_smaller1, Panel_vit, lag = lag, 
                                                      max.lead = max.lead, M = M, method = method)),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Pscore"|method == "CBPS") {
      
      return(list("treatment" = treatment, "qoi" = qoi, 
                  "dependent" = dependent, "covariate" = covariate, 
                  "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                  "covariate.only" = covariate.only, "lag" = lag, 
                  "covariate_names" = covariate_names,
                  "max.lead" = max.lead, "data" = d2, "method" = method,
                  "ATT_matches" = lapply(even_smaller1, Panel_vit, lag = lag,
                                         weighting = weighting,
                                         max.lead = max.lead, M = M, method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x[,1]))-1)))
    } else {
      cat("Please select either NULL or one of the following
          estimation methods: Synth, Maha, Pscore, SynthPscore and SynthCBPS")
    }
    } else {
      if (qoi == "atc") {
        dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
        # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
        
        ### finding matches using the cpp function ###
        # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
        ### finding matches using the cpp function ###
        
        ### cleaning the output from cpp ###
        # delete both higher level and lower level null entries
        if (naive == FALSE) {
          smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
        } else {
          smallerlist <- lapply(Filter(function (x) !is.null(x), findDDNaive(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
        }
        # further cleaning
        smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
        # use function dframelist.rb_dup to turn every list element into a data.frame
        even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
        
        if (length(even_smaller2) == 0)
          stop("There are no matches for ATC")
        
        if (method == "Pscore"|method == "CBPS"|method == "SynthPscore"|method == "SynthCBPS") {
          only.t0 <- lapply(even_smaller2, function (x) {
            x <- x[x$V1 == unique(x$V1)[lag+1], ]
            return(x)
          })
          
          # # take the forward periods from each subset:
          # # IMPORTANT
          # Fs <- lapply(even_smaller2, function(x) {
          #   x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
          #   return(x)
          # })
          # 
          # # to only include lag and the first treatment period
          # even_smaller2 <- lapply(even_smaller2, function(x) {
          #   x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
          #   return(x)
          # })
          
          # add varnames
          even_smaller2 <- lapply(even_smaller2, function (x) {
            colnames(x) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
            return(x)
          })
          
          # # add varnames to the forward-periods-subset
          # Fs <- lapply(Fs, function (x) {
          #   colnames(x) <- varnames
          #   return(x)
          # })
          
          pooled <- rbindlist(only.t0) # get a dataset for propensity score generation
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
          
          even_smaller2 <- lapply(even_smaller2, function (x) {
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
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate,
                      "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                      "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                      "data" = d2, "ATC_matches" = even_smaller2, 
                      "covariate_names" = covariate_names,
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Synth") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, 
                      "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                      "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                      "data" = d2, "method" = method, "ATC_matches" = lapply(even_smaller2, 
                                                                             Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                             method = method),
                      "covariate_names" = covariate_names,
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "SynthPscore"|method == "SynthCBPS") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                      "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                      "M" = M, "covariate.only" = FALSE, "lag" = lag, "max.lead" = max.lead, 
                      "data" = d2, "method" = method,
                      "covariate_names" = covariate_names,
                      "ATT_matches" = lapply(even_smaller2, 
                                             Panel_vit, lag = lag, 
                                             max.lead = max.lead, 
                                             M = M, method = method),
                      "NC_ATT" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Maha") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                      "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                      "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                      "max.lead" = max.lead, "method" = method,
                      "covariate_names" = covariate_names,
                      "data" = d2, "ATC_matches" = delete.NULLs(lapply(even_smaller2, Panel_vit, 
                                                                       lag = lag, max.lead = max.lead, 
                                                                       M = M,method = method)),
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Pscore"|method == "CBPS") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                      "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                      "covariate.only" = covariate.only, "lag" = lag, 
                      "max.lead" = max.lead, "data" = d2, "method" = method,
                      "covariate_names" = covariate_names,
                      "ATC_matches" = lapply(even_smaller2, weighting = weighting, 
                                             Panel_vit, lag = lag,
                                             max.lead = max.lead, M = M, method = method),
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x[,1]))-1)))
        } else {
          cat("Please select either NULL or one of the following
              estimation methods: Synth, Maha, Pscore, SynthPscore and SynthCBPS")
        }
        } else {
          if (qoi == "ate") {
            ### cleaning the output from cpp ###
            # delete both higher level and lower level null entries
            if (naive == FALSE) {
              smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
            } else {
              smallerlist <- lapply(Filter(function (x) !is.null(x), findDDNaive(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
            }
            # further cleaning
            smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
            # use function dframelist.rb_dup to turn every list element into a data.frame
            even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
            
            if (length(even_smaller1) == 0)
              stop("There are no matches for ATT")
            
            if (method == "Pscore"|method == "CBPS"|method == "SynthPscore"|method == "SynthCBPS") {
              only.t0 <- lapply(even_smaller1, function (x) {
                x <- x[x$V1 == unique(x$V1)[lag+1], ]
                return(x)
              })
              
              # # take the forward periods from each subset:
              # # IMPORTANT
              # Fs <- lapply(even_smaller1, function(x) {
              #   x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
              #   return(x)
              # })
              # 
              # # to only include lag and the first treatment period
              # even_smaller1 <- lapply(even_smaller1, function(x) {
              #   x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
              #   return(x)
              # })
              
              # add varnames
              even_smaller1 <- lapply(even_smaller1, function (x) {
                colnames(x) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
                return(x)
              })
              
              # # add varnames to the forward-periods-subset
              # Fs <- lapply(Fs, function (x) {
              #   colnames(x) <- varnames
              #   return(x)
              # })
              
              pooled <- rbindlist(only.t0) # get a dataset for propensity score generation
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
            
            # atc
            dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
            # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
            
            ### finding matches using the cpp function ###
            # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
            ### finding matches using the cpp function ###
            
            ### cleaning the output from cpp ###
            # delete both higher level and lower level null entries
            if (naive == FALSE) {
              smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
            } else {
              smallerlist <- lapply(Filter(function (x) !is.null(x), findDDNaive(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
            }
            # further cleaning
            smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
            # use function dframelist.rb_dup to turn every list element into a data.frame
            even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
            
            if (method == "Pscore"|method == "CBPS"|method == "SynthPscore"|method == "SynthCBPS") {
              only.t0 <- lapply(even_smaller2, function (x) {
                x <- x[x$V1 == unique(x$V1)[lag+1], ]
                return(x)
              })
              
              # add varnames
              even_smaller2 <- lapply(even_smaller2, function (x) {
                colnames(x) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
                return(x)
              })
              
              # # add varnames to the forward-periods-subset
              # Fs <- lapply(Fs, function (x) {
              #   colnames(x) <- varnames
              #   return(x)
              # })
              
              pooled <- rbindlist(only.t0) # get a dataset for propensity score generation
              if (length(pooled) > 0) {
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
                
                even_smaller2 <- lapply(even_smaller2, function (x) {
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
              
            }
            #########################
            if (length(even_smaller1) == 0 & length(even_smaller2) == 0) 
              stop("No matches found for ATT.")
            if (length(even_smaller2) == 0) {
              qoi <- "att"
            }
            #########################
            if (is.null(method)){
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                          "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                          "data" = d2, "ATT_matches" = even_smaller1, 
                          "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                          "ATC_matches" = even_smaller2, 
                          "covariate_names" = covariate_names,
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "Synth") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                          "max.lead" = max.lead, "data" = d2, "method" = method,
                          "ATT_matches" = lapply(even_smaller1,  Panel_vit, 
                                                 lag = lag, max.lead = max.lead, 
                                                 M = M,method = method), 
                          "covariate_names" = covariate_names,
                          "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                          "ATC_matches" = lapply(even_smaller2, 
                                                 Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                 method = method), 
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "SynthPscore"|method == "SynthCBPS") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                          "max.lead" = max.lead, "data" = d2, "method" = method,
                          "ATT_matches" = lapply(even_smaller1,  Panel_vit, 
                                                 lag = lag, max.lead = max.lead, 
                                                 M = M,method = method), 
                          "covariate_names" = covariate_names,
                          "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                          "ATC_matches" = lapply(even_smaller2, 
                                                 Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                 method = method), 
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "Maha") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                          "data" = d2, "method" = method,
                          "ATT_matches" = delete.NULLs(lapply(even_smaller1,Panel_vit, lag = lag, 
                                                              max.lead = max.lead, M = M,method = method)), 
                          "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                          "covariate_names" = covariate_names,
                          "ATC_matches" = delete.NULLs(lapply(even_smaller2, 
                                                              Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                              method = method)), 
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "Pscore"|method == "CBPS") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                          "max.lead" = max.lead, "data" = d2, method = method,
                          "covariate_names" = covariate_names,
                          "ATT_matches" = lapply(even_smaller1, Panel_vit, weighting = weighting, 
                                                 lag = lag, max.lead = max.lead, M = M,method = method), 
                          "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x[,1]))-1),
                          "ATC_matches" = lapply(even_smaller2, weighting = weighting,
                                                 Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                 method = method),
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x[,1]))-1)))
            } else {     cat("Please select either NULL or one of the following
          estimation methods: Synth, Maha, Pscore, SynthPscore and SynthCBPS")
            }
          } else {
            cat("Please supply one of the following quantity of interest:
                att, atc and ate")
          }
          }
      }
}












