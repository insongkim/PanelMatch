#' Illustration of PanelMatch
#' 
#' \code{PanelMatch} finds matched set for each treated observation  with the number of lags and forwards as well as covariates for matching and weighting refinement determined by the user
#' @importFrom DataCombine slide MoveFront
#' @importFrom CBPS CBPS
#' @importFrom lasso2 merge.formula
#' @importFrom data.table rbindlist
#' @importFrom plyr rbind.fill
#' @importFrom reshape2 dcast
PanelMatch <- function(lag, max.lead, time.id = "year", qoi = "ate",
                           unit.id = "ccode",
                           treatment, 
                           formula = y ~ treat, 
                           data,
                           weighting = FALSE,
                           M = 3, covariate.only = FALSE,
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
  
  formula <- suppressWarnings(merge.formula(reformulate(termlabels = c(time.id, unit.id), response = dependent),formula))
  
  
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
    smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
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
          fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = covariate), 
                       family = binomial(link = "logit"), data = pooled)
        } else {
          fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                      family = binomial(link = "logit"), data = pooled)
        }
        
      } else {
        if (method == "CBPS") {
          fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                       family = binomial(link = "logit"), data = pooled)
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
                  "max.lead" = max.lead, "data" = d2, "ATT_matches" = even_smaller1, 
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Synth") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                  "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                  "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                  "data" = d2, "method" = method,
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
                  "ATT_matches" = delete.NULLs(lapply(even_smaller1, Panel_vit, lag = lag, 
                                         max.lead = max.lead, M = M, method = method)),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Pscore"|method == "CBPS") {
      
      return(list("treatment" = treatment, "qoi" = qoi, 
                  "dependent" = dependent, "covariate" = covariate, 
                  "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                  "covariate.only" = covariate.only, "lag" = lag, 
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
        smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
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
          
          pooled <- data.table::rbindlist(only.t0) # get a dataset for propensity score generation
          colnames(pooled) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
          
          # get propensity scores
          if (covariate.only == TRUE|method == "SynthPscore"|method == "SynthCBPS") {
            if (method == "CBPS"|method == "SynthCBPS") {
              fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = covariate), 
                           family = binomial(link = "logit"), data = pooled)
            } else {
              fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                          family = binomial(link = "logit"), data = pooled)
            }
            
          } else {
            if (method == "CBPS") {
              fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                           family = binomial(link = "logit"), data = pooled)
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
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Synth") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, 
                      "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                      "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                      "data" = d2, "method" = method, "ATC_matches" = lapply(even_smaller2, 
                                                                             Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                             method = method),
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "SynthPscore"|method == "SynthCBPS") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent,
                      "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                      "M" = M, "covariate.only" = FALSE, "lag" = lag, "max.lead" = max.lead, 
                      "data" = d2, "method" = method,
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
                      "data" = d2, "ATC_matches" = delete.NULLs(lapply(even_smaller2, Panel_vit, 
                                                          lag = lag, max.lead = max.lead, 
                                                          M = M,method = method)),
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Pscore"|method == "CBPS") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                      "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                      "covariate.only" = covariate.only, "lag" = lag, 
                      "max.lead" = max.lead, "data" = d2, "method" = method,
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
            smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
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
                  fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = covariate), 
                               family = binomial(link = "logit"), data = pooled)
                } else {
                  fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                              family = binomial(link = "logit"), data = pooled)
                }
                
              } else {
                if (method == "CBPS") {
                  fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                               family = binomial(link = "logit"), data = pooled)
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
            smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
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
              
              pooled <- data.table::rbindlist(only.t0) # get a dataset for propensity score generation
              if (length(pooled) > 0) {
                colnames(pooled) <- c(time.id, unit.id, treatment, dependent, colnames(d2)[5:length(d2)])
                
                
                # get propensity scores
                if (covariate.only == TRUE|method == "SynthPscore"|method == "SynthCBPS") {
                  if (method == "CBPS"|method == "SynthCBPS") {
                    fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = covariate), 
                                 family = binomial(link = "logit"), data = pooled)
                  } else {
                    fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                                family = binomial(link = "logit"), data = pooled)
                  }
                  
                } else {
                  if (method == "CBPS") {
                    fit0 <- CBPS::CBPS(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                                 family = binomial(link = "logit"), data = pooled)
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
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "Synth") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                          "max.lead" = max.lead, "data" = d2, "method" = method,
                          "ATT_matches" = lapply(even_smaller1,  Panel_vit, 
                                                 lag = lag, max.lead = max.lead, 
                                                 M = M,method = method), 
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
                          "ATC_matches" = delete.NULLs(lapply(even_smaller2, 
                                                 Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                 method = method)), 
                          "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
            } else if (method == "Pscore"|method == "CBPS") {
              return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, 
                          "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, 
                          "M" = M, "covariate.only" = covariate.only, "lag" = lag, 
                          "max.lead" = max.lead, "data" = d2, method = method,
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






