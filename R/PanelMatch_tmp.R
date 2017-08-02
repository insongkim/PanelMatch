PanelMatch_tmp <- function(lag, max.lead, time.id = "year", qoi = "ate",
                       unit.id = "ccode",
                       treatment, covariate, dependent, data,
                       M = 3, covariate.only = FALSE,
                       method = NULL) {
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(data[varnames])
  
  if (method == "Pscore") {
    if (covariate.only == FALSE) {
      dlist <- lapply(1:lag, 
                      function (i) slide(data = d2, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                         NewVar = paste("dependent_l", i, sep="")))
      d2 <- Reduce(function(x, y) {merge(x, y)}, dlist)
      # to include ldvs in varnames
      varnames <- c(time.id, unit.id, treatment, dependent, c(covariate, colnames(d2)[(4 + length(covariate) + 1):length(d2)]))
    } else {
      varnames <- c(time.id, unit.id, treatment, dependent, covariate)
    }
    d2 <- na.omit(d2[varnames])
  }
  
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  if (method == "Pscore") {
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
    
    if (method == "Pscore") {
      # take the forward periods from each subset:
      # IMPORTANT
      Fs <- lapply(even_smaller1, function(x) {
        x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
        return(x)
      })
      
      # to only include lag and the first treatment period
      even_smaller1 <- lapply(even_smaller1, function(x) {
        x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
        return(x)
      })
      
      # add varnames
      even_smaller1 <- lapply(even_smaller1, function (x) {
        colnames(x) <- varnames
        return(x)
      })
      
      # add varnames to the forward-periods-subset
      Fs <- lapply(Fs, function (x) {
        colnames(x) <- varnames
        return(x)
      })
      
      pooled <- rbindlist(even_smaller1) # get a dataset for propensity score generation
      
      # get propensity scores
      if (covariate.only == TRUE) {
        fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                    family = binomial(link = "logit"), data = pooled)
      } else {
        fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                    family = binomial(link = "logit"), data = pooled)
      }
      
      pooled$ps <- fit0$fitted.values
      
      # aggregate to delete duplicates
      aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
      
      even_smaller1 <- lapply(even_smaller1, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
      even_smaller1 <- Map(rbind.fill, even_smaller1, Fs)
    }
    
    if (is.null(method)){
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = even_smaller1, 
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Synth") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                        Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                        method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Maha") {
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                        Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                        method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (method == "Pscore") {
      
      return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                                            Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                                            method = method),
                  "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x[,1]))-1)))
    } else {
      cat("Please either select NULL or chose one of the following three 
                       estimation methods: Synth, Maha and Pscore")
    }
  } else {
    if (qoi == "atc") {
      ### from 1 to zero ###
      dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
      
      ### finding matches using the cpp function ###
      # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
      ### finding matches using the cpp function ###
      
      ### cleaning the output from cpp ###
      # delete both higher level and lower level null entries
      smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
      # further cleaning
      smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
      # use function dframelist.rb_dup to turn every list element into a data.frame
      even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
      
      if (method == "Pscore") {
        # take the forward periods from each subset:
        # IMPORTANT
        Fs <- lapply(even_smaller2, function(x) {
          x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
          return(x)
        })
        
        # to only include lag and the first treatment period
        even_smaller2 <- lapply(even_smaller2, function(x) {
          x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
          return(x)
        })
        
        # add varnames
        even_smaller2 <- lapply(even_smaller2, function (x) {
          colnames(x) <- varnames
          return(x)
        })
        
        # add varnames to the forward-periods-subset
        Fs <- lapply(Fs, function (x) {
          colnames(x) <- varnames
          return(x)
        })
        
        pooled <- rbindlist(even_smaller2) # get a dataset for propensity score generation
        
        # get propensity scores
        if (covariate.only == TRUE) {
          fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                      family = binomial(link = "logit"), data = pooled)
        } else {
          fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                      family = binomial(link = "logit"), data = pooled)
        }
        
        pooled$ps <- fit0$fitted.values
        
        # aggregate to delete duplicates
        aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
        
        even_smaller2 <- lapply(even_smaller2, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
        even_smaller2 <- Map(rbind.fill, even_smaller2, Fs)
      }
      
      if (is.null(method)){
        return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate,
                    "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                    "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                    "data" = d2, "Matched sets for ATC" = even_smaller2, 
                    "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (method == "Synth") {
        return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, 
                    "unit.id" = unit.id, "time.id" = time.id, "M" = M, 
                    "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, 
                    "data" = d2, "Matched sets for ATC" = lapply(even_smaller2, 
                                                                                                                                                                                                                                                          Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                          method = method),
                    "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (method == "Maha") {
        return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATC" = lapply(even_smaller2, 
                                                                                                                                                                                                                                                          Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                          method = method),
                    "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (method == "Pscore") {
        return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATC" = lapply(even_smaller2, 
                                                                                                                                                                                                                                                                              Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                                              method = method),
                    "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x[,1]))-1)))
      } else {
        cat("Please either select NULL or chose one of the following three 
                       estimation methods: Synth, Maha and Pscore")
      }
    } else {
      if (qoi == "ate") {
        # ate
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
        
        if (method == "Pscore") {
          # take the forward periods from each subset:
          # IMPORTANT
          Fs <- lapply(even_smaller1, function(x) {
            x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
            return(x)
          })
          
          # to only include lag and the first treatment period
          even_smaller1 <- lapply(even_smaller1, function(x) {
            x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
            return(x)
          })
          
          # add varnames
          even_smaller1 <- lapply(even_smaller1, function (x) {
            colnames(x) <- varnames
            return(x)
          })
          
          # add varnames to the forward-periods-subset
          Fs <- lapply(Fs, function (x) {
            colnames(x) <- varnames
            return(x)
          })
          
          pooled <- rbindlist(even_smaller1) # get a dataset for propensity score generation
          
          # get propensity scores
          if (covariate.only == TRUE) {
            fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                        family = binomial(link = "logit"), data = pooled)
          } else {
            fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                        family = binomial(link = "logit"), data = pooled)
          }
          
          pooled$ps <- fit0$fitted.values
          
          # aggregate to delete duplicates
          aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
          
          even_smaller1 <- lapply(even_smaller1, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
          even_smaller1 <- Map(rbind.fill, even_smaller1, Fs)
        }
        
        # atc
        dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
        
        ### finding matches using the cpp function ###
        # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
        ### finding matches using the cpp function ###
        
        ### cleaning the output from cpp ###
        # delete both higher level and lower level null entries
        smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = max.lead, dmatrix)), delete.NULLs) 
        # further cleaning
        smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
        # use function dframelist.rb_dup to turn every list element into a data.frame
        even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
        
        if (method == "Pscore") {
          # take the forward periods from each subset:
          # IMPORTANT
          Fs <- lapply(even_smaller2, function(x) {
            x <- x[x$V1 %in% unique(x$V1)[(lag+2):(lag+1+max.lead)], ]
            return(x)
          })
          
          # to only include lag and the first treatment period
          even_smaller2 <- lapply(even_smaller2, function(x) {
            x <- x[x$V1 %in% sort(unique(x$V1))[1:(lag+1)], ]
            return(x)
          })
          
          # add varnames
          even_smaller2 <- lapply(even_smaller2, function (x) {
            colnames(x) <- varnames
            return(x)
          })
          
          # add varnames to the forward-periods-subset
          Fs <- lapply(Fs, function (x) {
            colnames(x) <- varnames
            return(x)
          })
          
          pooled <- rbindlist(even_smaller2) # get a dataset for propensity score generation
          
          # get propensity scores
          if (covariate.only == TRUE) {
            fit0 <- glm(reformulate(response = treatment, termlabels = covariate), 
                        family = binomial(link = "logit"), data = pooled)
          } else {
            fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                        family = binomial(link = "logit"), data = pooled)
          }
          
          pooled$ps <- fit0$fitted.values
          
          # aggregate to delete duplicates
          aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
          
          even_smaller2 <- lapply(even_smaller2, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
          even_smaller2 <- Map(rbind.fill, even_smaller2, Fs)
        }
        
        if (is.null(method)){
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = even_smaller1, 
                      "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = even_smaller2, 
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Synth") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                            Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                            method = method), 
                      "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                      Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                      method = method), 
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Maha") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                            Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                            method = method), 
                      "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                      Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                      method = method), 
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (method == "Pscore") {
          return(list("treatment" = treatment, "qoi" = qoi, "dependent" = dependent, "covariate" = covariate, "unit.id" = unit.id, "time.id" = time.id, "M" = M, "covariate.only" = covariate.only, "lag" = lag, "max.lead" = max.lead, "data" = d2, "Matched sets for ATT" = lapply(even_smaller1, 
                                                                                                                                                                                                                                                                                Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                                                                                                                                                                                                                                                method = method), 
                      "NC_ATT" = lapply(even_smaller1, function (x) length(unique(x[,1]))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                      Panel_vit, lag = lag, max.lead = max.lead, M = M,
                                                      method = method),
                      "NC_ATC" = lapply(even_smaller2, function (x) length(unique(x[,1]))-1)))
        } else { (cat("Please either select NULL or chose one of the following three 
                       estimation methods: Synth, Maha and Pscore")) 
        }
      } else {
        cat("Please supply one of the following quantity of interest:
               att, atc and ate")
      }
    }
  }
}






