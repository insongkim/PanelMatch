PanelMatch <- function(lag, lead, time.id = "year", qoi = "ate",
                       unit.id = "ccode",
                       treatment, covariate, dependent, data,
                       M = 3, 
                       scheme = NULL) {
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(data[varnames])

  if (scheme == "Pscore") {
    if (length(d2)  <= 4) {
      dlist <- lapply(1:lag, 
                      function (i) slide(data = d2, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                         NewVar = paste("dependent_l", i, sep="")))
      d2 <- Reduce(function(x, y) {merge(x, y)}, dlist)
      # to include ldvs in varnames
      varnames <- c(time.id, unit.id, treatment, dependent, c(covariate, colnames(data)[(4 + length(covariate) + 1):length(data)]))
    } else {
      varnames <- c(time.id, unit.id, treatment, dependent, covariate)
    }
    d2 <- na.omit(d2[varnames])
  }
  
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  if (scheme == "Pscore") {
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
    smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = lead, dmatrix)), delete.NULLs) 
    # further cleaning
    smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
    # use function dframelist.rb_dup to turn every list element into a data.frame
    even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
    if (is.null(scheme)){
      return(list("Matched sets for ATT" = even_smaller1, 
                  "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (scheme == "Synth") {
      return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                 Panel_vit, lag = lag, lead = lead, M = M,
                                                           scheme = scheme),
                  "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (scheme == "Maha") {
      return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                 Panel_vit, lag = lag, lead = lead, M = M,
                                                           scheme = scheme),
                  "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else if (scheme == "Pscore") {
      
      return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                 Panel_vit, lag = lag, lead = lead, M = M,
                                                           scheme = scheme),
                  "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
    } else {
      cat("Please either select NULL or chose one of the following three 
                       estimation schemes: Synth, Maha and Pscore")
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
      smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = lead, dmatrix)), delete.NULLs) 
      # further cleaning
      smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
      # use function dframelist.rb_dup to turn every list element into a data.frame
      even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
      if (is.null(scheme)){
        return(list("Matched sets for ATC" = even_smaller2, 
                    "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (scheme == "Synth") {
        return(list("Matched sets for ATC" = lapply(even_smaller2, 
                                                   Panel_vit, lag = lag, lead = lead, M = M,
                                                             scheme = scheme),
                    "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (scheme == "Maha") {
        return(list("Matched sets for ATC" = lapply(even_smaller2, 
                                                   Panel_vit, lag = lag, lead = lead, M = M,
                                                             scheme = scheme),
                    "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else if (scheme == "Pscore") {
        return(list("Matched sets for ATC" = lapply(even_smaller2, 
                                                   Panel_vit, lag = lag, lead = lead, M = M,
                                                             scheme = scheme),
                    "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
      } else {
        cat("Please either select NULL or chose one of the following three 
                       estimation schemes: Synth, Maha and Pscore")
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
        smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = lead, dmatrix)), delete.NULLs) 
        # further cleaning
        smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
        # use function dframelist.rb_dup to turn every list element into a data.frame
        even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
        
        # atc
        dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
        
        ### finding matches using the cpp function ###
        # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
        ### finding matches using the cpp function ###
        
        ### cleaning the output from cpp ###
        # delete both higher level and lower level null entries
        smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = lag, F = lead, dmatrix)), delete.NULLs) 
        # further cleaning
        smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
        # use function dframelist.rb_dup to turn every list element into a data.frame
        even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
        if (is.null(scheme)){
          return(list("Matched sets for ATT" = even_smaller1, 
                      "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = even_smaller2, 
                      "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (scheme == "Synth") {
          return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme), 
                      "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme), 
                      "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (scheme == "Maha") {
          return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme), 
                      "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme), 
                      "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else if (scheme == "Pscore") {
          return(list("Matched sets for ATT" = lapply(even_smaller1, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme), 
                      "# Controls for ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                      "Matched sets for ATC" = lapply(even_smaller2, 
                                                     Panel_vit, lag = lag, lead = lead, M = M,
                                                               scheme = scheme),
                      "# Controls for ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
        } else { (cat("Please either select NULL or chose one of the following three 
                       estimation schemes: Synth, Maha and Pscore")) 
        }
      } else {
        cat("Please supply one of the following quantity of interest:
               att, atc and ate")
      }
    }
  }
}
  
  

  
  
  
