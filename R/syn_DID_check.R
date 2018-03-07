syn_DID_check <- function(L, FORWARD, time.id = "year", qoi = "ATE",
                          unit.id = "ccode",
                          treatment, covariate, dependent, d) {

  
  varnames <- c(time.id, unit.id, treatment, covariate, dependent)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
  
  ### finding matches using the cpp function ###
  # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
  ### finding matches using the cpp function ###
  
  ### cleaning the output from cpp ###
  # delete both higher level and lower level null entries
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L, F = FORWARD, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  
  # brute forcing qoi
  all.diffs.weighted <- lapply(even_smaller1, cscwdid_tmp, L = L, FORWARD = FORWARD)
  ATT <- Reduce("+", all.diffs.weighted)/length(all.diffs.weighted)
  if (qoi == "ATT") {
    return(list("ATT" = ATT, 
                "No.controls" = lapply(even_smaller1, function (x) length(unique(x$V2))-1)))
  } else {
    ### from 1 to zero ###
    dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
    
    ### finding matches using the cpp function ###
    # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
    ### finding matches using the cpp function ###
    
    ### cleaning the output from cpp ###
    # delete both higher level and lower level null entries
    smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = L, F = FORWARD, dmatrix)), delete.NULLs) 
    # further cleaning
    smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
    # use function dframelist.rb_dup to turn every list element into a data.frame
    even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
    
    # brute forcing qoi
    all.diffs.weighted2 <- lapply(even_smaller2, cscwdid_tmp, L = L, FORWARD = FORWARD)
    ATC <- Reduce("+", all.diffs.weighted2)/length(all.diffs.weighted2)
    ATC <- -(ATC) # make the sign correction
    
    DID_ATE <- ifelse(length(ATC)>0,
                      (ATT * length(all.diffs.weighted) + ATC * length(all.diffs.weighted2))
                      /(length(all.diffs.weighted) + length(all.diffs.weighted2))
                      , ATT)
    
    return(list("ATE" = DID_ATE, 
                "No.controls_for_ATT" = lapply(even_smaller1, function (x) length(unique(x$V2))-1),
                "No.controls_for_ATC" = lapply(even_smaller2, function (x) length(unique(x$V2))-1)))
  }
  
}