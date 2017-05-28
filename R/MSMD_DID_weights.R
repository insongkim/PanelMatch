MSMD_DID_weights <- function(L, FORWARD, time.id = "year", qoi = "ate",
                             unit.id = "ccode", M = 3,
                             treatment, covariate, dependent, d) {
  
  d <- d # set dataset
  
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
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
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = L, F = FORWARD, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  smallerlist <- lapply(smallerlist, dframelist.rb_dup)
  # subset out any dataframe that have 2 or fewer than 2 units
  smallerlist <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), smallerlist)
  
  # only focus on ATT
  even_smaller1 <- Filter(function (x) x[x$V2 == unique(x$V2)[2] & x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
  
  # d.sum1 <- length(even_smaller1) 
  
  # calibrating weights
  weights_and_dits <- lapply(even_smaller1, MSMD_weight, L = L, FORWARD = FORWARD, M = M,
                             d2 = d2, unit.id = unit.id, time.id = time.id)
  all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
  all.dits <- (lapply(weights_and_dits, function (x) x$dit))
  d2$weights_att <- Reduce("+", all.weights)
  d2$dit_att <- Reduce("+", all.dits)
  if (qoi == "att") {
    return(d2)
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
    smallerlist <- lapply(smallerlist, dframelist.rb_dup)
    # subset out any dataframe that have 2 or fewer than 2 units
    smallerlist <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), smallerlist)
    
    # only focus on ATC
    even_smaller2 <- Filter(function (x) x[x$V2 == unique(x$V2)[2] & x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
    
    # d.sum2 <- length(even_smaller2) 
    
    # calibrating weights
    weights_and_dits <- lapply(even_smaller2, MSMD_weight, L = L, FORWARD = FORWARD, M = M,
                               d2 = d2, unit.id = unit.id, time.id = time.id)
    all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    d2$weights_atc <- Reduce("+", all.weights)
    d2$dit_atc <- Reduce("+", all.dits)
    
    return(d2)
  }
  
}