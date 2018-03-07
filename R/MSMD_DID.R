MSMD_DID <- function (L, FORWARD = 0, M = 1, time.id = "year", unit.id = "ccode", treatment, 
                          covariate, dependent, data, qoi = "ate") 
{
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
  d2 <- na.omit(data[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x)))
  dmatrix <- as.matrix(d2)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L = L, 
                                                                       F = FORWARD, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  # even_smaller1 <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), even_smaller1)
  # even_smaller1 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
  #                                         x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
  all.diffs.MSMD <- sapply(even_smaller1, MSMD_result, L = L, FORWARD = FORWARD, M = M)
  ATT <- mean(all.diffs.MSMD, na.rm = T)
  if (qoi == "att") {
    return(ATT)
  } else {
    dmatrix[, 3] <- ifelse(dmatrix[, 3] == 1, 0, 1)
    smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L = L, 
                                                                         F = FORWARD, dmatrix)), delete.NULLs)
    smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
    even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
    # even_smaller2 <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), even_smaller2)
    # even_smaller2 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
    #                                         x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
    all.diffs.MSMD2 <- sapply(even_smaller2, MSMD_result, L = L, FORWARD = FORWARD, M = M)
    ATC <- mean(all.diffs.MSMD2, na.rm = T)
    ATC <- -(ATC)
    DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(all.diffs.MSMD) + 
                                          ATC * length(all.diffs.MSMD2))/(length(all.diffs.MSMD) + 
                                                                            length(all.diffs.MSMD2)), ATT)
    return(DID_ATE)
  }
 
    
}
