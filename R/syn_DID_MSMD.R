syn_DID_MSMD <- function (L, F, time.id = "year", unit.id = "ccode", treatment, 
                          covariate, dependent, d) 
{
  L <- L
  F <- F
  FORWARD <- F
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
  
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x)))
  dmatrix <- as.matrix(d2)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L, 
                                                                       F, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  smallerlist <- lapply(smallerlist, dframelist.rb_dup)
  smallerlist <- Filter(function(x) nrow(x) > 2 * (L + F + 
                                                     1), smallerlist)
  even_smaller1 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
                                          x$V1 == unique(x$V1)[1], ]$V3 == 0, smallerlist)
  all.diffs.weighted <- lapply(even_smaller1, MSMD_result, L = L, FORWARD = FORWARD)
  ATT <- mean(unlist(all.diffs.weighted), na.rm = T)
  dmatrix[, 3] <- ifelse(dmatrix[, 3] == 1, 0, 1)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L = L, 
                                                                       F, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  smallerlist <- lapply(smallerlist, dframelist.rb_dup)
  smallerlist <- Filter(function(x) nrow(x) > 2 * (L + F + 
                                                     1), smallerlist)
  even_smaller2 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
                                          x$V1 == unique(x$V1)[1], ]$V3 == 0, smallerlist)
  all.diffs.weighted2 <- lapply(even_smaller2, MSMD_result, L = L, FORWARD = FORWARD)
  ATC <- mean(unlist(all.diffs.weighted2), na.rm = T)
  ATC <- -(ATC)
  DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(all.diffs.weighted) + 
                                        ATC * length(all.diffs.weighted2))/(length(all.diffs.weighted) + 
                                                                              length(all.diffs.weighted2)), ATT)
  return(DID_ATE)
}
