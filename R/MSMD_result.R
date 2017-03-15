MSMD_result <- function(x, L, FORWARD) {
  L <- L
  FORWARD <- FORWARD
  
  testid <- unique(x$V2)
  
  timeid <- unique(x$V1)[-length(unique(x$V1))]
  
  matched_set <- x[x$V1 %in% timeid, ]
  
  MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set)
  
  MSMD <- Reduce("+", MSMDlist)/length(MSMDlist)
  
  MSMD_clone <- MSMD
  
  MSMD_clone[2] <- NA
  
  matchid <- which.min(abs(MSMD_clone - MSMD[2]))
  
  first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + FORWARD] - 
    x[x$V2 == testid[2], ]$V4[L]
  
  second.diff <- x[x$V2 == testid[matchid], ]$V4[L + 1 + FORWARD] -
    x[x$V2 == testid[matchid], ]$V4[L]
  
  return(first.diff - second.diff)
  
}
