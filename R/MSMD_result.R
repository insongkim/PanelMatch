MSMD_result <- function(x, L, FORWARD) {
  L <- L
  FORWARD <- FORWARD
  
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-length(unique(x$V1))]
  
  matched_set <- x[x$V1 %in% timeid, ]
  MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
  MSMD <- append(Reduce("+", MSMDlist)/length(MSMDlist), 0, after = 1)
  matchid <- order(MSMD)[2]
  
  first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + F] - 
    x[x$V2 == testid[2], ]$V4[L]
  second.diff <- x[x$V2 == testid[matchid], ]$V4[L + 1 + F] -
    x[x$V2 == testid[matchid], ]$V4[L]
  
  return(first.diff - second.diff)
  
}
