MSMD_result <- function(x, L, FORWARD, M = 3) {
  L <- L
  FORWARD <- F
  M <- M
  
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-length(unique(x$V1))]
  
  matched_set <- x[x$V1 %in% timeid, ]
  MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
  MSMD <- append(Reduce("+", MSMDlist)/length(MSMDlist), 0, after = 1)
  
  first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + FORWARD] - x[x$V2 == testid[2], ]$V4[L]
  
  if (M < length(testid)) {
    matchid <- order(MSMD)[2:(M+1)]
  } else {
    matchid <- order(MSMD)
  }
  
  second.diff <- mean(x[x$V2 %in% testid[matchid] & 
                          x$V1 == unique(x$V1)[L + 1 + FORWARD], ]$V4 -
                        x[x$V2 %in% testid[matchid] & 
                            x$V1 == unique(x$V1)[L], ]$V4, na.rm = T) 
  return(first.diff - second.diff)
}