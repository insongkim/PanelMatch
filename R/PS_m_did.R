PS_m_did <- function(x, L, FORWARD, M = M) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  colnames(x)[1:5] <- c("V2", "V1", "ps", "V4", "V5")
  x <- x[!duplicated(x[c("V2", "V1")]),]
  treated.id <- x[x$V4 == 1 & x$V1 == (max(x$V1)-FORWARD), ]$V2
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-((length(unique(x$V1)) - 1):length(unique(x$V1)))]
  matched_set <- x[x$V1 %in% timeid, ]
  
  PS_distance <- abs(tapply(matched_set$ps, matched_set$V2, mean) - mean(matched_set$ps[which(matched_set$V2 == treated.id)]))
  
  first.diff <- x[x$V2 == treated.id & x$V1 == max(x$V1), ]$V5 - 
    x[x$V2 == treated.id & x$V1 == sort(unique(x$V1))[L], ]$V5
  
  if (M < length(testid)) {
    matchid <- as.numeric(names(sort(PS_distance))[2:(M+1)])
  } else {
    matchid <- as.numeric(names(sort(PS_distance)))
  }
  
  if (FORWARD > 0) {
    second.diff <- mean(x[x$V2 %in% matchid & 
                            x$V1 == sort(unique(x$V1))[L + 1 + 1], ]$V5 -
                          x[x$V2 %in% matchid & 
                              x$V1 == sort(unique(x$V1))[L], ]$V5, na.rm = T) 
  } else {
    second.diff <- mean(x[x$V2 %in% matchid & 
                            x$V1 == sort(unique(x$V1))[L + 1], ]$V5 -
                          x[x$V2 %in% matchid & 
                              x$V1 == sort(unique(x$V1))[L], ]$V5, na.rm = T) 
  }
  return(first.diff - second.diff)
}