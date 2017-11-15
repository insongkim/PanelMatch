match_brute <- function(x, lead, max.lead, lag) {
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x$V2[which(x$V3 == 1 & x$V1 == max(x$V1) - max.lead)]
  diff <- x$V4[which(x$V2 == treated.id & x$V1 == min(x$V1) + lag + lead)] -
   sum(x$V4[which(x$V2 != treated.id & x$V1 == min(x$V1) + lag + lead)] * 
          x$w.weight[which(x$V2 != treated.id & x$V1 == min(x$V1) + lag + lead)])
  return(diff)
}

brute_matching <- function(Matches, lead, max.lead, lag) {
  return(mean(sapply(Matches, match_brute, lead = lead, lag = lag,
              max.lead = max.lead)))
}


#      
# brute_matching(Matches = Matches_MahaL4_2$ATT_matches, 
#                lead = 1,
#                max.lead = 2,
#                lag = 4)
# 
# match_brute(x = Matches_MahaL4_2$ATT_matches[[2]], lead = 0, max.lead = 0,
#             lag = 4)
# 
# result <- rep(NA, length(Matches_MahaL4_2$ATT_matches))
# for (i in 1:length(Matches_MahaL4_2$ATT_matches)) {
#   result[i] <- match_brute(x = Matches_MahaL4_2$ATT_matches[[i]], lead = 2, max.lead = 2,
#               lag = 4)
# 
# }