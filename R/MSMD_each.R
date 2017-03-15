MSMD_each <- function(time_id, matched_set, testid) {
  sub_sub <- matched_set[matched_set$V1 == time_id, ]
  return(tryCatch(mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)], 
                              center = as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)]),
                              cov = cov(sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)])), 
                  error=function(e) NULL))
}