MSMD_each <- function(timeid, matched_set) {
  sub_sub <- matched_set[matched_set$V1 == timeid, ]
  # D2 <- mahalanobis(sub_sub[, 4:length(sub_sub)], colMeans(sub_sub[, 4:length(sub_sub)]), cov(sub_sub[, 4:length(sub_sub)]))
  return(tryCatch(mahalanobis(sub_sub[, 4:length(sub_sub)], colMeans(sub_sub[, 4:length(sub_sub)]), cov(sub_sub[, 4:length(sub_sub)])), error=function(e) NULL))
}