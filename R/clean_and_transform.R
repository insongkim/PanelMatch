delete.NULLs <- function(x.list) {
  x.list[unlist(lapply(x.list, length) != 0)]
}
dframelist.rb_dup <- function(x) {
  test <- rbindlist(lapply(x, as.data.frame))
  return(test[!duplicated(test), ])
}