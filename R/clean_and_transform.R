delete.NULLs <- function(x.list) {
  x.list[unlist(lapply(x.list, length) != 0)]
}
dframelist.rb_dup <- function(x) {
  test <- rbindlist(lapply(x, as.data.frame))
  return(test[!duplicated(test), ])
}

grab_index <- function(x) {
  return(list("unit.index" = unique(x$V2), 
              "time.index" = unique(x$V1)))
}
# indexes <- lapply(matched_sets$`Matched sets for ATT`, grab_index)

gen_sets <- function(x, data, unit.id, time.id) {
  return(data[data[unit.id][[1]] %in% x$unit.index &
                data[time.id][[1]] %in% x$time.index, ])
}

# new.matched_sets <- lapply(indexes, gen_sets, data = Data.obs, unit.id = "unit",
 #                          time.id = "time")