delete.NULLs <- function(x.list) {
  x.list[unlist(lapply(x.list, length) != 0)]
}
dframelist.rb_dup <- function(x) {
  test <- data.table::rbindlist(lapply(x, as.data.frame))
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


# 
# merge.formula <- function(form1, form2, ...){
#   
#   # get character strings of the names for the responses 
#   # (i.e. left hand sides, lhs)
#   lhs1 <- deparse(form1[[2]])
#   lhs2 <- deparse(form2[[2]])
#   if(lhs1 != lhs2) stop('both formulas must have the same response')
#   
#   # get character strings of the right hand sides
#   rhs1 <- strsplit(deparse(form1[[3]]), " \\+ ")[[1]]
#   rhs2 <- strsplit(deparse(form2[[3]]), " \\+ ")[[1]]
#   
#   # create the merged rhs and lhs in character string form
#   rhs <- c(rhs1, rhs2)
#   lhs <- lhs1
#   
#   # put the two sides together with the amazing 
#   # reformulate function
#   out <- reformulate(rhs, lhs)
#   
#   # set the environment of the formula (i.e. where should
#   # R look for variables when data aren't specified?)
#   environment(out) <- parent.frame()
#   
#   return(out)
# }

# this is how you get the addition operator working for formulas
Ops.formula <- function(e1, e2){
  FUN <- .Generic
  if(FUN == '+'){
    out <- merge(e1, e2)
    environment(out) <- parent.frame()
    return(out)
  }
  else stop('can not yet subtract formula objects')
}