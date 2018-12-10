#' @export
subset_expanded.data <- function(expanded_data, matched_set_indices)
{
  create_subset <- function(index, expanded.data)
  {
    return(expanded.data[index, ])
  }
  res <- lapply(matched_set_indices, FUN = create_subset, expanded.data = expanded_data)
  return(res)
}

expand.treated.ts <- function(lag, treated.ts)
{
  helper <- function(treated.t)
  {
    return(seq(from =  (treated.t - lag), to = treated.t, by = 1))
  }
  lapply(treated.ts, helper)
}

##OPTIMIZE THIS
build_maha_mats <- function(idx, ordered_expanded_data)
{
  subset.per.matchedset <- function(sub.idx)
  {
    ordered_expanded_data[sub.idx,]
  }
  unnest <- function(mset.idx)
  {
    lapply(mset.idx, subset.per.matchedset)
  }
  result <- lapply(idx, unnest)
  return(result)
}

parse_and_prep <- function(formula, data, unit.id)
{
  #browser()
  #check if formula empty or no covariates provided -- error checks
  terms <- attr(terms(formula),"term.labels")
  terms <- gsub(" ", "", terms) #remove whitespace
  lag.calls <- terms[grepl("lag(*)", terms)] #regex to get calls to lag function
  if(any(grepl("=", lag.calls))) stop("fix lag calls to use only unnamed arguments in the correct positions")
  data <- data.table::as.data.table(data) #check sorting
  # lag <- function(var.name, lag.window, data = data, unitid = unit.id)
  # {
  #   browser()
  #   tr <- data[, data.table::shift(.SD, n = lag.window, fill = NA, type = "lag", give.names = TRUE), .SDcols=var.name, by = unitid]
  #   return(tr[,-1])
  # }
  results.unmerged <- mapply(FUN = handle.calls, call.as.string = lag.calls, MoreArgs =  list(.data = data, .unitid = unit.id), SIMPLIFY = FALSE)
  names(results.unmerged) <- NULL
  full.data <- cbind(data, do.call("cbind", results.unmerged))
  return(full.data)
}

handle.calls <- function(call.as.string, .data, .unitid)
{
  lag <- function(var.name, lag.window, data = .data, unitid = .unitid) #want to make sure its being called appropriately
  {
    tr <- data[, data.table::shift(.SD, n = lag.window, fill = NA, type = "lag", give.names = TRUE), .SDcols=var.name, by = unitid]
    tr <- tr[, -1]
    colnames(tr) <- paste0(var.name, "_l", lag.window)
    return(tr)
  }
  return(eval(parse(text = call.as.string)))
}

#use col.index to determine which columns we want to "scan" for missing data
handle.missing.data <- function(data, col.index)
{
  new.names <- paste0(colnames(data)[col.index], "_NA")  
  missing.mat <- as.data.frame(apply(data[, col.index], 2, is.na)) * 1
  colnames(missing.mat) <- new.names
  data <- cbind(data, missing.mat)
  data[is.na(data)] <- 0
  return(data)
}

#OPTIMIZATION OPPORTUNITIES HERE
handle_mahalanobis_calculations <- function(mahal.nested.list)
{
  do.calcs <- function(year.df)
  { #might need to change that index number
    browser()
    cov.data <- year.df[1:(nrow(year.df) - 1), 5:ncol(year.df)]
    cov.matrix <- cov(cov.data)
    center.data <- year.df[nrow(year.df), 5:ncol(year.df)]
    #rolling without all of the error checking stuff for now? hopefully we don't need it anymore for whatever reason
    return(mahalanobis(x = cov.data, center = center.data, cov = cov.matrix))
  }
  handle_set <- function(sub.list)
  {
    return(lapply(sub.list, do.calcs))
  }
  return(lapply(mahal.nested.list, handle_set))
}








