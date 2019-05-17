# all functions taken and modified from plm package
make.pconsecutive.indexes <- function(x, index, balanced = FALSE, ...) {
  # make.pconsecutive.indexes: helper function, not exported
  # returns list with 3 elements:
  #   1 "consec_index":        consecutive data.frame to serve as the new index data.frame in other functions, 
  #   2 "NArows_former_index": information about dropped lines (logical vector with length of original data)
  #   3 "has_fancy_rownames":  logical whether fancy row.names were used in original data (can only be TRUE for pdata.frame or pseries)
  
  if (inherits(x, "pdata.frame") || inherits(x, "pseries")) { # this whole branch can go
    stop("data must be in data.frame object")
  }
  if (inherits(x, "data.frame") && !inherits(x, "pdata.frame")) {
    # x is a data.frame, but no pdata.frame
    pdataframe_or_pseries <- FALSE
    has_fancy_rownames    <- FALSE
    index_orig <- x[ , index]
    id_orig    <- index_orig[[1]]
    times_orig <- index_orig[[2]]
    id_orig_typeof    <- typeof(id_orig)
    times_orig_typeof <- typeof(times_orig)
    rownames_mode <- mode(attr(x, "row.names"))
    rownames_typeof <- typeof(attr(x, "row.names"))
    
  }
  
  df_index <- data.frame(id = id_orig, times = times_orig)
  
  # remove any rows with NA in id or time variable as it is impossible to
  # infer their values, thus: drop them
  is_NA <- is.na(id_orig) | is.na(times_orig)
  df_index <- df_index[!is_NA, ]
  
  n_id_orig <- length(unique(id_orig))
  
  if (!balanced) { 
    min_values <- by(df_index[ , "times"], df_index[ , "id"], min)
    max_values <- by(df_index[ , "times"], df_index[ , "id"], max)
    
    times_filled_list <- sapply(seq_len(n_id_orig), function(i) {
      seq(from = min_values[i], to = max_values[i], by = 1)
    }, simplify = FALSE)
    
  } else {
    min_value <- min(df_index[, "times"])
    max_value <- max(df_index[, "times"])
    
    times_filled_list <- sapply(seq_len(n_id_orig), function(i) {
      seq(from = min_value, to = max_value, by = 1)
    }, simplify = FALSE, USE.NAMES = FALSE)
  }
  
  times_filled_vector <- unlist(times_filled_list)
  id_times <- sapply(times_filled_list, length) # lengths (with an "s") would be more efficient, but requires R >= 3.2
  id_filled_vector <- unlist(mapply(rep, unique(id_orig), id_times, SIMPLIFY = FALSE))
  # SIMPLIFY = FALSE => always return list
  
  df_index_filled <- data.frame(id = id_filled_vector, times = times_filled_vector)
  names(df_index_filled)[1:2] <- names(index_orig)[1:2] # set original index names
  
  
  if (pdataframe_or_pseries) {
    df_index_filled[ , 1] <- as.factor(df_index_filled[ , 1])
    df_index_filled[ , 2] <- as.factor(df_index_filled[ , 2])
    class(df_index_filled) <- c("pindex", class(df_index_filled))
  } else {
    if (typeof(df_index_filled[ , 1]) != id_orig_typeof)    { mode(df_index_filled[ , 1]) <- id_orig_typeof    }
    if (typeof(df_index_filled[ , 2]) != times_orig_typeof) { mode(df_index_filled[ , 2]) <- times_orig_typeof }
  }
  
  # restore mode of row.names attribute
  # [was changed by above code due to some simplification by R's standard behaviour]
  mode(attr(df_index_filled, "row.names")) <- rownames_typeof
  
  res <- list(consec_index         = df_index_filled,
              NArows_former_index  = is_NA,
              has_fancy_rownames   = has_fancy_rownames)
  
  return(res)
} ### END: make.pconsecutive.indexes


make.pconsecutive.data.frame <- function(x, balanced = FALSE, index = NULL, ...){
  # if not NULL, index is must be character of length 2
  if (!is.null(index) & length(index) != 2)
    stop("if argument 'index' is not NULL, 'index' needs to specify
         'individual' and 'time' dimension for make.pconsecutive to work on a data.frame")
  
  # assume first two columns to be the index vars
  if (is.null(index)) index_orig_names <- names(x)[1:2]
  else index_orig_names <- index
  
  list_ret_make_index <- make.pconsecutive.indexes(x, index_orig_names, balanced = balanced, ...)
  index_df_filled    <- list_ret_make_index[["consec_index"]]
  NArows_old_index   <- list_ret_make_index[["NArows_former_index"]]
  has_fancy_rownames <- list_ret_make_index[["has_fancy_rownames"]]
  
  # silently drop rows with NA in either individual or time variable of original index
  x <- x[!NArows_old_index, ]
  
  index_df_filled_plus_x <- merge(index_df_filled, x, by.x = names(index_df_filled)[1:2],
                                  by.y = index_orig_names,
                                  all.x = TRUE)
  
  # restore mode of row.names attribute [was changed by above code due to some simplification as R's standard behaviour]
  mode(attr(index_df_filled_plus_x, "row.names")) <- typeof(attr(index_df_filled, "row.names"))
  
  # restore original order of columns, esp. place index vars at original position
  index_df_filled_plus_x <- index_df_filled_plus_x[ , names(x)]
  
  return(index_df_filled_plus_x)
} ### END: make.pconsecutive.data.frame


make.pconsecutive <- function(x, ...){
  UseMethod("make.pconsecutive")
}


############# make.pbalanced #############
## make.pbalanced.* methods make the input balanced (but not consecutive).
## It does so by either 
## balance.type = "fill": filling in only those missing time periods are 
##                        introduced that are present for at least one individual
##                        (union of time periods)
##
## balance.type = "shared.times": remove all observations with time periods
##                                not shared among all individuals
##                                (keep intersect of time periods)
##
##                "shared.individuals": drop individuals which don't have all time periods
##                                      (symmetric to "shared.times")

make.pbalanced.data.frame <- function(x, balance.type = c("fill", "shared.times", "shared.individuals"), index = NULL, ...) {
  # NB: for data.frame interface: the data is also sorted as stack time series
  
  if (length(balance.type) == 1 && balance.type == "shared") {
    # accept "shared" for backward compatibility
    balance.type <- "shared.times"
    warning("Use of balanced.type = 'shared' discouraged, set to 'shared.times'")
  }
  balance.type <- match.arg(balance.type)
  
  ## identify index of data.frame  
  # if not NULL, index is must be character of length 2
  if (!is.null(index) & length(index) != 2)
    stop("if argument 'index' is not NULL, 'index' needs to specify
             'individual' and 'time' dimension for make.pconsecutive to work on a data.frame")
  
  # assume first two columns to be the index vars
  if (is.null(index)) index_orig_names <- names(x)[1:2]
  else index_orig_names <- index
  
  index_df <- x[ , index_orig_names]
  
  switch(balance.type,
         "fill" = {
           x_consec_bal <- make.pconsecutive(x, index = index_orig_names, balanced = TRUE)
           
           # delete time periods that were not present for any individual, but introduced by
           # making data consecutive
           # result: no time periods are added that are not present for at least one individual
           times_present_orig <- x_consec_bal[ , index_orig_names[2]] %in% unique(index_df[[2]])
           result <- x_consec_bal[times_present_orig , ]},
         
         "shared.times" = {
           keep <- intersect_index(index_df, "time")
           result <- x[keep, ]},
         
         "shared.individuals" = {
           keep <- intersect_index(index_df, "individual")
           result <- x[keep, ]
         })
  return(result)
} ## END make.pbalanced.data.frame

make.pbalanced <- function(x, balance.type = c("fill", "shared.times", "shared.individuals"), ...) {
  UseMethod("make.pbalanced")
}

# helper function: returns logical vector which rows/entries to keep
#                  when balance.type = "shared.times" or "shared.individuals"
#                  (intersect of all time periods or individuals)
intersect_index <- function(index, by) {
  # intersect() is defined on vectors (not factors)
  #  -> convert respective index to character before
  
  switch(by,
         "time" = {
           id <- index[[1]]
           time <- as.character(index[[2]])
         },
         "individual" = {
           id <- index[[2]]
           time <- as.character(index[[1]])
         })
  
  times_by_ids <- split(time, id)
  common_times <- Reduce(intersect, times_by_ids) 
  keep_entries <- time %in% common_times
  return(keep_entries)
}

check_time_data <- function(data, time.id)
{
  if(class(data[, time.id]) != "integer") stop("time data is not integer")
  u.times <- unique(data[, time.id])
  increase.by.one <- all(seq(min(u.times), max(u.times), by = 1) %in% u.times)
  if(increase.by.one)
  {
    return(TRUE)
  }
  else
  {
    stop("integer representation of time data has problematic gaps, as it does not increase by one. Perhaps time data for observations is irregular/not uniform across units?")
  }
}


