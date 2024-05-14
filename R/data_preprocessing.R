#' check_time_data
#'
#' enforces the requirements for time data, with some reasonable defaults
#' @description Time data should be consecutive integers: When it is not, try to convert it as best we can or throw an error. If function does not fail, returns the data as data frame object, either processed or not as appropriately
#' @param data data.frame object.
#' @param time.id string specifying the time id variable.
#' @return data.frame object with the data. If function throws error, nothing is returned.
#' @keywords internal
check_time_data <- function(data, time.id)
{
  if (!class(data[, time.id]) %in% c("numeric", "integer"))
  { # if data is not numeric or integer, just throw an error. too hard to convert
    stop("Time data not consecutive integer")
  }
  is.not.int <- (!inherits(data[, time.id], "integer")) 
  u.times <- unique(data[, time.id])
  increase.by.one <- all(seq(min(u.times), max(u.times), by = 1) %in% u.times)
  if(is.not.int || !increase.by.one) 
  { # if we can reasonably perform some kind of conversion (e.g. numeric data), do so
    warning("Data is not consecutive integer: Attempting automatic conversion, which may cause undefined behavior")
    # Assuming sorted data here
    data[, time.id] <- as.integer(as.factor(data[,time.id]))
    return(data)
  }
  
  if ("numeric" %in% class(data[, time.id]))
  {
    warning("time data is numeric: attempting to convert to integer")
    data[, time.id] <- as.integer(data[,time.id])
    return(data)
  }
  
  # at this point, if not integer, it is unclear what the problem is. So, throw an error
  if(is.not.int)
  {
    stop("time data is not integer")
  } else
  {
    return(data)
  }
  
}

#' handle_missing_data
#'
#' Tags missing data
#' @description use col.index to determine which columns we want to "scan" for missing data. Note that in earlier points in the code, we rearrange the columns and prepare the data frame such that cols 1-4 are bookkeeping (unit id, time id, treated variable, unlagged outcome variable) and all remaining columns are used in the calculations after going through parse_and_prep function, so col.index should usually be 5:ncol(data). In practice, this function just looks over the data in the specified columns in the "data" data frame for missing data. Then it creates columns with indicator variables about the missingness of those variables: 1 for missing data, 0 for present
#' @param data data.frame object.
#' @param col.index numeric vector specifying which columns to inspect
#' @return data.frame object with the data and the missingness indicators described above.
#' @keywords internal
handle_missing_data <- function(data, col.index)
{
  
  new.names <- paste0(colnames(data)[col.index], "_NA")
  if(length(col.index) == 1)
  {
    new.col <- is.na(data[, col.index]) * 1
    data <- cbind(data, new.col)
    colnames(data)[ncol(data)] <- new.names
  }
  else
  {
    missing.mat <- as.data.frame(apply(data[, col.index], 2, is.na)) * 1
    colnames(missing.mat) <- new.names
    data <- cbind(data, missing.mat)
  }
  
  data[is.na(data)] <- 0
  return(data)
}
