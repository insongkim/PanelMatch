#' build_ps_data
#' 
#' @param idxlist 
#' @param data data.frame object with the data
#' @param lag see PanelMatch() documentation
#'
#' @return Returns a list of length equal to the number of matched sets. Each item is a data frame and each data frame contains information at time = t + 0 for each treated unit and their corresponding controls.
#' @keywords internal
build_ps_data <- function(idxlist, data, lag)
{
  
  obtain.t.rows <- function(idx)
  {
    return(idx[length(idx)])
  }
  unnest <- function(subidxlist,  lag)
  {
    temp <- sapply(subidxlist[[lag + 1]], obtain.t.rows)
    return(data.frame(data[temp, ]))
  }
  results <- lapply(idxlist, unnest, lag = lag)
  return(results)
}


#' find_ps
#'
#' @param sets matched sets
#' @param fitted.model Result of a fitted (CB) PS model call
#'
#' @return Returns a list of data frames with propensity score weights for each unit in a matched set. Each element in the list is a data frame which corresponds to a matched set of 1 treatment and all matched control units
#' @keywords internal
find_ps <- function(sets, fitted.model)
{
  
  apply_formula <- function (x, B)
  {
    xx <- cbind(1, as.matrix(x[, 4:ncol(x)]))
    x[, (ncol(x) + 1)] <- 1 - 1/(1+exp(xx %*% B))
    names(x)[ncol(x)] <- "ps"
    return(x[, c(1:3, ncol(x))])
  }
  sets_with_ps <- lapply(sets, apply_formula, B = fitted.model$coefficients)
  return(sets_with_ps)
}