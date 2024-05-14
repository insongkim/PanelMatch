#' parse_and_prep
#'
#' accepts formula object and data, creates the data used for refinement
#' @param formula formula object specifying how to construct the data used for refinement. This is likely to be some variation of the covs.formula argument.
#' @param data data.frame object to be used to create the data needed for refinement. data has unit, time, treatment columns in that order, followed by everything else 
#' @return data.frame object with the data prepared for refinement. Data will have unit, time, treatment columns in that order, followed by everything else. 
#' @keywords internal
parse_and_prep <- function(formula, data)
{
  internal.lag <- function (x, n = 1L, default = NA)
  {
    if (n == 0) return(x)
    xlen <- length(x)
    n <- pmin(n, xlen)
    out <- c(rep(default, n), x[seq_len(xlen - n)])
    attributes(out) <- attributes(x)
    out
  }
  
  lag <- function(y, lwindow)
  {
    sapply(lwindow, internal.lag, x = y)
  }
  
  
  apply_formula <- function(x, form)
  {
    attr(form, ".Environment") <- environment()
    tdf <- model.frame(form, x, na.action = NULL)
    # unit, time, treatment columns are safely in these positions
    cbind(x[, c(1, 2, 3)], model.matrix(form, tdf)[, -1])
  }
  
  t.data <- do.call(rbind, by(data, as.factor(data[, 1]), 
                              FUN = apply_formula, form = formula))
  t.data <- t.data[order(t.data[,1], t.data[,2]), ]
  rownames(t.data) <- NULL
  return(t.data)
}

#' merge_formula
#'
#' Simple helper function for merging formula objects
#' @param form1 formula object
#' @param form2 formula object
#'
#' @return Returns a formula object, which is the concatenation of two provided formula objects. 
#'
#' @keywords internal
merge_formula <- function(form1, form2)
{
  
  rhs1 <- trimws(unlist(strsplit(as.character(form1)[2], "\\+")))
  rhs2 <- strsplit(deparse(form2[[2]]), " \\+ ")[[1]]
  
  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  
  # put the two sides together with 
  # reformulate function
  out <- reformulate(rhs)
  
  # set the environment of the formula
  environment(out) <- environment(form1)
  
  return(out)
}