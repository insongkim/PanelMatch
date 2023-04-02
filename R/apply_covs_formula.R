#### helper functions for applying the formula supplied to the covs.formula argument

# accepts formula and data, creates the data used for refinement
# data has unit, time, treatment, everything else column order at this point
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
    cbind(x[, c(1, 2, 3)], model.matrix(form, tdf)[, -1])
  }
  
  t.data <- do.call(rbind, by(data, as.factor(data[, 1]), 
                              FUN = apply_formula, form = formula))
  t.data <- t.data[order(t.data[,1], t.data[,2]), ]
  rownames(t.data) <- NULL
  return(t.data)
}


merge_formula <- function(form1, form2)
{
  
  rhs1 <- trimws(unlist(strsplit(as.character(form1)[2], "\\+")))
  rhs2 <- strsplit(deparse(form2[[2]]), " \\+ ")[[1]]
  
  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  
  # put the two sides together with the amazing
  # reformulate function
  out <- reformulate(rhs)
  
  # set the environment of the formula (i.e. where should
  # R look for variables when data aren't specified?)
  environment(out) <- environment(form1)
  
  return(out)
}