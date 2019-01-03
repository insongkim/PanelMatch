## anything related to matched set objects
#' @export
matched_set <- function(matchedsets, id, t, L, t.var, id.var, treated.var)
{
  if(length(matchedsets) != length(id) | length(matchedsets) != length(t) | length(id) != length(t))
  {
    stop("Number of matched sets, length of t, length of id specifications do not match up")
  }
  names(matchedsets) <- paste0(id, ".", t)
  class(matchedsets) <- c("matched.set", "list")
  attr(matchedsets, 'refinement.method') <- NULL
  attr(matchedsets, "lag") <- L
  attr(matchedsets, "t.var") <- t.var
  attr(matchedsets, "id.var" ) <- id.var
  attr(matchedsets, "treated.var") <- treated.var
  return(matchedsets)
}

#' @export
summary.matched.set <- function(set, verbose = T)
{
  Lengthcol <- sapply(set, length)
  temp <-unlist(strsplit(names(set), split = ".", fixed = TRUE))
  ids <- temp[c(T,F)]
  ts <- temp[c(F,T)]
  df <- data.frame(i = ids, t = ts, matched.set.size = Lengthcol)
  colnames(df)[1:2] <- c(attr(set, "id.var"), attr(set, "t.var"))
  rownames(df) <- NULL
  #class(df) <- c("overview.matched.set", "data.frame")
  if(verbose)
  {
    summary.result <- list()
    summary.result$overview <- df
    summary.result$set.size.summary <- summary(Lengthcol)
    summary.result$number.of.treated.units <- length(set)
    summary.result$num.units.empty.set <- sum(Lengthcol == 0)
    summary.result$lag <- attr(set, "lag")
    #class(summary.result) <- "summary.matched.set"
    return(summary.result)
  }
  else
  {
    #class(df) <- "summary.matched.set"
    return(df)
  }
}
#' @export
plot.matched.set <- function(set, border = NA, col = "grey", xlim = NULL, ylim = NULL, ylab = "", xlab ="" , lwd = NULL,
                             main = "Distribution of matched set sizes", ...)
{
  
    lvec <- sapply(set, length)
    if(is.null(xlim))
    {
      xlim = c(0, length(lvec))
    }
    if(is.null(ylim))
    {
      ylim = c(0, max(lvec, max(sum(lvec == 0))))
    }
    if(is.null(lwd))
    {
      lwd = 4
    }
    hist(x = lvec, freq = TRUE, border = border, col = col, xlim = xlim, ylim = ylim, ylab = ylab, xlab = xlab, main = main, ...)
    lines(x = c(0,0), y = c(0, length(rep(0, sum(lvec == 0) )) ), col = "red", lwd = lwd)
  
  
}
#' @export
extract.set <- function(set, id, t)
{
  strid <- paste0(id, ".", t)
  subset <- set[names(set) == strid]
  if(length(subset) != 1) stop('t,id pair invalid')
  #if(class(subset) != 'matched.set') stop("Problem extracting individual matched set")
  class(subset) <- "matched.set"
  return(subset)
}

#' @export
print.matched.set <- function(set, verbose = F)
{
  if(verbose) 
  {
    class(set) <- "list"
    print(set)
  }
  
  else {
    print(summary(set, verbose = F))
  }
}

#' @export
`[.matched.set` <- function(x, i, j = NULL, drop = NULL)
{

  if(!is.null(j)) stop("matched.set object is a list.")
  class(x) <- "list"
  temp <- x[i]
  attr(temp, "lag") <- attr(x, "lag")
  attr(temp, "refinement.method") <- attr(x, "refinement.method")
  attr(temp, "t.var") <- attr(x, "t.var")
  attr(temp, "id.var" ) <- attr(x, "id.var" )
  attr(temp, "treated.var") <- attr(x, "treated.var")
  attr(temp, "distances") <- attr(x, "distances")
  attr(temp, "max.match.size") <- attr(x, "max.match.size")
  attr(temp, "covs.formula") <- attr(x, "covs.formula")
  attr(temp, "match.missing") <- attr(x, "match.missing")
  class(temp) <- "matched.set"
  
  return(temp)
}

#' #' @export
# `[.overview.matched.set` <- function(x, i, j = NULL, drop = NULL)
#' {
#'   
#'   # if(is.null(sets)) stop("please specify matched sets")
#'   # class(x) <- "data.frame"
#'   # ids <- x[i, c("i")]
#'   # ts <- x[i, c("t")]
#'   # 
#'   # attr(temp, "lag") <- attr(x, "lag")
#'   # attr(temp, "refinement.method") <- attr(x, "refinement.method")
#'   # class(temp) <- "matched.set"
#'   # return(temp)
#' }
#' 

# `[<-.matched.set` <- function(set)
# {
#   stop("not implemented yet, use a list")
# }
