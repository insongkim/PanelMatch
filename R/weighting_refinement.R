########################################################################
######################## calculates weights for units for weighting-based
######################## refinement methods ########################
########################################################################
handle_ps_weighted <- function(just.ps.sets, msets, refinement.method)
{
  handle_set <- function(set)
  {
    control.ps.set <- set[1:(nrow(set) - 1), ncol(set)]
    if(length(control.ps.set) == 1)
    {
      return(1)
    }
    vec.ratio <- control.ps.set / (1 - control.ps.set) #just for clarity
    if(sum(vec.ratio) == 0)
    {
      wts <- rep(1 / length(control.ps.set), length(control.ps.set))
    }
    else
    {
      wts <- ( vec.ratio ) / sum( vec.ratio )
    }
    return(as.vector(wts))
  }
  wts <- lapply(just.ps.sets, handle_set)
  for(i in 1:length(msets))
  {
    names(wts[[i]]) <- msets[[i]]
    attr(msets[[i]], "weights") <- wts[[i]]
  }
  attr(msets, "refinement.method") <- refinement.method
  return(msets)
}