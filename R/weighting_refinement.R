#' handle_ps_weighted
#'
#' @param just.ps.sets results of find_ps()
#' @param msets list of matched sets of treated and control observations
#' @param refinement.method string specifying the refinement method
#'
#' @return matched.set object with treated and matched control observations, with weights as determined by the specification
#' @keywords internal
handle_ps_weighted <- function(just.ps.sets, 
                               msets, 
                               refinement.method)
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


gather_msm_sets <- function(lead.data.list)
{
  number.of.sets <- sapply(lead.data.list, length)
  if (length(unique(number.of.sets)) != 1) stop("error with matched sets in msm calculations")
  number.of.sets <- unique(number.of.sets)
  long.data.lead.list <- unlist(lead.data.list, recursive = FALSE)
  
  long.weights.list <- lapply(long.data.lead.list, function(x){return(as.vector(x[, 4]))})
  multiplied.weights <-  multiply_weights_msm(long.weights.list, number.of.sets)
  reassembled.sets <- long.data.lead.list[1:number.of.sets]
  reassemble.weights <- function(set, weights)
  {
    set[, "ps"] <- weights #again this ps is misleading but for consistency with the other functions lets go with it
    return(set)
  }
  reassembled.sets <- mapply(FUN = reassemble.weights, set = reassembled.sets,
                             weights = multiplied.weights, SIMPLIFY = FALSE)
  
  return(reassembled.sets)
}

#check old Rcpp code if this does not work
multiply_weights_msm <- function(weights, number_of_sets) {
  final_weights <- vector("list", number_of_sets)
  
  for (i in seq_len(number_of_sets)) {
    base_mult <- weights[[i]]
    for (j in seq(i, length(weights), by = number_of_sets)) {
      if (j != i) {
        temp <- weights[[j]]
        base_mult <- base_mult * temp
      }
    }
    final_weights[[i]] <- base_mult
  }
  
  return(final_weights)
}