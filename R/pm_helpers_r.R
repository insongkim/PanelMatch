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