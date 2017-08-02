PanelHist <- function(matched_sets, ...) {
  if ("NC_ATT" %in% names(matched_sets)) {
    hist(unlist(matched_sets$NC_ATT), ...)
  }
  
  if ("NC_ATC" %in% names(matched_sets)) {
    hist(unlist(matched_sets$NC_ATC), ...)
  }
  
}

