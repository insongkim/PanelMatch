PanelHist <- function(matched_sets, qoi = "att",...) {
  if (qoi == "att") {
    hist(unlist(matched_sets$NC_ATT), ...)
  } else if (qoi == "atc") {
    hist(unlist(matched_sets$NC_ATC), ...)
  } else {
    stop("please select either att or atc for qoi")
  }
  
}

