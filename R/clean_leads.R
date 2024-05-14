#' clean_leads
#' Function to check the lead windows in treated and control units for missing outcome data. If data is missing, remove those units from matched sets.
#' @param matched_sets matched.set object contained pre-filtered matched sets
#' @param ordered.data data.frame object to be checked for missing data. This should have been passed through data preparation functions already.
#' @param max.lead Integer specifying the biggest value of the lead window.
#' @param t.var string specifying the time id variable
#' @param id.var string specifying the unit id variable
#' @param outcome.var string specifying the outcome variable.
#' @return a cleaned/filtered matched.set object
#' @keywords internal
clean_leads <- function(matched_sets, 
                        ordered.data, 
                        max.lead, 
                        t.var, 
                        id.var, 
                        outcome.var)
{
  old.attributes <- attributes(matched_sets)[names(attributes(matched_sets)) != "names"]
  
  ts <- as.numeric(sub(".*\\.", "", names(matched_sets)))
  tids <- as.numeric(sub("\\..*", "", names(matched_sets)))
  class(matched_sets) <- "list" # For Rcpp
  
  #check to make sure we have outcome data at for each period in the lead window
  idx <- 
    check_missing_data_treated_units(subset_data = 
                                       as.matrix(ordered.data[, c(id.var, 
                                                                  t.var, 
                                                                  outcome.var)]),
                                          sets = matched_sets,
                                          tid_pairs = 
                                       paste0(ordered.data[, id.var], 
                                              ".", ordered.data[, t.var]),
                                          treated_tid_pairs = names(matched_sets),
                                          treated_ids = tids, 
                                          lead =  max.lead)
  
  if(all(!idx)) stop("estimation not possible: All treated units are missing data necessary for the calculations to proceed")
  if(any(!idx))
  {
    #to get the matched.set subsetting with attributes
    class(matched_sets) <- c("matched.set", "list") 
    matched_sets <- matched_sets[idx]
    ts <- ts[idx]
    
  }
  class(matched_sets) <- "list" # for Rcpp reasons again
  
  if(any(idx)) 
  {
    create_control_maps <- function(matched_set, time)
    {
      return(paste0(matched_set, ".", time))
    }
    
    prepped_sets <- mapply(create_control_maps, 
                           matched_set = matched_sets, 
                           time = ts, 
                           SIMPLIFY = FALSE)
    
    
    #check to make sure we have outcome data at for each period in the lead window
    #checking control units
    tpx <- 
      check_missing_data_control_units(subset_data = 
                                         as.matrix(ordered.data[, c(id.var, 
                                                                    t.var,
                                                                    outcome.var)]),
                                            sets = matched_sets,
                                            prepared_sets = prepped_sets,
                                            tid_pairs = 
                                         paste0(ordered.data[, id.var], ".", 
                                                ordered.data[, t.var]),
                                            lead =  max.lead)
    
    create_new_sets <- function(set, index)
    {
      return(set[index])
    }
    sub_sets <- mapply(FUN = create_new_sets, 
                       matched_sets, 
                       tpx, 
                       SIMPLIFY = FALSE)
    
    if(all(sapply(sub_sets, length) == 0)) stop('estimation not possible: none of the matched sets have viable control units due to a lack of necessary data')
    
    matched_sets <- sub_sets[sapply(sub_sets, length) > 0]
    
    for(idx in names(old.attributes))
    {
      attr(matched_sets, idx) <- old.attributes[[idx]]
    }
    
  }
  class(matched_sets) <- c("matched.set")
  return(matched_sets)
  
}