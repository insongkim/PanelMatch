#' enforce_lead_restrictions
#' check treatment and control units for treatment reversion in the lead window. Treated units must stay treated and control units must stay in control (according to the specified qoi)
#' 
#' @param matched_sets matched.set object
#' @param ordered.data parsed data as data.frame object
#' @param max.lead The largest lead value (e.g. the biggest F)
#' @param t.var string specifying the time variable
#' @param id.var string specifying the unit id variable
#' @param treatment.var string specifying the treatment variable.
#'
#' @return matched.set object with the matched sets that meet the conditions
#' @keywords internal
enforce_lead_restrictions <- function(matched_sets, 
                                      ordered.data, 
                                      max.lead, 
                                      t.var, 
                                      id.var, 
                                      treatment.var)
{
  
  ordered.data <- ordered.data[order(ordered.data[,id.var], ordered.data[,t.var]), ]
  compmat <- data.table::dcast(data.table::as.data.table(ordered.data),
                               formula = paste0(id.var, "~", t.var), value.var = treatment.var)
  ts <- as.numeric(sub(".*\\.", "", names(matched_sets)))
  tids <- as.numeric(sub("\\..*", "", names(matched_sets)))
  class(matched_sets) <- "list" #so that Rcpp::List is accurate when we pass it into cpp functions
  compmat <- data.matrix(compmat)
  
  idx <- check_treated_units_for_treatment_reversion(compmat = compmat, 
                                                     compmat_row_units = as.numeric(compmat[, 1]),
                                                     compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]),
                                                     lead = max.lead, treated_ids = tids, treated_ts = ts)
  if (all(!idx)) stop("estimation not possible: All treated units are missing data necessary for the calculations to proceed")
  if (any(!idx))
  {
    class(matched_sets) <- c("matched.set", "list") 
    #to get the matched.set subsetting with attributes
    matched_sets <- matched_sets[idx]
    ts <- ts[idx]
    
  }
  
  class(matched_sets) <- "list" # for rcpp reasons again
  ll <- check_control_units_for_treatment_restriction(compmat = compmat,
                                                      compmat_row_units = as.numeric(compmat[, 1]),
                                                      compmat_cols = as.numeric(colnames(compmat)[2:ncol(compmat)]),
                                                      lead = max.lead, sets = matched_sets, control_start_years = ts)
  # check which units need to be removed
  idx <- needs_adjustment(ll)
  class(matched_sets) <- c("matched.set", "list")
  if (any(idx))
  {
    sub.index <- ll[idx]
    sub.set <- matched_sets[idx]
    create_new_sets <- function(set, index)
    {
      return(set[index])
    }
    sub.set.new <- mapply(FUN = create_new_sets, 
                          sub.set, 
                          sub.index,
                          SIMPLIFY = FALSE)
    attributes(sub.set.new) <- attributes(sub.set)
    all.gone.counter <- sapply(sub.set.new, function(x){sum(x)})
    if (sum(all.gone.counter == 0) > 0) #case in which all the controls in a particular group were dropped
    {
      idx[all.gone.counter == 0] <- FALSE
      sub.index <- ll[idx]
      sub.set <- matched_sets[idx]
      create_new_sets <- function(set, index)
      {
        return(set[index])
      }
      sub.set.new <- mapply(FUN = create_new_sets, 
                            sub.set, 
                            sub.index, 
                            SIMPLIFY = FALSE)
      attributes(sub.set.new) <- attributes(sub.set)
    }
    if (all(sapply(sub.set.new, length) == 0)) stop('estimation not possible: none of the matched sets have viable control units due to a lack of necessary data')
    matched_sets[idx] <- sub.set.new
    matched_sets <- matched_sets[sapply(matched_sets, length) > 0]
  }
  return(matched_sets)
  
}