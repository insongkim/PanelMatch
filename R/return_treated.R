#' return_treated
#'
#' \code{return_treated} obtains unit-time information for treated observations that 
#' do not have matched control units given 
#' a user-specified quantity of interest and number of lags. The user will 
#' first run \code{PanelMatch} with \code{naive} set to TRUE and \code{lag} 
#' set to 1. This will determine the total number of treated observations 
#' for a specific quantity of interest. The user will pass the returned object to 
#' argument \code{matched_sets1} in \code{return_treated}.
#' Then, the user will run \code{PanelMatch}
#' with preferred specifications, and pass the returned object to \code{matched_sets2} in \code{return_treated}
#'  \code{return_treated} will return the 
#' missing treated observations by comparing the \code{matched_sets1} and \code{matched_sets2}.
#' 
#' @param matched_sets1 a list returned by \code{PanelMatch}, 
#' with argument \code{naive} set to TRUE and \code{lag} set to 1
#' @param matched_sets2 a list returned by \code{PanelMatch}
#' @param lag number of lags specified by the user when using \code{PanelMatch}
#' @param unit.id A character string indicating the name of unit identifier
#' variable in the \code{data}. The variable must be numeric. 
#' @param time.id the unit.id specified by the user when using \code{PanelMatch}
#' @param unit.name the name of the unit 
#' @param data The data frame in data.frame class
#' @param qoi the quantity of interest
#' 
#' @return \code{return_treated} returns a matrix in which the first row contains the unit names
#' and second row contains the time.id of the treated observations that do not have matched
#' control units
#' 
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#' 
#' @export
return_treated <- function(matched_sets1, 
                           matched_sets2,
                           lag, unit.id,
                           qoi = "att", 
                           time.id,
                           unit.name,
                           data) {
  matched_sets1 <- matched_sets1
  matched_sets2 <- matched_sets2
  if (qoi == "att") {
    missing_tre <- setdiff(lapply(matched_sets1$ATT_matches, get_treated2, lag = 1), 
                           lapply(matched_sets2$ATT_matches, get_treated, lag = lag))
  } else if (qoi == "atc") {
    missing_tre <- setdiff(lapply(matched_sets1$ATC_matches, get_treated2, lag = 1), 
                           lapply(matched_sets2$ATC_matches, get_treated2, lag = lag))
  }
  
  missing_tre <- sapply(missing_tre, get_names2, data = data,
         unit.id = unit.id, time.id = time.id, 
         unit.name = unit.name)
  
  return(missing_tre)
}


#'
#' @examples \dontrun{
#' CBPSW_matchesL1_naive <- PanelMatch(lag = 1, max.lead = 0, time.id = "year",
#' unit.id = "wbcode2",
#' treatment = "dem",
#' formula = y ~ dem,
#' naive = TRUE, restricted = FALSE,
#' method = "CBPS", weighting = TRUE,
#' qoi = "ate",  M = 1, data = dem)
#' 
#' CBPSW_matchesL4 <- PanelMatch(lag = 4, max.lead = 0, time.id = "year",
#'                               unit.id = "wbcode2",
#'                               treatment = "dem",
#'                               formula = y ~ dem,
#'                               method = "Pscore", weighting = TRUE,
#'                               naive = FALSE,
#'                               qoi = "ate",  M = 1, data = dem)
#' 
#' return_treated(matched_sets1 = CBPSW_matchesL1_naive,
#'                matched_sets2 = CBPSW_matchesL4,
#'                lag = 4, unit.id = "wbcode2", qoi = "att",
#'                unit.name = "wbcode2",
#'                data = dem,
#'                time.id = "year")

#' 
#' }




