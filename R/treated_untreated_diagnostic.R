#' Compare pre-treatment covariate trajectories of matched and unmatched
#' treated observations
#'
#' For a fitted \code{PanelMatch} object, identifies treated unit-times that
#' did and did not receive at least one matched control unit, and compares
#' their average pre-treatment covariate values over the matching lag window.
#' This serves as a diagnostic for whether treated observations that
#' successfully matched are systematically different from those that did
#' not, which has implications for the generalizability of the estimated
#' treatment effect to the full set of treated observations.
#'
#' For each covariate, the comparison is computed within each treatment
#' cohort (i.e., treated unit-times sharing the same treatment time) that
#' contains both matched and unmatched treated units. Within a cohort, the
#' covariate is averaged over the full lag window for matched units and,
#' separately, for unmatched units; the difference between these two
#' averages is computed; and the resulting per-cohort differences are then
#' aggregated across all valid cohorts two ways: an unweighted average that
#' gives equal weight to each cohort regardless of size, and a weighted
#' average that weights each cohort by its number of treated units.
#'
#' This function returns a plain list (no class) -- see
#' \code{\link{diagnostic_summary}} for a formatted report built on top of
#' it.
#'
#' @param pm.object A \code{PanelMatch} object, as returned by
#'   \code{\link[PanelMatch]{PanelMatch}}.
#' @param panel.data A \code{PanelData} object corresponding to
#'   \code{pm.object}, containing the unit identifier, time identifier,
#'   and all variables named in \code{covariates}. This should be the same
#'   \code{PanelData} object passed to \code{PanelMatch()}.
#' @param covariates A character vector of covariate names (columns in
#'   \code{panel.data}) to compare between matched and unmatched treated
#'   observations.
#'
#' @return A plain list (no class) with the following components:
#'   \describe{
#'     \item{n_has_match}{Number of treated unit-times with at least one
#'       matched control.}
#'     \item{n_no_match}{Number of treated unit-times with no matched
#'       control.}
#'     \item{pct_no_match}{Percentage of all treated unit-times that had no
#'       matched control (\code{n_no_match} divided by the total number of
#'       treated unit-times).}
#'     \item{n_no_viable}{Number of unmatched treated unit-times whose
#'       treatment cohort contains no matched treated units at all, i.e.
#'       no viable basis for comparison.}
#'     \item{pct_no_viable}{Percentage of unmatched treated unit-times
#'       (\code{n_no_match}) that had no viable comparison cohort
#'       (\code{n_no_viable} divided by \code{n_no_match}).}
#'     \item{covariate_diffs}{A \code{data.frame} with one row per
#'       covariate, giving the unweighted average difference (matched
#'       minus unmatched) in pre-treatment covariate values across
#'       treatment cohorts (\code{mean_diff}), and the same difference
#'       weighted by cohort size (\code{weighted_mean_diff}). \code{NULL}
#'       if no cohort contained both matched and unmatched treated units.}
#'   }
#'
#' @examples
#' \dontrun{
#' PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4,
#'                           refinement.method = "ps.match",
#'                           match.missing = TRUE,
#'                           covs.formula = ~ tradewb,
#'                           size.match = 5, qoi = "att",
#'                           lead = 0:4,
#'                           forbid.treatment.reversal = FALSE)
#'
#' result <- compare_treated_observations(PM.results, dem.sub.panel,
#'                            covariates = c("tradewb", "y"))
#' str(result)
#' }
#'
#' @export
compare_treated_observations <- function(pm.object, panel.data, covariates) {
  
  if (!inherits(panel.data, "PanelData")) {
    stop("Please provide a PanelData object.")
  }
  
  qoi          <- attr(pm.object, "qoi")
  matched_sets <- pm.object[[qoi]]
  
  id.var <- attr(matched_sets, "id.var")
  t.var  <- attr(matched_sets, "t.var")
  lag    <- attr(matched_sets, "lag")
  
  # Parse treated unit-times from matched set names
  set_names <- names(matched_sets)
  parsed    <- do.call(rbind, strsplit(set_names, "\\."))
  
  treated_df <- data.frame(
    unit      = type.convert(parsed[, 1], as.is = TRUE),
    time      = type.convert(parsed[, 2], as.is = TRUE),
    has_match = sapply(matched_sets, length) > 0,
    stringsAsFactors = FALSE
  )
  
  # Counts of matched and unmatched treated unit-times
  n_has_match <- sum(treated_df$has_match)
  n_no_match  <- sum(!treated_df$has_match)
  
  # Percentage of all treated unit-times with no matched control
  pct_no_match <- 100 * n_no_match / nrow(treated_df)
  
  # Unmatched treated unit-times with no viable comparison cohort
  times_with_matches <- unique(treated_df$time[treated_df$has_match])
  n_no_viable        <- sum(!treated_df$has_match & !(treated_df$time %in% times_with_matches))
  
  # Percentage of unmatched treated unit-times with no viable comparison
  pct_no_viable <- if (n_no_match > 0) 100 * n_no_viable / n_no_match else NA_real_
  
  # For each treated unit-time, expand to lag window: t-1, t-2, ..., t-lag
  lag_rows <- do.call(rbind, lapply(seq_len(nrow(treated_df)), function(i) {
    data.frame(
      unit       = treated_df$unit[i],
      time_treat = treated_df$time[i],
      time_lag   = treated_df$time[i] - seq_len(lag),
      has_match  = treated_df$has_match[i]
    )
  }))
  
  # Merge in covariate values at the lag time periods
  panel_sub             <- panel.data[, c(id.var, t.var, covariates)]
  names(panel_sub)[1:2] <- c("unit", "time_lag")
  merged                <- merge(lag_rows, panel_sub, by = c("unit", "time_lag"), all.x = TRUE)
  
  # Only retain treatment cohorts that have both matched and unmatched treated units
  cohort_counts <- tapply(treated_df$has_match, treated_df$time, function(x) length(unique(x)))
  valid_cohorts <- names(cohort_counts[cohort_counts > 1])
  merged        <- merged[merged$time_treat %in% valid_cohorts, ]
  
  if (nrow(merged) == 0) {
    warning("No treatment cohorts contain both matched and unmatched treated units. Cannot compute comparison.")
    out <- list(
      n_has_match     = n_has_match,
      n_no_match      = n_no_match,
      pct_no_match    = pct_no_match,
      n_no_viable     = n_no_viable,
      pct_no_viable   = pct_no_viable,
      covariate_diffs = NULL
    )
    return(out)
  }
  
  # For each covariate: within each cohort, average over the lag window
  # separately for matched and unmatched treated units, take the difference,
  # then aggregate that difference across cohorts both unweighted (equal
  # weight per cohort) and weighted by cohort size (number of treated
  # units in the cohort).
  results <- lapply(covariates, function(cov) {
    
    cohort_diffs <- lapply(valid_cohorts, function(t) {
      t   <- type.convert(t, as.is = TRUE)
      sub <- merged[merged$time_treat == t, ]
      
      vals_with    <- sub[sub$has_match == TRUE,  cov]
      vals_without <- sub[sub$has_match == FALSE, cov]
      
      if (length(vals_with) == 0 || length(vals_without) == 0) return(NULL)
      
      # Cohort size = number of treated units (matched + unmatched),
      # recovered from unit-lag row counts divided by the lag length.
      n_units_with    <- length(vals_with)    / lag
      n_units_without <- length(vals_without) / lag
      cohort_size     <- n_units_with + n_units_without
      
      data.frame(
        diff = mean(vals_with, na.rm = TRUE) - mean(vals_without, na.rm = TRUE),
        size = cohort_size
      )
    })
    
    cohort_diffs <- do.call(rbind, cohort_diffs)
    
    data.frame(
      covariate          = cov,
      mean_diff          = mean(cohort_diffs$diff, na.rm = TRUE),
      weighted_mean_diff = stats::weighted.mean(cohort_diffs$diff, w = cohort_diffs$size, na.rm = TRUE)
    )
  })
  
  out <- list(
    n_has_match     = n_has_match,
    n_no_match      = n_no_match,
    pct_no_match    = pct_no_match,
    n_no_viable     = n_no_viable,
    pct_no_viable   = pct_no_viable,
    covariate_diffs = do.call(rbind, results)
  )
  out
}