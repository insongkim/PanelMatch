#' get_set_treatment_effects
#'
#' Calculates the treatment effect size at the matched set level
#'
#'
#' Calculate the size of treatment effects for each matched set.
#' @param pm.obj an object of class \code{PanelMatch}
#' @param data data.frame with the time series cross sectional data used for matching, refinement, and estimation
#' @param lead integer (or integer vector) indicating the time period(s) in the future for which the treatment effect size will be calculated. Calculations will be made for the period t + lead, where t is the time of treatment. If more than one lead value is provided, then calculations will be performed for each value.
#' @return a list equal in length to the number of lead periods specified to the \code{lead} argument. Each element in the list is a vector of the matched set level effects.
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' # create subset of data for simplicity
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "ps.match",
#'                          data = dem.sub, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
#'                          placebo.test = FALSE)
#' set.effects <- get_set_treatment_effects(pm.obj = PM.results, data = dem.sub, lead = 0)
#'
#'
#' @export

get_set_treatment_effects <- function(pm.obj, data, lead)
{
  return(lapply(lead, calculate_set_effects, 
                pm.obj = pm.obj, data.in = data))
  
}

# Helper functions for calculating set level effects
calculate_set_effects <- function(pm.obj, data.in, lead)
{
  if (identical(attr(pm.obj, "qoi"), "att"))
  {
    msets <- pm.obj[["att"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
    
  } else if (identical(attr(pm.obj, "qoi"), "atc"))
  {
    msets <- pm.obj[["atc"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
  } else if (identical(attr(pm.obj, "qoi"), "ate"))
  {
    msets <- pm.obj[["att"]]
    msets.atc <- pm.obj[["atc"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
    
  } else if (identical(attr(pm.obj, "qoi"), "art"))
  {
    msets <- pm.obj[["art"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
  } else {
    stop("invalid qoi")
  }
  
  rownames(data.in) <- paste0(data.in[, id.var], ".", data.in[, t.var])
  # unlike in PanelEstimate(), we calculate the set level effects using brute force estimator
  # get_ind_effects implements brute force approach. Point estimates from both methods should match
  get_ind_effects <- function(mset, data_in, 
                              lead.val,
                              mset.name,
                              outcome, 
                              use.abs.value = FALSE,
                              is.atc = FALSE)
  {
    
    if ( identical(length(mset), 0L)) return(NA)
    t.val <- as.numeric(sub(".*\\.", "", mset.name))
    id.val <- as.numeric(sub("\\..*", "", mset.name))
    
    
    past.lookups <- paste0(mset, ".", (t.val - 1))
    future.lookups <- paste0(mset, ".", (t.val + lead.val))
    
    t.past.lookup <- paste0(id.val, ".", (t.val - 1))
    t.future.lookup <- paste0(id.val, ".", (t.val + lead.val))
    
    control.diffs <- 
      data_in[future.lookups, outcome] - data_in[past.lookups, outcome]
    
    treat.diff <- 
      data_in[t.future.lookup, outcome] - data_in[t.past.lookup, outcome]
    
    if (is.atc)
    {
      ind.effects <- sum(attr(mset, "weights") * control.diffs) - treat.diff
    } else {
      ind.effects <- treat.diff - sum(attr(mset, "weights") * control.diffs)
    }
    denom <- attr(mset, "treatment.change")
    if (use.abs.value) denom <- abs(denom)
    if (is.atc) denom <- 1
    ind.effects <- ind.effects / denom
    return(ind.effects)
    
  }
  
  
  if ( identical(attributes(pm.obj)[["qoi"]], "att"))
  { #using simplify = TRUE because we should always expect a vector, so nothing unexpected should happen
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]]),
                      SIMPLIFY = TRUE)
    
    return(effects)
  } else if (identical(attributes(pm.obj)[["qoi"]], "atc"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]],
                                      is.atc = TRUE),
                      SIMPLIFY = TRUE)
    
    return(effects)
  } else if (identical(attr(pm.obj, "qoi"), "ate"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]]),
                      SIMPLIFY = TRUE)
    
    
    effects.atc <- mapply(FUN = get_ind_effects,
                          mset = msets.atc,
                          mset.name = names(msets.atc),
                          MoreArgs = list(lead.val = lead,
                                          data_in = data.in,
                                          outcome = attributes(pm.obj)[["outcome.var"]]),
                          SIMPLIFY = TRUE)
    
    
    return(list(att = effects,
                atc = effects.atc))
    
  } else if (identical(attr(pm.obj, "qoi"), "art"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]],
                                      use.abs.value = TRUE),
                      SIMPLIFY = TRUE)
    return(effects)
  } else
  {
    stop("invalid qoi")
  }
}
