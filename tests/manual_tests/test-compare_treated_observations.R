## Tests for compare_treated_observations()
##
## Covers: a fully deterministic, hand-computable case (built from a
## synthetic PanelMatch-like object rather than real matching, so every
## expected number can be verified by inspection); the "no valid cohorts"
## edge case (all treated units matched, or all unmatched); and a smoke test
## against real PanelMatch()/dem data for structural sanity.

test_that("compare_treated_observations() produces hand-computable results on a synthetic example", {
  ## --- Build a minimal, fully controlled PanelMatch-like object ---
  ## Two treated units, both treated at time 5:
  ##   unit 1: matched (to control unit 3, whose data isn't actually used --
  ##           the comparison is between treated units' own pre-treatment
  ##           trajectories, matched vs unmatched, not treated vs control)
  ##   unit 2: unmatched
  fake.matched.sets <- list(
    "1.5" = 3L,          # matched
    "2.5" = integer(0)   # unmatched
  )
  attr(fake.matched.sets, "id.var") <- "unit"
  attr(fake.matched.sets, "t.var")  <- "time"
  attr(fake.matched.sets, "lag")    <- 2
  
  fake.pm <- list(att = fake.matched.sets)
  attr(fake.pm, "qoi") <- "att"
  class(fake.pm) <- "PanelMatch"
  
  ## Pre-treatment covariate "x" at t-1 (time 4) and t-2 (time 3):
  ##   unit 1 (matched):   x = 12 at t-1, 10 at t-2
  ##   unit 2 (unmatched): x = 6  at t-1, 5  at t-2
  panel.data <- data.frame(
    unit = c(1, 1, 2, 2),
    time = c(4, 3, 4, 3),
    x    = c(12, 10, 6, 5)
  )
  class(panel.data) <- c("PanelData", "data.frame")
  
  result <- compare_treated_observations(fake.pm, panel.data, covariates = "x")
  
  # --- counts: hand-computable ---
  expect_equal(result$n_has_match, 1)
  expect_equal(result$n_no_match, 1)
  expect_equal(result$pct_no_match, 50)   # 1 of 2 treated unit-times unmatched
  expect_equal(result$n_no_viable, 0)     # the unmatched unit's cohort (time 5)
  # does contain a matched unit, so it's
  # not "no viable comparison"
  expect_equal(result$pct_no_viable, 0)
  
  # --- covariate_diffs: hand-computable ---
  # matched unit's mean x over the lag window:   mean(12, 10) = 11
  # unmatched unit's mean x over the lag window: mean(6, 5)   = 5.5
  # diff (matched - unmatched) = 11 - 5.5 = 5.5
  # only one treatment cohort (time 5), so weighted and unweighted means
  # are identical and both equal to that single cohort's diff.
  expect_equal(nrow(result$covariate_diffs), 1)
  expect_equal(result$covariate_diffs$covariate, "x")
  expect_equal(result$covariate_diffs$mean_diff, 5.5)
  expect_equal(result$covariate_diffs$weighted_mean_diff, 5.5)
})

test_that("compare_treated_observations() correctly pools multiple units within a group (hand-computable, multi-unit)", {
  ## Same idea as the single-unit test above, but with two treated units
  ## in the matched group and two in the unmatched group, all in the same
  ## treatment cohort (time 5). This is the case that distinguishes
  ## "average of per-unit means" from "mean pooled across all unit-periods
  ## in the group" -- compare_treated_observations() does the latter, so
  ## the expected numbers below are computed the same way (pooling all
  ## individual lag-period observations within each group before taking
  ## the mean), not by first averaging each unit's own two periods.
  fake.matched.sets <- list(
    "1.5" = 5L,          # matched
    "2.5" = 6L,          # matched
    "3.5" = integer(0),  # unmatched
    "4.5" = integer(0)   # unmatched
  )
  attr(fake.matched.sets, "id.var") <- "unit"
  attr(fake.matched.sets, "t.var")  <- "time"
  attr(fake.matched.sets, "lag")    <- 2
  
  fake.pm <- list(att = fake.matched.sets)
  attr(fake.pm, "qoi") <- "att"
  class(fake.pm) <- "PanelMatch"
  
  ## Pre-treatment covariate "x" at t-1 (time 4) and t-2 (time 3):
  ##   unit 1 (matched):   x = 12 at t-1, 10 at t-2
  ##   unit 2 (matched):   x = 14 at t-1, 12 at t-2
  ##   unit 3 (unmatched): x = 6  at t-1, 5  at t-2
  ##   unit 4 (unmatched): x = 8  at t-1, 7  at t-2
  panel.data <- data.frame(
    unit = c(1, 1, 2, 2, 3, 3, 4, 4),
    time = c(4, 3, 4, 3, 4, 3, 4, 3),
    x    = c(12, 10, 14, 12, 6, 5, 8, 7)
  )
  class(panel.data) <- c("PanelData", "data.frame")
  
  result <- compare_treated_observations(fake.pm, panel.data, covariates = "x")
  
  # --- counts: hand-computable ---
  expect_equal(result$n_has_match, 2)
  expect_equal(result$n_no_match, 2)
  expect_equal(result$pct_no_match, 50)   # 2 of 4 treated unit-times unmatched
  expect_equal(result$n_no_viable, 0)     # both unmatched units' cohort (time 5)
  # does contain matched units
  expect_equal(result$pct_no_viable, 0)
  
  # --- covariate_diffs: hand-computable ---
  # matched group's pooled values:   12, 10, 14, 12 -> mean = 48 / 4 = 12
  # unmatched group's pooled values: 6, 5, 8, 7      -> mean = 26 / 4 = 6.5
  # diff (matched - unmatched) = 12 - 6.5 = 5.5
  #
  # note this is NOT the same as averaging each unit's own mean first
  # (unit 1 mean = 11, unit 2 mean = 13 -> average = 12; that happens to
  # give the same matched-group mean here only because both matched units
  # have the same number of lag periods -- the pooling, not per-unit
  # averaging, is what the function actually does).
  #
  # only one treatment cohort (time 5), so weighted and unweighted means
  # are again identical here; this test does not distinguish the two --
  # a separate multi-cohort test would be needed for that.
  expect_equal(nrow(result$covariate_diffs), 1)
  expect_equal(result$covariate_diffs$covariate, "x")
  expect_equal(result$covariate_diffs$mean_diff, 5.5)
  expect_equal(result$covariate_diffs$weighted_mean_diff, 5.5)
})

test_that("compare_treated_observations() weighted and unweighted mean_diff differ across cohorts of different sizes", {
  ## Two treatment cohorts (different treatment times), deliberately of
  ## different sizes and with different per-cohort diffs, so the
  ## unweighted (equal-weight-per-cohort) and weighted (weighted by
  ## cohort size) aggregates diverge:
  ##
  ##   cohort at time 5: 1 matched unit, 1 unmatched unit -> size 2, diff = 10
  ##   cohort at time 6: 2 matched units, 2 unmatched units -> size 4, diff = 2
  ##
  ## unweighted mean_diff          = mean(10, 2)                      = 6
  ## weighted_mean_diff            = weighted.mean(c(10, 2), c(2, 4)) = 28/6 ~ 4.667
  fake.matched.sets <- list(
    "1.5" = 10L,          # matched,   cohort (time) 5
    "2.5" = integer(0),   # unmatched, cohort (time) 5
    "3.6" = 11L,          # matched,   cohort (time) 6
    "4.6" = 12L,          # matched,   cohort (time) 6
    "5.6" = integer(0),   # unmatched, cohort (time) 6
    "6.6" = integer(0)    # unmatched, cohort (time) 6
  )
  attr(fake.matched.sets, "id.var") <- "unit"
  attr(fake.matched.sets, "t.var")  <- "time"
  attr(fake.matched.sets, "lag")    <- 2
  
  fake.pm <- list(att = fake.matched.sets)
  attr(fake.pm, "qoi") <- "att"
  class(fake.pm) <- "PanelMatch"
  
  ## Pre-treatment covariate "x":
  ##   cohort at time 5 (t-1 = time 4, t-2 = time 3):
  ##     unit 1 (matched):   15, 13 -> pooled mean = 14
  ##     unit 2 (unmatched):  5,  3 -> pooled mean = 4
  ##     diff = 14 - 4 = 10
  ##   cohort at time 6 (t-1 = time 5, t-2 = time 4):
  ##     unit 3 (matched):    9,  7
  ##     unit 4 (matched):    9,  7  -> pooled mean = (9+7+9+7)/4 = 8
  ##     unit 5 (unmatched):  7,  5
  ##     unit 6 (unmatched):  7,  5  -> pooled mean = (7+5+7+5)/4 = 6
  ##     diff = 8 - 6 = 2
  panel.data <- data.frame(
    unit = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6),
    time = c(4, 3, 4, 3, 5, 4, 5, 4, 5, 4, 5, 4),
    x    = c(15, 13, 5, 3, 9, 7, 9, 7, 7, 5, 7, 5)
  )
  class(panel.data) <- c("PanelData", "data.frame")
  
  result <- compare_treated_observations(fake.pm, panel.data, covariates = "x")
  
  # --- counts: hand-computable ---
  expect_equal(result$n_has_match, 3)   # units 1, 3, 4
  expect_equal(result$n_no_match, 3)    # units 2, 5, 6
  expect_equal(result$pct_no_match, 50)
  expect_equal(result$n_no_viable, 0)   # both cohorts contain matched units
  expect_equal(result$pct_no_viable, 0)
  
  # --- covariate_diffs: this is the point of the test ---
  expect_equal(nrow(result$covariate_diffs), 1)
  expect_equal(result$covariate_diffs$mean_diff, 6)
  expect_equal(result$covariate_diffs$weighted_mean_diff, 28 / 6)
  # the two should NOT be equal here, unlike the single-cohort tests above
  expect_false(isTRUE(all.equal(
    result$covariate_diffs$mean_diff,
    result$covariate_diffs$weighted_mean_diff
  )))
})


test_that("compare_treated_observations() handles the case where no cohort has both matched and unmatched treated units", {
  ## Two treated units, in different time cohorts, both matched --
  ## no cohort ever contains an unmatched unit, so there's no valid basis
  ## for any within-cohort comparison.
  fake.matched.sets <- list(
    "1.5" = 3L,   # matched, cohort (time) 5
    "2.6" = 4L    # matched, cohort (time) 6
  )
  attr(fake.matched.sets, "id.var") <- "unit"
  attr(fake.matched.sets, "t.var")  <- "time"
  attr(fake.matched.sets, "lag")    <- 2
  
  fake.pm <- list(att = fake.matched.sets)
  attr(fake.pm, "qoi") <- "att"
  class(fake.pm) <- "PanelMatch"
  
  panel.data <- data.frame(
    unit = c(1, 1, 2, 2),
    time = c(4, 3, 5, 4),
    x    = c(12, 10, 6, 5)
  )
  class(panel.data) <- c("PanelData", "data.frame")
  
  expect_warning(
    result <- compare_treated_observations(fake.pm, panel.data, covariates = "x"),
    "No treatment cohorts contain both matched and unmatched treated units"
  )
  
  expect_equal(result$n_has_match, 2)
  expect_equal(result$n_no_match, 0)
  expect_equal(result$pct_no_match, 0)
  expect_equal(result$n_no_viable, 0)
  expect_true(is.na(result$pct_no_viable))  # NA when n_no_match == 0
  expect_null(result$covariate_diffs)
})

test_that("compare_treated_observations() runs cleanly on real PanelMatch/dem data and is internally consistent", {
  dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
  dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
  
  PM.results <- PanelMatch(
    panel.data                = dem.sub.panel,
    lag                        = 4,
    refinement.method          = "mahalanobis",
    match.missing              = TRUE,
    covs.formula               = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
    size.match                 = 5,
    qoi                        = "att",
    lead                       = 0:4,
    forbid.treatment.reversal  = FALSE
  )
  
  result <- suppressWarnings(
    compare_treated_observations(PM.results, dem.sub.panel, covariates = c("tradewb", "y"))
  )
  
  # returns a plain, unclassed list
  expect_true(is.list(result))
  expect_equal(class(result), "list")
  
  # counts should be internally consistent with the total number of
  # treated unit-times in the matched.set object
  n.treated.total <- length(PM.results[["att"]])
  expect_equal(result$n_has_match + result$n_no_match, n.treated.total)
  
  # if any comparison was possible, covariate_diffs should have one row
  # per covariate requested, in the order requested
  if (!is.null(result$covariate_diffs)) {
    expect_equal(nrow(result$covariate_diffs), 2)
    expect_equal(result$covariate_diffs$covariate, c("tradewb", "y"))
    expect_true(all(c("mean_diff", "weighted_mean_diff") %in% names(result$covariate_diffs)))
  }
})

test_that("compare_treated_observations() requires panel.data to be a PanelData object", {
  fake.matched.sets <- list("1.5" = 3L, "2.5" = integer(0))
  attr(fake.matched.sets, "id.var") <- "unit"
  attr(fake.matched.sets, "t.var")  <- "time"
  attr(fake.matched.sets, "lag")    <- 2
  
  fake.pm <- list(att = fake.matched.sets)
  attr(fake.pm, "qoi") <- "att"
  class(fake.pm) <- "PanelMatch"
  
  raw.data.frame <- data.frame(
    unit = c(1, 1, 2, 2),
    time = c(4, 3, 4, 3),
    x    = c(12, 10, 6, 5)
  )
  
  expect_error(
    compare_treated_observations(fake.pm, raw.data.frame, covariates = "x"),
    "Please provide a PanelData object"
  )
})
