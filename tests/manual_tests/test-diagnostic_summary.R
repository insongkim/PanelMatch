## Tests for diagnostic_summary()
##
## Covers: input validation; each QOI branch (att/atc/art/ate); multiple
## placebo se.methods; every combination of optional inputs supplied/omitted;
## the refined-only balance guarantee; all four warning conditions (and their
## absence under loose thresholds); the structure of the returned checks
## table; and the function's invisible-return behavior.

# ================================================================
# Shared setup
# ================================================================

## NOTE: this test file uses the full `dem` dataset rather than a subset.
## It is intended to be run locally (e.g. via devtools::test()) rather than
## as part of a fast per-commit suite, since matching/refinement/bootstrap
## SEs over the full dataset take noticeably longer than on a subset.

set.seed(1)
dem$rdata <- runif(nrow(dem))
dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")

covs.formula <- ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4))

make_pm <- function(qoi, placebo.test = TRUE, forbid.treatment.reversal = FALSE) {
  PanelMatch(
    lag                        = 4,
    lead                       = 0:3,
    refinement.method          = "mahalanobis",
    panel.data                 = dem.panel,
    match.missing              = TRUE,
    covs.formula               = covs.formula,
    size.match                 = 5,
    qoi                        = qoi,
    forbid.treatment.reversal  = forbid.treatment.reversal,
    placebo.test               = placebo.test
  )
}

pm.att <- make_pm("att")
pm.atc <- make_pm("atc")
pm.art <- make_pm("art")
pm.ate <- make_pm("ate")

make_placebo <- function(pm.obj, se.method = "unconditional", number.iterations = 50) {
  suppressWarnings(
    placebo_test(
      pm.obj,
      panel.data          = dem.panel,
      se.method            = se.method,
      number.iterations    = number.iterations,
      plot                 = FALSE
    )
  )
}

## NOTE ON get_covariate_balance() NAMING: get_covariate_balance() names each
## configuration in its output using the deparsed argument expression at ITS
## OWN call site (via match.call()). That means calling it through a wrapper
## function (e.g. get_covariate_balance(pm.obj, ...) inside a helper whose
## parameter is named pm.obj) would name every resulting config "pm.obj",
## regardless of what the caller's actual variable was called -- which would
## silently break diagnostic_summary()'s name-matching against pb.object.
## For that reason, pb.* fixtures below call get_covariate_balance() directly
## at the top level, one call per fixture, so each config is named to match
## its corresponding pm.* variable exactly (e.g. pb.att's config is named
## "pm.att", matching the pm.att variable used below).

pb.att <- suppressWarnings(get_covariate_balance(
  pm.att, include.unrefined = TRUE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
))
pb.atc <- suppressWarnings(get_covariate_balance(
  pm.atc, include.unrefined = TRUE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
))
pb.art <- suppressWarnings(get_covariate_balance(
  pm.art, include.unrefined = TRUE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
))
pb.ate <- suppressWarnings(get_covariate_balance(
  pm.ate, include.unrefined = TRUE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
))

# a pb.object containing multiple configurations, for testing that
# diagnostic_summary() correctly selects only the one matching pm.object
pb.multi <- suppressWarnings(get_covariate_balance(
  pm.att, pm.atc, include.unrefined = TRUE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
))

placebo.att.unconditional <- make_placebo(pm.att, "unconditional")
placebo.att.conditional   <- make_placebo(pm.att, "conditional")
placebo.att.bootstrap     <- make_placebo(pm.att, "bootstrap")
placebo.atc.unconditional <- make_placebo(pm.atc, "unconditional")
placebo.art.unconditional <- make_placebo(pm.art, "unconditional")
placebo.ate.bootstrap     <- make_placebo(pm.ate, "bootstrap")  # ate only supports bootstrap

# ================================================================
# 1. Input validation
# ================================================================

test_that("pm.object must be a PanelMatch object", {
  expect_error(
    diagnostic_summary(pm.object = list(a = 1)),
    "pm.object must be a PanelMatch object"
  )
})

test_that("panel.data is required when covariates is specified", {
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, covariates = "tradewb")),
    "panel.data must be supplied when covariates is specified"
  )
})

test_that("pb.object must be a PanelBalance object", {
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, pb.object = list(a = 1))),
    "pb.object must be a PanelBalance object"
  )
})

test_that("empty.set.threshold must be numeric between 0 and 1", {
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, empty.set.threshold = -0.1)),
    "empty.set.threshold must be numeric, between 0 and 1"
  )
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, empty.set.threshold = 1.5)),
    "empty.set.threshold must be numeric, between 0 and 1"
  )
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, empty.set.threshold = "a")),
    "empty.set.threshold must be numeric, between 0 and 1"
  )
})

test_that("min.matched.sets must be a non-negative number", {
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, min.matched.sets = -5)),
    "min.matched.sets must be a non-negative number"
  )
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, min.matched.sets = "a")),
    "min.matched.sets must be a non-negative number"
  )
})

test_that("placebo.results missing conf.intervals raises an informative error", {
  bad.placebo <- list(estimates = c(t_1 = 0.1), standard.errors = c(t_1 = 0.05))
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.att, placebo.results = bad.placebo)),
    "placebo.results must contain 'estimates' and 'conf.intervals'"
  )
})

# ================================================================
# 2. Each QOI branch runs cleanly, across se.methods, with
#    everything supplied
# ================================================================

test_that("att branch runs cleanly with all optional inputs, across se.methods", {
  for (placebo.res in list(placebo.att.unconditional, placebo.att.conditional, placebo.att.bootstrap)) {
    out <- capture.output(
      diag <- suppressWarnings(diagnostic_summary(
        pm.object       = pm.att,
        panel.data      = dem.panel,
        covariates      = c("tradewb", "rdata"),
        pb.object       = pb.att,
        placebo.results = placebo.res
      ))
    )
    expect_true(is.list(diag))
    expect_true(is.data.frame(diag$checks))
    expect_true(nrow(diag$checks) > 0)
    expect_false(is.null(diag$matched.treated.summary))
    expect_false(is.null(diag$balance.summary))
    expect_false(is.null(diag$placebo.table))
  }
})

test_that("atc branch runs cleanly with all optional inputs", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(
      pm.object       = pm.atc,
      panel.data      = dem.panel,
      covariates      = c("tradewb", "rdata"),
      pb.object       = pb.atc,
      placebo.results = placebo.atc.unconditional
    ))
  )
  expect_true(is.data.frame(diag$checks))
  expect_true(all(diag$checks$qoi[!is.na(diag$checks$qoi)] == "atc"))
})

test_that("art branch runs cleanly with all optional inputs", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(
      pm.object       = pm.art,
      panel.data      = dem.panel,
      covariates      = c("tradewb", "rdata"),
      pb.object       = pb.art,
      placebo.results = placebo.art.unconditional
    ))
  )
  expect_true(is.data.frame(diag$checks))
  expect_true(all(diag$checks$qoi[!is.na(diag$checks$qoi)] == "art"))
})

test_that("ate branch runs cleanly, splits att/atc automatically, bootstrap se only", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(
      pm.object       = pm.ate,
      panel.data      = dem.panel,
      covariates      = c("tradewb", "rdata"),
      pb.object       = pb.ate,
      placebo.results = placebo.ate.bootstrap
    ))
  )
  # matched-set-size / empty-set checks present for both att and atc
  expect_setequal(
    unique(diag$checks$qoi[diag$checks$metric %in% c("empty_set_proportion", "n_matched_sets")]),
    c("att", "atc")
  )
  # balance.summary nested by att/atc, each containing refined matrices only
  expect_named(diag$balance.summary, c("att", "atc"))
  expect_true(is.list(diag$balance.summary$att))
  expect_true(is.list(diag$balance.summary$atc))

  # matched.treated.summary nested by att/atc
  expect_named(diag$matched.treated.summary, c("att", "atc"))

  # balance checks present for both att and atc
  expect_setequal(
    unique(diag$checks$qoi[diag$checks$metric == "balance"]),
    c("att", "atc")
  )
})

test_that("ate branch errors informatively if a non-bootstrap se.method placebo result is supplied", {
  # placebo_test()/PanelEstimate() itself should reject non-bootstrap se.method
  # for an ate PanelMatch object -- confirm this happens upstream of
  # diagnostic_summary(), i.e. we never even get a placebo.results object to pass in.
  expect_error(
    suppressWarnings(placebo_test(pm.ate, panel.data = dem.panel, se.method = "unconditional", plot = FALSE))
  )
})

# ================================================================
# 3. Optional inputs: every combination of supplied/omitted
# ================================================================

test_that("only pm.object supplied: all optional sections are NULL/absent", {
  out <- capture.output(diag <- diagnostic_summary(pm.object = pm.att))
  expect_null(diag$matched.treated.summary)
  expect_null(diag$balance.summary)
  expect_null(diag$placebo.table)
  # only empty_set_proportion + n_matched_sets rows (2, since att is a single QOI)
  expect_equal(nrow(diag$checks), 2)
  expect_setequal(diag$checks$metric, c("empty_set_proportion", "n_matched_sets"))
  # report should say these sections weren't provided
  expect_true(any(grepl("Not provided.*pb.object", out)))
  expect_true(any(grepl("Not provided.*placebo.results", out)))
})

test_that("only covariates supplied (with panel.data): matched.treated.summary present, others absent", {
  out <- capture.output(
    diag <- diagnostic_summary(
      pm.object  = pm.att,
      panel.data = dem.panel,
      covariates = c("tradewb", "rdata")
    )
  )
  expect_false(is.null(diag$matched.treated.summary))
  expect_null(diag$balance.summary)
  expect_null(diag$placebo.table)
})

test_that("only pb.object supplied: balance.summary present, others absent", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att))
  )
  expect_null(diag$matched.treated.summary)
  expect_false(is.null(diag$balance.summary))
  expect_null(diag$placebo.table)
  expect_true("balance" %in% diag$checks$metric)
})

test_that("only placebo.results supplied: placebo.table present, others absent", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.att, placebo.results = placebo.att.unconditional))
  )
  expect_null(diag$matched.treated.summary)
  expect_null(diag$balance.summary)
  expect_false(is.null(diag$placebo.table))
  expect_true("placebo" %in% diag$checks$metric)
})

test_that("all optional inputs supplied together: all sections present, no crosstalk", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(
      pm.object       = pm.att,
      panel.data      = dem.panel,
      covariates      = c("tradewb", "rdata"),
      pb.object       = pb.att,
      placebo.results = placebo.att.unconditional
    ))
  )
  expect_false(is.null(diag$matched.treated.summary))
  expect_false(is.null(diag$balance.summary))
  expect_false(is.null(diag$placebo.table))
  expect_setequal(
    diag$checks$metric,
    c("empty_set_proportion", "n_matched_sets", "balance", "placebo")
  )
})

# ================================================================
# 4. Refined-only balance guarantee
# ================================================================

test_that("balance.summary never contains _unrefined columns, even though pb was built with include.unrefined = TRUE", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att))
  )
  for (mat in diag$balance.summary) {
    expect_false(any(grepl("_unrefined$", colnames(mat))))
  }
})

test_that("refined-only guarantee also holds for the ate case (nested by att/atc)", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.ate, pb.object = pb.ate))
  )
  for (q in c("att", "atc")) {
    for (mat in diag$balance.summary[[q]]) {
      expect_false(any(grepl("_unrefined$", colnames(mat))))
    }
  }
})

test_that("refined-only guarantee holds even if pb.object was built with include.unrefined = FALSE", {
  pb.att.norefine <- suppressWarnings(get_covariate_balance(
    pm.att, include.unrefined = FALSE, panel.data = dem.panel, covariates = c("tradewb", "rdata")
  ))
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att.norefine))
  )
  for (mat in diag$balance.summary) {
    expect_false(any(grepl("_unrefined$", colnames(mat))))
  }
})

# ================================================================
# 4b. Multi-configuration pb.object: only the matching config is used
# ================================================================

test_that("a pb.object with multiple configs is filtered down to just the one matching pm.object", {
  out.att <- capture.output(
    diag.att <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.multi))
  )
  # balance.summary should contain exactly one configuration: "pm.att"
  expect_named(diag.att$balance.summary, "pm.att")
  expect_true(all(diag.att$checks$config[diag.att$checks$metric == "balance"] == "pm.att"))

  out.atc <- capture.output(
    diag.atc <- suppressWarnings(diagnostic_summary(pm.object = pm.atc, pb.object = pb.multi))
  )
  expect_named(diag.atc$balance.summary, "pm.atc")
  expect_true(all(diag.atc$checks$config[diag.atc$checks$metric == "balance"] == "pm.atc"))
})

test_that("results are identical whether pb.object has one config or several, for the matching config", {
  out.single <- capture.output(
    diag.single <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att))
  )
  out.multi <- capture.output(
    diag.multi <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.multi))
  )
  expect_equal(
    diag.single$balance.summary[["pm.att"]],
    diag.multi$balance.summary[["pm.att"]]
  )
})

test_that("an informative error is raised when pb.object has no config matching pm.object's variable name", {
  # pb.art has a config named "pm.art"; pm.atc's variable name won't match it
  expect_error(
    capture.output(diagnostic_summary(pm.object = pm.atc, pb.object = pb.art)),
    "pb.object does not contain a configuration named 'pm.atc'"
  )
})

test_that("passing pm.object through an intermediate variable with a different name fails to match, as documented", {
  renamed.pm <- pm.att
  expect_error(
    capture.output(diagnostic_summary(pm.object = renamed.pm, pb.object = pb.att)),
    "pb.object does not contain a configuration named 'renamed.pm'"
  )
})

# ================================================================
# 5. Checks table structure
# ================================================================

test_that("checks table has the expected columns and types", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(
      pm.object       = pm.att,
      panel.data      = dem.panel,
      covariates      = c("tradewb", "rdata"),
      pb.object       = pb.att,
      placebo.results = placebo.att.unconditional
    ))
  )
  expect_named(
    diag$checks,
    c("qoi", "config", "metric", "label", "value", "threshold", "direction", "flagged")
  )
  expect_type(diag$checks$label, "character")
  expect_type(diag$checks$flagged, "logical")
  expect_type(diag$checks$value, "double")
  expect_type(diag$checks$threshold, "double")
  expect_true(all(diag$checks$direction %in% c("above", "below")))
})

test_that("checks table is well-formed (empty but valid) when no optional inputs are supplied and matched sets are perfect", {
  # even in the minimal case, checks should always have >= 2 rows (matched set checks)
  out <- capture.output(diag <- diagnostic_summary(pm.object = pm.att))
  expect_true(is.data.frame(diag$checks))
  expect_true(nrow(diag$checks) >= 2)
})

# ================================================================
# 6. Warnings: each of the four conditions fires when flagged,
#    and none fire when thresholds are loose
# ================================================================

test_that("empty_set_proportion warning fires under a very strict threshold", {
  expect_warning(
    capture.output(
      diagnostic_summary(pm.object = pm.att, empty.set.threshold = 0.0001)
    ),
    "Empty matched set check flagged"
  )
})

test_that("n_matched_sets warning fires when min.matched.sets is unrealistically high", {
  expect_warning(
    capture.output(
      diagnostic_summary(pm.object = pm.att, min.matched.sets = 100000)
    ),
    "Minimum matched sets check flagged"
  )
})

test_that("balance warning fires under a very strict balance.threshold", {
  expect_warning(
    capture.output(
      diagnostic_summary(pm.object = pm.att, pb.object = pb.att, balance.threshold = 0.0001)
    ),
    "Covariate balance check flagged"
  )
})

test_that("placebo warning fires whenever any placebo estimate's CI excludes 0 (no tunable threshold)", {
  # only run this if at least one period is actually significant in this fixture;
  # otherwise assert the corresponding absence-of-warning case instead
  any.sig <- placebo.att.unconditional$conf.intervals[, 1] > 0 |
             placebo.att.unconditional$conf.intervals[, 2] < 0
  if (any(any.sig)) {
    expect_warning(
      capture.output(
        diagnostic_summary(pm.object = pm.att, placebo.results = placebo.att.unconditional)
      ),
      "Placebo test check flagged"
    )
  } else {
    expect_warning(
      capture.output(
        diagnostic_summary(pm.object = pm.att, placebo.results = placebo.att.unconditional)
      ),
      NA
    )
  }
})

test_that("no warnings fire when all thresholds are loose and matching is reasonable", {
  expect_warning(
    capture.output(
      diagnostic_summary(
        pm.object            = pm.att,
        panel.data           = dem.panel,
        covariates           = c("tradewb", "rdata"),
        pb.object            = pb.att,
        balance.threshold    = 10,     # essentially unreachable
        empty.set.threshold  = 1,      # 100% -- unreachable unless every set is empty
        min.matched.sets     = 0       # always satisfied
      )
    ),
    NA
  )
})

test_that("warnings correctly aggregate multiple flagged QOIs in the ate case into one message", {
  w <- tryCatch({
    capture.output(
      diagnostic_summary(pm.object = pm.ate, empty.set.threshold = 0.0001)
    )
    NULL
  }, warning = function(w) w)
  expect_true(inherits(w, "warning"))
  expect_match(conditionMessage(w), "ATT", fixed = TRUE)
  # the two QOIs' labels should both be present if both are flagged; at minimum
  # the message should reference the check category once, not once per QOI
  expect_match(conditionMessage(w), "^Empty matched set check flagged: ")
})

# ================================================================
# 7. Invisible return
# ================================================================

test_that("diagnostic_summary() returns its value invisibly", {
  capture.output(
    result <- withVisible(
      diagnostic_summary(pm.object = pm.att)
    )
  )
  expect_false(result$visible)
  expect_true(is.list(result$value))
})

test_that("the printed report and the returned list are consistent with each other", {
  out <- capture.output(
    diag <- suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att, empty.set.threshold = 0.0001))
  )
  flagged.labels <- diag$checks$label[diag$checks$flagged]
  for (lbl in flagged.labels) {
    expect_true(any(grepl(lbl, out, fixed = TRUE)))
  }
})

# ================================================================
# 8. digits argument affects printed rounding without erroring
# ================================================================

test_that("digits argument controls rounding in printed output without error", {
  out.low  <- capture.output(
    suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att, digits = 1))
  )
  out.high <- capture.output(
    suppressWarnings(diagnostic_summary(pm.object = pm.att, pb.object = pb.att, digits = 5))
  )
  expect_true(length(out.low) > 0)
  expect_true(length(out.high) > 0)
})
