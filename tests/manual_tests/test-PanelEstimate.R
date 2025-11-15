test_that("(ATT) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  comp.results <- c(-0.9036306, -0.3174628, 0.5733430, 1.6002885)
  names(comp.results) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, comp.results, tolerance = 1e-6)
})

test_that("(ATC) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  comp.results <- c(-0.7690320, -0.2760459, -0.6871066, -0.3738689)
  names(comp.results) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, comp.results, tolerance = 1e-6)
})

test_that("(ART) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  comp.results <- c(-5.1717610, -7.8607830, -8.4489830, -7.8208880)
  names(comp.results) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, comp.results, tolerance = 1e-6)
})

test_that("(ATE) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "ate"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "bootstrap"
  )
  
  comp.results <- c(-0.7745270, -0.2777368, -0.6356488, -0.2932741)
  names(comp.results) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, comp.results, tolerance = 1e-6)
})

test_that("(ATT) bootstrap SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100
  )
  
  est <- c(-0.9036306, -0.3174628, 0.5733430, 1.6002885)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.7508209, 1.1889910, 1.5156986, 1.7993030)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATC) bootstrap SEs", {
  qoi_ <- "atc"
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  set.seed(1)
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100
  )
  
  est <- c(-0.7690320, -0.2760459, -0.6871066, -0.3738689)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.6917734, 1.3272229, 1.7502791, 2.1681358)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ART) bootstrap SEs", {
  qoi_ <- "art"
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  set.seed(1)
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100
  )
  
  est <- c(-5.1717610, -7.8607830, -8.4489830, -7.8208880)
  names(est) <- paste0("t+", 0:3)
  se <- c(1.3650310, 2.1535930, 2.8435180, 3.1568320)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATE) bootstrap SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  set.seed(1)
  qoi_ <- "ate"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100
  )
  
  est <- c(-0.7745270, -0.2777368, -0.6356488, -0.2932741)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.6824677, 1.2995159, 1.7102891, 2.1120251)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("bootstrap SEs - Pooled ", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  ## ATT
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100,
    pooled = TRUE
  )
  expect_equal(pe.results$estimate, 0.2381345, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, 1.200568, tolerance = 1e-6)
  
  ## ATC
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  set.seed(1)
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100,
    pooled = TRUE
  )
  expect_equal(pe.results$estimate, -0.5265133, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, 1.425316, tolerance = 1e-6)
  
  ## ART
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  set.seed(1)
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100,
    pooled = TRUE
  )
  expect_equal(pe.results$estimate, -7.325604, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, 2.323442, tolerance = 1e-6)
  
  ## ATE
  qoi_ <- "ate"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  set.seed(1)
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    number.iterations = 100,
    pooled = TRUE
  )
  expect_equal(pe.results$estimate, -0.4952967, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, 1.391672, tolerance = 1e-6)
})

test_that("(ATT) PanelEstimate Runs: analytical SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  est <- c(-0.9036306, -0.3174628, 0.5733430, 1.6002885)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.6394206, 1.0766697, 1.4164110, 1.7346946)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATC) PanelEstimate Runs: analytical SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  est <- c(-0.7690320, -0.2760459, -0.6871066, -0.3738689)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.5929128, 1.1390635, 1.5049599, 1.8968186)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ART) PanelEstimate Runs: analytical SEs ", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "conditional"
  )
  
  est <- c(-5.1717610, -7.8607830, -8.4489830, -7.8208880)
  names(est) <- paste0("t+", 0:3)
  se <- c(1.0227660, 1.4848030, 1.9151960, 2.1476240)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATT) PanelEstimate Runs: unconditional analytical SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "unconditional"
  )
  
  est <- c(-0.9036306, -0.3174628, 0.5733430, 1.6002885)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.7929164, 1.3280250, 1.7474357, 2.1441677)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATC) PanelEstimate Runs: unconditional analytical SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "unconditional"
  )
  
  est <- c(-0.7690320, -0.2760459, -0.6871066, -0.3738689)
  names(est) <- paste0("t+", 0:3)
  se <- c(0.7102239, 1.3605176, 1.7980118, 2.2655081)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ART) PanelEstimate Runs: unconditional analytical SEs ", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets = pm1,
    panel.data = dem.panel,
    se.method = "unconditional"
  )
  
  est <- c(-5.1717610, -7.8607830, -8.4489830, -7.8208880)
  names(est) <- paste0("t+", 0:3)
  se <- c(1.5748440, 2.3154610, 2.8570590, 3.0942460)
  names(se) <- paste0("t+", 0:3)
  
  expect_equal(pe.results$estimate, est, tolerance = 1e-6)
  expect_equal(pe.results$standard.error, se, tolerance = 1e-6)
})

test_that("(ATE) PanelEstimate fails: analytical SEs ", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  qoi_ <- "ate"
  pm1 <- PanelMatch(
    lag = 4,
    refinement.method = "mahalanobis",
    panel.data = dem.panel,
    match.missing = FALSE,
    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match = 5,
    qoi = qoi_,
    lead = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  expect_error(
    PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  )
  expect_error(
    PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  )
})