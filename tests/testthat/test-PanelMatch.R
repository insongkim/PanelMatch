# This file contains only a bare bones set of tests.
# Please add future tests to manual_tests/test-*.R
# Similarly, please test package against that set of tests before submitting.
# For this, use testthat::test_dir("/Users/adamrauh/Code/PanelMatch/PanelMatch/tests/manual_tests")

test_that("(ATT) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.9036306, -0.3174628,  0.5733430,  1.6002885)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
})


test_that("(ATC) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.7690320, -0.2760459, -0.6871066, -0.3738689)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  
})

test_that("(ART) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-5.171761, -7.860783, -8.448983, -7.820888)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  
})


test_that("(ATE) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, panel.data = dem.panel, se.method = "bootstrap", number.iterations = 300)
  comp.results <-  c(-0.7745270, -0.2777368, -0.6356488, -0.2932741)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .0000001)
  
})
