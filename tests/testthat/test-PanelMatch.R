# This file contains only a bare bones set of tests.
# Please add future tests to manual_tests/test-*.R
# Similarly, please test package against that set of tests before submitting.
# For this, use testthat::test_dir()

test_that("(ATT) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimate, comp.results)
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
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
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
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimate, comp.results)
  
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
  comp.results <-  c(-0.73400424, -0.14920097, -0.45276678, -0.06580353)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .0000001)
  
})
