test_that("testing MSM - failure", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  expect_error(PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.msm.weight",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE))
  
  expect_error(PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "CBPS.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "att",
                          lead = 0:4,
                          forbid.treatment.reversal = FALSE))
  
  
  expect_error(PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "ps.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "ate",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE))
  
  expect_error(PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "ps.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "atc",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE))
  
})


test_that("testing MSM - success", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.result <- PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "ps.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "att",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE)
  
  expect_no_error(PanelEstimate(sets = pm.result, 
                panel.data = dem.panel, 
                se.method = "unconditional"))
  
  pm.result <- PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "CBPS.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "att",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE)
  
  expect_no_error(PanelEstimate(sets = pm.result, 
                               panel.data = dem.panel, 
                               se.method = "unconditional"))
  
  
  pm.result <- PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "ps.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "art",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE)
  
  expect_no_error(PanelEstimate(sets = pm.result, 
                               panel.data = dem.panel, 
                               se.method = "unconditional"))
  
  
  pm.result <- PanelMatch(panel.data = dem.panel, lag = 4,
                          refinement.method = "CBPS.msm.weight",
                          match.missing = TRUE,
                          covs.formula = ~ tradewb,
                          size.match = 5, qoi = "art",
                          lead = 0:4,
                          forbid.treatment.reversal = TRUE)
  
  expect_no_error(PanelEstimate(sets = pm.result, 
                               panel.data = dem.panel, 
                               se.method = "unconditional"))
})
