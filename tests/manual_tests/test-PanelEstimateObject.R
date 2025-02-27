test_that("print.PanelEstimate tests", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  expect_output(print(pe.results, regexp = "Point estimates:"))
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  expect_output(print(pe.results, regexp = "Point estimates:"))
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  expect_output(print(pe.results, regexp = "Point estimates:"))
  
  
  ##### trying ate
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  expect_output(print(pe.results, regexp = "Point estimates:"))
  
})

test_that("summary.PanelEstimate tests", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  trt <- summary(pe.results)
  expect_true(all(dim(trt) == c(4, 4)))
  expect_equal(trt$estimate, c(-0.5933998, -0.3212602,  0.4563113, 1.7318216), tolerance = .000001)
  # check lower bounds and that one can specify a new confidence level
  expect_equal(trt[,3], c(-2.367325, -3.223384 ,-3.282366 ,-2.653756), tolerance = .000001)
  expect_equal(summary(pe.results, confidence.level = .9)[,3], 
               c(-2.082125, -2.756800, -2.681285 ,-1.948671), tolerance = .000001)
  
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  trt <- summary(pe.results)
  expect_true(all(dim(trt) == c(4, 4)))
  
  
  
  
})

# no good way to test plotting results directly. 
test_that("plot.PanelEstimate tests", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  plot(pe.results)
  
  plot(pe.results, confidence.level = .9)
  plot(pe.results, confidence.level = .99)

  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  plot(pe.results)

})