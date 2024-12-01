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
  expect_output(print(pe.results))
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  expect_output(print(pe.results))
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  expect_output(print(pe.results))
  
  
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
  expect_output(print(pe.results))
  
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
# test_that("plot.PanelEstimate tests", {
#   dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#   
#   qoi_ <- "att"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
#   plot(pe.results)
#   
#   qoi_ <- "ate"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
#   plot(pe.results)
#   
# })