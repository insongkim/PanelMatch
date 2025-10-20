test_that("test matched set getter function", {
  dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
  dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
  # create subset of data for simplicity
  PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  expect_true(inherits(PM.results, "PanelMatch"))
  att.set <- extract(PM.results, qoi = "att")
  expect_true(inherits(att.set, "matched.set"))
  att.set2 <- extract(PM.results)
  expect_true(inherits(att.set2, "matched.set"))
  expect_true(identical(att.set, att.set2))
  
  PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  expect_true(inherits(PM.results, "PanelMatch"))
  att.set <- extract(PM.results, qoi = "att")
  att.set2 <- extract(PM.results, qoi = "atc")
  expect_true(inherits(att.set, "matched.set"))
  expect_true(inherits(att.set2, "matched.set"))
  expect_error(extract(PM.results, qoi = "ate"))
  
  
  PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  expect_true(inherits(PM.results, "PanelMatch"))
  att.set <- extract(PM.results, qoi = "art")
  expect_true(inherits(att.set, "matched.set"))
  
  
})

test_that("test getter results", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "att")
  expect_true(identical(r1, PM.results[["att"]]))
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "art")
  expect_true(identical(r1, PM.results[["art"]]))
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "atc",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "atc")
  expect_true(identical(r1, PM.results[["atc"]]))
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  r1 <- extract(PM.results, qoi = "att")
  expect_true(identical(r1, PM.results[["att"]]))
  r1 <- extract(PM.results, qoi = "atc")
  expect_true(identical(r1, PM.results[["atc"]]))
})

test_that("PanelMatch summary method", {
  
  dem.sub.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "att",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  #l1 <- summary(PM.results, verbose = TRUE)
  l2 <- summary(PM.results)
  #expect_true(length(l1) == 1)
  expect_true(length(l2) == 1)
  expect_true(nrow(l2[["att"]]) == 8)
  expect_true(inherits(l2[['att']], "data.frame"))
  
  
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "art",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  #l1 <- summary(PM.results, verbose = TRUE)
  l2 <- summary(PM.results)
  expect_true(length(l2) == 1)
  expect_true(nrow(l2[["art"]]) == 8)
  expect_true(inherits(l2[['art']], "data.frame"))
  
  
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "ate",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  #l1 <- summary(PM.results, verbose = TRUE)
  l2 <- summary(PM.results)
  expect_true(length(l2) == 2)
  expect_true(nrow(l2[["att"]]) == 8)
  expect_true(inherits(l2[['att']], "data.frame"))
  expect_true(nrow(l2[["atc"]]) == 8)
  expect_true(inherits(l2[['atc']], "data.frame"))
  
})


test_that("print.PanelMatch", {
  
  dem.sub.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "att",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  
  expect_output(print(PM.results))
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "ate",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_output(print(PM.results))  
})


test_that("plot.PanelMatch", {
  dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
  dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
  PM.results <- PanelMatch(panel.data = dem.sub.panel,
                           lag = 4, 
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
                           size.match = 5, qoi = "att",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  plot(PM.results)
  plot(PM.results, include.empty.sets = TRUE)
})