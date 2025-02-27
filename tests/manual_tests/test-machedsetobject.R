test_that("test weights function", {
  
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "att")
  lt <- weights(r1)
  expect_true(length(lt) == 105)
  expect_true(sum(lt[[1]]) == 1) #just verify that the first matched set has the properties it should
  
  
  r1 <- extract(PM.results) #check default behavior
  lt <- weights(r1)
  expect_true(length(lt) == 105)
  expect_true(sum(lt[[1]]) == 1) #just verify that the first matched set has the properties it should
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "art")
  lt <- weights(r1)
  expect_true(length(lt) == 54)
  expect_true(sum(lt[[1]]) == 1) #just verify that the first matched set has the properties it should
  
  
  # PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
  #                          refinement.method = "ps.match",
  #                          match.missing = TRUE,
  #                          covs.formula = ~ tradewb,
  #                          size.match = 5, qoi = "atc",
  #                          lead = 0:4,
  #                          forbid.treatment.reversal = FALSE)
  # r1 <- extract(PM.results, qoi = "atc")
  
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  r1 <- extract(PM.results, qoi = "att")
  lt <- weights(r1)
  expect_true(length(lt) == 105)
  expect_true(sum(lt[[1]]) == 1) #just verify that the first matched set has the properties it should
  r1 <- extract(PM.results, qoi = "atc")
  lt <- weights(r1)
  expect_true(length(lt) == 3287)
  expect_true(sum(lt[[1]]) == 1) #just verify that the first matched set has the properties it should
  
  
})


test_that("test distances function", {
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "mahalanobis",
                           verbose = TRUE,
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "att")
  lt <- distances(r1)
  expect_true(length(lt) == 105)
  
  
  
  r1 <- extract(PM.results) #check default behavior
  lt <- distances(r1)
  expect_true(length(lt) == 105)
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "mahalanobis",
                           verbose = TRUE,
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "art")
  lt <- distances(r1)
  expect_true(length(lt) == 54)

  
  
  # PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
  #                          refinement.method = "ps.match",
  #                          match.missing = TRUE,
  #                          covs.formula = ~ tradewb,
  #                          size.match = 5, qoi = "atc",
  #                          lead = 0:4,
  #                          forbid.treatment.reversal = FALSE)
  # r1 <- extract(PM.results, qoi = "atc")
  
  
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "mahalanobis",
                           match.missing = TRUE,
                           verbose = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  r1 <- extract(PM.results, qoi = "att")
  lt <- distances(r1)
  expect_true(length(lt) == 105)
  
  r1 <- extract(PM.results, qoi = "atc")
  lt <- distances(r1)
  expect_true(length(lt) == 3287)
  
  
  
  
  #check that distances fails for non matching based method
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.weight",
                           verbose= TRUE,
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  r1 <- extract(PM.results, qoi = "art")
  expect_error(lt <- distances(r1))
  
})

test_that("test print method", {
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
 
  mso <- extract(PM.results)
  expect_output(print(mso))
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  expect_output(print(mso))
  
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                         refinement.method = "ps.match",
                         match.missing = TRUE,
                         covs.formula = ~ tradewb,
                         size.match = 5, qoi = "atc",
                         lead = 0:4,
                         forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  expect_output(print(mso))
  
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  expect_error(mso <- extract(PM.results))
  expect_error(mso <- extract(PM.results, qoi = "ate"))
  mso <- extract(PM.results, qoi = "att")
  expect_output(print(mso))
  
  mso <- extract(PM.results, qoi = "atc")
  expect_output(print(mso))
})

test_that("test summary method", {
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  smso <- summary(mso)
  expect_true(length(smso) == 5L)
  expect_true(all(names(smso) == c("overview", "set.size.summary", "number.of.treated.units", "num.units.empty.set", "lag")))
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "atc",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  smso <- summary(mso)
  expect_true(length(smso) == 5L)
  expect_true(all(names(smso) == c("overview", "set.size.summary", "number.of.treated.units", "num.units.empty.set", "lag")))
  
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  smso <- summary(mso)
  expect_true(length(smso) == 5L)
  expect_true(all(names(smso) == c("overview", "set.size.summary", "number.of.treated.units", "num.units.empty.set", "lag")))
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results, qoi = "att")
  smso <- summary(mso)
  expect_true(length(smso) == 5L)
  expect_true(all(names(smso) == c("overview", "set.size.summary", "number.of.treated.units", "num.units.empty.set", "lag")))
  
  mso <- extract(PM.results, qoi = "atc")
  smso <- summary(mso)
  expect_true(length(smso) == 5L)
  expect_true(all(names(smso) == c("overview", "set.size.summary", "number.of.treated.units", "num.units.empty.set", "lag")))
})

test_that("test plot method (manual)", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  plot(mso, panel.data = dem.panel)
  plot(mso, panel.data = dem.panel, include.missing = FALSE)
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.weight",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  mso <- extract(PM.results)
  plot(mso, panel.data = dem.panel)
  plot(mso, panel.data = dem.panel, include.missing = FALSE)
  
})
