
test_that("(ATT) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimate, comp.results)
})



test_that("(ATC) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  
})

test_that("(ART) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimate, comp.results)
  
})


test_that("(ATE) PanelEstimate Runs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  comp.results <-  c(-0.73400424, -0.14920097, -0.45276678, -0.06580353)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .0000001)
  
})


test_that("(ATT) bootstrap SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100)
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimate, comp.results)
  expect_equivalent(pe.results$standard.error, c(0.8701336, 1.3611202, 1.7512925, 2.0019833), 
                    tolerance = .0000001)
})



test_that("(ATC) bootstrap SEs", {
  qoi_ <- "atc"
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100)
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .0000001)
  expect_equivalent(pe.results$standard.error, c(0.7027994, 1.3144841, 1.7200521, 2.1352594 ), tolerance = .0000001)
})

test_that("(ART) bootstrap SEs", {
  qoi_ <- "art"
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100)
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimate, comp.results)
  expect_equivalent(pe.results$standard.error, c(1.342158, 2.159589, 2.869549, 3.231702), tolerance = .000001)
})


test_that("(ATE) bootstrap SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  set.seed(1)
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100)
  comp.results <-  c(-0.73400424, -0.14920097, -0.45276678, -0.06580353)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  expect_equivalent(pe.results$standard.error, c(0.6953794, 1.2931504, 1.6894504, 2.0902084), tolerance = .0000001)
})


test_that("bootstrap SEs - Pooled ", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100, pooled = TRUE)
  
  expect_equal(pe.results$estimate, 0.3183682, tolerance = .000001)
  expect_equal(pe.results$standard.error, 1.404445, tolerance = .000001)
  
  
  qoi_ <- "atc"
  
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  set.seed(1)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100, pooled = TRUE)
  expect_equal(pe.results$estimate, -0.3789102, tolerance = .000001)
  expect_equal(pe.results$standard.error, 1.40745, tolerance = .000001)
  
  
  qoi_ <- "art"
  
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  set.seed(1)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100, pooled = TRUE)
  expect_equal(pe.results$estimate, -7.529892, tolerance = .000001)
  expect_equal(pe.results$standard.error, 2.341633, tolerance = .000001)
  
  
  
  qoi_ <- "ate"
  
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  set.seed(1)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, number.iterations = 100, pooled = TRUE)
  expect_equal(pe.results$estimate, -0.3504439, tolerance = .000001)
  expect_equal(pe.results$standard.error, 1.381604, tolerance = .000001)
})



test_that("(ATT) PanelEstimate Runs: analytical SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimate, comp.results)
  comp.results <- c(0.7386351, 1.2103820 ,1.5592321 ,1.8248642)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000001)
})



test_that("(ATC) PanelEstimate Runs: analytical SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  comp.results <- c(0.5917567, 1.1136940, 1.4693785, 1.8571621)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})

test_that("(ART) PanelEstimate Runs: analytical SEs ", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimate, comp.results)
  comp.results <- c(1.026000, 1.483106 ,1.919551 ,2.150267)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})


test_that("(ATT) PanelEstimate Runs: unconditional analytical SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimate, comp.results)
  comp.results <- c(0.905080467048022,1.48070275292364,1.90752332999506,2.23758074403404)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000001)
})



test_that("(ATC) PanelEstimate Runs: unconditional analytical SEs", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimate, comp.results, tolerance = .000001)
  comp.results <- c(0.7086883, 1.3301189, 1.7552049, 2.2180256)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})

test_that("(ART) PanelEstimate Runs: unconditional analytical SEs ", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimate, comp.results)
  comp.results <- c(1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})


test_that("(ATE) PanelEstimate fails: analytical SEs ", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  expect_error(PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional"))
  expect_error(PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional"))
  
})