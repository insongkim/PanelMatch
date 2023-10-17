test_that("ensuring data.frame check", {
  
  d2 <- as.matrix(dem)
  expect_error(PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
             treatment = "dem", refinement.method = "none",
             data = d2, match.missing = FALSE,
             size.match = 5, qoi = "att",
             outcome.var = "y",
             lead = 0:4, forbid.treatment.reversal = FALSE))
  PM.object <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
             treatment = "dem", refinement.method = "none",
             data = dem, match.missing = FALSE,
             size.match = 5, qoi = "att",
             outcome.var = "y",
             lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_error(PanelEstimate(PM.object, data = as.matrix(dem), se.method = "bootstrap"))
  expect_error(DisplayTreatment(unit.id = "wbcode2",
                    time.id = "year", legend.position = "none",
                    xlab = "year", ylab = "Country Code",
                    treatment = "dem", data = d2))
  
})



test_that("edge case matching checks", {
  
  PM.object <- PanelMatch(lag = 1, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "att",
                          outcome.var = "y",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$att), "matched.set")
  set.seed(1)
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "bootstrap")
  expect_equivalent(pe.results$estimates, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 0.776390146867648, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "conditional")
  expect_equivalent(pe.results$estimates, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "conditional")
  expect_equivalent(pe.results$standard.error, 0.618327233714141, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "unconditional")
  expect_equivalent(pe.results$estimates, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "unconditional")
  expect_equivalent(pe.results$standard.error, 0.745063157005526, tolerance = .000001)
  
  
  #art
  PM.object <- PanelMatch(lag = 1, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "art",
                          outcome.var = "y",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$art), "matched.set")
  set.seed(1)
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "bootstrap")
  expect_equivalent(pe.results$estimates, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 1.11077012122371, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "conditional")
  expect_equivalent(pe.results$estimates, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "conditional")
  expect_equivalent(pe.results$standard.error, 0.90258389202935, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "unconditional")
  expect_equivalent(pe.results$estimates, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "unconditional")
  expect_equivalent(pe.results$standard.error, 1.21738691542365, tolerance = .000001)
  
  
})





test_that("ATT: refinement method does not affect unrefined matched sets", {
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "mahalanobis",
                          covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "att",
                          outcome.var = "y",
                          lead = 0:4, forbid.treatment.reversal = FALSE)
  
  pm2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "ps.weight",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm3 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "ps.match", covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    data = dem, match.missing = FALSE,
                    size.match = 5, qoi = "att",
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm4 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "CBPS.weight",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm5 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "CBPS.match",
                    data = dem, match.missing = FALSE,covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_equal(length(pm1$att), 105)
  expect_equivalent(pm1$att, pm2$att)
  expect_equivalent(pm2$att, pm3$att)
  expect_equivalent(pm3$att, pm4$att)
  expect_equivalent(pm4$att, pm5$att)

  
})


test_that("ATC: refinement method does not affect unrefined matched sets", {
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  pm2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "ps.weight",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm3 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "ps.match",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm4 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "CBPS.weight",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm5 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "CBPS.match",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:4, forbid.treatment.reversal = FALSE)

  expect_equal(length(pm1$atc), 3287)
  expect_equivalent(pm1$atc, pm2$atc)
  expect_equivalent(pm2$atc, pm3$atc)
  expect_equivalent(pm3$atc, pm4$atc)
  expect_equivalent(pm4$atc, pm5$atc)
  
  
})


test_that("matched.set methods work", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  sets <- pm1$att  
  expect_equal(class(sets), "matched.set")
  expect_equal(length(summary(sets)), 5)
  expect_true("set.size.summary" %in% names(summary(sets)))
  expect_equal(names(sets)[1], "4.1992")
  hardcoded.set <- c(3,13,16,19,28,29,31,35,36,37,43,45,47,51,53,57,62,64,65,67,70,71,81,84,87,93,95,
  96,97,103,104,105,109,110,112,114,115,116,118,123,124,125,128,129,134,140,142,150,155,
  156,157,159,161,163,168,171,172,173,175,176,178,179,180,182,184,186,187,190,193,196,197,199,200,202)
  expect_equal(as.numeric(sets[[1]]), hardcoded.set)
  subset <- sets[1:3]
  expect_equal(class(subset), "matched.set") 
})


test_that("balance checking functions are sensible", {
  set.seed(1)
  dem$rdata <- runif(runif(nrow(dem)))
  
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, time.id = "year", unit.id = "wbcode2", treatment = "dem",
                       outcome.var = "y", refinement.method = "mahalanobis", 
                       data = dem, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj$att, dem, covariates = c("tradewb", "rdata"), 
                            plot = FALSE, ylim = c(-2,2))
  compmat <- matrix(data = c(0.0601621099966881,-0.00601042376056592,
                             0.00212772294115459,0.104946475743932,
                             0.0252020903167403,0.051861118295163,-0.125475177405918,
                             0.0268787670968836), 
                    ncol = 2, nrow = 4)
  expect_equal(nrow(balmat), 4)
  expect_equal(ncol(balmat), 2)
  expect_equivalent(balmat, compmat)
  
  
  balmat <- get_covariate_balance(pm.obj$att, dem, covariates = c("tradewb", "rdata"), 
                                  plot = FALSE, ylim = c(-2,2), use.equal.weights = TRUE)
  
  
  expect_false(isTRUE(all.equal(balmat, compmat, check.attributes = FALSE)))
  
})


test_that("(ATT) PanelEstimate Runs", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimates, comp.results)
})



test_that("(ATC) PanelEstimate Runs", {
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .000001)
  
})

test_that("(ART) PanelEstimate Runs", {
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  
})


test_that("(ATE) PanelEstimate Runs", {
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, number.iterations = 300)
  comp.results <-  c(-0.73400424, -0.14920097, -0.45276678, -0.06580353)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .0000001)
  
})


test_that("(ATT) bootstrap SEs", {
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, number.iterations = 100)
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(0.8701336, 1.3611202, 1.7512925, 2.0019833), tolerance = .0000001)
})



test_that("(ATC) bootstrap SEs", {
  qoi_ <- "atc"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, number.iterations = 100)
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .0000001)
  expect_equivalent(pe.results$standard.error, c(0.7027994, 1.3144841, 1.7200521, 2.1352594 ), tolerance = .0000001)
})

test_that("(ART) bootstrap SEs", {
  qoi_ <- "art"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, number.iterations = 100)
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(1.342158, 2.159589, 2.869549, 3.231702), tolerance = .000001)
})


test_that("(ATE) bootstrap SEs", {
  set.seed(1)
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, number.iterations = 100)
  comp.results <-  c(-0.73400424, -0.14920097, -0.45276678, -0.06580353)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .000001)
  expect_equivalent(pe.results$standard.error, c(0.6953794, 1.2931504, 1.6894504, 2.0902084), tolerance = .0000001)
})



test_that("(ATT) PanelEstimate Runs: analytical SEs", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(0.7386351, 1.2103820 ,1.5592321 ,1.8248642)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000001)
})



test_that("(ATC) PanelEstimate Runs: analytical SEs", {
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .000001)
  comp.results <- c(0.5917567, 1.1136940, 1.4693785, 1.8571621)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})

test_that("(ART) PanelEstimate Runs: analytical SEs ", {
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(1.026000, 1.483106 ,1.919551 ,2.150267)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equal(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})


test_that("(ATT) PanelEstimate Runs: unconditional analytical SEs", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(0.905080467048022,1.48070275292364,1.90752332999506,2.23758074403404)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000001)
})



test_that("(ATC) PanelEstimate Runs: unconditional analytical SEs", {
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  comp.results <-  c(-0.7399887, -0.1418777, -0.4914594, -0.1423150)
  expect_equivalent(pe.results$estimates, comp.results, tolerance = .000001)
  comp.results <- c(0.7086883, 1.3301189, 1.7552049, 2.2180256)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})

test_that("(ART) PanelEstimate Runs: unconditional analytical SEs ", {
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423)
  names(comp.results) <- paste0("t+", 0:3)
  expect_equivalent(pe.results$standard.error, comp.results, tolerance = .0000002)
  
})


test_that("(ATE) PanelEstimate fails: analytical SEs ", {
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  expect_error(PanelEstimate(pm1, data = dem, se.method = "conditional"))
  expect_error(PanelEstimate(pm1, data = dem, se.method = "unconditional"))
  
})

test_that("summary.PanelEstimate (conditional)", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.738635065855285,1.21038195434255,1.55923210008201,1.82486418813223,-2.04109789825896,-2.69356525042577,-2.59972747285187,-1.84484646286253,0.854298355330497,2.05104482567144,3.51235004654712,5.30848970796966)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "att")
  
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.591756652995951,1.11369402788093,1.46937851945315,1.85716210624063,-1.89981045123073,-2.32467787575185,-3.37138841788079,-3.78228584150992,0.419833003737329,2.04092249313608,2.38846953768912,3.49765584185843)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "atc")
  
  
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "conditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.02600037922772,1.48310644776257,1.91955116019568,2.15026718759732,-7.22864265630046,-10.9282208645128,-12.5187199096139,-12.3384409599025,-3.20679507347894,-5.11455041880524,-4.99421762868267,-3.90954847024455)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "art")
})

test_that("summary.PanelEstimate (unconditional)", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.905080467048022,1.48070275292364,1.90752332999506,2.23758074403404,-2.36732488998905,-3.22338427991681,-3.2823657396126,-2.65375604825349,1.18052534706058,2.58086385516248,4.19498831330785,6.11739929336062)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "att")
  
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.708688258862262,1.33011893273585,1.75520486661622,2.21802563150778,-2.12899218738313,-2.74886289462501,-3.93159776415305,-4.4895653543677,0.649014739889732,2.46510751200924,2.94867888396138,4.20493535471621)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "atc")
  
  
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem, se.method = "unconditional")
  expect_output(summary(pe.results)) 
  expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423,-8.29736530389636,-12.5517704314982,-14.369074947577,-14.1843445061431,-2.13807242588304,-3.49100085181983,-3.14386259071958,-2.06364492400392)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "art")
  
  
})


test_that("summary.PanelEstimate (bootstrap)", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  set.seed(1)

  pe.results <- PanelEstimate(pm1, data = dem)
  
  expect_output(summary(pe.results)) 
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL
  cmat <- unlist(cmat)
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.906705857062448,1.46876653085356,1.87876179592695,2.21608840568008,-2.47319633748927,-3.25290871525436,-3.34863750081872,-2.68949846129008,1.12752178115124,2.43516177434591,4.04890571881315,5.89171531990835)
  expect_equal(cmat, comp.vec, tolerance = .0000001)
  
})

test_that("set level treatment effects", {
  
  PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                           treatment = "dem", refinement.method = "ps.match",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "att",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, data = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem, se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                           treatment = "dem", refinement.method = "ps.match",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "art",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, data = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem, se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)
  
  PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                           treatment = "dem", refinement.method = "ps.match",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "atc",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, data = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem, se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)

  
  
})

test_that("test time data checker", {
  
  good <- data.frame(id = as.integer(unlist(lapply(1:5, function(x) rep(x, 20)))),
                     time.in = rep(as.integer(seq(1,20,by=1)), 5))
  
  ok <- data.frame(id = unlist(lapply(1:5, function(x) rep(x, 20))),
                   time.in = rep(seq(1,20,by=1), 5))
  
  
  bad <- data.frame(id = unlist(lapply(1:5, function(x) rep(x, 10))),
                    time.in = rep(seq(1,20,by=2), 5))
  
  
  fail <- data.frame(id = as.integer(unlist(lapply(1:5, function(x) rep(x, 20)))),
                     time.in = "A")
  
  r1 <- PanelMatch:::check_time_data(good, "time.in")
  r1 <- expect_warning(PanelMatch:::check_time_data(ok, "time.in"))
  r1 <- expect_warning(PanelMatch:::check_time_data(bad, "time.in"))
  r1 <- expect_error(PanelMatch:::check_time_data(fail, "time.in"))
})

test_that("test placebo test", {
  
  PM.results <- PanelMatch(lag = 4, 
                           time.id = "year",
                           unit.id = "wbcode2",
                           treatment = "dem", 
                           refinement.method = "mahalanobis",
                           data = dem,
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "art",
                           outcome.var = "y", lead = 0:4, 
                           forbid.treatment.reversal = FALSE,
                           placebo.test = TRUE)
  set.seed(1)
  pt.res <- placebo_test(PM.results, 
               data = dem, 
               number.iterations = 100, 
               se.method = "bootstrap",
               plot = FALSE)
  
  est.comps <- c(20.438656,
                 4.110167,
                 -21.562320)
  
  st.comps <- c(28.54986,
                23.67567,
                17.30092)
  expect_equivalent(pt.res$estimates, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  PM.results <- PanelMatch(lag = 4, 
                           time.id = "year",
                           unit.id = "wbcode2",
                           treatment = "dem", 
                           refinement.method = "mahalanobis",
                           data = dem,
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "att",
                           outcome.var = "y", lead = 0:4, 
                           forbid.treatment.reversal = FALSE,
                           placebo.test = TRUE)
  set.seed(1)
  pt.res <- placebo_test(PM.results, 
               data = dem, 
               number.iterations = 100, 
               se.method = "bootstrap",
               plot = FALSE)
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  
  st.comps <- c(18.33094, 
                17.76234,
                13.37736)
  
  expect_equivalent(pt.res$estimates, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  
  pt.res <- placebo_test(PM.results, 
                         data = dem, 
                         se.method = "conditional",
                         plot = FALSE)
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  
  st.comps <- c(16.30824,
                15.11819,
                11.29278)
  
  expect_equivalent(pt.res$estimates, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  pt.res <- placebo_test(PM.results, 
                         data = dem, 
                         se.method = "unconditional",
                         plot = FALSE)
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  st.comps <- c(19.87698,
                18.35898,
                13.72441)
  expect_equivalent(pt.res$estimates, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
})

test_that("max, numeric, adjusting caliper", {
  
  #####Next few tests is just changing the caliper up to see if the refinement responds appropriately
  ## using max method, numeric data
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  ###How does caliper work? 
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"max", 1, "numeric", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3)))
  expect_true(all(PM.results$att[[2]] == c(3,5)))
  expect_true(all(PM.results$att[[3]] == c(5, 7)))
  
  
  ####*****#####******
  
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  ###How does caliper work? 
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"max", 2, "numeric", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3)))
  expect_true(all(PM.results$att[[2]] == c(3,5)))
  expect_true(all(PM.results$att[[3]] == c(5,7,8)))
  
  ####*****#####******
  
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  ###How does caliper work? 
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"max", 3, "numeric", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3, 5)))
  expect_true(all(PM.results$att[[2]] == c(1, 3, 5, 7)))
  expect_true(all(PM.results$att[[3]] == c(3, 5, 7, 8, 9)))
  
  
})


test_that("average, numerical", {
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1
  input.data$cal.data <- input.data$id
  input.data[input.data$id %in% c(1), 'cal.data'] <- c(rep(1, 5), rep(0, 5))
  input.data[input.data$id %in% c(3), 'cal.data'] <- c(rep(3, 5), rep(2, 5))
  
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"average", 1, "numeric", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(3)))
  expect_true(all(PM.results$att[[2]] == c(5)))
  expect_true(all(PM.results$att[[3]] == c(5, 7)))
})



test_that("categorical, max", {
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 #created treated units
  input.data$cal.data <- input.data$id
  input.data[input.data$id %in% c(1,3), 'cal.data'] <- 2
  input.data[input.data$id %in% c(5), 'cal.data'] <- 4
  input.data[input.data$id %in% c(7,8,9,10), 'cal.data'] <- 6
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"max", 1, "categorical", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3)))
  expect_true(all(PM.results$att[[2]] == c(5)))
  expect_true(all(PM.results$att[[3]] == c(7,8,9, 10) ))
})



test_that("categorical, average", {
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 #created treated units
  input.data$cal.data <- input.data$id
  input.data[input.data$id %in% c(1,3), 'cal.data'] <- 2
  input.data[input.data$id %in% c(1,3) & input.data$time %in% 4:5, 'cal.data'] <- 3 #not 2
  
  input.data[input.data$id %in% c(5), 'cal.data'] <- 4
  input.data[input.data$id %in% c(7,8,9,10), 'cal.data'] <- 6
  input.data[input.data$id %in% c(7,8) & input.data$time %in% 3:5, 'cal.data'] <- 7 #not six, these should drop out
  
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           caliper.formula = ~ I(caliper(cal.data,"average", .5, "categorical", "raw")) )
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3)))
  expect_true(all(PM.results$att[[2]] == c(5)))
  expect_true(all(PM.results$att[[3]] == c(9, 10) ))
})


test_that("Testing Continuous Matching: basic, att", {
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7, 8, 9, 10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7, 8,9,10)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8,9,10)))
})

test_that("Continuous Matching, att, unmatchable controls", {
  
  # making two control units too far away to be matched with anything
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(9, 10) & input.data$time < 5, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8)))
})


test_that("continuous matching, att, various exceptions", {
  
  # making one unit have one time period that will violate the caliper
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time == 3, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8,9)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8,9)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8,9)))
  
  
  
  
  ####*****#####*****
  # making one unit have one time period that will violate the caliper
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time == 3, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8,9)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8,9)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8,9)))
  
  
  
  
  # making one unit have one time period that will be different but still be included in the matching
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time == 3, 'treatment'] <- .2
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8,9, 10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8,9, 10)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8,9, 10)))
  
  
  ####*****#####*****
  # making one unit have one time period that will violate the caliper, one that has different values but still will be matched
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time %in% 1:3, 'treatment'] <- .5
  input.data[input.data$id %in% c(9) & input.data$time %in% 2:3, 'treatment'] <- 52
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8,10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8,10)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8,10)))
  
})
test_that("checking new definition of ATT (continuous) explicitly", {
  
  ####*****#####*****
  # CHECKING NEW DEFINITION OF ATT
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- 1 + .43
  
  input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- -1.43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time %in% 1:3, 'treatment'] <- .5
  input.data[input.data$id %in% c(9) & input.data$time %in% 2:3, 'treatment'] <- 52
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5, 
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 2) 
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8,10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8,10)))
  
  
  
  
  ####*****#####*****
  # changing some numbers around to make sure the thresholds work as intended
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1.2
  input.data[input.data$id %in% c(2,4,6) & input.data$time <= 5, 'treatment'] <- 1.1
  #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  
  input.data[input.data$id %in% c(1,3), 'treatment'] <- .5
  input.data[input.data$id %in% c(5,7), 'treatment'] <- 2
  input.data[input.data$id %in% c(9), 'treatment'] <- 52
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .05,  control.threshold = .05,
                                    units = "raw", matching.threshold = 1) #small threshold for att matching 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7))) ## 9 has giant fluctuation, 8 and 10 are all zeroes which are too far away by the threshold. 
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7)))
  
  
  
  #########*************************************************
  ## including some control units in "both directions"
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 10
  
  input.data[input.data$id %in% c(1,3,5,7) & input.data$time %in% 1:10, 'treatment'] <- 1
  input.data[input.data$id %in% c(8) & input.data$time %in% 1:10, 'treatment'] <- -1
  input.data[input.data$id %in% c(9, 10) & input.data$time %in% 1:10, 'treatment'] <- 4
  
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5,  
                                    units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 3)
  expect_true(all(PM.results$att[[1]] == c(1,3,5,7,8))) # 9 and 10 are now too far away
  expect_true(all(PM.results$att[[2]] == c(1,3,5,7,8)))
  expect_true(all(PM.results$att[[3]] == c(1,3,5,7,8)))
  
  
  ##################################################################################################################
  ######################################testing the matching threshold
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- .6
  input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- .4
  #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  input.data[input.data$id %in% c(1,3,5, 6, 7, 8, 9, 10), 'treatment']  <- 100
  
  continuous.treatment.info <- list(treatment.threshold = .5,  control.threshold = .5,  
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 2) 
  expect_true(all(c("2.6","4.6") %in% names(PM.results$att)))
  
  expect_true(all(PM.results$att[[1]] == c(1,3,5, 6, 7, 8, 9, 10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,6, 7, 8,9,10))) #six is no longer a treated unit so it will be matched like everything else  
})




# test_that("empty sets exist, continuous treatment", {
#   
#   ####*****#####*****
#   input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
#   input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
#   input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- .6
#   input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- .4
#   #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
#   
#   
#   input.data$cal.data <- input.data$id
#   input.data$outcome <- rnorm(nrow(input.data))
#   
#   input.data[input.data$id %in% c(1,3,5, 6, 7, 8, 9, 10), 'treatment']  <- 100
#   
#   continuous.treatment.info <- list(treatment.threshold = .5, control.threshold = .5, 
#                                     method = "max", units = "raw", matching.threshold = 2) #include everything 
#   
#   PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
#                            treatment = "treatment", refinement.method = "none", # should be none for all of them
#                            data = input.data, match.missing = TRUE, 
#                            size.match = 5, qoi = "att" , outcome.var = "outcome",
#                            lead = 0:4, forbid.treatment.reversal = FALSE,
#                            continuous.treatment.info = continuous.treatment.info)
#   
#   
# })





test_that("testing directionality (ATT)", {
  ## testing directionality
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,6) & input.data$time > 5, 'treatment'] <- (1 + .43)
  input.data[input.data$id %in% c(4) & input.data$time > 5, 'treatment'] <- (1 + .43) * -1 
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time == 3, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, 
                                    units = "raw", control.threshold = .5, 
                                    matching.threshold = 0) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(all(names(PM.results$att) == c("2.6", "6.6")))
  
  continuous.treatment.info <- list(treatment.threshold = .5, 
                                    units = "raw", control.threshold = .5, 
                                    matching.threshold = 0) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(all(names(PM.results$att) == c("2.6", "6.6")))
  
  
})
# 
# commenting out until we are sure about definitions
test_that("simple brute force match? (continuous)", {

  brute_force_did_continuous <- function(matched.set, data, treatment.variable,
                                         outcome.variable, F_lead, time, id)
  {
    all.data <- as.numeric(unlist(strsplit(names(matched.set), split = "[.]")))
    tids <- all.data[seq(from = 1, to = length(all.data), by = 2)]
    ts <- all.data[seq(from = 2, to = length(all.data), by = 2)]
    set.dids <- list()
    for (i in 1:length(tids)) {

      t0 <- ts[i] + F_lead
      true.t0 <- ts[i] + 0
      controls <- matched.set[[i]]
      controls <- controls[attr(controls, "weights") > 0]

      checkdf1 <- data[data[, id] %in% controls & data[, time] == (ts[i] - 1), ]
      checkdf2 <- data[data[, id] %in% controls & data[, time] == t0, ]
      treated1 <- data[data[, id] == tids[i] & data[, time] == (ts[i] - 1), ]
      treated2 <- data[data[, id] == tids[i] & data[, time] == t0, ]

      std1 <- data[data[, id] == tids[i] & data[, time] == (ts[i] - 1), ]
      std2 <- data[data[, id] == tids[i] & data[, time] == true.t0, ]

      d1 <- sum((checkdf2[,outcome.variable] - checkdf1[, outcome.variable]) * attr(matched.set[[i]], "weights")[attr(matched.set[[i]], "weights") > 0])

      standardized.denom <- std2[,treatment.variable] - std1[, treatment.variable]

      d2 <- treated2[,outcome.variable] - treated1[, outcome.variable]

      set.dids[i] <- (d2 - d1) / standardized.denom

      #set.dids[i] <- treated[,outcome.variable] - (mean(treated[,outcome.variable] - checkdf[, outcome.variable])) ## should always be of length one
    }

    return(mean(unlist(set.dids), na.rm = T))
  }


  continuous.treatment.info <- list(treatment.threshold = 1,  control.threshold = 1,
                                    units = "raw", matching.threshold = 20)


  sdt <- data.frame(time = rep(1:10, 10),
                    id = unlist(sapply(1:10, function(x) rep(x, 10), simplify = FALSE)),
                    treatment = 0, outcome = 4)

  sdt$treatment[5:10] <- 2
  sdt$outcome[5:10] <- 5

  sdt$treatment[16:20] <- 3
  sdt$outcome[16:20] <- 6

  sdt$treatment[27:30] <- 4
  sdt$outcome[27:30] <- 7


  s <- PanelMatch(lag = 4, time.id = "time", unit.id = "id",
                  treatment = "treatment", refinement.method = "none", # should be none for all of them
                  data = sdt, match.missing = TRUE,
                  size.match = 5, qoi = "att" , outcome.var = "outcome",
                  lead = 0, forbid.treatment.reversal = FALSE,
                  continuous.treatment.info = continuous.treatment.info)

  brute.force.result <- brute_force_did_continuous(matched.set = s$att, data = sdt,
                                                   treatment.variable = 'treatment',
                                                   outcome.variable = "outcome",
                                                   F_lead = 0, time = "time", id = "id")


  pe <- PanelEstimate(sets = s, data = sdt)

  expect_equivalent(brute.force.result, pe$estimates)

  sdt$refine <- rnorm(n = nrow(sdt))


  s <- PanelMatch(lag = 2, time.id = "time", unit.id = "id",
                  treatment = "treatment", refinement.method = "ps.weight", # should be none for all of them
                  data = sdt, match.missing = TRUE, covs.formula = ~ I(lag(refine, 1:2)),
                  size.match = 3, qoi = "att" , outcome.var = "outcome",
                  lead = 0, forbid.treatment.reversal = FALSE,
                  continuous.treatment.info = continuous.treatment.info)

  brute.force.result <- brute_force_did_continuous(matched.set = s$att, data = sdt,
                                                   treatment.variable = 'treatment',
                                                   outcome.variable = "outcome",
                                                   F_lead = 0, time = "time", id = "id")


  pe <- PanelEstimate(sets = s, data = sdt)

  expect_equivalent(brute.force.result, pe$estimates)


})