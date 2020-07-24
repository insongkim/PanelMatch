test_that("ATT (no refinement) runs", {
  PM.object <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
             treatment = "dem", refinement.method = "none",
             data = dem, match.missing = FALSE,
             size.match = 5, qoi = "att",
             outcome.var = "y",
             lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$att), "matched.set")
  expect_equal(length(PM.object), 1)
  expect_equal(length(PM.object$att), 105)
  expect_equal(names(PM.object)[1], "att")
})


test_that("ATC (no refinement) runs", {
  PM.object <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "atc",
                          outcome.var = "y",
                          lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$atc), "matched.set")
  expect_equal(length(PM.object), 1)
  expect_equal(length(PM.object$atc), 55)
  expect_equal(names(PM.object)[1], "atc")
})


test_that("ATE (no refinement) runs", {
  PM.object <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "ate",
                          outcome.var = "y",
                          lead = 0:4, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$att), "matched.set")
  expect_equal(class(PM.object$atc), "matched.set")
  expect_equal(length(PM.object), 2)
  expect_equal(length(PM.object$att), 105)
  expect_equal(length(PM.object$atc), 55)
  expect_equal(names(PM.object)[1], "att")
  expect_equal(names(PM.object)[1], "att")
})


test_that("ATT: refinement method (excluding MSM) does not affect unrefined matched sets", {
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


test_that("ATC: refinement method (excluding MSM) does not affect unrefined matched sets", {
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

  expect_equal(length(pm1$atc), 55)
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
                       outcome.var ="y", refinement.method = "mahalanobis", 
                       data = dem, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj$att, dem, covariates = c("tradewb", "rdata"), 
                            plot = FALSE, ylim = c(-2,2))
  compmat <- matrix(data = c(0.0601621099966881,-0.00601042376056592,
                             0.00212772294115459,0.104946475743932,0.133482367739366,
                             0.0252020903167403,0.051861118295163,-0.125475177405918,
                             0.0268787670968836,-0.0117184403355875)
                    , ncol = 2, nrow = 5)
  expect_equal(nrow(balmat), 5)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  c(1.40908029917124,2.55357045354071,3.31650068953231,3.93452991794896)
  expect_equivalent(pe.results$estimates, comp.results)
  
})


test_that("summary.PanelEstimate", {
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem)
  expect_output(summary(pe.results)) 
  
  
})

















