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
  
  
  
  #atc
  PM.object <- PanelMatch(lag = 1, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "atc",
                          outcome.var = "y",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$atc), "matched.set")
  set.seed(1)
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "bootstrap")
  expect_equivalent(pe.results$estimates, 4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 1.11077012122371, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "conditional")
  expect_equivalent(pe.results$estimates, 4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "conditional")
  expect_equivalent(pe.results$standard.error, 0.90258389202935, tolerance = .000001)
  
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "unconditional")
  expect_equivalent(pe.results$estimates, 4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "unconditional")
  expect_equivalent(pe.results$standard.error, 1.21738691542365, tolerance = .000001)
  
  #ate
  PM.object <- PanelMatch(lag = 1, time.id = "year", unit.id = "wbcode2",
                          treatment = "dem", refinement.method = "none",
                          data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "ate",
                          outcome.var = "y",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$att), "matched.set")
  expect_equal(class(PM.object$atc), "matched.set")
  set.seed(1)
  pe.results <- PanelEstimate(PM.object, data = dem, se.method = "bootstrap")
  expect_equivalent(pe.results$estimates, 0.766157910045996, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 0.605100986438279, tolerance = .000001)
  
  
  
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
                       outcome.var = "y", refinement.method = "mahalanobis", 
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

test_that("(ART) PanelEstimate Runs", {
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem)
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

test_that("(ATT) bootstrap SEs", {
  qoi_ <- "att"
  set.seed(1)
  pm1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                    treatment = "dem", refinement.method = "mahalanobis",
                    data = dem, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    outcome.var = "y",
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(0.9067059,1.4687665,1.8787618,2.2160884), tolerance = .0000001)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(1.436841, 2.279820, 3.030895, 3.451612), tolerance = .0000001)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  -c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(1.436841, 2.279820, 3.030895, 3.451612), tolerance = .0000001)
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
  
  pe.results <- PanelEstimate(pm1, data = dem)
  comp.results <-  c(1.40908029917124,2.55357045354071,3.31650068953231,3.93452991794896)
  expect_equivalent(pe.results$estimates, comp.results)
  expect_equivalent(pe.results$standard.error, c(0.767560957486421,1.17704405447648,
                                                 1.50237292151594,1.76640885892815), tolerance = .0000001)
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
  comp.results <-  c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(1.026000, 1.483106 ,1.919551 ,2.150267)
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
  comp.results <-  c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353)
  expect_equivalent(pe.results$estimates, comp.results)
  comp.results <- c(1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423)
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
  comp.vec <- c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353,1.02600037922772,1.48310644776257,1.91955116019568,2.15026718759732,3.20679507347894,5.11455041880524,4.99421762868267,3.90954847024455,7.22864265630046,10.9282208645128,12.5187199096139,12.3384409599025)
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
  comp.vec <- c(5.2177188648897,8.02138564165901,8.75646876914828,8.12399471507353,1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423,2.13807242588304,3.49100085181983,3.14386259071958,2.06364492400392,8.29736530389636,12.5517704314982,14.369074947577,14.1843445061431)
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
                           treatment = "dem", refinement.method = "mahalanobis",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "att",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- getSetTreatmentEffects(pm.obj = PM.results, data.in = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem)
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)
  
  
  PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                           treatment = "dem", refinement.method = "mahalanobis",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "art",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- getSetTreatmentEffects(pm.obj = PM.results, data.in = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem)
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)
  
  PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
                           treatment = "dem", refinement.method = "mahalanobis",
                           data = dem, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "atc",
                           outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- getSetTreatmentEffects(pm.obj = PM.results, data.in = dem, lead = 0:1)
  #mean(set.effects[[1]], na.rm = TRUE)
  #mean(set.effects[[2]], na.rm = TRUE)
  pe.results <- PanelEstimate(PM.results, data = dem)
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimates[1:2], tolerance = .000001)
  
  
  # PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
  #                          treatment = "dem", refinement.method = "mahalanobis",
  #                          data = dem, match.missing = TRUE,
  #                          covs.formula = ~ I(lag(tradewb, 1:4)),
  #                          size.match = 5, qoi = "ate",
  #                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
  #                          placebo.test = FALSE)
  # set.effects <- getSetTreatmentEffects(pm.obj = PM.results, data.in = dem, lead = 0:1)
  # #mean(set.effects[[1]], na.rm = TRUE)
  # #mean(set.effects[[2]], na.rm = TRUE)
  # pe.results <- PanelEstimate(PM.results, data = dem)
  # expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
  #                     mean(set.effects[[2]], na.rm = TRUE)), 
  #                   pe.results$estimates[1:2], tolerance = .000001)
  
  
})


