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
  compmat <- matrix(data = c(0.152972260454625,0.0923632722068149,0.109795148587497,0.249427241400355,0.295472804606837,0.166692980558207,
                             0.103662692289609,0.00896485955293765,0.048526876268204,0.0604550292871631), ncol = 2, nrow = 5)
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
  comp.results <-  c(-0.903871831205692,-0.497466128634409,0.403842635007254,1.24893191622705)
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
  comp.results <-  c(5.21729365330116,7.97913255878523,8.6319948682598,7.93186896829045)
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
  comp.results <-  c(1.20544870737437,2.42352395959804,3.23921941808752,3.55183590038403)
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

















