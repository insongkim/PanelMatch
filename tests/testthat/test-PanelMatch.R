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
  compmat <- matrix(data = c(0.0601621099966881,-0.00582951116008258,0.00211544759632536,0.124463506760345,
                             0.168548575854201,0.0252020903167403,0.0540683319809353,
                             -0.117257087671973,0.0260535302880658,-0.0130941237795542)
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

test_that("Basic Network test", {
  #each odd unit is neighboring the even unit after it, each even unit is treated
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8, 10), 'treatment'] <- 1
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  ## expect proportions to all be 100% and 1 treated neighbor
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9), 'neighborhood_t_count.1'] == 1))
  
})


test_that("Variation in treated neighbors", {
  #each odd unit is neighboring the even unit after it, not everything is treated
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6), 'treatment'] <- 1
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  ## expect some to be 100 percent and 1, others to have no trated neighbors, 0%
  expect_true(all(d2[d2$id %in% c(1,3,5), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(1,3,5), 'neighborhood_t_count.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(7,9), 'neighborhood_t_prop.1'] == 0))
  
  expect_true(all(d2[d2$id %in% c(7,9), 'neighborhood_t_count.1'] == 0))
  
})

test_that("Multiple Neighbors, mix of treated", {
  ##mix of treated and untreated neighbors
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  ejmat = rbind(ejmat, data.frame(v1 = c(1,5), v2 = 3, e =1))
  ### we give 1 and 5 a second neighbor: 3, which is untreated, so they have a mix of treated/untreated neighbors
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8, 10), 'treatment'] <- 1
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  
  expect_true(all(d2[d2$id %in% c(1,5), 'neighborhood_t_prop.1'] == .5))
  
  expect_true(all(d2[d2$id %in% c(1,5), 'neighborhood_t_count.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(7,9), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(7,9), 'neighborhood_t_count.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(3), 'neighborhood_t_prop.1'] == (1/3)))
  
  expect_true(all(d2[d2$id %in% c(3), 'neighborhood_t_count.1'] == 1))
  
})



test_that("time variation", {
  ##adding some time variation across units for treated neighbor count and status
  
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8,10) & input.data$time %in% c(2,4,6,8,10), 'treatment'] <- 1
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  
  
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9) & input.data$time %in% c(2,4,6,8,10), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9) & input.data$time %in% c(2,4,6,8,10), 'neighborhood_t_count.1'] == 1))
  
  
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9) & input.data$time %in% c(1,3,5,7,9), 'neighborhood_t_prop.1'] == 0))
  
  expect_true(all(d2[d2$id %in% c(1,3,5,7,9) & input.data$time %in% c(1,3,5,7,9), 'neighborhood_t_count.1'] == 0))
  
})

test_that("more neighbor variation", {
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  ejmat = rbind(ejmat, data.frame(v1 = c(1,5), v2 = 4, e =1))
  
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8, 10), 'treatment'] <- 1
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  
  expect_true(all(d2[d2$id %in% c(1,5), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(1,5), 'neighborhood_t_count.1'] == 2))
  
  expect_true(all(d2[d2$id %in% c(3, 7,9), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(3, 7,9), 'neighborhood_t_count.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(4), 'neighborhood_t_prop.1'] == 0))
  
  expect_true(all(d2[d2$id %in% c(4), 'neighborhood_t_count.1'] == 0))
})


test_that("complex variation", {
  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)
  ejmat = rbind(ejmat, data.frame(v1 = c(1), v2 = 10, e =1))
  
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8,10) & input.data$time %in% c(2,4,6,8,10), 'treatment'] <- 1
  input.data[input.data$id %in% c(2) & input.data$time %in% c(1,3), 'treatment'] <- 1
  
  
  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')
  
  
  expect_true(all(d2[d2$id %in% c(1) & d2$time %in% c(1,3) , 'neighborhood_t_prop.1'] == .5))
  expect_true(all(d2[d2$id %in% c(1) & d2$time %in% c(1,3) , 'neighborhood_t_count.1'] == 1))
  expect_true(all(d2[d2$id %in% c(1) & d2$time %in% c(2) , 'neighborhood_t_prop.1'] == 1))
  expect_true(all(d2[d2$id %in% c(1) & d2$time %in% c(2) , 'neighborhood_t_count.1'] == 2))
  
  expect_true(all(d2[d2$id %in% c(1) & d2$time %in% c(4,6,8,10) , 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(3,5,7,9) & input.data$time %in% c(2,4,6,8,10), 'neighborhood_t_prop.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(3,5,7,9) & input.data$time %in% c(2,4,6,8,10), 'neighborhood_t_count.1'] == 1))
  
  expect_true(all(d2[d2$id %in% c(3,5,7,9) & input.data$time %in% c(1,3,5,7,9), 'neighborhood_t_prop.1'] == 0))
  
  expect_true(all(d2[d2$id %in% c(3,5,7,9) & input.data$time %in% c(1,3,5,7,9), 'neighborhood_t_count.1'] == 0))
})



