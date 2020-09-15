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


test_that("network caliper and refinement tests", {
  ### verify that formulas are updated:

  network.refinement.info = list(use.proportion.data = TRUE,
                                 proportion.lags = 0:4,
                                 use.count.data = TRUE,
                                 count.lags = 0:4)
  network.caliper.info = list(use.proportion.data = TRUE,
                              proportion.caliper.method = "average",
                              proportion.caliper.threshold = .5,
                              prop.unit.type = "raw",
                              use.count.data = TRUE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 1, 
                              count.unit.type = "raw")


  ejmat = data.frame(v1 = seq(from = 1, to = 10, by = 2), v2 = seq(from = 2, to = 10, by = 2), e = 1)

  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6,8, 10), 'treatment'] <- 1
  input.data$outcome <- rnorm(nrow(input.data))
  input.data$cal.data <- input.data$id

  d2 <- PanelMatch:::calculate_neighbor_treatment(input.data, ejmat, 1, 'id', 'time', 'treatment')

  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           NULL, NULL)
  
  covs.formula <- ll[[1]]
  # expect_true(identical(covs.formula, ~I(lag(neighborhood_t_prop.1, 0:4)) + I(lag(neighborhood_t_count.1, 
  #                                                                                 0:4))))
  expect_true(covs.formula == ~I(lag(neighborhood_t_prop.1, 0:4)) + I(lag(neighborhood_t_count.1, 0:4)))
  
  caliper.formula <- ll[[2]]
  expect_true(caliper.formula == ~ I(caliper(neighborhood_t_prop.1, "average", 0.5, "numeric", "raw")) + I(caliper(neighborhood_t_count.1, "max", 1, "numeric", "raw")))

  ##Unclear why identical fails but == works...guessing because of environment which is not as important right now.
  network.refinement.info = list(use.proportion.data = FALSE,
                                 proportion.lags = 0,
                                 use.count.data = TRUE,
                                 count.lags = 0:2)
  network.caliper.info = list(use.proportion.data = TRUE,
                              proportion.caliper.method = "max",
                              proportion.caliper.threshold = .5,
                              prop.unit.type = "raw",
                              use.count.data = TRUE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 1,
                              count.unit.type = "raw")

  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           ~ outcome, NULL)
  covs.formula <- ll[[1]]
  caliper.formula <- ll[[2]]


  covs.formula <- ll[[1]]
  expect_true(identical(covs.formula, ~ outcome + I(lag(neighborhood_t_count.1, 0:2))))
  caliper.formula <- ll[[2]]
  expect_true(identical(caliper.formula, ~I(caliper(neighborhood_t_prop.1, "max", 0.5, "numeric", "raw")) + I(caliper(neighborhood_t_count.1, "max", 1, "numeric",
                                                                                                                          "raw"))))




  network.refinement.info = list(use.proportion.data = FALSE,
                                 proportion.lags = 0,
                                 use.count.data = FALSE,
                                 count.lags = 2)
  network.caliper.info = list(use.proportion.data = TRUE,
                              proportion.caliper.method = "average",
                              proportion.caliper.threshold = .2,
                              prop.unit.type = "raw",
                              use.count.data = TRUE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 4,
                              count.unit.type = "sd")

  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           ~ outcome,
                                                           ~ I(caliper(cal.data,"average", .5, "categorical", "raw")))


  covs.formula <- ll[[1]]
  expect_true(identical(covs.formula, ~ outcome))
  caliper.formula <- ll[[2]]
  expect_true(identical(caliper.formula, ~I(caliper(cal.data,"average", .5, "categorical", "raw")) + I(caliper(neighborhood_t_prop.1, "average", 0.2, "numeric", "raw")) + I(caliper(neighborhood_t_count.1, "max", 4, "numeric",
                                                                                                                      "sd"))))


  network.refinement.info = list(use.proportion.data = TRUE,
                                 proportion.lags = 0:4,
                                 use.count.data = TRUE,
                                 count.lags = 0:4)
  network.caliper.info = list(use.proportion.data = FALSE,
                              proportion.caliper.method = "average",
                              proportion.caliper.threshold = .2,
                              prop.unit.type = "raw",
                              use.count.data = FALSE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 4,
                              count.unit.type = "sd")
  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           NULL,
                                                           ~ I(caliper(cal.data,"average", .5, "categorical", "raw")))

  covs.formula <- ll[[1]]
  expect_true(covs.formula == ~I(lag(neighborhood_t_prop.1, 0:4)) + I(lag(neighborhood_t_count.1,
                                                                                  0:4)))
  caliper.formula <- ll[[2]]
  expect_true(identical(caliper.formula, ~ I(caliper(cal.data,"average", .5, "categorical", "raw"))))



  network.refinement.info = list(use.proportion.data = TRUE,
                                 proportion.lags = 0:2,
                                 use.count.data = TRUE,
                                 count.lags = 0:2)
  network.caliper.info = list(use.proportion.data = FALSE,
                              proportion.caliper.method = "average",
                              proportion.caliper.threshold = .2,
                              prop.unit.type = "raw",
                              use.count.data = TRUE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 4,
                              count.unit.type = "sd")
  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           ~ outcome,
                                                           ~ I(caliper(cal.data,"max", .2, "numeric", "sd")))
  covs.formula <- ll[[1]]
  expect_true(identical(covs.formula, ~outcome + I(lag(neighborhood_t_prop.1, 0:2)) + I(lag(neighborhood_t_count.1,
                                                                                  0:2))))
  caliper.formula <- ll[[2]]
  expect_true(identical(caliper.formula, ~ I(caliper(cal.data,"max", .2, "numeric", "sd"))+ I(caliper(neighborhood_t_count.1, "max", 4, "numeric", "sd"))))



  network.refinement.info = list(use.proportion.data = TRUE,
                                 proportion.lags = 0:4,
                                 use.count.data = TRUE,
                                 count.lags = 0:4)
  network.caliper.info = list(use.proportion.data = TRUE,
                              proportion.caliper.method = "average",
                              proportion.caliper.threshold = .2,
                              prop.unit.type = "raw",
                              use.count.data = FALSE,
                              count.caliper.method = "max",
                              count.caliper.threshold = 4,
                              count.unit.type = "sd")
  ll <- PanelMatch:::handle_network_caliper_and_refinement(network.caliper.info, network.refinement.info, d2,
                                                           ejmat, 1, 'id', 'time', 'treatment',
                                                           ~ outcome,
                                                           ~ I(caliper(cal.data, "max", .2, "categorical", "raw")))
  covs.formula <- ll[[1]]
  expect_true(identical(covs.formula, ~outcome + I(lag(neighborhood_t_prop.1, 0:4)) + I(lag(neighborhood_t_count.1,
                                                                                            0:4))))
  caliper.formula <- ll[[2]]
  expect_true(identical(caliper.formula, ~ I(caliper(cal.data,"max", .2, "categorical", "raw"))+ I(caliper(neighborhood_t_prop.1, "average", .2, "numeric", "raw"))))
})


test_that("Testing Continuous Matching", {
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  ####*****#####******
  
  # making two control units too far away to be matched with anything
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(9, 10) & input.data$time < 5, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  
  ####*****#####******
  
  # making one unit have one time period that will violate the caliper
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1 + .43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time == 3, 'treatment'] <- 5
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  
  
  ####*****#####*****
  # some treated units are only treated under atc conditions
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- 1 + .43
  
  input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- -1.43
  input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data[input.data$id %in% c(10) & input.data$time %in% 1:3, 'treatment'] <- .5
  input.data[input.data$id %in% c(9) & input.data$time %in% 2:3, 'treatment'] <- 52
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "att" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$att) == 2) # six is now included in the matched sets because we are only looking from 2,3,4,5, t = 6
  expect_true(all(PM.results$att[[1]] == c(1,3,5,6, 7,8,10)))
  expect_true(all(PM.results$att[[2]] == c(1,3,5,6, 7,8,10)))
  
  
  
  
  ####*****#####*****
  # changing some numbers around to make sure the thresholds work as intended
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  
  input.data[input.data$id %in% c(2,4,6) & input.data$time > 5, 'treatment'] <- 1.2
  input.data[input.data$id %in% c(2,4,6) & input.data$time <= 5, 'treatment'] <- 1.1
  #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  vals <- runif(n = 4, min = .1, max = 2.1) #these actual values should not matter because of the distribution
  input.data[input.data$id %in% c(1,3,5,7) & input.data$time %in% 2:5, 'treatment'] <- unlist(sapply(vals, rep, times = 4, simplify = FALSE))  #none of these should ever violate matching restriction
  input.data[input.data$id %in% c(9) & input.data$time %in% 2:3, 'treatment'] <- 52
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  
  continuous.treatment.info <- list(treatment.threshold = .05, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 1) #small threshold for att matching 
  
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
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
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
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
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
  
  #TODO: remove the extra list items, add safeguards into the code to check for this...
  
  
  ##################################################################################################################
  ######################################testing ATC functionality
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- -.6
  input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- .4
  #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  continuous.treatment.info <- list(treatment.threshold = .5, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "atc" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  
  expect_true(length(PM.results$atc) == 2) 
  expect_true(all(c("2.6","4.6") %in% names(PM.results$atc)))
  
  expect_true(all(PM.results$atc[[1]] == c(1,3,5, 6, 7, 8, 9, 10)))
  expect_true(all(PM.results$atc[[2]] == c(1,3,5,6, 7, 8,9,10))) #six is no longer a treated unit so it will be matched like everything else
  
  
  
  ##################################################################################################################
  ######################################testing ATE functionality
  
  input.data = data.frame(id = rep(1:10, 10), time = unlist(lapply(1:10, FUN = function(x) rep(x, 10))), treatment = 0)
  input.data <- input.data[order(input.data[,'id'], input.data[,'time']), ]
  input.data[input.data$id %in% c(2,4) & input.data$time > 5, 'treatment'] <- -.6
  input.data[input.data$id %in% c(6) & input.data$time > 5, 'treatment'] <- .4
  #input.data[input.data[, 'treatment'] == 0, 'treatment'] <- .035
  
  input.data$cal.data <- input.data$id
  input.data$outcome <- rnorm(nrow(input.data))
  
  continuous.treatment.info <- list(treatment.threshold = .2, type = "numeric", 
                                    method = "max", units = "raw", matching.threshold = 2) #include everything 
  
  PM.results <- PanelMatch(lag = 4, time.id = "time", unit.id = "id", 
                           treatment = "treatment", refinement.method = "none", # should be none for all of them
                           data = input.data, match.missing = TRUE, 
                           size.match = 5, qoi = "ate" , outcome.var = "outcome",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           continuous.treatment.info = continuous.treatment.info)
  
  expect_true(length(PM.results$atc) == 2) 
  expect_true(all(c("2.6","4.6") %in% names(PM.results$atc)))
  
  expect_true(all(PM.results$atc[[1]] == c(1,3,5, 6, 7, 8, 9, 10)))
  expect_true(all(PM.results$atc[[2]] == c(1,3,5,6, 7, 8,9,10))) 
  
  expect_true(length(PM.results$att) == 1) 
  expect_true(all(c("6.6") %in% names(PM.results$att)))
  
  expect_true(all(PM.results$att[[1]] == c(1,2,3,4,5,7,8,9,10)))
  
  
})

