test_that("Checking that non-PanelData objects return errors at various levels", {
  
  d2 <- as.matrix(dem)
  #passing any non-PanelData object directly should fail now
  expect_error(PanelMatch(lag = 4, refinement.method = "none",
             panel.data = d2, match.missing = FALSE,
             size.match = 5, qoi = "att",
             lead = 0:4, forbid.treatment.reversal = FALSE))
  
  #passing any non-PanelData object directly should fail now
  expect_error(PanelMatch(lag = 4,
                          refinement.method = "none",
                          panel.data = dem, match.missing = FALSE,
                          size.match = 5, qoi = "att",
                          lead = 0:4, forbid.treatment.reversal = FALSE))
  
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(panel.data = dem.panel, lag = 4,
                           refinement.method = "ps.match",
                           match.missing = TRUE,
                           covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4,
                           forbid.treatment.reversal = FALSE)
  
  expect_error(PanelEstimate(sets = PM.results, 
                             panel.data = as.matrix(dem), se.method = "bootstrap"))
  
  expect_error(PanelEstimate(sets = PM.results, 
                             panel.data = dem, se.method = "bootstrap"))
  
  expect_silent(PanelEstimate(sets = PM.results, 
                             panel.data = dem.panel, se.method = "bootstrap"))
  
  expect_error(DisplayTreatment(panel.data = dem))
  expect_silent(DisplayTreatment(panel.data = dem.panel))
})

test_that("edge case matching checks", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  
  PM.object <- PanelMatch(lag = 1, 
                          refinement.method = "none",
                          panel.data = dem.panel, match.missing = FALSE,
                          size.match = 5, qoi = "att",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$att), "matched.set")
  
  set.seed(1)
  pe.results <- PanelEstimate(sets = PM.object, 
                              panel.data = dem.panel, 
                              se.method = "bootstrap")
  expect_equivalent(pe.results$estimate, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 0.776390146867648, tolerance = .000001)
  
  pe.results <- PanelEstimate(sets = PM.object, 
                              panel.data = dem.panel, se.method = "conditional")
  expect_equivalent(pe.results$estimate, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "conditional")
  expect_equivalent(pe.results$standard.error, 0.618327233714141, tolerance = .000001)
  
  pe.results <- PanelEstimate(sets = PM.object, 
                              panel.data = dem.panel, se.method = "unconditional")
  expect_equivalent(pe.results$estimate, -1.20423408195041, tolerance = .000001)
  expect_identical(pe.results$se.method, "unconditional")
  expect_equivalent(pe.results$standard.error, 0.745063157005526, tolerance = .000001)
  
  
  #art
  PM.object <- PanelMatch(lag = 1, 
                          refinement.method = "none",
                          panel.data = dem.panel, match.missing = FALSE,
                          size.match = 5, qoi = "art",
                          lead = 0, forbid.treatment.reversal = FALSE)
  expect_equal(class(PM.object), "PanelMatch")
  expect_equal(class(PM.object$art), "matched.set")
  set.seed(1)
  pe.results <- PanelEstimate(sets = PM.object, panel.data = dem.panel, se.method = "bootstrap")
  expect_equivalent(pe.results$estimate, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "bootstrap")
  expect_equivalent(pe.results$standard.error, 1.11077012122371, tolerance = .000001)
  
  pe.results <- PanelEstimate(sets = PM.object, panel.data = dem.panel, se.method = "conditional")
  expect_equivalent(pe.results$estimate, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "conditional")
  expect_equivalent(pe.results$standard.error, 0.90258389202935, tolerance = .000001)
  
  pe.results <- PanelEstimate(sets = PM.object, panel.data = dem.panel, se.method = "unconditional")
  expect_equivalent(pe.results$estimate, -4.22192232677815, tolerance = .000001)
  expect_identical(pe.results$se.method, "unconditional")
  expect_equivalent(pe.results$standard.error, 1.21738691542365, tolerance = .000001)
  
  
})


test_that("ATT: refinement method does not affect unrefined matched sets", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    panel.data = dem.panel, match.missing = FALSE,
                    size.match = 5, qoi = "att",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  pm2 <- PanelMatch(lag = 4, 
                    refinement.method = "ps.weight",
                    panel.data = dem.panel, match.missing = FALSE, 
                    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm3 <- PanelMatch(lag = 4, 
                    refinement.method = "ps.match", 
                    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    panel.data = dem.panel, match.missing = FALSE,
                    size.match = 5, qoi = "att",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm4 <- PanelMatch(lag = 4, 
                    refinement.method = "CBPS.weight",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm5 <- PanelMatch(lag = 4, 
                    refinement.method = "CBPS.match",
                    panel.data = dem.panel, match.missing = FALSE,
                    covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = "att",
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  expect_equal(length(pm1$att), 105)
  expect_equivalent(pm1$att, pm2$att)
  expect_equivalent(pm2$att, pm3$att)
  expect_equivalent(pm3$att, pm4$att)
  expect_equivalent(pm4$att, pm5$att)

  
})


test_that("ATC: refinement method does not affect unrefined matched sets", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  pm2 <- PanelMatch(lag = 4,
                    refinement.method = "ps.weight",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm3 <- PanelMatch(lag = 4, 
                    refinement.method = "ps.match",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm4 <- PanelMatch(lag = 4, 
                    refinement.method = "CBPS.weight",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:4, forbid.treatment.reversal = FALSE)
  
  
  pm5 <- PanelMatch(lag = 4, 
                    refinement.method = "CBPS.match",
                    panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:4, forbid.treatment.reversal = FALSE)

  expect_equal(length(pm1$atc), 3287)
  expect_equivalent(pm1$atc, pm2$atc)
  expect_equivalent(pm2$atc, pm3$atc)
  expect_equivalent(pm3$atc, pm4$atc)
  expect_equivalent(pm4$atc, pm5$atc)
  
  
})

test_that("matched.set methods work", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
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


test_that("set level treatment effects", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(lag = 4,
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "att",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, 
                                           panel.data = dem.panel, lead = 0:1)
  
  pe.results <- PanelEstimate(sets = PM.results, panel.data = dem.panel, 
                              se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimate[1:2], tolerance = .000001)
  
  
  PM.results <- PanelMatch(lag = 4,
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "art",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, 
                                           panel.data = dem.panel, lead = 0:1)
  
  pe.results <- PanelEstimate(sets = PM.results, panel.data = dem.panel, 
                              se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimate[1:2], tolerance = .000001)
  
  PM.results <- PanelMatch(lag = 4,
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "atc",
                           lead = 0:4, forbid.treatment.reversal = FALSE,
                           placebo.test = FALSE)
  set.effects <- get_set_treatment_effects(pm.obj = PM.results, 
                                           panel.data = dem.panel, lead = 0:1)
  pe.results <- PanelEstimate(sets = PM.results, panel.data = dem.panel, 
                              se.method = "conditional")
  expect_equivalent(c(mean(set.effects[[1]], na.rm = TRUE), 
                      mean(set.effects[[2]], na.rm = TRUE)), 
                    pe.results$estimate[1:2], tolerance = .000001)

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

test_that("ensure parallel bootstrap runs without warning, error", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(lag = 4, 
                           refinement.method = "ps.match",
                           panel.data = dem.panel, 
                           match.missing = TRUE, covs.formula = ~ tradewb,
                           size.match = 5, qoi = "att",
                           lead = 0:4, forbid.treatment.reversal = FALSE)

  
  expect_no_error(PanelEstimate(sets = PM.results,
                                number.iterations = 1000,
                                panel.data = dem.panel,
                                se.method = "bootstrap",
                                parallel = TRUE,
                                num.cores = 4))
  
  expect_no_warning(PanelEstimate(sets = PM.results,
                                number.iterations = 1000,
                                panel.data = dem.panel,
                                se.method = "bootstrap",
                                parallel = TRUE,
                                num.cores = 4))
  
  PM.results <- PanelMatch(lag = 4,
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE, covs.formula = ~ tradewb,
                           size.match = 5, qoi = "art",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  
  expect_no_error(PanelEstimate(sets = PM.results,
                                number.iterations = 1000,
                                panel.data = dem.panel,
                                se.method = "bootstrap",
                                parallel = TRUE,
                                num.cores = 4))
  
  expect_no_warning(PanelEstimate(sets = PM.results,
                                  number.iterations = 1000,
                                  panel.data = dem.panel,
                                  se.method = "bootstrap",
                                  parallel = TRUE,
                                  num.cores = 4))
  
  
  PM.results <- PanelMatch(lag = 4,
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE, covs.formula = ~ tradewb,
                           size.match = 5, qoi = "atc",
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  
  expect_no_error(PanelEstimate(sets = PM.results,
                                number.iterations = 1000,
                                panel.data = dem.panel,
                                se.method = "bootstrap",
                                parallel = TRUE,
                                num.cores = 4))
  
  expect_no_warning(PanelEstimate(sets = PM.results,
                                  number.iterations = 1000,
                                  panel.data = dem.panel,
                                  se.method = "bootstrap",
                                  parallel = TRUE,
                                  num.cores = 4))
  
  PM.results <- PanelMatch(lag = 4, 
                           refinement.method = "ps.match",
                           panel.data = dem.panel, match.missing = TRUE, covs.formula = ~ tradewb,
                           size.match = 5, qoi = "ate", 
                           lead = 0:4, forbid.treatment.reversal = FALSE)
  
  expect_no_error(PanelEstimate(sets = PM.results,
                                number.iterations = 1000,
                                panel.data = dem.panel,
                                se.method = "bootstrap",
                                parallel = TRUE,
                                num.cores = 4))
  
  expect_no_warning(PanelEstimate(sets = PM.results,
                                  number.iterations = 1000,
                                  panel.data = dem.panel,
                                  se.method = "bootstrap",
                                  parallel = TRUE,
                                  num.cores = 4))
})


test_that("forbid treatment reversal restrictions are triggered correctly", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  
  expect_no_error(PanelMatch(lag = 4, 
                               refinement.method = "mahalanobis",
                               panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                               size.match = 5, qoi = "att",
                               lead = 0:3, forbid.treatment.reversal = TRUE))
  
  expect_no_error(PanelMatch(lag = 4, 
                             refinement.method = "mahalanobis",
                             panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                             size.match = 5, qoi = "art",
                             lead = 0:3, forbid.treatment.reversal = TRUE))
  
  expect_error(PanelMatch(lag = 4, 
                               refinement.method = "mahalanobis",
                               panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                               size.match = 5, qoi = "atc",
                               lead = 0:3, forbid.treatment.reversal = TRUE))
  
  expect_error(PanelMatch(lag = 4, 
                             refinement.method = "mahalanobis",
                             panel.data = dem.panel, match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                             size.match = 5, qoi = "ate",
                             lead = 0:3, forbid.treatment.reversal = TRUE))
  
})
