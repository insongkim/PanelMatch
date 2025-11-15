# test_that("summary.PanelEstimate (conditional)", {
#   dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#   qoi_ <- "att"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "conditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.738635065855285,1.21038195434255,1.55923210008201,1.82486418813223,-2.04109789825896,-2.69356525042577,-2.59972747285187,-1.84484646286253,0.854298355330497,2.05104482567144,3.51235004654712,5.30848970796966)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "att")
#   
#   qoi_ <- "atc"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "conditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.591756652995951,1.11369402788093,1.46937851945315,1.85716210624063,-1.89981045123073,-2.32467787575185,-3.37138841788079,-3.78228584150992,0.419833003737329,2.04092249313608,2.38846953768912,3.49765584185843)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "atc")
#   
#   
#   qoi_ <- "art"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "conditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.02600037922772,1.48310644776257,1.91955116019568,2.15026718759732,-7.22864265630046,-10.9282208645128,-12.5187199096139,-12.3384409599025,-3.20679507347894,-5.11455041880524,-4.99421762868267,-3.90954847024455)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "art")
# })
# 
# 
# test_that("summary.PanelEstimate (unconditional)", {
#   dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#   qoi_ <- "att"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "unconditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.905080467048022,1.48070275292364,1.90752332999506,2.23758074403404,-2.36732488998905,-3.22338427991681,-3.2823657396126,-2.65375604825349,1.18052534706058,2.58086385516248,4.19498831330785,6.11739929336062)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "att")
#   
#   qoi_ <- "atc"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "unconditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.708688258862262,1.33011893273585,1.75520486661622,2.21802563150778,-2.12899218738313,-2.74886289462501,-3.93159776415305,-4.4895653543677,0.649014739889732,2.46510751200924,2.94867888396138,4.20493535471621)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "atc")
#   
#   
#   qoi_ <- "art"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   
#   pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
#                               se.method = "unconditional")
#   expect_no_error(summary(pe.results)) 
#   # expect_true(all(length(summary(pe.results)) == 3))
#   comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423,-8.29736530389636,-12.5517704314982,-14.369074947577,-14.1843445061431,-2.13807242588304,-3.49100085181983,-3.14386259071958,-2.06364492400392)
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL #just compare the numbers
#   cmat <- unlist(cmat)
#   expect_equal(cmat, comp.vec, tolerance = .000001)
#   expect_identical(pe.results$qoi, "art")
#   
# })
# 
# 
# test_that("summary.PanelEstimate (bootstrap)", {
#   dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#   qoi_ <- "att"
#   pm1 <- PanelMatch(lag = 4, 
#                     refinement.method = "mahalanobis",
#                     panel.data = dem.panel, 
#                     match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#                     size.match = 5, qoi = qoi_,
#                     lead = 0:3, forbid.treatment.reversal = FALSE)
#   set.seed(1)
#   
#   pe.results <- PanelEstimate(sets = pm1, 
#                               panel.data = dem.panel,
#                               se.method = "bootstrap")
#   
#   expect_no_error(summary(pe.results)) 
#   cmat <- summary(pe.results, verbose = FALSE)
#   attributes(cmat) <- NULL
#   cmat <- unlist(cmat)
#   comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.906705857062448,1.46876653085356,1.87876179592695,2.21608840568008,-2.47319633748927,-3.25290871525436,-3.34863750081872,-2.68949846129008,1.12752178115124,2.43516177434591,4.04890571881315,5.89171531990835)
#   expect_equal(cmat, comp.vec, tolerance = .0000001)
#   
# })
test_that("summary.PanelEstimate (conditional)", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  ## ATT, conditional SEs
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "conditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_att_cond <- rbind(
    "t+0" = c(estimate = -0.9036306, std.error = 0.6394206, "2.5%" = -2.156872,  "97.5%" = 0.3496108),
    "t+1" = c(estimate = -0.3174628, std.error = 1.0766697, "2.5%" = -2.427697,  "97.5%" = 1.7927710),
    "t+2" = c(estimate =  0.5733430, std.error = 1.4164110, "2.5%" = -2.202771,  "97.5%" = 3.3494575),
    "t+3" = c(estimate =  1.6002885, std.error = 1.7346946, "2.5%" = -1.799650,  "97.5%" = 5.0002274)
  )
  
  expect_equal(cmat, expected_att_cond, tolerance = 1e-6)
  
  ## ATC, conditional SEs
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "conditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_atc_cond <- rbind(
    "t+0" = c(estimate = -0.7690320, std.error = 0.5929128, "2.5%" = -1.931120,  "97.5%" = 0.3930557),
    "t+1" = c(estimate = -0.2760459, std.error = 1.1390635, "2.5%" = -2.508569,  "97.5%" = 1.9564775),
    "t+2" = c(estimate = -0.6871066, std.error = 1.5049599, "2.5%" = -3.636774,  "97.5%" = 2.2625606),
    "t+3" = c(estimate = -0.3738689, std.error = 1.8968186, "2.5%" = -4.091565,  "97.5%" = 3.3438272)
  )
  
  expect_equal(cmat, expected_atc_cond, tolerance = 1e-6)
  
  ## ART, conditional SEs
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "conditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_art_cond <- rbind(
    "t+0" = c(estimate = -5.171761, std.error = 1.022766, "2.5%" = -7.176345,  "97.5%" = -3.167178),
    "t+1" = c(estimate = -7.860783, std.error = 1.484803, "2.5%" = -10.770945, "97.5%" = -4.950622),
    "t+2" = c(estimate = -8.448983, std.error = 1.915196, "2.5%" = -12.202698, "97.5%" = -4.695268),
    "t+3" = c(estimate = -7.820888, std.error = 2.147624, "2.5%" = -12.030154, "97.5%" = -3.611622)
  )
  
  expect_equal(cmat, expected_art_cond, tolerance = 1e-6)
})

test_that("summary.PanelEstimate (unconditional)", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  ## ATT, unconditional SEs
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "unconditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_att_uncond <- rbind(
    "t+0" = c(estimate = -0.9036306, std.error = 0.7929164, "2.5%" = -2.457718, "97.5%" = 0.6504569),
    "t+1" = c(estimate = -0.3174628, std.error = 1.3280250, "2.5%" = -2.920344, "97.5%" = 2.2854185),
    "t+2" = c(estimate =  0.5733430, std.error = 1.7474357, "2.5%" = -2.851568, "97.5%" = 3.9982540),
    "t+3" = c(estimate =  1.6002885, std.error = 2.1441677, "2.5%" = -2.602203, "97.5%" = 5.8027799)
  )
  
  expect_equal(cmat, expected_att_uncond, tolerance = 1e-6)
  
  ## ATC, unconditional SEs
  qoi_ <- "atc"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "unconditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_atc_uncond <- rbind(
    "t+0" = c(estimate = -0.7690320, std.error = 0.7102239, "2.5%" = -2.161045, "97.5%" = 0.6229814),
    "t+1" = c(estimate = -0.2760459, std.error = 1.3605176, "2.5%" = -2.942611, "97.5%" = 2.3905196),
    "t+2" = c(estimate = -0.6871066, std.error = 1.7980118, "2.5%" = -4.211145, "97.5%" = 2.8369319),
    "t+3" = c(estimate = -0.3738689, std.error = 2.2655081, "2.5%" = -4.814183, "97.5%" = 4.0664453)
  )
  
  expect_equal(cmat, expected_atc_uncond, tolerance = 1e-6)
  
  ## ART, unconditional SEs
  qoi_ <- "art"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "unconditional"
  )
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_art_uncond <- rbind(
    "t+0" = c(estimate = -5.171761, std.error = 1.574844, "2.5%" = -8.258399,  "97.5%" = -2.085124),
    "t+1" = c(estimate = -7.860783, std.error = 2.315461, "2.5%" = -12.399004, "97.5%" = -3.322563),
    "t+2" = c(estimate = -8.448983, std.error = 2.857059, "2.5%" = -14.048716, "97.5%" = -2.849250),
    "t+3" = c(estimate = -7.820888, std.error = 3.094246, "2.5%" = -13.885500, "97.5%" = -1.756277)
  )
  
  expect_equal(cmat, expected_art_uncond, tolerance = 1e-6)
})

test_that("summary.PanelEstimate (bootstrap)", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  qoi_ <- "att"
  pm1 <- PanelMatch(
    lag               = 4,
    refinement.method = "mahalanobis",
    panel.data        = dem.panel,
    match.missing     = FALSE,
    covs.formula      = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
    size.match        = 5,
    qoi               = qoi_,
    lead              = 0:3,
    forbid.treatment.reversal = FALSE
  )
  
  set.seed(1)
  pe.results <- PanelEstimate(
    sets       = pm1,
    panel.data = dem.panel,
    se.method  = "bootstrap"
  )
  
  expect_no_error(summary(pe.results))
  
  cmat <- as.matrix(summary(pe.results, verbose = FALSE))
  
  expected_att_boot <- rbind(
    "t+0" = c(estimate = -0.9036306, std.error = 0.7895802, "2.5%" = -2.513913, "97.5%" = 0.575451),
    "t+1" = c(estimate = -0.3174628, std.error = 1.3293180, "2.5%" = -3.003532, "97.5%" = 2.107150),
    "t+2" = c(estimate =  0.5733430, std.error = 1.7296406, "2.5%" = -3.068349, "97.5%" = 3.720226),
    "t+3" = c(estimate =  1.6002885, std.error = 2.1206464, "2.5%" = -2.712990, "97.5%" = 5.825178)
  )
  
  expect_equal(cmat, expected_att_boot, tolerance = 1e-6)
})