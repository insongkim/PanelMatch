test_that("summary.PanelEstimate (conditional)", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "conditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.738635065855285,1.21038195434255,1.55923210008201,1.82486418813223,-2.04109789825896,-2.69356525042577,-2.59972747285187,-1.84484646286253,0.854298355330497,2.05104482567144,3.51235004654712,5.30848970796966)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "att")
  
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "conditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.591756652995951,1.11369402788093,1.46937851945315,1.85716210624063,-1.89981045123073,-2.32467787575185,-3.37138841788079,-3.78228584150992,0.419833003737329,2.04092249313608,2.38846953768912,3.49765584185843)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "atc")
  
  
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "conditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.02600037922772,1.48310644776257,1.91955116019568,2.15026718759732,-7.22864265630046,-10.9282208645128,-12.5187199096139,-12.3384409599025,-3.20679507347894,-5.11455041880524,-4.99421762868267,-3.90954847024455)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "art")
})


test_that("summary.PanelEstimate (unconditional)", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "unconditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.905080467048022,1.48070275292364,1.90752332999506,2.23758074403404,-2.36732488998905,-3.22338427991681,-3.2823657396126,-2.65375604825349,1.18052534706058,2.58086385516248,4.19498831330785,6.11739929336062)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "att")
  
  qoi_ <- "atc"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "unconditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-0.7399887237467,-0.141877691307886,-0.491459440095836,-0.142314999825745,0.708688258862262,1.33011893273585,1.75520486661622,2.21802563150778,-2.12899218738313,-2.74886289462501,-3.93159776415305,-4.4895653543677,0.649014739889732,2.46510751200924,2.94867888396138,4.20493535471621)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "atc")
  
  
  qoi_ <- "art"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel,
                              se.method = "unconditional")
  expect_no_error(summary(pe.results)) 
  # expect_true(all(length(summary(pe.results)) == 3))
  comp.vec <- c(-5.2177188648897,-8.02138564165901,-8.75646876914828,-8.12399471507353,1.57127705575129,2.31146328482272,2.86362720065278,3.09207201707423,-8.29736530389636,-12.5517704314982,-14.369074947577,-14.1843445061431,-2.13807242588304,-3.49100085181983,-3.14386259071958,-2.06364492400392)
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL #just compare the numbers
  cmat <- unlist(cmat)
  expect_equal(cmat, comp.vec, tolerance = .000001)
  expect_identical(pe.results$qoi, "art")
  
})


test_that("summary.PanelEstimate (bootstrap)", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4, 
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel, 
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  set.seed(1)
  
  pe.results <- PanelEstimate(sets = pm1, 
                              panel.data = dem.panel,
                              se.method = "bootstrap")
  
  expect_no_error(summary(pe.results)) 
  cmat <- summary(pe.results, verbose = FALSE)
  attributes(cmat) <- NULL
  cmat <- unlist(cmat)
  comp.vec <- c(-0.593399771464233,-0.321260212377162,0.456311286847623,1.73182162255356,0.906705857062448,1.46876653085356,1.87876179592695,2.21608840568008,-2.47319633748927,-3.25290871525436,-3.34863750081872,-2.68949846129008,1.12752178115124,2.43516177434591,4.04890571881315,5.89171531990835)
  expect_equal(cmat, comp.vec, tolerance = .0000001)
  
})