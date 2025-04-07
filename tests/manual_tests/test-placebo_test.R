test_that("test placebo test", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  PM.results <- PanelMatch(lag = 4, 
                           refinement.method = "mahalanobis",
                           panel.data = dem.panel,
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "art",
                           lead = 0:4, 
                           forbid.treatment.reversal = FALSE,
                           placebo.test = TRUE)
  set.seed(1)
  pt.res <- placebo_test(pm.obj = PM.results, 
                         panel.data = dem.panel, 
                         number.iterations = 100, 
                         se.method = "bootstrap",
                         plot = FALSE)
  
  est.comps <- c(-3.740490,
                 -2.136306,
                 -1.702141)
  
  st.comps <- c(2.1743547, 1.5384570, 0.7407867)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  
  PM.results <- PanelMatch(lag = 4, 
                           refinement.method = "mahalanobis",
                           panel.data = dem.panel,
                           match.missing = TRUE,
                           covs.formula = ~ I(lag(tradewb, 1:4)),
                           size.match = 5, qoi = "att",
                           lead = 0:4, 
                           forbid.treatment.reversal = FALSE,
                           placebo.test = TRUE)
  set.seed(1)
  pt.res <- placebo_test(PM.results, 
                         panel.data = dem.panel, 
                         number.iterations = 100, 
                         se.method = "bootstrap",
                         plot = FALSE)
  
  est.comps <- c(-7.378622, -5.835059, -2.673118)
  
  st.comps <- c(2.0099124, 1.3899878, 0.9683683)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  
  pt.res <- placebo_test(PM.results, 
                         panel.data = dem.panel, 
                         se.method = "conditional",
                         plot = FALSE)
  
  est.comps <- c(-7.378622, -5.835059, -2.673118)
  
  st.comps <- c(1.6581540, 1.2276371, 0.8310945)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  pt.res <- placebo_test(PM.results, 
                         panel.data = dem.panel, 
                         se.method = "unconditional",
                         plot = FALSE)
  
  est.comps <- c(-7.378622, -5.835059, -2.673118)
  st.comps <- c(2.187916, 1.632308, 1.066801)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
})