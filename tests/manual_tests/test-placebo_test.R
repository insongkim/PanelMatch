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
  
  est.comps <- c(20.438656,
                 4.110167,
                 -21.562320)
  
  st.comps <- c(28.54986,
                23.67567,
                17.30092)
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
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  
  st.comps <- c(18.33094, 
                17.76234,
                13.37736)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  
  pt.res <- placebo_test(PM.results, 
                         panel.data = dem.panel, 
                         se.method = "conditional",
                         plot = FALSE)
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  
  st.comps <- c(16.30824,
                15.11819,
                11.29278)
  
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
  pt.res <- placebo_test(PM.results, 
                         panel.data = dem.panel, 
                         se.method = "unconditional",
                         plot = FALSE)
  
  est.comps <- c(20.490151,
                 8.484286,  
                 8.649806)
  st.comps <- c(19.87698,
                18.35898,
                13.72441)
  expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
  expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
  
})