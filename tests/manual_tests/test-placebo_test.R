# This is the old version of the placebo test tests before changing the mahalanobis distance calculation changes
# test_that("test placebo test", {
#   dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#   PM.results <- PanelMatch(lag = 4, 
#                            refinement.method = "mahalanobis",
#                            panel.data = dem.panel,
#                            match.missing = TRUE,
#                            covs.formula = ~ I(lag(tradewb, 1:4)),
#                            size.match = 5, qoi = "art",
#                            lead = 0:4, 
#                            forbid.treatment.reversal = FALSE,
#                            placebo.test = TRUE)
#   set.seed(1)
#   pt.res <- placebo_test(pm.obj = PM.results, 
#                          panel.data = dem.panel, 
#                          number.iterations = 100, 
#                          se.method = "bootstrap",
#                          plot = FALSE)
#   
#   est.comps <- c(-3.740490,
#                  -2.136306,
#                  -1.702141)
#   
#   st.comps <- c(2.1743547, 1.5384570, 0.7407867)
#   
#   expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
#   expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
#   
#   
#   PM.results <- PanelMatch(lag = 4, 
#                            refinement.method = "mahalanobis",
#                            panel.data = dem.panel,
#                            match.missing = TRUE,
#                            covs.formula = ~ I(lag(tradewb, 1:4)),
#                            size.match = 5, qoi = "att",
#                            lead = 0:4, 
#                            forbid.treatment.reversal = FALSE,
#                            placebo.test = TRUE)
#   set.seed(1)
#   pt.res <- placebo_test(PM.results, 
#                          panel.data = dem.panel, 
#                          number.iterations = 100, 
#                          se.method = "bootstrap",
#                          plot = FALSE)
#   
#   est.comps <- c(-7.378622, -5.835059, -2.673118)
#   
#   st.comps <- c(2.0099124, 1.3899878, 0.9683683)
#   
#   expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
#   expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
#   
#   
#   pt.res <- placebo_test(PM.results, 
#                          panel.data = dem.panel, 
#                          se.method = "conditional",
#                          plot = FALSE)
#   
#   est.comps <- c(-7.378622, -5.835059, -2.673118)
#   
#   st.comps <- c(1.6581540, 1.2276371, 0.8310945)
#   
#   expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
#   expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
#   
#   pt.res <- placebo_test(PM.results, 
#                          panel.data = dem.panel, 
#                          se.method = "unconditional",
#                          plot = FALSE)
#   
#   est.comps <- c(-7.378622, -5.835059, -2.673118)
#   st.comps <- c(2.187916, 1.632308, 1.066801)
#   
#   expect_equivalent(pt.res$estimate, est.comps, tolerance = .00001)
#   expect_equivalent(pt.res$standard.errors, st.comps, tolerance = .00001)
#   
# })

#These are the new version of the tests
test_that("placebo_test returns stable estimates and SEs", {
  dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
  
  ## -------------------------
  ## ART, bootstrap SEs
  ## -------------------------
  PM.results <- PanelMatch(
    lag                       = 4,
    refinement.method         = "mahalanobis",
    panel.data                = dem.panel,
    match.missing             = TRUE,
    covs.formula              = ~ I(lag(tradewb, 1:4)),
    size.match                = 5,
    qoi                       = "art",
    lead                      = 0:4,
    forbid.treatment.reversal = FALSE,
    placebo.test              = TRUE
  )
  
  set.seed(1)
  pt.res <- placebo_test(
    pm.obj            = PM.results,
    panel.data        = dem.panel,
    number.iterations = 100,
    se.method         = "bootstrap",
    plot              = FALSE
  )
  
  est_art_boot <- c(-3.559363, -2.399573, -1.802547)
  se_art_boot  <- c(2.1572803, 1.5549341, 0.7578422)
  names(est_art_boot) <- paste0("t-", 4:2)
  names(se_art_boot)  <- paste0("t-", 4:2)
  
  expect_equal(pt.res$estimate,         est_art_boot, tolerance = 1e-6)
  expect_equal(pt.res$standard.errors,  se_art_boot,  tolerance = 1e-6)
  
  ## -------------------------
  ## ATT, bootstrap SEs
  ## -------------------------
  PM.results <- PanelMatch(
    lag                       = 4,
    refinement.method         = "mahalanobis",
    panel.data                = dem.panel,
    match.missing             = TRUE,
    covs.formula              = ~ I(lag(tradewb, 1:4)),
    size.match                = 5,
    qoi                       = "att",
    lead                      = 0:4,
    forbid.treatment.reversal = FALSE,
    placebo.test              = TRUE
  )
  
  set.seed(1)
  pt.res <- placebo_test(
    pm.obj            = PM.results,
    panel.data        = dem.panel,
    number.iterations = 100,
    se.method         = "bootstrap",
    plot              = FALSE
  )
  
  est_att_boot <- c(-7.295516, -5.458616, -2.363859)
  se_att_boot  <- c(1.8858259, 1.3524194, 0.9659462)
  names(est_att_boot) <- paste0("t-", 4:2)
  names(se_att_boot)  <- paste0("t-", 4:2)
  
  expect_equal(pt.res$estimate,        est_att_boot, tolerance = 1e-6)
  expect_equal(pt.res$standard.errors, se_att_boot,  tolerance = 1e-6)
  
  ## -------------------------
  ## ATT, conditional SEs
  ## -------------------------
  pt.res <- placebo_test(
    pm.obj     = PM.results,
    panel.data = dem.panel,
    se.method  = "conditional",
    plot       = FALSE
  )
  
  est_att_cond <- c(-7.295516, -5.458616, -2.363859)
  se_att_cond  <- c(1.5778832, 1.1573531, 0.7964844)
  names(est_att_cond) <- paste0("t-", 4:2)
  names(se_att_cond)  <- paste0("t-", 4:2)
  
  expect_equal(pt.res$estimate,        est_att_cond, tolerance = 1e-6)
  expect_equal(pt.res$standard.errors, se_att_cond,  tolerance = 1e-6)
  
  ## -------------------------
  ## ATT, unconditional SEs
  ## -------------------------
  pt.res <- placebo_test(
    pm.obj     = PM.results,
    panel.data = dem.panel,
    se.method  = "unconditional",
    plot       = FALSE
  )
  
  est_att_uncond <- c(-7.295516, -5.458616, -2.363859)
  se_att_uncond  <- c(2.139582, 1.572801, 1.042982)
  names(est_att_uncond) <- paste0("t-", 4:2)
  names(se_att_uncond)  <- paste0("t-", 4:2)
  
  expect_equal(pt.res$estimate,        est_att_uncond, tolerance = 1e-6)
  expect_equal(pt.res$standard.errors, se_att_uncond,  tolerance = 1e-6)
})