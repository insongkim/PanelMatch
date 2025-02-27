test_that("balance checking functions are sensible", {
  set.seed(1)
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj, 
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  compmat <- matrix(data = c(0.0601621099966881,-0.00601042376056592,
                             0.00212772294115459,0.104946475743932, 0.133482367739366,
                             0.0252020903167403,0.051861118295163,-0.125475177405918,
                             0.0268787670968836,-0.0117184403355875), 
                    ncol = 2, nrow = 5)
  just.mat <- summary(balmat, qoi = 'att', include.unrefined = F)[[1]]
  expect_equal(nrow(just.mat), 5)
  expect_equal(ncol(just.mat), 2)
  expect_equivalent(just.mat, compmat)
  

  expect_warning(just.mat <- summary(balmat, qoi = 'att', unrefined.only = TRUE)[[1]])
  
  expect_false(isTRUE(all.equal(balmat, compmat, check.attributes = FALSE)))
  
})


test_that("Handles different QOIs, especially ATE", {
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj, 
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 1)
  expect_true(length(balmat[[1]]) == 1)
  expect_true(names(balmat[[1]]) == "att")
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "ate")
  
  balmat <- get_covariate_balance(pm.obj, 
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 1)
  expect_true(length(balmat[[1]]) == 2)
  expect_true(all(c("att", "atc") %in% names(balmat[[1]])))
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "art")
  
  balmat <- get_covariate_balance(pm.obj, 
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 1)
  expect_true(length(balmat[[1]]) == 1)
  expect_true(names(balmat[[1]]) == "art")
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "atc")
  
  balmat <- get_covariate_balance(pm.obj, 
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 1)
  expect_true(length(balmat[[1]]) == 1)
  expect_true(names(balmat[[1]]) == "atc")
})


test_that("Handles multiple PanelMatch configurations", {
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  pm.obj1 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.match", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  
  pm.obj.art <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "art")
  
  
  pm.obj.ate <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                           size.match = 5, qoi = "ate")
  
  # balmat <- get_covariate_balance(pm.obj, pm.obj1, pm.obj.art, pm.obj.ate,
  #                                 use.equal.weights = FALSE,
  #                                 panel.data = dem.panel, 
  #                                 covariates = c("tradewb", "rdata"))
  
  balmat <- get_covariate_balance(pm.obj, pm.obj1, pm.obj.art, pm.obj.ate,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 4)
  expect_true(length(balmat[1]) == 1)
  expect_true(length(balmat[2]) == 1)
  expect_true(length(balmat[3]) == 1)
  expect_true(length(balmat[4]) == 1)
  expect_true(length(balmat[[1]]) == 1)
  expect_true(length(balmat[[2]]) == 1)
  expect_true(length(balmat[[3]]) == 1)
  expect_true(length(balmat[[4]]) == 2)
  
  jk <- balmat[[1]]
  
  expect_true(identical(dim(jk[["att"]]), c(5L, 2L)))
  
  jk <- balmat[[2]]
  
  expect_true(identical(dim(jk[["att"]]), c(5L, 2L)))
  
  
  jk <- balmat[[3]]
  
  expect_true(identical(dim(jk[["art"]]), c(5L, 2L)))
  
  jk <- balmat[[4]]
  
  expect_true(identical(dim(jk[["att"]]), c(5L, 2L)))
  expect_true(identical(dim(jk[["atc"]]), c(5L, 2L)))
})


test_that("summary function tests", {
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  # balmat <- get_covariate_balance(pm.obj,
  #                                 use.equal.weights = FALSE,
  #                                 panel.data = dem.panel, 
  #                                 covariates = c("tradewb", "rdata"))
  
  balmat <- get_covariate_balance(pm.obj,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  
  expect_error(summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = TRUE))
  a2 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_error(a3 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = TRUE))
  expect_error(a4 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = FALSE))
  expect_error(summary(balmat))
  #expect_true(length(a1) == 1)
  expect_true(length(a2) == 1)
  # expect_true(length(a3) == 1)
  # expect_true(length(a4) == 1)
  
  # expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  # expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  # expect_false(identical(a1[[1]], a2[[1]]))
  # expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  # expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  # 
  pm.obj1 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.weight", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj, pm.obj1,
                                  include.unrefined = TRUE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  
  a1 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = TRUE)
  a2 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_warning(a3 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = TRUE))
  a4 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = FALSE)
  #expect_error(summary(balmat))
  
  expect_true(length(a1) == 2)
  expect_true(length(a2) == 2)
  expect_true(length(a3) == 2)
  expect_true(length(a4) == 2)
  
  expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  expect_false(identical(a1[[1]], a2[[1]]))
  expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  
  expect_true(all(dim(a1[[2]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[2]]) == c(5L, 2L)))
  expect_false(identical(a1[[2]], a2[[2]]))
  expect_true(all(dim(a3[[2]]) == c(5L, 2L)))
  expect_true(all(dim(a4[[2]]) == c(5L, 4L)))
  # results should not be the same
  expect_false(identical(a4[[1]], a4[[2]]))
  
  
  ### do again for ATE
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "ate")
  
  balmat <- get_covariate_balance(pm.obj,
                                  #use.equal.weights = FALSE,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  
  expect_error(a1 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = TRUE))
  a2 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_error(a3 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = TRUE))
  expect_error(a4 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = FALSE))
  expect_error(summary(balmat))
  #expect_true(length(a1) == 1)
  expect_true(length(a2) == 1)
  #expect_true(length(a3) == 1)
  #expect_true(length(a4) == 1)
  
  #expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  #expect_false(identical(a1[[1]], a2[[1]]))
  #expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  #expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  
  expect_error(a1 <- summary(balmat, qoi = "atc", include.unrefined = FALSE, unrefined.only = TRUE))
  a2 <- summary(balmat, qoi = "atc", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_error(a3 <- summary(balmat, qoi = "atc", include.unrefined = TRUE, unrefined.only = TRUE))
  expect_error(a4 <- summary(balmat, qoi = "atc", include.unrefined = TRUE, unrefined.only = FALSE))
  expect_error(summary(balmat))
  #expect_true(length(a1) == 1)
  expect_true(length(a2) == 1)
  #expect_true(length(a3) == 1)
  #expect_true(length(a4) == 1)
  
  #expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  #expect_false(identical(a1[[1]], a2[[1]]))
  #expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  #expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  
  
  pm.obj1 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.weight", 
                        panel.data = dem.panel, match.missing = TRUE,
                        covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                        size.match = 5, qoi = "ate")
  
  balmat <- get_covariate_balance(pm.obj, pm.obj1,
                                  #use.equal.weights = FALSE,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  
  expect_error(a1 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = TRUE))
  a2 <- summary(balmat, qoi = "att", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_error(a3 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = TRUE))
  expect_error(a4 <- summary(balmat, qoi = "att", include.unrefined = TRUE, unrefined.only = FALSE))
  expect_error(summary(balmat))
  
  #expect_true(length(a1) == 2)
  expect_true(length(a2) == 2)
  #expect_true(length(a3) == 2)
  #expect_true(length(a4) == 2)
  
  #expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  # expect_false(identical(a1[[1]], a2[[1]]))
  # expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  # expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  
  # expect_true(all(dim(a1[[2]]) == c(5L, 2L)))
  # expect_true(all(dim(a2[[2]]) == c(5L, 2L)))
  # expect_false(identical(a1[[2]], a2[[2]]))
  # expect_true(all(dim(a3[[2]]) == c(5L, 2L)))
  # expect_true(all(dim(a4[[2]]) == c(5L, 4L)))
  # # results should not be the same
  # expect_false(identical(a4[[1]], a4[[2]]))
  
  
  expect_error(a1 <- summary(balmat, qoi = "atc", include.unrefined = FALSE, unrefined.only = TRUE))
  a2 <- summary(balmat, qoi = "atc", include.unrefined = FALSE, unrefined.only = FALSE)
  expect_error(a3 <- summary(balmat, qoi = "atc", include.unrefined = TRUE, unrefined.only = TRUE))
  expect_error(a4 <- summary(balmat, qoi = "atc", include.unrefined = TRUE, unrefined.only = FALSE))
  expect_error(summary(balmat))
  #expect_true(length(a1) == 2)
  expect_true(length(a2) == 2)
  #expect_true(length(a3) == 2)
  #expect_true(length(a4) == 2)
  
  #expect_true(all(dim(a1[[1]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[1]]) == c(5L, 2L)))
  #expect_false(identical(a1[[1]], a2[[1]]))
  #expect_true(all(dim(a3[[1]]) == c(5L, 2L)))
  #expect_true(all(dim(a4[[1]]) == c(5L, 4L)))
  
  #expect_true(all(dim(a1[[2]]) == c(5L, 2L)))
  expect_true(all(dim(a2[[2]]) == c(5L, 2L)))
  #expect_false(identical(a1[[2]], a2[[2]]))
  #expect_true(all(dim(a3[[2]]) == c(5L, 2L)))
  #expect_true(all(dim(a4[[2]]) == c(5L, 4L)))
  # results should not be the same
  #expect_false(identical(a4[[1]], a4[[2]]))
  
})

test_that("unrefined parameters", {
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj,
                                  include.unrefined = TRUE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_false(is.null(attr(balmat, "unrefined.balance.results")))
  balmat <- get_covariate_balance(pm.obj,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  
  expect_true(is.null(attr(balmat, "unrefined.balance.results")))
  # expect_identical(dim(balmat[[1]][["att"]]), 
  #                  dim(attr(balmat, "unrefined.balance.results")[[1]][["att"]]))
})
  



test_that("print.PanelBalance", {
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  balmat <- get_covariate_balance(pm.obj,
                                  #use.equal.weights = TRUE,
                                  include.unrefined = TRUE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_output(print(balmat))
})



test_that("plot.PanelBalance", {

  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  pb <- get_covariate_balance(pm.obj,
                              include.unrefined = TRUE,
                              panel.data = dem.panel, 
                              covariates = c("tradewb", "rdata"))
  plot(pb, type = "panel")
  plot(pb, type = "panel", include.unrefined.panel = FALSE)
  plot(pb, type = "scatter")
  
    
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  
  pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata, 
                       size.match = 5, qoi = "att")
  
  pb <- get_covariate_balance(pm.obj, pm2,
                                  #use.equal.weights = TRUE,
                                  include.unrefined = TRUE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  plot(pb, type = "panel")
  plot(pb, type = "scatter")
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "art")
  
  
  pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                    panel.data = dem.panel, match.missing = TRUE,
                    covs.formula = ~ tradewb + rdata, 
                    size.match = 5, qoi = "art")
  
  pb <- get_covariate_balance(pm.obj, pm2,
                              #use.equal.weights = TRUE,
                              include.unrefined = TRUE,
                              panel.data = dem.panel, 
                              covariates = c("tradewb", "rdata"))
  plot(pb, type = "panel")
  plot(pb, type = "scatter")
  
  
  
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "ate")
  
  
  pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                    panel.data = dem.panel, match.missing = TRUE,
                    covs.formula = ~ tradewb + rdata, 
                    size.match = 5, qoi = "ate")
  
  pb <- get_covariate_balance(pm.obj, pm2,
                              #use.equal.weights = TRUE,
                              include.unrefined = TRUE,
                              panel.data = dem.panel, 
                              covariates = c("tradewb", "rdata"))
  plot(pb, type = "panel")
  plot(pb, type = "scatter")
  
})


test_that("Subsetting works", {
  dem$rdata <- runif(runif(nrow(dem)))
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                       panel.data = dem.panel, match.missing = TRUE,
                       covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                       size.match = 5, qoi = "att")
  
  pm.obj1 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.match", 
                        panel.data = dem.panel, match.missing = TRUE,
                        covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                        size.match = 5, qoi = "att")
  
  
  pm.obj.art <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                           size.match = 5, qoi = "art")
  
  
  pm.obj.ate <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
                           panel.data = dem.panel, match.missing = TRUE,
                           covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
                           size.match = 5, qoi = "ate")
  
  balmat <- get_covariate_balance(pm.obj, pm.obj1, pm.obj.art, pm.obj.ate,
                                  include.unrefined = FALSE,
                                  panel.data = dem.panel, 
                                  covariates = c("tradewb", "rdata"))
  expect_true(length(balmat) == 4)
  expect_true(length(balmat[1]) == 1)
  expect_true(length(balmat[4]) == 1)
  expect_true(all(class(balmat[1]) == c("PanelBalance", "list")))
  expect_true(class(balmat[[1]]) == c("list"))
  expect_true(1 == length(balmat[[1]]))
  expect_true(all(dim(balmat[[1]][[1]]) == c(5, 2)))
  expect_true(all(class(balmat[[1]][[1]]) == c("matrix", "array")))
  
  expect_true(all(class(balmat[4]) == c("PanelBalance", "list")))
  expect_true(class(balmat[[4]]) == c("list"))
  expect_true(2 == length(balmat[[4]]))
  expect_true(all(dim(balmat[[4]][[1]]) == c(5, 2)))
  expect_true(all(dim(balmat[[4]][[2]]) == c(5, 2)))
  expect_true(all(class(balmat[[4]][[1]]) == c("matrix", "array")))
  expect_true(all(class(balmat[[4]][[2]]) == c("matrix", "array")))
  
  summary(balmat[1], include.unrefined = FALSE)
  summary(balmat[4], qoi = "att", include.unrefined = FALSE)
  summary(balmat[4], qoi = "atc", include.unrefined = FALSE)
  
})