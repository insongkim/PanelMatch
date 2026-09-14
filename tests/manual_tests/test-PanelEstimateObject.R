test_that("print.PanelEstimate tests", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)

  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "conditional")
  expect_output(print(pe.results, regexp = "Point estimates:"))

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)

  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  expect_output(print(pe.results, regexp = "Point estimates:"))

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)

  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  expect_output(print(pe.results, regexp = "Point estimates:"))


  ##### trying ate
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)

  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  expect_output(print(pe.results, regexp = "Point estimates:"))

})

test_that("summary.PanelEstimate (object tests)", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  trt <- summary(pe.results)
 
  # check lower bounds and that one can specify a new confidence level
  expect_equal(trt[,3], c(-2.457718, -2.920344, -2.851568, -2.602203),
               tolerance = .000001)

  expect_equal(summary(pe.results, confidence.level = .9)[,3],
               c(-2.207862, -2.501870, -2.300933, -1.926554),
               tolerance = .000001)

  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  trt <- summary(pe.results)
  expect_true(all(dim(trt) == c(4, 4)))

})

# no good way to test plotting results directly.
test_that("plot.PanelEstimate tests", {
  dem.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')

  qoi_ <- "att"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "unconditional")
  plot(pe.results)

  plot(pe.results, confidence.level = .9)
  plot(pe.results, confidence.level = .99)

  qoi_ <- "ate"
  pm1 <- PanelMatch(lag = 4,
                    refinement.method = "mahalanobis",
                    panel.data = dem.panel,
                    match.missing = FALSE, covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
                    size.match = 5, qoi = qoi_,
                    lead = 0:3, forbid.treatment.reversal = FALSE)
  pe.results <- PanelEstimate(sets = pm1, panel.data = dem.panel, se.method = "bootstrap")
  plot(pe.results)

})


make_test_panel_match <- function()
{
  msets <- matched_set(
    matchedsets = list(c(2, 3, 4), c(1, 3)),
    id = c(1, 2),
    t = c(5, 6),
    L = 2,
    t.var = "time",
    id.var = "unit",
    treatment.var = "treatment"
  )
  
  attr(msets, "refinement.method") <- "ps.weight"
  attr(msets, "covs.formula") <- ~ x1 + x2
  attr(msets, "match.missing") <- TRUE
  attr(msets, "max.match.size") <- 5
  
  attr(msets[[1]], "weights") <- c(0.2, 0.3, 0.5)
  attr(msets[[1]], "distances") <- c(0.4, 0.2, 0.1)
  attr(msets[[1]], "treatment.change") <- 1
  attr(msets[[1]], "control.change") <- c(0, 0, 0)
  
  attr(msets[[2]], "weights") <- c(0.6, 0.4)
  attr(msets[[2]], "distances") <- c(0.3, 0.5)
  attr(msets[[2]], "treatment.change") <- 1
  attr(msets[[2]], "control.change") <- c(0, 0)
  
  res <- list(att = msets)
  class(res) <- "PanelMatch"
  attr(res, "qoi") <- "att"
  attr(res, "outcome.var") <- "y"
  attr(res, "lead") <- 0:2
  attr(res, "forbid.treatment.reversal") <- FALSE
  attr(res, "placebo.test") <- FALSE
  
  res
}


test_that("equal PanelMatch objects compare equal", {
  x <- make_test_panel_match()
  y <- x
  
  expect_true(isTRUE(all.equal(x, y)))
})


test_that("differences in matched control units are detected", {
  x <- make_test_panel_match()
  y <- x
  
  y$att[[1]][2] <- 10
  
  expect_false(isTRUE(all.equal(x, y)))
})


test_that("differences in control-unit weights are detected", {
  x <- make_test_panel_match()
  y <- x
  
  weights <- attr(y$att[[1]], "weights")
  weights[2] <- 0.4
  attr(y$att[[1]], "weights") <- weights
  
  expect_false(isTRUE(all.equal(x, y)))
})


test_that("differences in matched.set metadata are detected", {
  x <- make_test_panel_match()
  y <- x
  
  attr(y$att, "lag") <- 3
  
  expect_false(isTRUE(all.equal(x, y)))
})


test_that("differences in PanelMatch metadata are detected", {
  x <- make_test_panel_match()
  y <- x
  
  attr(y, "lead") <- 0:3
  
  expect_false(isTRUE(all.equal(x, y)))
})


test_that("comparison is sensitive to matched-control ordering", {
  x <- make_test_panel_match()
  y <- x
  
  controls <- y$att[[1]]
  attrs <- attributes(controls)
  
  controls <- controls[c(3, 2, 1)]
  attr(controls, "weights") <- attrs$weights[c(3, 2, 1)]
  attr(controls, "distances") <- attrs$distances[c(3, 2, 1)]
  attr(controls, "treatment.change") <- attrs$treatment.change
  attr(controls, "control.change") <- attrs$control.change[c(3, 2, 1)]
  y$att[[1]] <- controls
  
  expect_false(isTRUE(all.equal(x, y)))
})


test_that("additional all.equal arguments are respected", {
  x <- make_test_panel_match()
  y <- x
  
  weights <- attr(y$att[[1]], "weights")
  weights[1] <- weights[1] + 1e-8
  attr(y$att[[1]], "weights") <- weights
  
  expect_true(isTRUE(all.equal(x, y, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(x, y, tolerance = 1e-10)))
})


test_that("PanelMatch objects are not equal to objects of another class", {
  x <- make_test_panel_match()
  y <- unclass(x)
  
  comparison <- all.equal(x, y)
  
  expect_false(isTRUE(comparison))
  expect_match(comparison, "not a PanelMatch object", fixed = TRUE)
})