test_that("ps.match keeps a control whose propensity score equals the treated unit's", {
  set.seed(3)
  n.units <- 12
  n.periods <- 8
  d <- data.frame(id = rep(1:n.units, each = n.periods),
                  t = rep(1:n.periods, times = n.units))
  d$d <- 0
  d$d[d$id %in% 1:3 & d$t >= 6] <- 1
  d$x <- rnorm(nrow(d))
  d$y <- rnorm(nrow(d))
  # unit 12 has the covariate and outcome values of treated unit 1, but is never treated
  d[d$id == 12, c("x", "y")] <- d[d$id == 1, c("x", "y")]
  d.panel <- PanelData(d, "id", "t", "d", "y")
  
  for (method in c("ps.match", "mahalanobis"))
  {
    PM.results <- PanelMatch(lag = 4, refinement.method = method,
                             panel.data = d.panel, covs.formula = ~ x,
                             size.match = 3, qoi = "att", lead = 0)
    s <- PM.results$att[["1.6"]]
    w <- attr(s, "weights")
    expect_equal(unname(w[as.numeric(s) == 12]), 1 / 3)
    expect_equal(sum(w > 0), 3)
  }
})
