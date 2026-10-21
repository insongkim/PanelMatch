test_that("ATC ps.weight weights are the odds of being untreated at t", {
  set.seed(4)
  n.units <- 60
  n.periods <- 8
  d <- data.frame(id = rep(1:n.units, each = n.periods),
                  t = rep(1:n.periods, times = n.units))
  d$x <- rnorm(nrow(d))
  # units with a larger x at t = 6 are more likely to become treated at t = 6
  x6 <- d$x[d$t == 6]
  switchers <- which(runif(n.units) < plogis(-0.5 + 1.5 * x6))
  d$d <- as.numeric(d$id %in% switchers & d$t >= 6)
  d$y <- rnorm(nrow(d))
  d.panel <- PanelData(d, "id", "t", "d", "y")
  
  PM.results <- PanelMatch(lag = 4, refinement.method = "ps.weight",
                           panel.data = d.panel, covs.formula = ~ x,
                           qoi = "atc", lead = 0)
  sets <- PM.results$atc[sapply(PM.results$atc, length) > 0]
  ts <- as.numeric(sub(".*\\.", "", names(sets)))
  ids <- as.numeric(sub("\\..*", "", names(sets)))
  
  # propensity score model fitted to the matched observations and their controls at t
  rows <- unique(unlist(lapply(seq_along(sets), function(k) {
    which(d$id %in% c(as.numeric(sets[[k]]), ids[k]) & d$t == ts[k])
  })))
  fit <- glm(d ~ x, family = binomial(link = "logit"), data = d[rows, ])
  for (k in seq_along(sets))
  {
    controls <- as.numeric(sets[[k]])
    control.rows <- match(paste(controls, ts[k]), paste(d$id, d$t))
    e <- predict(fit, newdata = d[control.rows, ], type = "response")
    odds.untreated <- (1 - e) / e
    expect_equal(unname(attr(sets[[k]], "weights")),
                 unname(odds.untreated / sum(odds.untreated)), tolerance = 1e-6)
  }
  
  # the ATC half of the ATE uses the same weights
  PM.ate <- PanelMatch(lag = 4, refinement.method = "ps.weight",
                       panel.data = d.panel, covs.formula = ~ x,
                       qoi = "ate", lead = 0)
  expect_equal(lapply(PM.ate$atc[names(sets)], attr, "weights"),
               lapply(sets, attr, "weights"))
})
