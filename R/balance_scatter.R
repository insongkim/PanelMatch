#' balance_scatter
#'
#' Visualizing the standardized mean differences for covariates via a scatter plot.
#'
#' \code{balance_scatter} visualizes the standardized mean differences for each covariate.
#' Although users can use the scatter plot in a variety of ways, it is recommended that
#' the x-axis refers to balance for covariates before refinement, and y-axis
#' refers to balance after refinement. Users can utilize parameters powered by \code{plot()}
#' in base R to further customize the figure.
#' @param matched_set_list a list of one or more \code{matched.set} objects
#' @param xlim xlim of the scatter plot. This is the same as the \code{xlim} argument in \code{plot()}
#' @param ylim ylim of the scatter plot. This is the same as the \code{ylim} argument in \code{plot()}
#' @param main title of the scatter plot. This is the same as the \code{main} argument in \code{plot()}
#' @param x.axis.label x axis label
#' @param y.axis.label y axis label
#' @param pchs one or more pch indicators for the symbols on the scatter plot. You should specify a pch symbol for each matched.set you specify in matched_set_list. See \code{plot()} for more information
#' @param covariates variables for which balance is displayed
#' @param data the same time series cross sectional data set used to create the matched sets.
#' @param ... optional arguments to be passed to \code{plot()}
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <amrauh@umich.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' # create subset of data for simplicity
#' # get a matched set without refinement
#' sets0 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "none",
#'                     data = dem.sub, match.missing = FALSE,
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#'
#' # get a matched set with refinement using propensity score matching, setting the
#' # size of matched set to 5
#' sets1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "ps.match",
#'                     data = dem.sub, match.missing = FALSE,
#'                     covs.formula = ~ tradewb,
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#'
#' # get another matched set with refinement using propensity score weighting
#' sets2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "ps.weight",
#'                     data = dem.sub, match.missing = FALSE,
#'                     covs.formula = ~ tradewb,
#'                     size.match = 10, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#'
#'
#' # use the function to produce the scatter plot
#' balance_scatter(matched_set_list = list(sets0$att, sets1$att, sets2$att),
#'               data = dem.sub,
#'               covariates = c("y", "tradewb"))
#' # add legend
#' legend(x = 0, y = 0.8,
#' legend = c("mahalanobis",
#'            "PS weighting"),
#' y.intersp = 0.65,
#' x.intersp = 0.3,
#' xjust = 0,
#' pch = c(1, 3), pt.cex = 1,
#' bty = "n", ncol = 1, cex = 1, bg = "white")
#'
#'
#'
#' @export
balance_scatter <- function(matched_set_list,
                            xlim = c(0, .8),
                            ylim = c(0, .8),
                            main = "Standardized Mean Difference of Covariates",
                            pchs = c(2,3),
                            covariates, data,
                            x.axis.label = "Before refinement",
                            y.axis.label = "After refinement",
                            ...) {
  if (length(matched_set_list) < 1) stop("Please provide at least one matched.set object")
  # first, get balance for non-refined set
  # use first matched set and generate unrefined balance using get_covariate_balance
  non_refined_balance <- get_covariate_balance(matched.sets = matched_set_list[[1]],
                                               data = data,
                                               covariates = covariates,
                                               use.equal.weights = TRUE)
  
  # second, get balance for refined sets
  refined_balance <- list()
  for (i in 1:length(matched_set_list)) {
    refined_balance[[i]] <-
      get_covariate_balance(matched.sets = matched_set_list[[i]],
                            data = data,
                            covariates = covariates)
  }
  
  # extract values for x-axis from the non-refined sets
  benchmark <- as.vector(non_refined_balance)
  
  compared <- sapply(refined_balance, 
                     function(x) x <- x[1:(nrow(x)),])
  ## create the plot
  graphics::plot(abs(as.numeric(benchmark)),
                 abs(as.numeric(compared[,1])), 
                 pch = 1,
                 xlab = x.axis.label,
                 ylab = y.axis.label,
                 xlim = xlim,
                 ylim = ylim,
                 main = main,
                 font.main = 1, 
                 ...)
  
  if (length(refined_balance) > 1) {
    for (j in 2:length(refined_balance)) {
      graphics::points(abs(as.numeric(benchmark)),
                       abs(as.numeric(compared[,j])),
                       pch = pchs[j])
    }
  }
  
  graphics::abline(h = 0, 
                   lty = "dashed")
  graphics::abline(0, 1, 
                   lty = 2,
                   col = "red")
  
  
}
