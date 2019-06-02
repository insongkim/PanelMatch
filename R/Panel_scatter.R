#' PanelScatter
#' 
#' \code{DisplayTreatment} visualizes the treatment distribution
#' across unit and time in a panel dataset
#'
#' @param non_refined_set a matched.set object produced by setting `refinement.method` to "none" in `PanelMatch` 
#' @param refined_list a list of one or two matched.set objects
#' @param xlim xlim of the scatter plot
#' @param ylim ylim of the scatter plot
#' @param main title of the scatter plot
#' @param pchs one or two pch for the symbols on the scatter plot
#' @param covariates variables for which balance is displayed
#' @param data the dataset
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <adamrauh@mit.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples 
#' \dontrun{
#' sets1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "CBPS.match",
#'                     data = dem, match.missing = F,
#'                     covs.formula = ~ lag("y", 1:4) + lag("tradewb", 1:4),
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 4, forbid.treatment.reversal = FALSE)
#' 
#' sets0 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "none",
#'                     data = dem, match.missing = F,
#'                     covs.formula = ~ lag("y", 1:4) + lag("tradewb", 1:4),
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 4, forbid.treatment.reversal = FALSE)
#' 
#' sets2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "CBPS.weight",
#'                     data = dem, match.missing = F,
#'                     covs.formula = ~ lag("y", 1:4) + lag("tradewb", 1:4),
#'                     size.match = 10, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 4, forbid.treatment.reversal = FALSE)
#' 
#' 
#' 
#' Panel_scatter(non_refined_set = sets0$att,
#'               refined_lists = list(sets1$att, sets2$att),
#'               
#'               data = dem,
#'               covariates = c("y", "tradewb"))

#' }
#'
#' @export
PanelScatter <- function(non_refined_set, refined_list,
                          xlim = c(0, .8),
                          ylim = c(0, .8),
                          main = "",
                          #legend_text = legend_text,
                          pchs = c(2,3),
                          covariates, data) {
  # first, get balance for non-refined set
  non_refined_balance <- get_covariate_balance(matched.sets = non_refined_set, 
                        data = data,
                        covariates = covariates)
  
  # second, get balance for refined sets
  refined_balance <- list()
  for (i in 1:length(refined_list)) {
    refined_balance[[i]] <-
      get_covariate_balance(matched.sets = refined_list[[i]], 
                            data = data,
                            covariates = covariates)
  }
  
  # extract values for x-axis from the non-refined sets
  benchmark <- non_refined_balance
  benchmark <- as.vector(benchmark[1:(nrow(benchmark)-1),]) # delete balance results after t-1
  
  # extract values for y-axis from refined sets and delete balance results after t-1
  compared <- sapply(refined_balance, function(x) x <- x[1:(nrow(x)-1),])
  
  plot(abs(as.numeric(benchmark)), 
       abs(as.numeric(compared[,1])), pch = 1,
       xlab = "",
       ylab = "",
       xlim = xlim,
       ylim = ylim,
       main = main,
       font.main = 1)
  # logical statement for the length of the refined_balance
  if (length(refined_balance) > 1) {
    for (j in 2:length(refined_balance)) {
      points(abs(as.numeric(benchmark)), 
             abs(as.numeric(compared[,j])),
             pch = pchs[j])  
    }
  } 
  # else {
  #   points(abs(as.numeric(benchmark)), 
  #          abs(as.numeric(compared[,1])),
  #          pch = pchs[1])  
  # }
  abline(h = 0, lty = "dashed")
  abline(0, 1, lty = 2, col = "red")
  
  # legend(x = 0, y = 1.1,  legend = legend_text,
  #        y.intersp = 0.65,
  #        x.intersp = 0.3,
  #        xjust = 0,
  #        pch = c(j-1, j), pt.cex = 1,
  #        bty = "n", ncol = 1, cex = 1, bg = "white")

  
} 



