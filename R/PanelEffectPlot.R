#' PanelEffectPlot
#'
#'\code{PanelEffectPlot} visualizes the results returned by \code{PanelEstimate} 
#' when the option \code{inference} in \code{PanelEstimate} is set to 
#' "bootstrap" via 
#' plotting point estimates along with their confidence intervals by time periods
#' specified by the user. 
#'
#' @param result An object returned by PanelEstimate
#' @param CI A numeric value indicating the confidence interval set by the user. 
#' By default, it is .95  
#' @param atc.reverse A logical value indicating whether the user wants to reverse the signs 
#' for effects when the quantity of interest is "atc" in PanelEstimate
#' @param bias.correction A logical value indicating whether bias correction standard errors
#' are used.
#' @param xlab A character string for x-axis label
#' @param ylab A character string for y-axis label
#' @param x.size A numeric value indicating the font size for the x-axis label
#' @param y.size A numeric value indicating the font size for the y-axis label
#' @param title A character string for the title of the plot
#' @param point_size A numeric value indicating the area of the circle that is 
#' going to appear in the plot to represent the point estimate for a given time
#' period. 
#'
#' @return It returns a plot
#' @export
#'
#' @examples \dontrun{
#' matches.cbps <- PanelMatch(lag = 4, max.lead = 4, time.id = "year",
#' unit.id = "wbcode2", treatment = "dem", formula = y ~ dem, method =
#' "CBPS", weighting = FALSE, qoi = "ate", M = 5, data = dem)
#'
#' ## bootstrap
#' 
#' mod.bootSE <- PanelEstimate(lead = 0:4, inference =
#' "bootstrap", matched_sets = matches.cbps, qoi = "att", CI = .95,
#' ITER = 500) 
#' 
#' PanelEffectPlot(mod.bootSE)
#' 
#' }
PanelEffectPlot <- function(result, CI = .95, atc.reverse = FALSE,
                            bias.correction = TRUE,
                            xlab = "Effect by Time Period",
                            ylab = "Effect Size",x.size = 13, y.size = 15,
                            title = "",
                            point_size = 2) {
  x <- NULL
  if (atc.reverse == TRUE) {
    if (result$qoi == "atc") {
      result$o.coef <- -(result$o.coef)
      result$boots <- -(result$boots)
    }
  }
  if (bias.correction == TRUE) {
    P <- 2*result$o.coef - colMeans(result$boots, na.rm = T)
    L <- t(colQuantiles(2*matrix(nrow = result$ITER, ncol = length(result$o.coef), 
                                 result$o.coef, byrow = TRUE) - result$boots,
                        probs = c((1-CI)/2, CI+(1-CI)/2), 
                        na.rm = T, drop = FALSE))[1,] # bc percentile CI
    U <- t(colQuantiles(2*matrix(nrow = result$ITER, ncol = length(result$o.coef), 
                                 result$o.coef, byrow = TRUE) - result$boots,
                        probs = c((1-CI)/2, CI+(1-CI)/2), 
                        na.rm = T, drop = FALSE))[2,] # bc percentile CI
  } else {
    P <- result$o.coef
    
    L <- t(colQuantiles(result$boots,
                        probs = c((1-CI)/2, CI+(1-CI)/2), 
                        na.rm = T, drop = FALSE))[1,]
    U <- t(colQuantiles(result$boots,
                        probs = c((1-CI)/2, CI+(1-CI)/2), 
                        na.rm = T, drop = FALSE))[2,]
  }
  
  df_effect <- data.frame(x = result$lead, P = P, L = L, U = U)
  
  ggplot(df_effect, aes(x = x, y = P)) + 
    scale_x_continuous(breaks = seq(min(df_effect$x), max(df_effect$x), 
                                    by = 1)) +
    geom_point(size = point_size) + theme_bw() +
    geom_pointrange(aes(ymax = U, ymin = L)) +
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = x.size),
          axis.title.y = element_text(size = y.size))                     
  
  
}