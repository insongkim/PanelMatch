PanelEffectPlot <- function(result, CI = .95, atc.reverse = FALSE,
                            bias.correction = TRUE,
                            xlab = "Effect by Time Period",
                            ylab = "Effect Size",x.size = 13, y.size = 15,
                            title = "") {
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
    geom_point(size = 4) + theme_bw() +
    geom_pointrange(aes(ymax = U, ymin = L)) +
    xlab(xlab) +
    ylab(ylab) +
    geom_hline(yintercept = 0, linetype = 2, colour = "red") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(size = x.size),
          axis.title.y = element_text(size = y.size))                     
  
  
}

# P <- rep(NA, length(result$lead))
# for (i in 1:length(result$lead)) {
#   if (bias.correction == TRUE) {
#     P[i] <- 2*result$o.coef[i] - mean(result$boots[,i],na.rm = T)
#   } else {
#     P[i] <- result$o.coef[i]
#   }
#  
# }
# 
# L <- rep(NA, length(result$lead)) {
#   for (i in 1:length(result$lead))
# }
# 
# get_plot <- function(result, bias.correction = TRUE) {
#   if (bias.correction == TRUE) {
#     P <- 2*result$o.coef - colMeans(result$boots, na.rm = T)
#     L <- t(colQuantiles(2*matrix(nrow = result$ITER, ncol = length(result$o.coef), 
#                                  result$o.coef, byrow = TRUE) - result$boots,
#                         probs = c((1-CI)/2, CI+(1-CI)/2), 
#                         na.rm = T, drop = FALSE))[1,] # bc percentile CI
#     U <- t(colQuantiles(2*matrix(nrow = result$ITER, ncol = length(result$o.coef), 
#                                  result$o.coef, byrow = TRUE) - result$boots,
#                         probs = c((1-CI)/2, CI+(1-CI)/2), 
#                         na.rm = T, drop = FALSE))[2,] # bc percentile CI
#   }
#  
# }
# 
# 
