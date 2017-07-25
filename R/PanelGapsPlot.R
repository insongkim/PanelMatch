PanelGapsPlot <- function(matched_sets, lag, lead, treatment,
                          show.covariate = NULL, 
                          post.treatment = TRUE, 
                          vline = TRUE, xlab = "Time periods", 
                          ylab = "Differences between the treated units 
                          and their synthetic control units",
                          legend.position = "none",
                          dependent, theme_bw = FALSE,
                          linetype = "dashed", colour = "blue"
                          ) {
  # matched_sets <- PanelObject$matched_sets
  if (length(matched_sets) !=3) {
    stop("Gaps plot is only used for att or atc")
  }
  plot.materials <- lapply(matched_sets[[names(matched_sets)[2]]], 
                           gaps_plot, lag = lag, lead = lead,
                           show.covariate = show.covariate)
  
  df <- data.frame(x=rep(-lag:lead, length(plot.materials)), 
                   val=unlist(lapply(plot.materials, function(x) x$gap)), 
                   variable=rep(paste0("unit", 
                                       unlist(lapply(plot.materials, function(x) x$unit))), 
                                each=lag+1+lead))
  
  if(is.null(show.covariate)) {
    df$val <- df$val/sd(matched_sets$data[matched_sets$data[treatment] == 1, dependent], 
                        na.rm = T)
  } else {
    df$val <- df$val/sd(matched_sets$data[matched_sets$data[treatment] == 1, show.covariate], 
                        na.rm = T)
  }
  
  t <- tapply(df$val, df$x, quantile)
  
  # aggregate(val ~ x, data = df, FUN = quantile)
  y2 <- apply(do.call(rbind,t),2,  FUN = unlist)[,1]
  for (i in 2:5) {
    y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  }
  
  df2 <- data.frame(x2=rep(-lag:lead,5), 
                    y2=y2,
                    g2=gl(5,(lag+1+lead), 
                          labels=c("Min", "1st Qu.", "Median",
                                   "3rd Qu.", "Max.")))
  
  if(post.treatment == FALSE){
    df <- df[df$x < 0, ]
    df2 <- df2[df2$x2 < 0, ]
  }
  
  p <- ggplot() + 
    geom_line(aes(x=factor(x), y=val, group = variable,
                  colour = variable), df) +
    labs(x = xlab, y = ylab) + 
    if (post.treatment == FALSE) {
      scale_x_discrete(breaks = -lag:-1,
                       labels=paste("t", -lag:-1, sep = "")) 
    } else {
      scale_x_discrete(breaks = -lag:lead,
                       labels=c(paste("t", -lag:-1, sep = ""),
                                paste("t+", 0:lead, sep = ""))) 
    } 
    p <- p + geom_line(aes(x=factor(x2), y=y2, group = g2), 
              colour = "black",
              df2) + 
      if (theme_bw == TRUE) {
        theme_bw()
      } 
    p + theme(legend.position = legend.position) +
    if(vline == TRUE) {
      geom_vline(xintercept = lag+1, linetype = linetype, 
                 colour = colour)
    } else {
      geom_vline()
    }  

}
