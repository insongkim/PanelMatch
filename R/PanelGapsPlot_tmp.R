PanelGapsPlot_tmp <- function(matched_sets,
                          covariate = NULL,
                          qoi = NULL,
                          post.treatment = TRUE, 
                          vline = TRUE, xlab = "Time periods", 
                          ylab = "Differences between the treated units 
                          and their synthetic control units",
                          legend.position = "none",
                          theme_bw = FALSE,
                          linetype = "dashed", colour = "blue"
                          ) {
  lag = matched_sets$lag;lead = matched_sets$max.lead;
  treatment = matched_sets$treatment; dependent = matched_sets$dependent
  
  if (is.null(qoi)) {
    if (matched_sets$qoi == "att") {
      plot.materials <- lapply(matched_sets$`Matched sets for ATT`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "att",
                               data = matched_sets$data)
    } else if (matched_sets$qoi == "atc") {
      plot.materials <- lapply(matched_sets$`Matched sets for ATC`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "atc",
                               data = matched_sets$data)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  } else {
    if (qoi == "att") {
      plot.materials <- lapply(matched_sets$`Matched sets for ATT`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "att",
                               data = matched_sets$data)
    } else if (qoi == "atc") {
      plot.materials <- lapply(matched_sets$`Matched sets for ATC`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "atc",
                               data = matched_sets$data)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  }

  
  df <- data.frame(x=rep(-lag:lead, length(plot.materials)), 
                   val=unlist(lapply(plot.materials, function(x) x$gap)), 
                   variable=rep(paste0("unit", 
                                       unlist(lapply(plot.materials, function(x) x$unit))), 
                                each=lag+1+lead))
  
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
