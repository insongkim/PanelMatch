PanelCaliper <- function(matched_sets,
                              covariate = NULL,
                              qoi = NULL,
                              post.treatment = TRUE, number = .1) {
  lag = matched_sets$lag;lead = matched_sets$max.lead;
  treatment = matched_sets$treatment; dependent = matched_sets$dependent
  
  if (is.null(qoi)) {
    if (matched_sets$qoi == "att") {
      plot.materials <- lapply(matched_sets$`ATT_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "att",
                               data = matched_sets$data)
    } else if (matched_sets$qoi == "atc") {
      plot.materials <- lapply(matched_sets$`ATC_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "atc",
                               data = matched_sets$data)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  } else {
    if (qoi == "att") {
      plot.materials <- lapply(matched_sets$`ATT_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               qoi = "att",
                               data = matched_sets$data)
    } else if (qoi == "atc") {
      plot.materials <- lapply(matched_sets$`ATC_matches`, 
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
  
  # to_drop <- sort(abs(tapply(df$val[df$x < 0], df$variable[df$x <0], mean)), decreasing = T)[1:number]
  
  if (is.null(qoi)) {
    if (matched_sets$qoi == "att") {
      ind <- as.logical(sapply(matched_sets$`ATT_matches`, qoi = matched_sets$qoi,
                               gaps_caliper, lag = lag, covariate = covariate,
                               data = matched_sets$data) < number)
      matched_sets$`ATT_matches` <- matched_sets$`ATT_matches`[which(ind)]
    } else if (matched_sets$qoi == "atc") {
      ind <- as.logical(sapply(matched_sets$`ATC_matches`, qoi = matched_sets$qoi,
                               gaps_caliper, lag = lag, covariate = covariate,
                               data = matched_sets$data) < number)
      matched_sets$`ATC_matches` <- matched_sets$`ATC_matches`[which(ind)]
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  } else {
    if (qoi == "att") {
      ind <- as.logical(sapply(matched_sets$`ATT_matches`, qoi = qoi,
                               gaps_caliper, lag = lag, covariate = covariate,
                               data = matched_sets$data) < number)
      matched_sets$`ATT_matches` <- matched_sets$`ATT_matches`[which(ind)]
    } else if (qoi == "atc") {
      ind <- as.logical(sapply(matched_sets$`ATC_matches`, qoi = qoi,
                               gaps_caliper, lag = lag, covariate = covariate,
                               data = matched_sets$data) < number)
      matched_sets$`ATC_matches` <- matched_sets$`ATC_matches`[which(ind)]
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  }
  

  return(matched_sets)

}  

