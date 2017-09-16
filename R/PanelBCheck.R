PanelBCheck <- function(matched_sets,
                              qoi = NULL,
                            #  vline = TRUE, xlab = "Time periods", 
                             # ylab = "Average Outcome across 
                            #  \n Treatment and Control Groups",
                              legend.position = "right",
                             # theme_bw = TRUE,
                              adjustment = FALSE
                             # linetype = "dashed", colour = "blue"
) {
  # if (theme_bw == TRUE) {
  #   white <- theme_bw()
  # } else {
  #   white <- NULL
  # }
  lag = matched_sets$lag;lead = matched_sets$max.lead;
  # treatment = matched_sets$treatment; dependent = matched_sets$dependent
  
  if (is.null(qoi)) {
    if (matched_sets$qoi == "att") {
      outcomes <- lapply(matched_sets$ATT_matches, parallel_trends, 
                         adjustment = adjustment,
                         lead = matched_sets$max.lead,
                         lag = matched_sets$lag)
    } else if (matched_sets$qoi == "atc") {
      outcomes <- lapply(matched_sets$ATC_matches, parallel_trends, 
                         adjustment = adjustment,
                         lead = matched_sets$max.lead,
                         lag = matched_sets$lag)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  } else {
    if (qoi == "att") {
      outcomes <- lapply(matched_sets$ATT_matches, parallel_trends, 
                         adjustment = adjustment,
                         lead = matched_sets$max.lead,
                         lag = matched_sets$lag)
    } else if (qoi == "atc") {
      outcomes <- lapply(matched_sets$ATC_matches, parallel_trends, 
                         adjustment = adjustment,
                         lead = matched_sets$max.lead,
                         lag = matched_sets$lag)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  }
  
  treated.outcome <- rowMeans(sapply(outcomes, function (x) x$treated.outcome))
  control.outcome <- rowMeans(sapply(outcomes, function (x) x$control.outcome))
  return(list("treated.outcome" = treated.outcome,
         "control.outcome" = control.outcome))
  
  # 
  # df <- data.frame(x=rep(-lag:lead, 2), 
  #                  val=c(treated.outcome, control.outcome), 
  #                  group=c(rep("treated",
  #                              each=lag+1+lead), rep("control", 
  #                                                    each = lag + 1 + lead))) 
  # 
  # p <- ggplot() + 
  #   geom_line(aes(x=factor(x), y=val, group = group,
  #                 colour = group), df) +
  #   labs(x = xlab, y = ylab) + 
  #   scale_x_discrete(breaks = -lag:lead,
  #                    labels=c(paste("t", -lag:-1, sep = ""),
  #                             paste("t+", 0:lead, sep = ""))) 
  # 
  # p <- p + white
  # 
  # p + theme(legend.position = legend.position) +
  #   if(vline == TRUE) {
  #     geom_vline(xintercept = lag, linetype = linetype, 
  #                colour = colour)
  #   } else {
  #     geom_vline()
  #   }  
  
}
