#' PanelBCheck
#'
#' @param matched_sets 
#' @param covariate 
#' @param qoi 
#' @param post.treatment 
#' @param vline 
#' @param xlab 
#' @param ylab 
#' @param legend.position 
#' @param theme_bw 
#' @param adjustment 
#' @param int_size 
#' @param plot 
#' @param linetype 
#' @param colour 
#'
#' @return
#'
#' @examples
#' @export
PanelBCheck <- function(matched_sets,
                        covariate = NULL,
                        qoi = NULL,
                        post.treatment = TRUE, 
                        vline = TRUE, xlab = "Time periods", 
                        ylab = "Differences between the treated units 
                        and their synthetic control units",
                        legend.position = "none",
                        theme_bw = FALSE,
                        adjustment = TRUE,
                        int_size = .7,
                        plot = TRUE,
                        linetype = "dashed", colour = "blue"
) {
  lag = matched_sets$lag;lead = matched_sets$max.lead;
  treatment = matched_sets$treatment; dependent = matched_sets$dependent;
  method = matched_sets$method; covariate_names = matched_sets$covariate_names
  
  if (is.null(qoi)) {
    if (matched_sets$qoi == "att") {
      treated_set <- as.data.frame(data.table::rbindlist(lapply(matched_sets$`ATT_matches`, 
                                                                function(x) {
                                                                  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
                                                                  return(x[x$V2 == treated.id,])
                                                                  
                                                                })))
      treated_set$big_L <- rep(1:(lag+1+lead), nrow(treated_set)/(lag+1+lead))
      plot.materials <- lapply(matched_sets$`ATT_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               treated_set = treated_set,
                               method = method,
                               covariate_names = covariate_names,
                               qoi = "att", adjustment = adjustment,
                               data = matched_sets$data)
    } else if (matched_sets$qoi == "atc") {
      treated_set <- as.data.frame(data.table::rbindlist(lapply(matched_sets$`ATC_matches`, 
                                                                function(x) {
                                                                  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
                                                                  return(x[x$V2 == treated.id,])
                                                                  
                                                                })))
      treated_set$big_L <- rep(1:(lag+1+lead), nrow(treated_set)/(lag+1+lead))
      plot.materials <- lapply(matched_sets$`ATC_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               treated_set = treated_set,
                               method = method,
                               covariate_names = covariate_names,
                               qoi = "atc", adjustment = adjustment, 
                               data = matched_sets$data)
    } else {
      stop("Please specify either att or atc for `qoi`.")
    }
  } else {
    if (qoi == "att") {
      treated_set <- as.data.frame(data.table::rbindlist(lapply(matched_sets$`ATT_matches`, 
                                                                function(x) {
                                                                  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
                                                                  return(x[x$V2 == treated.id,])
                                                                  
                                                                })))
      treated_set$big_L <- rep(1:(lag+1+lead), nrow(treated_set)/(lag+1+lead))
      plot.materials <- lapply(matched_sets$`ATT_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               treated_set = treated_set,
                               method = method,
                               covariate_names = covariate_names,
                               qoi = "att", adjustment = adjustment,
                               data = matched_sets$data)
    } else if (qoi == "atc") {
      treated_set <- as.data.frame(data.table::rbindlist(lapply(matched_sets$`ATC_matches`, 
                                                                function(x) {
                                                                  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
                                                                  return(x[x$V2 == treated.id,])
                                                                  
                                                                })))
      treated_set$big_L <- rep(1:(lag+1+lead), nrow(treated_set)/(lag+1+lead))
      plot.materials <- lapply(matched_sets$`ATC_matches`, 
                               gaps_plot_tmp, lag = lag, lead = lead,
                               covariate = covariate,
                               treated_set = treated_set,
                               method = method,
                               covariate_names = covariate_names,
                               qoi = "atc", adjustment = adjustment, 
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
  
  # t <- tapply(df$val, df$x, quantile) # quantile hashed out
  box <- boxplot(df$val ~ df$x, plot = F, range = 0)
  
  
  # # aggregate(val ~ x, data = df, FUN = quantile)
  # y2 <- apply(do.call(rbind,t),2,  FUN = unlist)[,1]
  # for (i in 2:5) {
  #   y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  # }
  
  df2 <- data.frame(x2=rep(-lag:lead,5), 
                    y2=unlist(as.list(t(box$stats))),
                    g2=gl(5,(lag+1+lead), 
                          labels=c("Lowest", "1st Qu.", "Median",
                                   "3rd Qu.", "Highest")))
  
  if(post.treatment == FALSE){
    df <- df[df$x < 0, ]
    df2 <- df2[df2$x2 < 0, ]
  }
  
  if (plot == FALSE){
    colnames(df2) <- c("Time to Treatment",
                       "Balance",
                       "Interquartiles")
    output <- tapply(df$val, df$x, mean, na.rm = T)
    output <- data.frame("time to treatment" = as.numeric(names(output)), 
                         "balance" = output)
    return(list("Mean_Balance" = output,
                "Interquartile_Balance" = df2,
                "All_Balance" = df))
  } else {
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
    p <- p + geom_line(aes(x=factor(x2), y=y2, group = g2), size = int_size,
                       colour = "black",
                       df2) + 
      if (theme_bw == TRUE) {
        theme_bw()
      } 
    p + theme(legend.position = legend.position) +
      if(vline == TRUE) {
        geom_vline(xintercept = lag, linetype = linetype, 
                   colour = colour)
      } else {
        geom_vline()
      }  
    
  }
  
  
  
}