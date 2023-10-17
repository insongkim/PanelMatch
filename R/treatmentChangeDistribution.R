#' treatmentChangeDistribution
#'
#' Calculates results for a placebo test
#'
#'
#' Plots the levels of treatment at time = t-1 and time = t
#' 
#' @export
#' 
#' 
treatmentChangeDistribution <- function(pm.object, 
                                        data.in,
                                        type, 
                                        x.label = NULL, 
                                        y.label = NULL,
                                        title = NULL,
                                        bins = NULL,
                                        jitter = FALSE,
                                        ...)
{
  if (identical(pm.object[["qoi"]], "ate"))
  {
    stop("Diagnostic not available for ate")
  }
  qoi.in <-  attr(pm.object, "qoi")
  matched.sets <- pm.object[[qoi.in]]
  
  
  tvar <- attr(matched.sets, "t.var")
  id.var <- attr(matched.sets, "id.var")
  treatment.var <- attr(matched.sets, "treatment.var")
  rownames(data.in) <- paste0(data.in[, id.var], ".", data.in[, tvar])
  pretreatment.ts <- as.integer(sub(".*\\.", "", names(matched.sets))) - 1
  treated.ids <- as.integer(sub("\\..*", "", names(matched.sets)))
  t.treatment <- data.in[names(matched.sets), treatment.var]
  t1.treatment <- data.in[paste0(treated.ids, ".", pretreatment.ts), treatment.var]
  
  args <- list(...)
  if (all(!"xlab" %in% args))
  {
    xlab <- "Treatment level, time = t-1"
  }
  
  if (all(!"ylab" %in% args))
  {
    ylab <- "Treatment level, time = t-1"
  }
  
  if (all(!"main" %in% args))
  {
    main.t <- "Treatment levels at t-1 vs t for Treated Units"
  }

  
  plot(x = t1.treatment,
       y = t.treatment,
       xlab = xlab,
       ylab = ylab,
       main = main.t,
       pch = 16,
       col=rgb(red=0, green=0, blue=0, alpha=0.2),
       ...)
  
}