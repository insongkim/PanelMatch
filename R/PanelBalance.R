#' plot covariate balance results
#'
#' @param x PanelBalance object
#' @param ... additional parameters to be passed to \code{base::plot()}
#' @param include.unrefined.panel logical indicating whether or not unrefined balance plots should be returned for panel plot. NEEDS TO BE IMPLEMENTED. 
#' @return
#' @export
#'
#' @examples
plot.PanelBalance <- function(x, 
                              ..., 
                              type = "panel",
                              reference.line = TRUE,
                              legend = TRUE, 
                              ylab = "SD",
                              include.treatment.period = TRUE,
                              include.unrefined.panel = TRUE,
                              legend.position = "topleft")
{
  if (!type %in% c("panel", "scatter"))
  {
    stop("type must be 'panel' or 'scatter'")
  }
  attr(x, 'treatment') -> treatment
  if (type == "panel")
  {
    for (pb in x) {
      for (pointmat in pb) {
        
        panel_plot_helper(treatment, pointmat,
                    reference.line,
                    legend,
                    ylab,
                    include.treatment.period,
                    legend.position,
                    ...)
        
      }
    }
    if (include.unrefined.panel)
    {
      unrefined <- attr(x, "unrefined.balance.results")
      if (is.null(unrefined)) stop("unrefined balance results not found.")
      unrefined.mat <- unrefined
    } else {
      unrefined.mat <- NULL
    }
    if (!is.null(unrefined.mat))
    {
      plot(unrefined.mat, type = "panel", reference.line = reference.line,
           legend = legend,
           ylab = ylab, 
           include.treatment.period = include.treatment.period,
           include.unrefined.panel = FALSE,
           legend.position = legend.position, 
           ...)
    }
  } else { # should be scatter
    process_scatter_balance(x, ...)
  }
  
}

panel_plot_helper <- function(treatment,
                        pointmatrix,
                        reference.line,
                        legend, 
                        ylab,
                        include.treatment.period,
                        legend.position,
                        ...)
{
  
  if (!include.treatment.period)
  {
    pointmatrix <- pointmatrix[-nrow(pointmatrix), ,drop = FALSE]
    stop.val <- 1
    start.val <- nrow(pointmatrix)
  } else {
    stop.val <- 0
    start.val <- nrow(pointmatrix) - 1
  }
  
  treated.included <- treatment %in% colnames(pointmatrix)
  
  if (treated.included)
  {
    treated.data <- pointmatrix[,which(colnames(pointmatrix) == treatment)] # treated data
    pointmatrix <- pointmatrix[,-which(colnames(pointmatrix) == treatment)] #all non-treatment variable data
    graphics::matplot(pointmatrix, pch = 19,
                      type = "b", 
                      col = 1:ncol(pointmatrix), 
                      lty = 1, ylab = ylab, xaxt = "n", ...)
    graphics::lines(x = 1:nrow(pointmatrix), 
                    y = as.numeric(treated.data), 
                    type = "b",
                    lty = 2, lwd = 3)
    graphics::axis(side = 1, labels = paste0("t-", start.val:stop.val), 
                   at = 1:nrow(pointmatrix), ...)  
  } else
  {
    
    graphics::matplot(pointmatrix, type = "b",
                      pch = 19,
                      col = 1:ncol(pointmatrix), 
                      lty = 1, ylab = ylab, xaxt = "n", ...)
    graphics::axis(side = 1, labels = paste0("t-", start.val:stop.val), 
                   at = 1:nrow(pointmatrix), ...)  
  }
  
  if (legend) {
    if (treated.included)
    {
      legend(legend.position, 
             legend = c(colnames(pointmatrix), treatment), 
             col = c(1:ncol(pointmatrix), "black"), 
             lty = c(rep(1, ncol(pointmatrix)), 2))  
    } else {
      legend(legend.position, 
             legend = colnames(pointmatrix), 
             col = 1:ncol(pointmatrix), lty = 1)
    }
    
  }
  if(reference.line) graphics::abline(h = 0, lty = "dashed")
  
}

#' Summarize covariate balance over time (placeholder for now)
#'
#' @param PanelBalance object 
#' @param qoi Character string "att", "art", or "atc" -- specifying which qoi to extract.
#' @param refined logical. return refined or unrefined balance results. Default is TRUE.
#' @param ... Not used
#'
#' @return
#' @export
#'
#' @examples
summary.PanelBalance <- function(object, qoi = NULL, 
                                 include.unrefined = TRUE, 
                                 unrefined.only = FALSE,
                                 ...)
{
  #if (missingArg(qoi)) stop("Please specify a qoi")
  if (is.null(qoi))
  {
    qois <- unique(names(unlist(object, recursive = FALSE)))
    if (length(qois) > 1) 
    {
      stop("More than one QOI present in PanelBalance object. Please specify a QOI.")
    } else {
      qoi <- qois
    }  
  } else {
    if (length(qoi) > 1)
    {
      stop("please specify one qoi.")
    }
    if (!qoi %in% c("att", "art", "atc"))
    {
      stop("please specify a valid QOI: att, art, or atc.")
    }
  } 
  
  if (include.unrefined || unrefined.only)
  {
    if (is.null(attr(object, "unrefined.balance.results")))
    {
      stop("Unrefined balance results not found. Set include.unrefined = TRUE when using get_covariate_balance function or set include.unrefined = FALSE.")
    }
  }
  
  if (unrefined.only && include.unrefined)
  {
    warning("Note: only unrefined balance results are shown. Consider setting unrefined.only = FALSE and include.unrefined = TRUE?")
  }
  u.list <- NULL
  if (!is.null(attr(object, "unrefined.balance.results")) && 
      include.unrefined)
  {
    unrefined <-attr(object, "unrefined.balance.results")
    class(unrefined) <- c("list")
    u.list <- lapply(unrefined, 
                       function(x) x[[qoi]])
    
    u.list <- Filter(Negate(is.null), u.list)
  }
  
  if (unrefined.only)
  {
    object <- attr(object, "unrefined.balance.results")
    class(object) <- c("list")
  }
  
  class(object) <- c("list")
  ret.list <- lapply(object, 
         function(x) x[[qoi]])
  ret.list <- Filter(Negate(is.null), ret.list)
  if (!unrefined.only)
  {
    if (!is.null(u.list))
    {
      cbind_refined_results <- function(unrefined, refined)
      {
        colnames(unrefined) <- paste0(colnames(unrefined), "_unrefined")
        res <- cbind(unrefined, refined)
        return(res)
      }
      ret.list <- mapply(cbind_refined_results, 
                         unrefined = u.list, 
                         refined = ret.list, SIMPLIFY = FALSE)
    }
  }
  
  
  return(ret.list)
}


#' Print basic information about PanelBalance objects
#'
#' @param x PanelBalance object
#' @param ... Not used
#'
#' @return
#' @export
#'
#' @examples
print.PanelBalance <- function(x, ...)
{
  # print("PanelBalance object...")
  # print(paste("Contains balance statistics for", length(x), "configuration(s)"))
  # cov.str <- paste0(attr(x, "covariates"), collapse = ", ")
  # print(paste("Contains calculations for following covariates: ", cov.str))
  
  for (j in x) {
    k <- as.list(j)
    print(k)
  }
  
}

#' Title
#'
#' @param x 
#' @param i 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
`[.PanelBalance` <- function(x, i, ...) {
  # Subset the list while keeping the class and attributes
  #subset <- x[i]
  tmp <- x
  class(tmp) <- "list"
  subset <- tmp[i]
  attributes(subset) <- attributes(x)
  class(subset) <- class(x)
  return(subset)
}

#' #' Title
#' #'
#' #' @param x 
#' #' @param i 
#' #' @param ... 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' `[[.PanelBalance` <- function(x, i, ...) {
#'   # Subset the list while keeping the class and attributes
#'   #subset <- x[[i]]
#'   tmp <- x
#'   class(tmp) <- "list"
#'   subset <- base::`[[`(tmp, i, ...)
#'   return(subset)
#' }