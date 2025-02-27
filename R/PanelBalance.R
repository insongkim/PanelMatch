#' Plot covariate balance results
#' Create figures displaying covariate balance results for one or more \code{PanelMatch} configurations. Users can customize these visualizations. 
#' @param x \code{PanelBalance} object
#' @param ... additional parameters to be passed to \code{base::plot()}
#' @param type character specifying which type of plot to produce. Can be either "panel" or "scatter". When "panel," covariate balance results for covariates are shown over the lag period. When "scatter," the figure has the following characteristics. Each point on the plot represents a specific covariate at a particular time period in the lag window from t-L to t-1. The horizontal axis represents the covariate balance for this particular variable and time period before refinement is applied, while the vertical axis represents the post-refinement balance value.
#' @param reference.line logical. Include a reference line at y = 0? Only applicable to the panel plot.
#' @param legend logical. Describes whether or not to include a legend.
#' @param ylab character. Y-axis label. 
#' @param include.treatment.period Logical. Describes whether or not the treatment period should be included on the panel plot. Default is TRUE.
#' @param include.unrefined.panel logical indicating whether or not unrefined balance plots should be returned for panel plot. Only applicable to panel plot. Default is TRUE.
#' @param legend.position character. Describes where the legend should be placed on the figure. Uses base R syntax. 
#' @return returns a set of base R plots, depending on the specification of "panel" or "scatter" above. When \code{type = "panel"} and \code{include.unrefined.panel = TRUE}, two sets of plots are returned. The first set shows covariate balance levels for the specified \code{PanelMatch} configurations. The second set shows covariate balance levels for the same \code{PanelMatch} configurations, but with all control units receiving equal weight (i.e., balance levels prior to refinement). If \code{include.unrefined.panel = FALSE}, only the first set of figures are returned. The sets of figures are both returned in the same order as the \code{PanelMatch} configurations specified to \code{get_covariate_balance()} that compose the \code{PanelBalance} object. When \code{type = "scatter"}, the visualization described above is produced, with all configurations shown on the same plot with different symbols.
#' @export
#'
#' @examples
#' dem$rdata <- runif(runif(nrow(dem)))
#' dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
#'                      panel.data = dem.panel, match.missing = TRUE,
#'                      covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                      size.match = 5, qoi = "att")
#' 
#' # create multiple configurations to compare
#' pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.match", 
#'                   panel.data = dem.panel, match.missing = TRUE,
#'                   covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                   size.match = 5, qoi = "att")
#' 
#' pb <- get_covariate_balance(pm.obj, pm2,
#'                             include.unrefined = TRUE,
#'                             panel.data = dem.panel, 
#'                             covariates = c("tradewb", "rdata"))
#' plot(pb, type = "panel", include.unrefined.panel = TRUE)
#' plot(pb, type = "scatter")
#' # only show refined balance figures
#' plot(pb, type = "panel", include.unrefined.panel = FALSE)
plot.PanelBalance <- function(x, 
                              ..., 
                              type = "panel",
                              reference.line = TRUE,
                              legend = TRUE, 
                              ylab = NULL,
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
        if (is.null(ylab))
        {
          ylab <- "SD"
        }
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
      if (is.null(ylab))
      {
        ylab <- "SD"
      }
      plot(unrefined.mat, type = "panel", reference.line = reference.line,
           legend = legend,
           ylab = ylab, 
           include.treatment.period = include.treatment.period,
           include.unrefined.panel = FALSE,
           legend.position = legend.position, 
           ...)
    }
  } else { # should be scatter
    process_scatter_balance(x, ylab = ylab, ...)
  }
  
}

#' Summarize covariate balance over time
#'
#' @param object \code{PanelBalance} object 
#' @param qoi Character. Valid values include "att", "art", or "atc". Specifying which QOI information to extract and summarize.
#' @param include.unrefined logical. Indicates whether or not unrefined balance results should be included in the summary.
#' @param unrefined.only logical. Indicates whether or not only unrefined balance results should be included in the summary.
#' @param ... Not used
#'
#' @return returns a list of matrices with covariate balance results calculated. Each element in the list corresponds to a \code{PanelMatch} configuration given to \code{get_covariate_balance()} and are returned in order. Note that if a configuration has \code{qoi = "ate"}, the corresponding element in the returned list will also be a list, containing balance results corresponding to the ATT and ATC. Otherwise, each element in the returned list will be a matrix. Each matrix entry corresponds to balance results for a particular covariate in a particular period. When unrefined balance results are included, users will see additional columns with "_unrefined" appended to covariate names. These correspond to the unrefined balance results for a particular covariate-period. 
#' @export
#'
#' @examples
#' dem$rdata <- runif(runif(nrow(dem)))
#' dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
#'                      panel.data = dem.panel, match.missing = TRUE,
#'                      covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                      size.match = 5, qoi = "att")
#' pb <- get_covariate_balance(pm.obj,
#'                             include.unrefined = TRUE,
#'                             panel.data = dem.panel, 
#'                             covariates = c("tradewb", "rdata"))
#' summary(pb)
summary.PanelBalance <- function(object, qoi = NULL, 
                                 include.unrefined = TRUE, 
                                 unrefined.only = FALSE,
                                 ...)
{
  
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
#' @param x \code{PanelBalance} object 
#' @param ... Not used
#'
#' @return Nothing
#' @export
#'
#' @examples
#' dem$rdata <- runif(runif(nrow(dem)))
#' dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
#'                      panel.data = dem.panel, match.missing = TRUE,
#'                      covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                      size.match = 5, qoi = "att")
#' pb <- get_covariate_balance(pm.obj,
#'                             include.unrefined = TRUE,
#'                             panel.data = dem.panel, 
#'                             covariates = c("tradewb", "rdata"))
#' print(pb)
print.PanelBalance <- function(x, ...)
{

  for (j in x) {
    k <- as.list(j)
    print(k)
  }
  
}

#' Subset PanelBalance objects
#'
#' @param x \code{PanelBalance} object 
#' @param i numeric. Specifies which element to extract. Substantively, it specifies which \code{PanelMatch} configuration data to extract.
#' @param ... Not used
#'
#' @return Returns balance information for specified \code{PanelMatch} configuration. Note that results are still returned as a \code{PanelBalance} object. In order to return a list, use the [[ operator
#' @export
#' @examples
#' dem$rdata <- runif(runif(nrow(dem)))
#' dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
#'                      panel.data = dem.panel, match.missing = TRUE,
#'                      covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                      size.match = 5, qoi = "att")
#' 
#' # create multiple configurations to compare
#' pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.match", 
#'                   panel.data = dem.panel, match.missing = TRUE,
#'                   covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                   size.match = 5, qoi = "att")
#' 
#' pb <- get_covariate_balance(pm.obj, pm2,
#'                             include.unrefined = TRUE,
#'                             panel.data = dem.panel, 
#'                             covariates = c("tradewb", "rdata"))
#' bal.maha <- pb[1]
#' bal.ps <- pb[2]                             
`[.PanelBalance` <- function(x, i, ...) {
  # Subset the list while keeping the class and attributes
  tmp <- x
  tru <- attr(x, "unrefined.balance.results")
  if (!is.null(tru))
  {
    ubr <- tru[i]  
  } else {
    ubr <- NULL
  }
  class(tmp) <- "list"
  subset <- tmp[i]
  attributes(subset) <- attributes(x)
  attr(subset, "unrefined.balance.results") <- ubr
  class(subset) <- class(x)
  return(subset)
}

#' Extract just the unrefined covariate balance results, if they exist
#' @param pb.object \code{PanelBalance} object
#' 
#' @export
get_unrefined_balance <- function(pb.object) {
  UseMethod("get_unrefined_balance", pb.object)
}

#' Extract unrefined covariate balance results, if they exist
#'
#' @param pb.object \code{PanelBalance} object
#'
#' @return A \code{PanelBalance} object, with just the unrefined balance results
#' @export
#' @examples
#' dem$rdata <- runif(runif(nrow(dem)))
#' dem.panel <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "mahalanobis", 
#'                      panel.data = dem.panel, match.missing = TRUE,
#'                      covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                      size.match = 5, qoi = "att")
#' 
#' # create multiple configurations to compare
#' pm2 <- PanelMatch(lead = 0:3, lag = 4, refinement.method = "ps.match", 
#'                   panel.data = dem.panel, match.missing = TRUE,
#'                   covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                   size.match = 5, qoi = "att")
#' 
#' pb <- get_covariate_balance(pm.obj, pm2,
#'                             include.unrefined = TRUE,
#'                             panel.data = dem.panel, 
#'                             covariates = c("tradewb", "rdata"))
#' get_unrefined_balance(pb)
get_unrefined_balance.PanelBalance <- function(pb.object) {
  return(attr(pb.object, "unrefined.balance.results"))
}


# Helper functions for plotting

process_scatter_balance <- function(PanelBalance, 
                                    xlim = c(0, .8),
                                    ylim = c(0, .8),
                                    main = "Standardized Mean Difference of Covariates",
                                    #pchs = c(2,3),
                                    #covariates, data,
                                    xlab = "Before refinement",
                                    ylab = "After refinement", 
                                    ...)
{
  unrefined.balance <- attr(PanelBalance, "unrefined.balance.results")
  if (is.null(unrefined.balance))
  {
    stop("Unrefined balance results required. Please re-run covariate balance calculations with include.unrefined argument = TRUE")  
  }

  
  unrefined.vectors <- unlist(lapply(unrefined.balance, 
                                     function(x) lapply(x, as.vector)), recursive = FALSE)
  refined.vectors <- unlist(lapply(PanelBalance, 
                                   function(x) lapply(x, as.vector)), recursive = FALSE)
  

  
  for (i in 1:length(unrefined.vectors)) {
    
    unrefined.pm.vector <- unrefined.vectors[[i]]
    refined.pm.vector <- refined.vectors[[i]]
    if (i == 1)
    {
      graphics::plot(abs(as.numeric(unrefined.pm.vector)),
                     abs(as.numeric(refined.pm.vector)),
                     pch = 1,
                     xlab = xlab,
                     ylab = ylab,
                     xlim = xlim,
                     ylim = ylim,
                     main = main,
                     font.main = 1,
                     ...)
    } else
    {
      graphics::points(abs(unrefined.pm.vector),
                       abs(refined.pm.vector),
                       pch = i)
    }
  }
  graphics::abline(h = 0, 
                   lty = "dashed")
  graphics::abline(0, 1, 
                   lty = 2,
                   col = "red")
  
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