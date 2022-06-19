#' DisplayTreatment
#' 
#' \code{DisplayTreatment} visualizes the treatment distribution
#' across units and time in a panel dataset
#'
#' @param unit.id Name of the unit identifier variable as a character string
#' @param time.id Name of the time identifier variable as a character string
#' @param treatment Name of the treatment variable as a character string
#' @param data data.frame that contains the time series cross sectional data used for matching and estimation. Unit and time data must be integers. Time data must also be formatted as sequential integers that increase by one.
#' @param color.of.treated Color of the treated observations provided as a character string (this includes hex values). Default is red.
#' @param color.of.untreated Color of the untreated observations provided as a character string (this includes hex values). Default is blue.
#' @param title Title of the plot provided as character string
#' @param xlab Character label of the x-axis
#' @param ylab Character label of the y-axis
#' @param x.size Numeric size of the text for xlab or x axis tick labels. Assign x.size = NULL to use built in ggplot2 method of determining label size. 
#' When the length of the time period is long, consider setting to NULL and adjusting size and ratio of the plot.
#' @param y.size Numeric size of the text for ylab or y axis tick labels. Assign y.size = NULL to use built in ggplot2 method of determining label size. 
#' When the number of units is large, consider setting to NULL and adjusting size and ratio of the plot.
#' @param x.angle Angle (in degrees) of the tick labels for x-axis
#' @param y.angle Angle (in degrees) of the tick labels for y-axis
#' @param legend.position Position of the legend. Provide this according to ggplot2 standards. 
#' @param legend.labels Character vector of length two describing the
#' labels of the legend to be shown in the plot. ggplot2 standards are used.
#' @param decreasing Logical. Determines if display order should be increasing or decreasing by the amount of treatment received. Default is \code{decreasing} = FALSE.
#' @param matched.set a matched.set object (optional) containing a single treated unit and a set of matched controls. If provided, this set will be highlighted on the resulting plot.
#' @param show.set.only logical. If TRUE, only the treated unit and control units contained in the provided \code{matched.set} object will be shown on the plot. 
#' Default is FALSE. If no \code{matched.set} is provided, then this argument will have no effect.
#' @param gradient.weights logical. If TRUE, the "darkness"/shade of units in the provided \code{matched.set} object will be displayed according to their weight. Control units with higher weights will appear darker on the resulting plot. Control units with lower weights will appear lighter. This argument has no effect unless a \code{matched.set} is provided.
#' @param hide.x.tick.label logical. If TRUE, x axis tick labels are not shown. Default is FALSE. 
#' @param hide.y.tick.label logical. If TRUE, y axis tick labels are not shown. Default is FALSE.
#' @param dense.plot logical. if TRUE, lines between tiles are removed on resulting plot. This is useful for producing more readable plots in situations where the number of units and/or time periods is very high.
#' @return \code{DisplayTreatment} returns a treatment variation plot (using ggplot2),
#' which visualizes the variation of treatment across unit and time.
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <amrauh@umich.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples 
#' 
#' DisplayTreatment(unit.id = "wbcode2",
#'                  time.id = "year", legend.position = "none",
#'                  xlab = "year", ylab = "Country Code",
#'                  treatment = "dem", data = dem)
#' 
#'
#' @export
DisplayTreatment <- function(unit.id, time.id, treatment, data, 
                              color.of.treated = "red",
                              color.of.untreated = "blue", 
                              title = "Treatment Distribution \n Across Units and Time",
                              xlab = "Time", ylab = "Unit",
                              x.size = NULL, y.size = NULL,
                              legend.position= "none",
                              x.angle = NULL,
                              y.angle = NULL,
                              legend.labels = c("not treated", "treated"),
                              decreasing = FALSE,
                              matched.set = NULL,
                              show.set.only = FALSE,
                              hide.x.tick.label = FALSE,
                              hide.y.tick.label = FALSE,
                              gradient.weights = FALSE,
                              dense.plot = FALSE)

{
  
  alphaweight <- NULL #for some reason --as-cran checks need this
  
  
  if (!inherits(data, "data.frame")) stop("please convert data to data.frame class")
  
  if (any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  # if(!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if (!inherits(data[, time.id], "integer")) stop("please convert time id to consecutive integers")
  #if (class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  
  if (gradient.weights && is.null(matched.set)) stop("gradient.weights cannot be TRUE without a provided matched set")
  if (!gradient.weights && !show.set.only && !is.null(matched.set)) {
    warning("gradient.weights, show.set.only set to FALSE, but matched set provided. Ignoring matched set")
    matched.set <- NULL
  }
  
  data <- na.omit(data[c(unit.id, time.id, treatment)])
  # rename variables to match with the object names in the loop below
  colnames(data) <- c("unit.id", "time.id", "treatment")  
  
  data$trintens <- as.numeric(tapply(data$treatment, data$unit.id, mean, na.rm = TRUE)[as.character(data$unit.id)])
  data <- data[order(data$trintens, decreasing = decreasing), ] 
  
  
  data$unit.id <- factor(data$unit.id, 
                         levels = unique(as.character(data$unit.id)))
  data$time.id <- factor(x = data$time.id, 
                         levels = sort(unique(data$time.id)), 
                         ordered = TRUE)
  
  if (!is.null(matched.set))
  {
     #should only be one
    t.id <- as.numeric(unlist(strsplit(names(matched.set), 
                                     split = "[.]"))[c(T,F)])
    t.t <- as.numeric(unlist(strsplit(names(matched.set),
                                    split = "[.]"))[c(F,T)])
    lag <- attr(matched.set, "lag")
    t.range <- (t.t - lag):t.t
    control.ids <- unlist(matched.set)
    
    wts <- attr(matched.set[[1]], "weights")
    
    wtd.ids <- as.numeric(names(wts[wts > 0]))
    
  }
  
  if (show.set.only)
  {
    data <- data[data[, "unit.id"] %in% c(t.id, wtd.ids), ]
  }
  
  
  if (!gradient.weights)
  {
    if (!dense.plot)
    {
      p <- ggplot(data, aes(y = unit.id, x = time.id)) +
        geom_tile(aes(fill = treatment),colour = "white")
      
    } else
    {
      p <- ggplot(data, aes(y = unit.id, x= time.id)) +
        geom_raster(aes(fill = treatment), hjust = 0, vjust = .5)
      
    }
  }

  
  
  if (is.null(matched.set))
  {
    clrs <- NULL
    
    p <- p + scale_fill_gradient(low = color.of.untreated,
                        high = color.of.treated, guide = "legend",
                        breaks = c(0,1), labels = legend.labels) +
      theme_bw() +
      labs(title = title, x = xlab, y = ylab, fill = "") +
      theme(axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), panel.border = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
            axis.text.y = element_text(size = y.size, angle = y.angle),
            plot.title = element_text(hjust = 0.5))
  }
  
  
  if (!is.null(matched.set) & 
      length(matched.set) == 1 & 
      inherits(matched.set, "matched.set"))
  {
    
    
    if(gradient.weights)
    {
      
      wts <- attr(matched.set[[1]], "weights")
      max.wt <- max(wts)
      min.wt <- min(wts[wts > 0])
      #min.wt <- min(wts)
      low.wt <- min.wt * .5
      wts.z <- wts[wts > 0]
      wt <- wts.z[match(data$unit.id, names(wts.z))]
      wt[is.na(wt) | wt == 0] <- low.wt
      data$alphaweight <- wt
      
      data[!(data$time %in% t.range & data$unit.id %in% wtd.ids), "alphaweight"] <- low.wt
      data[data$unit.id == t.id & data$time %in% t.range, "alphaweight"] <- max.wt
      
      data$alphaweight <- (data$alphaweight - min(data$alphaweight)) / (max(data$alphaweight) - min(data$alphaweight))
      data$alphaweight[data$alphaweight == min(data$alphaweight)] <- min(data$alphaweight[data$alphaweight > 0]) * .5
      
      
      
    } else #if(!gradient.weights) 
    {
      low.wt <- 1 * .5
      max.wt <- 1
      data$alphaweight <- low.wt
      wtd.controls <- names(wts[wts > 0])
      data[data$unit.id %in% c(t.id, wtd.controls),
                       "alphaweight"] <- max.wt
      data[!data$time %in% t.range,
                       "alphaweight"] <- low.wt
        
    }
    
    if (!dense.plot)
    {
      p <- ggplot(data, aes(y = unit.id, x = time.id)) +
        geom_tile(aes(fill = treatment, alpha = alphaweight),colour = "white")
      
    } else
    {
      p <- ggplot(data, aes(y = unit.id, x= time.id)) +
        geom_raster(aes(fill = treatment, alpha = alphaweight), hjust = 0, vjust = .5)
      
    }
    
    
    p <- p + scale_fill_gradient(low = color.of.untreated,
                                 high = color.of.treated, guide = "legend",
                                 breaks = c(0,1), labels = legend.labels) +
      theme_bw() +
      labs(title = title, x = xlab, y = ylab, fill = "") +
      theme(axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), panel.border = element_blank(),
            legend.position = legend.position,
            panel.background = element_blank(),
            axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
            axis.text.y = element_text(size = y.size, angle = y.angle),
            plot.title = element_text(hjust = 0.5))
    
  }
    
  if(hide.x.tick.label)
  {
    p <- p + theme(axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank())
  }
  if(hide.y.tick.label)
  {
    p <- p + theme(axis.text.y = element_blank(), 
                     axis.ticks.y = element_blank())
  }
  
  if (hide.x.tick.label ||
      hide.y.tick.label)
  {
    pjp <- p
  } else {
    pjp <- p + scale_y_discrete(expand = c(0, 0), labels = unique(as.character(data$unit.id))) + 
      ggtitle(title) + xlab(xlab) + ylab(ylab)
  }

  return(pjp) # return the plot
}
