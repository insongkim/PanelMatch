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
#' @param x.size Numeric size of the text for xlab or x axis label. Default is 10. Assign x.size = NULL to use built in ggplot2 method of determining label size. 
#' When the length of the time period is long, consider setting to NULL and adjusting size and ratio of the plot.
#' @param y.size Numeric size of the text for ylab or y axis label. Default is 5. Assign y.size = NULL to use built in ggplot2 method of determining label size. 
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
#' @param hide.x.axis.label logical. If TRUE, x axis labels are not shown. Default is FALSE. 
#' @param hide.y.axis.label logical. If TRUE, y axis labels are not shown. Default is FALSE.
#' @param dense.plot logical. if TRUE, lines between tiles are removed on resulting plot. This is useful for producing more readable plots in situations where the number of units and/or time periods is very high.
#' @return \code{DisplayTreatment} returns a treatment variation plot (using ggplot2),
#' which visualizes the variation of treatment across unit and time.
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <adamrauh@mit.edu>, and Kosuke Imai <imai@harvard.edu>
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
                             x.size = 10, y.size = 5,
                             legend.position= "none",
                             x.angle = 45,
                             y.angle = NULL,
                             legend.labels = c("not treated", "treated"),
                             decreasing = FALSE,
                             matched.set = NULL,
                             show.set.only = FALSE,
                             hide.x.axis.label = FALSE,
                             hide.y.axis.label = FALSE,
                             gradient.weights = FALSE,
                             dense.plot = FALSE)
    
{
  #note: removing group_on and sort_by features for right now:
  group_on <- NULL
  sort_by <- NULL
  colref <- NULL #for some reason --as-cran checks need this
  color.manual <- NULL
  # @param sort_by Character name of a column containing data that will be used to help determine the order in which units are displayed on the plot. 
  # The average value for this variable will be calculated per unit. Units will then be displayed in ascending or descending order, according to that calculated value. 
  # Default value is the amount of treatment that units receive.
  # @param group_on Character name of column with a categorical variable that is not time dependent. 
  # If provided, units with shared values will be grouped together and each group will be highlighted on the resulting plot.
  ###############
  if(class(data) != "data.frame") stop("please convert data to data.frame class")
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  if(!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if(class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  if(show.set.only & !is.null(matched.set) & length(matched.set) == 1 & class(matched.set) == "matched.set")
  {
    info <- unlist(strsplit(names(matched.set)[1], split = ".", fixed = TRUE))
    id <- info[1]
    t <- info[2]
    .in.set <- function(id_comp)
    {
      if(id_comp == id | id_comp %in% unlist(matched.set))
      {
        return(TRUE)
      }
      else
      {
        return(FALSE)
      }
    }
    .has.weight <- function(.id)
    {
      if( .id == id | (.in.set(.id) & attr(matched.set[[1]], "weights")[as.character(.id)] > 0))
      {
        return(TRUE)
      }
      else
      {
        return(FALSE)
      }
    }
    data <- data[sapply(data[, unit.id], .in.set), ]
    if(gradient.weights)
    {
      data <- data[sapply(data[, unit.id], .has.weight), ]  
    }
    
  }
  #x.size = NULL
  #y.size = NULL
  if(!is.null(group_on))
  {
    if(!is.null(sort_by))
    {
      data <- na.omit(data[c(unit.id, time.id, treatment, group_on, sort_by)])  
      colnames(data) <- c("unit.id", "time.id", "treatment", group_on, sort_by)
    }
    else
    {
      data <- na.omit(data[c(unit.id, time.id, treatment, group_on)])  
      colnames(data) <- c("unit.id", "time.id", "treatment", group_on)  
    }
    
  }
  else
  {
    if(!is.null(sort_by))
    {
      data <- na.omit(data[c(unit.id, time.id, treatment, sort_by)])
      # rename variables to match with the object names in the loop below
      colnames(data) <- c("unit.id", "time.id", "treatment", sort_by)
    }
    else
    {
      data <- na.omit(data[c(unit.id, time.id, treatment)])
      # rename variables to match with the object names in the loop below
      colnames(data) <- c("unit.id", "time.id", "treatment")  
    }
      
  }
  
  # make unit.id a character: this is useful when the unit the user
  # passes to the function is numeric (e.g. dyad id)
  # data$unit.id <- as.character(data$unit.id)
  # Sorting units by treatment intensity -- default behavior 
  
  if(is.null(sort_by) & !decreasing)
  {
    data$trintens <- tapply(data$treatment, data$unit.id, mean, na.rm = T)[as.character(data$unit.id)]
    data <- data[order(data$trintens), ]  
  }
  else if(is.null(sort_by) & decreasing)
  {
    data$trintens <- tapply(data$treatment, data$unit.id, mean, na.rm = T)[as.character(data$unit.id)]
    data <- data[order(data$trintens, decreasing = TRUE), ] 
  }
  else if(!is.null(sort_by))
  {
    if(!decreasing)
    {
      data$sortvar <- tapply(data[, sort_by], data$unit.id, mean, na.rm = T)[as.character(data$unit.id)]
      data <- data[order(data$sortvar), ]  
    }
    else 
    {
      data$sortvar <- tapply(data[, sort_by], data$unit.id, mean, na.rm = T)[as.character(data$unit.id)]
      data <- data[order(data$sortvar, decreasing = TRUE), ]
    }
  }
  
  
  # then sort again by group for plotting, if specified by the user
  if(!is.null(group_on))
  {
    if(group_on %in% colnames(data))
    {
      if(length(unique(data[, group_on])) == 1)
      {
        group_on <- NULL
        warning("only one group possible for group provided")
      }
      else
      {
        data[, group_on] <- as.factor(data[, group_on])
        data <- data[order(data[, group_on]), ]
        #change the shades of treated/untreated colors to something more subdued to make group distinctions more clear
        color.of.untreated <- "#0b45a3"
        color.of.treated <- "#aa2f2f"  
      }
      
    }
    else
    {
      stop("group_on not valid column name")
    }
  }
  data$old.index <- data$unit.id
  data$unit.id <- match(data$unit.id, unique(data$unit.id) ) 
  data$unit.id <- factor(data$unit.id, levels=unique(as.character(data$unit.id)))
  data$time.id <- factor(x = data$time.id, levels = sort(unique(data$time.id)), ordered = T)
  if(!is.null(group_on))
  {
    lvls <- levels(data[, group_on])
    .get_y.max <- function(level)
    {
      data$unit.id[max(which(data[, group_on] == level))]
    }
    .get_y.min <- function(level)
    {
      data$unit.id[min(which(data[, group_on] == level))]
    }
    y.mins = sapply(lvls, .get_y.min)
    y.maxs = sapply(lvls, .get_y.max)
  }
  
  
  if(!is.null(matched.set) & length(matched.set) == 1 & class(matched.set) == "matched.set")
  {
    lag <- attr(matched.set, "lag")
    .in.set <- function(id_comp)
    {
      
      if(id_comp == id | id_comp %in% unlist(matched.set))
      {
        return(TRUE)
      }
      else
      {
        return(FALSE)
      }
    }
    info <- unlist(strsplit(names(matched.set)[1], split = ".", fixed = TRUE))
    id <- info[1]
    t <- info[2]
    
    
    # "#adc5ff" #control, not in set -- this will correspond to 0
    # "#ffadad" #treated, not in set -- corresponds to 1
    # "red" #treated, in set -- 3
    # "blue" #control, in set -- 2
    # #data$colref <- NA
    
    
    time <- as.numeric(t)
    .set.colref <- function(id_, treatment, time_)
    {
      if(.in.set(id_comp = id_))
      {
        if(treatment)
        {
          if(time_ %in% ((time-lag):time))
          {
            return(3)  
          }
          else
          {
            return(1)
          }
        }
        else
        {
          if(time_ %in% ((time-lag):time))
          {
            return(2)  
          }
          else
          {
            return(0)
          }
        }
      }
      else
      {
        if(treatment)
        {
          return(1)
        }
        else
        {
          return(0)
        }
      }
    }
    data[, "colref"] <- mapply(FUN = .set.colref, id_ = data$old.index, treatment = data$treatment, time_ = data$time.id)
    if(!gradient.weights)
    {
      clrs <- sapply(unique(data$old.index), FUN = .in.set)
      clrs <- ifelse(clrs, color.of.untreated, "#eaeaea")
      clrs[which(unique(data$old.index) == id)] <- color.of.treated  
    }
    
    title = paste0(title, "\n", "highlighted matched set for unit id: ", id, " at time t = ", t, " and lag = ", lag)
    if(gradient.weights)
    {
      info <- unlist(strsplit(names(matched.set)[1], split = ".", fixed = TRUE))
      id <- info[1]
      t <- info[2]
      time <- as.numeric(as.character(t))
      idx <- (data[, "time.id"] %in% ((time-lag):t)) & (data[, "old.index"] %in% unlist(matched.set)) #control units -- highlighting the periods in the lag window
      data$weights <- 0 
      min.non.zero.weight <- min(attr(matched.set[[1]], "weights")[attr(matched.set[[1]], "weights") > 0])
      
      top.weight <- max(attr(matched.set[[1]], "weights"))
      lowest.intensity <- top.weight - (top.weight * .9)
      data[idx, "weights"] <- attr(matched.set[[1]], "weights")[as.character(data[idx, "old.index"])]
      idx <- (data[, "time.id"] %in% ((time-lag):t)) & data[, "old.index"] == id
      data[idx, "weights"] <- top.weight
      data[data[, 'weights'] == 0, "weights"] <- lowest.intensity
      .has.weight <- function(.id)
      {
        if( .id == id | (.in.set(.id) & attr(matched.set[[1]], "weights")[as.character(.id)] > 0))
        {
          return(1)
        }
        else
        {
          return(0)
        }
      }
      clrs <- sapply(unique(data$old.index), FUN = .has.weight)
      clrs <- ifelse(clrs, color.of.untreated, "#eaeaea")
      clrs[which(unique(data$old.index) == id)] <- color.of.treated  
      
      setmanualcolor <- function(max.w, intensity.val, treated.val, color.of.t, color.of.unt)
      {
        if(treated.val)
        {
          return(color.factor(color.of.t, intensity.val, max.w))
        }
        else
        {
          return(color.factor(color.of.unt, intensity.val, max.w))
        }
      }
      data[, "color.manual"] <- mapply(FUN = setmanualcolor, intensity.val = data$weights, treated.val = data$treatment, 
                                  MoreArgs = list(max.w = top.weight, color.of.t = color.of.treated, color.of.unt = color.of.untreated))
     
      p <- ggplot(data, aes(unit.id, time.id)) + geom_tile(aes(fill = color.manual),
                                                           colour = "white") +
        scale_fill_identity()+
        theme_bw() +
        labs(list(title = title, x = ylab, y = xlab, fill = "")) +
        theme(axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.position = legend.position,
              panel.background = element_blank(), 
              axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
              axis.text.y = element_text(size = y.size, angle = y.angle),
              plot.title = element_text(hjust = 0.5))
      
    }
    if(!gradient.weights) 
    {
      colorvals <- c(color.factor(color.of.untreated, .3, 1), color.factor(color.of.treated, .3, 1), color.of.untreated, color.of.treated)
      if(!dense.plot)
      {
        p <- ggplot(data, aes(unit.id, time.id)) + geom_tile(aes(fill = as.factor(colref)),
                                                             colour = "white") +
          scale_fill_manual(values = colorvals)+
          theme_bw() +
          labs(list(title = title, x = ylab, y = xlab, fill = "")) +
          theme(axis.ticks.x=element_blank(),
                panel.grid.major = element_blank(), panel.border = element_blank(),
                legend.position = legend.position,
                panel.background = element_blank(), 
                axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
                axis.text.y = element_text(size = y.size, angle = y.angle),
                plot.title = element_text(hjust = 0.5))  
      } else
      {
        p <- ggplot(data, aes(unit.id, time.id)) + geom_raster(aes(fill = as.factor(colref)), hjust = 0, vjust = .5) +
          scale_fill_manual(values = colorvals)+
          theme_bw() +
          labs(list(title = title, x = ylab, y = xlab, fill = "")) +
          theme(axis.ticks.x=element_blank(),
                panel.grid.major = element_blank(), panel.border = element_blank(),
                legend.position = legend.position,
                panel.background = element_blank(), 
                axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
                axis.text.y = element_text(size = y.size, angle = y.angle),
                plot.title = element_text(hjust = 0.5))
      }
      
    }
  }
  else
  {
    clrs <- NULL
    if(!dense.plot)
    {
      p <- ggplot(data, aes(unit.id, time.id)) + geom_tile(aes(fill = treatment),
                                                           colour = "white") +
        scale_fill_gradient(low = color.of.untreated,
                            high = color.of.treated, guide = "legend",
                            breaks = c(0,1), labels = legend.labels) +
        theme_bw() +
        labs(list(title = title, x = ylab, y = xlab, fill = "")) +
        theme(axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.position = legend.position,
              panel.background = element_blank(),
              axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
              axis.text.y = element_text(size = y.size, angle = y.angle),
              plot.title = element_text(hjust = 0.5))
    } else
    {
      p <- ggplot(data, aes(unit.id, time.id)) + geom_raster(aes(fill = treatment), hjust = 0, vjust = .5) +
        scale_fill_gradient(low = color.of.untreated,
                            high = color.of.treated, guide = "legend",
                            breaks = c(0,1), labels = legend.labels) +
        theme_bw() +
        labs(list(title = title, x = ylab, y = xlab, fill = "")) +
        theme(axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.position = legend.position,
              panel.background = element_blank(),
              axis.text.x = element_text(angle=x.angle, size = x.size, vjust=0.5),
              axis.text.y = element_text(size = y.size, angle = y.angle),
              plot.title = element_text(hjust = 0.5))
    }
    
    
    
  }
  
  
  if(!is.null(group_on))
  {
    dat <- ggplot_build(p)$data[[1]] 
    
    min.vals.x <- (dat$xmin[match(y.mins, dat$x)])
    max.vals.x <- (dat$xmax[match(y.maxs, dat$x)])
    
    min.vals.y <- rep(min(dat$ymin), length(min.vals.x))
    max.vals.y <- rep(max(dat$ymax), length(max.vals.x))
    # adding in some extra space between the lines for clarity-- might be a better way to do this systematically, or just remove
    max.vals.x[1:length(max.vals.x) - 1] <- max.vals.x[1:length(max.vals.x) - 1] - .02
    min.vals.x[2:length(min.vals.x)] <- min.vals.x[2:length(min.vals.x)] + .02
    pj <-  p + annotate("rect", xmin = min.vals.x, xmax = max.vals.x, ymin = min.vals.y, ymax = max.vals.y, alpha = 0, color = factor(lvls), size = 1.25)  
  }
  else
  {
    pj <- p
  }
  pjp <- pj + scale_x_discrete(expand = c(0, 0), labels = unique(as.character(data$old.index))) +
    theme(axis.text.y = element_text(color = clrs)) + coord_flip() + ggtitle(title) + xlab(ylab) + ylab(xlab) #flip them because of how the plot is generated, should probably fix this
  if(hide.x.axis.label & !hide.y.axis.label)
  {
    return(pjp + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  }
  if(hide.y.axis.label & !hide.x.axis.label)
  {
    return(pjp + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()))
  }
  if(hide.x.axis.label & hide.y.axis.label)
  {
    return(pjp + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
  }
  return(pjp) # return the plot
}

#assuming integer t
find.window <- function(lag,t)
{
  return((lag - t):lag)
}

#helper function to set intensity of colors using the weights
color.factor<-function(color, value, max){
  l.color<-length(color)
  t.rgb=rep(grDevices::col2rgb(color),length(value))
  t.rgb[is.na(t.rgb)]<-255
  dim(t.rgb)<-c(3,l.color*length(value))
  value[is.na(value)]<-0
  t.rgb<-255-((255-t.rgb)*rep(value, each=3)/max)
  return(grDevices::rgb(red=t.rgb[1,], green=t.rgb[2,], blue=t.rgb[3,], maxColorValue=255))
  
}


