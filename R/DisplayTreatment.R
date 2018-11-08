#' DisplayTreatment
#' 
#' \code{DisplayTreatment} visualizes the treatment distribution
#' across unit and time in a panel dataset
#'
#' @param unit.id A numeric vector of unit identifiers
#' @param time.id A numeric vector of time identifiers
#' @param treatment Name of the treatment variable in character class
#' @param data The data frame in data.frame class
#' @param color.of.treated Color of the treated observations in
#' character. Default is red, but shade will depend on whether or not \code{group_on} is provided
#' @param color.of.untreated Color of the untreated observations in
#' character. Default is blue, but shade will depend on whether or not \code{group_on} is provided
#' @param title Title of the plot in character
#' set as in ggplot2
#' @param xlab Character label of the x-axis
#' @param ylab Character label of the y-axis
#' @param x.size Numeric size of the text for xlab. Default is 10. Assign x.size = NULL to use built in ggplot2 method of determining label size. 
#' When the length of the time period is long, consider setting to NULL and adjusting size and ratio of the plot.
#' @param y.size Numeric size of the text for ylab. Default is 5. Assign y.size = NULL to use built in ggplot2 method of determining label size. 
#' When the number of units is large, consider setting to NULL and adjusting size and ratio of the plot.
#' @param x.angle Angle (in degrees) of the tick labels for x-axis
#' @param y.angle Angle (in degrees) of the tick labels for y-axis
#' @param legend.position Position of the legend with the same choice
#' @param legend.labels Character vector of length two describing the
#' labels of the legend to be shown in the plot
#' @param group_on Character name of column with a categorical variable that is not time dependent. 
#' If provided, units with shared values will be grouped together and each group will be highlighted on the resulting plot.
#' @param sort_by Character name of a column containing data that will be used to help determine the order in which units are displayed on the plot. 
#' The average value for this variable will be calculated per unit. Units will then be displayed in ascending or descending order, according to that calculated value. 
#' Default value is the amount of treatment that units receive.
#' @param decreasing Logical. Determines if display order-- according to either the column specified in the \code{sort_by} argument or amount of treatment that units receive-- should be increasing or decreasing
#' @return \code{DisplayTreatment} returns a treatment variation plot,
#' which visualizes the variation of treatment across unit and time.
#' 
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples 
#' \dontrun{
#' DisplayTreatment(unit.id = "wbcode2",
#'                  time.id = "year", legend.position = "none",
#'                  xlab = "year", ylab = "Country Code",
#'                  treatment = "dem", data = dem)
#' }
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
                             # color.for.missing = "grey",
                             legend.labels = c("not treated", "treated"),
                             group_on = NULL,
                             sort_by = NULL,
                             decreasing = FALSE,
                             matched.set = NULL,
                             show.set.only = FALSE)
    
{
  # load the dataframe that the user specifies
  #TODO: fix layout here

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
    data <- data[sapply(data[, unit.id], .in.set), ]
  }
  x.size = NULL
  y.size = NULL
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
    # "#adc5ff" #control, not in set -- this will correspond to 0
    # "#ffadad" #treated, not in set -- corresponds to 1
    # "red" #treated, in set -- 3
    # "blue" #control, in set -- 2
    # #data$colref <- NA
    
    
    lag <- attr(matched.set, "lag")
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
    data$colref <- mapply(FUN = .set.colref, id_ = data$old.index, treatment = data$treatment, time_ = data$time.id)
    clrs <- sapply(unique(data$old.index), FUN = .in.set)
    clrs <- ifelse(clrs, "blue", "#eaeaea")
    clrs[which(unique(data$old.index) == id)] <- "red"
    title = paste0(title, "\n", "highlighted matched set for unit id: ", id, " at time t = ", t, " and lag = ", lag)
    p <- ggplot(data, aes(unit.id, time.id)) + geom_tile(aes(fill = as.factor(colref)),
                                                         colour = "white") +
      scale_fill_manual(values = c("#adc5ff", "#ffadad", "blue", "red"))+
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
  else
  {
    clrs <- NULL
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
    scale_y_continuous(limits = c(min(data$time.id)-1,max(data$time.id) + 1), expand = c(0,0)) + theme(axis.text.y = element_text(color = clrs)) +
    coord_flip()
  
  return(pjp) # return the plot
}

#assuming integer t
find.window <- function(lag,t)
{
  return((lag - t):lag)
}

