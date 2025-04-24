#' Pre-process and balance panel data
#'
#' @param panel.data A \code{data.frame} object containing time series cross sectional data. 
#' Time data should be sequential integers that increase by 1. Unit identifiers must be integers. Treatment data must be binary. If time data is non-integer, the package will attempt to sensibly convert it by converting the data to factor, then to integer. If a conversion is performed, a mapping will be returned as an attribute called "time.data.map"
#' @param unit.id A character string indicating the name of unit identifier in the data. This data must be integer.
#' @param time.id A character string indicating the name of the time variable in the data. 
#' @param treatment A character string indicating the name of the treatment variable.
#' The treatment must be a binary indicator variable (integer with 0 for the control group and 1 for the treatment group).
#' @param outcome A character string identifying the outcome variable
#' @return \code{PanelData()} returns an object of class \code{PanelData}. This takes the form of a \code{data.frame} object with the following properties and attributes. First, the data has been balanced and sorted. These properties are noted in the "is.balanced" and "is.sorted" attributes, respectively. So, each unit appears the same number of times in the resulting \code{PanelData} object, with NAs filling out missing data. Second, the data has been sorted to appear in order for each unit. Next, the \code{PanelData} object has the following attributes: "unit.id", "time.id", "treatment", and "outcome" reflecting the variables provided in the specification. If the function attempts to automatically convert time data to be consecutive integers, the mapping between the original time data and the "new" converted time data is provided as a \code{data.frame} object and stored as the "time.data.map" attribute. 
#'
#' @export
#' @examples
#' d <- PanelData(panel.data = dem, 
#'                unit.id = "wbcode2", 
#'                time.id = "year", 
#'                treatment = "dem", 
#'                outcome = "y")
PanelData <- function(panel.data,
                      unit.id,
                      time.id,
                      treatment,
                      outcome)
{
  if (!inherits(panel.data, "data.frame")) stop("please convert data to data.frame class")
  if (any(duplicated(panel.data[, c(unit.id, time.id)]))) stop("Time, unit combinations should uniquely identify rows. Please remove duplicates")
  if (!inherits(panel.data[, unit.id], "integer") && !inherits(panel.data[, unit.id], "numeric")) stop("please convert unit id column to integer or numeric")
  if ( !all(c(time.id, unit.id, treatment, outcome)  %in% colnames(panel.data)) ) stop("time id, unit id, outcome, or treatment column name invalid. Please ensure all variables are specified correctly.")
  
  if (!inherits(panel.data[, treatment], "integer") && 
      !inherits(panel.data[, treatment], "numeric")) 
  {
    stop("Please convert treatment data to integer or numeric")
  }
   
  if (!inherits(panel.data[, outcome], "integer") && 
      !inherits(panel.data[, outcome], "numeric")) 
  {
    stop("Please convert outcome data to integer or numeric")
  } 
  
  if (any(is.na(panel.data[, unit.id]))) stop("Cannot have NA unit ids")
  
  ## balance the panel 
  if (any(table(panel.data[, unit.id]) != max(table(panel.data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(panel.data), 
                                 formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table::data.table(testmat), 
                          id = unit.id, 
                          variable = time.id, 
                          value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    panel.data <- merge(data.table::data.table(d), 
                  data.table::data.table(panel.data), 
                  all.x = TRUE, 
                  by = c(unit.id, time.id))
    panel.data <- as.data.frame(panel.data)
    
  }
  
  ordered.data <- panel.data[order(panel.data[,unit.id], panel.data[,time.id]), ]
  ordered.data <- check_time_data(ordered.data, time.id)
  mapping.if.exists <- attr(ordered.data, "time.data.map")
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  
  class(ordered.data) <- c("PanelData", "data.frame")
  
  attr(ordered.data, "unit.id") <- unit.id
  attr(ordered.data, "time.id") <- time.id
  attr(ordered.data, "treatment") <- treatment
  attr(ordered.data, "outcome") <- outcome
  
  attr(ordered.data, "is.sorted") <- TRUE
  attr(ordered.data,  "is.balanced") <- TRUE
  attr(ordered.data, "time.data.map") <- mapping.if.exists
  
  return(ordered.data)
}

#' Create basic plots of PanelData objects
#' @param x \code{PanelData} object
#' @param ... Not used
#' @param plotting.variable character string specifying which variable to plot in the resulting figure. The values of this variable will be used to fill in cells on the resulting heatmap. Defaults to whatever is specified as the treatment variable. 
#' @return Returns a ggplot2 object created by \code{geom_tile()}. The basic figure shows units along the y-axis and time along the x-axis. The figure takes the form of a heatmap. The value of the plotting.variable argument is used to fill in the color of the cells. 
#'
#' @export
#' @examples
#' dem$rdata <- rnorm(nrow(dem))
#' d <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' plot(d)
#' plot(d, plotting.variable = "rdata")
plot.PanelData <- function(x, ..., plotting.variable = NA)
{
  if (!inherits(x, "PanelData")) {
    stop("Please provide a PanelData object")
  }
  if (is.na(plotting.variable)) {
      plotting.variable <- attr(x, "treatment")
  }
  time_ID <- attr(x, "time.id")
  unit_ID <- attr(x, "unit.id")
  # Convert unit_ID and time_ID to factors
  x[[time_ID]] <- as.factor(x[[time_ID]])
  x[[unit_ID]] <- as.factor(x[[unit_ID]])
  
  # Create the heat map
  heatmap_plot <- ggplot(x, aes(x = !!sym(time_ID), 
                                y = !!sym(unit_ID), 
                                fill = !!sym(plotting.variable))) +
    geom_tile() +
    scale_fill_gradient() + 
    theme_minimal() +
    scale_x_discrete() +
    scale_y_discrete() +
    labs(title = paste("Heatmap of", plotting.variable),
         x = time_ID,
         y = unit_ID,
         fill = plotting.variable)
  return(heatmap_plot)
  
}

#' Print PanelData objects and basic metadata
#' @param x \code{PanelData} object
#' @param ... additional arguments to be passed to \code{print.data.frame()}
#' @param n Integer. Number of rows to print by default for previewing data. Default is 5.
#' @param verbose Logical. Print the entire data frame, rather than just a preview. Default is FALSE.
#'
#' @return Returns nothing but prints \code{PanelData} object. This is a data frame that has been balanced, sorted, and tagged with important metadata to facilitate the use of other functions.
#' @export
#' @examples
#' d <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' print(d)
#' 
print.PanelData <- function(x, ..., n = 5, verbose = FALSE)
{
  if (!inherits(x, "PanelData")) {
    stop("Please provide a PanelData object")
  }
  
  unit.id <- attr(x, "unit.id")
  time.id <- attr(x, "time.id")
  treatment <- attr(x, "treatment")
  outcome <- attr(x, "outcome")
  
  print.str <- paste0(
    "Unit identifier: ", unit.id, "\n",
    "Time identifier: ", time.id, "\n",
    "Treatment variable: ", treatment, "\n",
    "Outcome variable: ", outcome, "\n"
  )
  cat(print.str)
  
  total.rows <- nrow(x)
  
  if (verbose) {
    print.data.frame(x, ...)
  } else {
    string.out <- paste0(total.rows, " x ", ncol(x))
    cat(paste0("Dimensions: ", string.out, "\n"))
    print.data.frame(utils::head(x, n = n))
    
    if (n < total.rows) {
      cat(paste0("... [", total.rows - n, " more row(s) not printed]\n"))
    }
  }
}

#' Summarize information about variable names, and unit, time, and treatment data in a PanelData object
#' @param object \code{PanelData} object
#' @param ... Not used
#' @return Returns a \code{data.frame} object, with columns "quantity" and "value." Within the data frame the following information is returned: The name of the unit id variable, the name of the time id variable, the name of the treatment variable, the name of the outcome variable, the number of unique units found in the data, the number of unique time periods found in the data and the percentage of treated periods that are missing treatment data. 
#' @export
#' @examples
#' d <- PanelData(dem, "wbcode2", "year", "dem", "y")
#' summary(d)
#' 
summary.PanelData <- function(object, ...)
{
  attr(object, "unit.id") -> unit.id
  attr(object, "time.id") -> time.id
  attr(object, "treatment") -> treatment
  attr(object, "outcome") -> outcome
  num.units <- length(unique(object[, unit.id]))
  num.periods <- length(unique(object[, time.id]))
  pct.missing <- round(mean(is.na(object[, treatment])) * 100, 4)
  quantities <- c("Unit ID", "Time ID", "Treatment", "Outcome",
    "# Unique Units", "# Unique Time Periods",
    "% of Treatment Periods Missing Data")
  values <- c(unit.id, time.id, treatment, outcome, 
    num.units, num.periods, pct.missing)
  dt <- data.frame(quantity = quantities,
             value = values)
  return(dt)
}