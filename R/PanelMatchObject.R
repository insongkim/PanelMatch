#' Extract matched.set objects from PanelMatch results
#' @param pm.object \code{PanelMatch} object
#' @param qoi character, specifying the qoi. Valid inputs include "att", "atc", "art", and NULL. If NULL, function extracts att, art, or atc results if possible. Otherwise, throws an error if ate is specified.
#'
#' @export
extract <- function(pm.object, qoi) {
  UseMethod("extract", pm.object)
}

#' Extract matched.set objects from PanelMatch results
#'
#' @param pm.object \code{PanelMatch} obect
#' @param qoi character, specifying the qoi. Valid inputs include "att", "atc", "art", and NULL. If NULL, function extracts att, art, or atc results if possible. Otherwise, throws an error if ate is specified.
#' @return a \code{matched.set} object
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
#' # create subset of data for simplicity
#' PM.results <- PanelMatch(panel.data = dem.sub.panel,
#'                          lag = 4, 
#'                          refinement.method = "mahalanobis",
#'                          match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, forbid.treatment.reversal = FALSE)
#' extract(PM.results, qoi = "att")
#' extract(PM.results) # valid since att is specified
#' @method extract PanelMatch
#' @export
extract.PanelMatch <- function(pm.object, 
                           qoi = NULL)
{
  
  
  if(is.null(qoi) && identical(attr(pm.object, "qoi"), "ate"))
  {
    stop("Please specify qoi = NULL (if the qoi is not ate), att, art, or atc. ATE is not a valid specification for extraction.")
  }
  
  if (!is.null(qoi) && identical(qoi, "ate"))
  {
    stop("ATE is not a valid specification.")
  }
  msets <- handle_pm_qoi(pm.object, qoi)
  return(msets)
  
}

#' Summarize information about a PanelMatch object and the matched sets contained within them.
#'
#'
#' A method for viewing summary data about the sizes of matched sets, the number of treated units, and the number of empty matched sets. If the quantity of interest is ate, then a summary will be provided for the matched sets associated with the att and the atc. 
#'
#' @param object a \code{PanelMatch} object
#' @param ... Not used
#' 
#' @return A list of data frame(s) containing information about matched sets associated with the specified qoi. If the qoi is "att", "art", or "atc", then the returned list contains one data frame and the element is named for the specified qoi. If the qoi is "ate", then a list of two elements is returned, with one data frame corresponding to the "att" and the other to the "atc". The data frame contains summary information about the sizes of matched sets, along with information about the number of treated observations and the number of empty sets. Specifically, it contains the minimum, 1st quartile, median, mean, 3rd quartile, and maximum matched set size. It also contains the number of treated units total and the number of empty matched sets.
#' 
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
#' PM.results <- PanelMatch(panel.data = dem.sub.panel,
#'                          lag = 4, 
#'                          refinement.method = "mahalanobis",
#'                          match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, forbid.treatment.reversal = FALSE)
#' summary(PM.results)
#'
#' @method summary PanelMatch
#' @export
summary.PanelMatch <- function(object, ...)
{
  qoi.in <- attr(object, "qoi")
  ll <- list()
  if (qoi.in %in% c("att", "atc", "art"))
  {
    msets <- object[[qoi.in]]
    
    dlt <- summary(msets, ..., verbose = TRUE)
    
    named_stats <- c(
      dlt$set.size.summary,
      "Number of treated units" = dlt$number.of.treated.units,
      "Number of empty matched sets" = dlt$num.units.empty.set
    )
    
    ll[[qoi.in]] <- data.frame(
      quantity = names(named_stats),
      value = as.numeric(named_stats),
      row.names = NULL
    )
    
    
  } else if (qoi.in == "ate")
  {
    mset1 <- object[["att"]]
    mset2 <- object[["atc"]]
    
    dlt <- summary(mset1, ..., verbose = TRUE)
    
    named_stats <- c(
      dlt$set.size.summary,
      "Number of treated units" = dlt$number.of.treated.units,
      "Number of empty matched sets" = dlt$num.units.empty.set
    )
    
    ll[["att"]] <- data.frame(
      quantity = names(named_stats),
      value = as.numeric(named_stats),
      row.names = NULL
    )
    
    dlt <- summary(mset2, ..., verbose = TRUE)
    
    named_stats <- c(
      dlt$set.size.summary,
      "Number of treated units" = dlt$number.of.treated.units,
      "Number of empty matched sets" = dlt$num.units.empty.set
    )
    
    ll[["atc"]] <- data.frame(
      quantity = names(named_stats),
      value = as.numeric(named_stats),
      row.names = NULL
    )
  } else {
    stop("QOI is misspecified.")
  }
  return(ll)  
}

#' Plot the distribution of the sizes of matched sets.
#'
#'
#' A plot method for creating a histogram of the distribution of the sizes of matched sets.
#' This method accepts all standard optional \code{hist} arguments via the \code{...} argument.
#' By default, empty matched sets (treated units that could not be
#' matched with any control units) are noted as a vertical bar at x = 0 and not included in the
#' regular histogram. See the \code{include.empty.sets} argument for more information about this. If the quantity of interest is ATE, a plot will be returned for the matched sets associated with the att and the atc. 
#'
#' @param x a \code{PanelMatch} object
#' @param ... optional arguments to be passed to \code{hist()}
#' @param border default is NA. This is the same argument as the standard argument for \code{hist()}
#' @param col default is "grey". This is the same argument as the standard argument for \code{hist()}
#' @param ylab default is "Frequency of Size". This is the same argument as the standard argument for \code{hist()}
#' @param xlab default is "Matched Set Size". This is the same argument as the standard argument for \code{hist()}
#' @param lwd default is NULL. This is the same argument as the standard argument for \code{hist()}
#' @param main default is "Distribution of Matched Set Sizes". This is the same argument as the standard argument for \code{hist}
#' @param freq default is TRUE. See \code{freq} argument in \code{hist()} function for more.
#' @param include.empty.sets logical value indicating whether or not empty sets should be included in the histogram. default is FALSE. If FALSE, then empty sets will be noted as a separate vertical bar at x = 0. If TRUE, empty sets will be included as normal sets.
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
#' PM.results <- PanelMatch(panel.data = dem.sub.panel,
#'                          lag = 4, 
#'                          refinement.method = "mahalanobis",
#'                          match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, forbid.treatment.reversal = FALSE)
#' plot(PM.results)
#' plot(PM.results, include.empty.sets = TRUE)
#'
#' @method plot PanelMatch
#' @export
plot.PanelMatch <- function(x, ..., border = NA, col = "grey", ylab = "Frequency of Size",
                             xlab ="Matched Set Size" , lwd = NULL,
                             main = "Distribution of Matched Set Sizes",
                             freq = TRUE, include.empty.sets = FALSE)
{
  qoi.in <- attr(x, "qoi")
  
  if (qoi.in == "ate")
  {
    qoi.in <- c("att", "atc")
  }
  
  for (q.in in qoi.in) {
    mset <- x[[q.in]]
    plot_matched_set(mset, border = border, col = col, ylab = ylab,
         xlab = xlab, lwd = lwd, main = main, 
         freq = freq, include.empty.sets, ...)
  }
  
}

#' Print PanelMatch objects.
#'
#' @param x a \code{PanelMatch} object
#' @param ... additional arguments to be passed to \code{print.matched.set}
#' @param verbose logical indicating whether or not underlying data should be printed in expanded/raw list form.
#' The verbose form is not recommended unless the data set is small. Default is FALSE
#' @param n Integer. Number of matched sets to display information about as a preview. Default is 5.
#' @param show.all Logical. By default (`show.all = FALSE`), the print method only shows a small preview of the sizes of matched sets. When set to TRUE, a full summary description of matched set sizes is shown. 
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem, 'wbcode2', 'year', 'dem', 'y')
#' PM.results <- PanelMatch(panel.data = dem.sub.panel,
#'                          lag = 4, 
#'                          refinement.method = "mahalanobis",
#'                          match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, forbid.treatment.reversal = FALSE)
#' print(PM.results)
#' 
#'
#' @method print PanelMatch
#' @export
print.PanelMatch <- function(x, ..., verbose = FALSE, n = 5, show.all = FALSE) {
  # Extract required attributes
  qoi.in <- attr(x, "qoi")
  unit.id <- attr(x, "unit.id") 
  time.id <- attr(x, "time.id") 
  outcome <- attr(x, "outcome") 
  treatment.var <- attr(x, "treatment") 
  
  # Extract optional attributes if present
  lag <- attr(x, "lag")
  lead <- attr(x, "lead")
  refinement.method <- attr(x, "refinement.method")
  
  # Header
  cat("PanelMatch Object Summary\n")
  cat(strrep("-", 40), "\n")
  cat(sprintf("Unit ID             : %s\n", unit.id))
  cat(sprintf("Time ID             : %s\n", time.id))
  cat(sprintf("Outcome Variable    : %s\n", outcome))
  cat(sprintf("Treatment Variable  : %s\n", treatment.var))
  
  # Print optional attributes if they exist
  if (!is.null(lag)) {
    cat(sprintf("Lag                 : %s\n", lag))
  }
  if (!is.null(lead)) {
    cat(sprintf("Max Lead            : %s\n", as.character(max(lead))))
  }
  if (!is.null(refinement.method)) {
    cat(sprintf("Refinement Method   : %s\n", refinement.method))
  }
  
  cat(strrep("-", 40), "\n")
  
  # Handle QOI expansion if "ate"
  if (qoi.in == "ate") {
    qoi.in <- c("att", "atc")
  }
  
  # Print each matched set
  for (q.in in qoi.in) {
    mset <- x[[q.in]]
    nsets <- length(mset)
    cat(sprintf("QOI: %s\n", toupper(q.in)))
    #cat(sprintf("Number of treated observations: %d\n", nsets))
    #cat("Summary of matched set sizes:\n")
    #print(summary(sapply(mset, length)))
    cat(strrep("-", 40), "\n")
    
    
    print(mset, ..., verbose = verbose, n = n, show.all = show.all)
    cat("\n")
  }
}


#' Helper function for plotting the distribution of matched set sizes
#'
#'
#'
#' @param x a \code{matched.set} object
#' @param ... optional arguments to be passed to \code{hist()}
#' @param border default is NA. This is the same argument as the standard argument for \code{hist()}
#' @param col default is "grey". This is the same argument as the standard argument for \code{hist()}
#' @param ylab default is "Frequency of Size". This is the same argument as the standard argument for \code{hist()}
#' @param xlab default is "Matched Set Size". This is the same argument as the standard argument for \code{hist()}
#' @param lwd default is NULL. This is the same argument as the standard argument for \code{hist()}
#' @param main default is "Distribution of Matched Set Sizes". This is the same argument as the standard argument for \code{hist}
#' @param freq default is TRUE. See \code{freq} argument in \code{hist()} function for more.
#' @param include.empty.sets logical value indicating whether or not empty sets should be included in the histogram. default is FALSE. If FALSE, then empty sets will be noted as a separate vertical bar at x = 0. If TRUE, empty sets will be included as normal sets.
#' @keywords internal
plot_matched_set <- function(x, border = NA, col = "grey", ylab = "Frequency of Size",
                             xlab ="Matched Set Size" , lwd = NULL,
                             main = "Distribution of Matched Set Sizes",
                             freq = TRUE, include.empty.sets = FALSE, ...)
{
  set <- x
  lvec <- sapply(set, length)
  
  if(include.empty.sets)
  {
    graphics::hist(x = lvec, freq = freq, 
                   border = border, col = col, 
                   ylab = ylab, xlab = xlab, 
                   main = main, ...)
  }
  else
  {
    lvec.nonempty <- lvec[lvec > 0]
    
    if(sum(lvec == 0) > 0)
    {
      num.empties <- as.character(sum(lvec == 0))
      graphics::hist(x = lvec.nonempty, freq = freq, 
                     border = border, col = col, ylab = ylab,
                     xlab = xlab, main = main, ...)
      graphics::lines(x = c(0,0),
                      y = c(0, num.empties),
                      lwd = 4,
                      col = "#ffc6c4", ...)
    }
    else
    {
      graphics::hist(x = lvec.nonempty, 
                     freq = freq, 
                     border = border, 
                     col = col, ylab = ylab,
                     xlab = xlab, main = main, ...)
    }
  }
  
}


#' Return the refinement formula used in a PanelMatch specification
#'
#' @param x A PanelMatch Object 
#' @param ... not used
#'
#' @return One sided formula object containing the variables/specification used in refinement. This corresponds to what was provided to the \code{covs.formula} argument.
#' @export
formula.PanelMatch <- function(x, ...) {
  if (!inherits(x, "PanelMatch")) {
    stop("This method is only for PanelMatch objects.")
  }
  
  covs_formula <- attr(x, "covs.formula")
  
  if (is.null(covs_formula)) {
    return(NULL)
  }
  
  return(covs_formula)
}