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
#' A method for viewing summary data about the sizes of matched sets and metadata about how they were created. This method
#' accepts all standard \code{summary} arguments. If the quantity of interest is ate, then a summary will be provided for the matched sets associated with the att and the atc. 
#'
#' @param object a \code{PanelMatch} object
#' @param ... Optional additional arguments to be passed to the \code{summary} function
#' @param verbose Logical value specifying whether or not a longer, more verbose summary should be calculated and returned. Default is \code{FALSE}.
#'
#' @return A list of lists containing a summary of the matched sets associated with the specified qoi. Each sublist object will either have 5 or 1 element(s), depending on whether or not \code{verbose} is set to \code{TRUE} or not.
#' \item{overview}{A \code{data.frame} object containing information about the treated units (unit id, time of treatment), and the number of matched control units with weights zero and above.}
#' \item{set.size.summary}{a \code{summary} object summarizing the minimum, maximum, and IQR of matched set sizes}
#' \item{number.of.treated.units}{The number of unit, time pairs that are considered to be "treated" units}
#' \item{num.units.empty.set}{The number of units treated at a particular time that were not able to be matched to any control units}
#' \item{lag}{The size of the lag window used for matching on treatment history. This affects which treated and control units are matched.}
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
summary.PanelMatch <- function(object, ..., verbose = FALSE)
{
  qoi.in <- attr(object, "qoi")
  ll <- list()
  if (qoi.in %in% c("att", "atc", "art"))
  {
    msets <- object[[qoi.in]]
    ll[[qoi.in]] <- summary(msets, ..., verbose = verbose)
    
  } else if (qoi.in == "ate")
  {
    mset1 <- object[["att"]]
    mset2 <- object[["atc"]]
    ll[["att"]] <- summary(mset1, ..., verbose = verbose)
    ll[["atc"]] <- summary(mset2, ..., verbose = verbose)
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
#' @param verbose logical indicating whether or not underlying data should be printed in expanded/raw list form.
#' The verbose form is not recommended unless the data set is small. Default is FALSE
#' @param ... additional arguments to be passed to \code{print}
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
print.PanelMatch <- function(x, ..., verbose)
{
  qoi.in <- attr(x, "qoi")
  unit.id <- attr(x, "unit.id") 
  time.id <- attr(x, "time.id") 
  outcome <- attr(x, "outcome") 
  treatment.var <- attr(x, "treatment") 
  cat(paste0(paste0("PanelMatch Object\n"), 
    paste0("Unit id: ", unit.id, "\n"),
    paste0("Time id: ", time.id, "\n"),
    paste0("Outcome variable: ", outcome, "\n"),
    paste0("Treatment variable: ", treatment.var, "\n")))
  if (qoi.in == "ate")
  {
    qoi.in <- c("att", "atc")
  }
  
  for (q.in in qoi.in) {
    mset <- x[[q.in]]
    cat(paste0("QOI: ", q.in), "\n")
    print(mset, ..., verbose)
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
