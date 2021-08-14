#' matched_set
#'
#' \code{matched_set} is a constructor for the \code{matched.set} class.
#'
#'
#' Users should never need to use this function by itself. See below for more about \code{matched.set} objects.
#'
#' @param matchedsets a list of treated units and matched control units. Each element in the list should be a vector of control unit ids.
#' @param id A vector containing the ids of treated units
#' @param t A vector containing the times of treatment for treated units.
#' @param L integer specifying the length of the lag window used in matching
#' @param t.var string specifying the time variable
#' @param id.var string specifying the unit id variable
#' @param treatment.var string specifying the treatment variable.
#'
#' The constructor function returns a \code{matched.set} object.
#' \code{matched.set} objects are a modified lists. Each element in the list is a vector of ids
#' corresponding to the control unit ids in a matched set.
#' Additionally, these vectors might have additional attributes -- "weights". These correspond to the
#' weights assigned to each control unit,
#' as determined by the specified refinement method.
#' Each element in the list also has a name, which corresponds to the unit id of the treated unit and time of treatment,
#' concatenated together and separated by a period. \code{matched.set} objects also have a number of
#' methods defined: \code{summary}, \code{plot}, and \code{`[`}. \code{matched.set} objects can be modified manually
#' as long as these conventions (and conventions about other attributes) are maintained. It is important to note that \code{matched.set} objects
#' are distinct from \code{PanelMatch} objects. \code{matched.set} objects are often contained within \code{PanelMatch} objects.
#' @return \code{matched.set} objects have additional attributes. These reflect the specified parameters when using the \code{PanelMatch} function:
#' \item{lag}{an integer value indicating the length of treatment history to be used for matching. Treated and control units are matched based on whether or not they have exactly matching treatment histories in the lag window.}
#' \item{t.var}{time variable name, represented as a character/string}
#' \item{id.var}{unit id variable name, represented as a character/string}
#' \item{treatment.var}{treatment variable name, represented as a character/string}
#' \item{class}{class of the object: should always be "matched.set"}
#' \item{refinement.method}{method used to refine and/or weight the control units in each set.}
#' \item{covs.formula}{One sided formula indicating which variables should be used for matching and refinement}
#' \item{match.missing}{Logical variable indicating whether or not units should be matched on the patterns of missingness in their treatment histories}
#' \item{max.match.size}{Maximum size of the matched sets after refinement. This argument only affects results when using a matching method}
#' @author Adam Rauh <adamrauh@mit.edu>, In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <imai@harvard.edu>
#' @export
matched_set <- function(matchedsets, id, t, L, t.var, id.var, treatment.var)
{
  if(length(matchedsets) != length(id) | length(matchedsets) != length(t) | length(id) != length(t))
  {
    stop("Number of matched sets, length of t, length of id specifications do not match up")
  }
  names(matchedsets) <- paste0(id, ".", t)
  class(matchedsets) <- c("matched.set", "list")
  attr(matchedsets, 'refinement.method') <- NULL
  attr(matchedsets, "lag") <- L
  attr(matchedsets, "t.var") <- t.var
  attr(matchedsets, "id.var" ) <- id.var
  attr(matchedsets, "treatment.var") <- treatment.var
  return(matchedsets)
}

#' Summarize information about a \code{matched.set} object and the matched sets contained within them.
#'
#'
#' A method for viewing summary data about the sizes of matched sets and metadata about how they were created. This method
#' accepts all standard \code{summary} arguments.
#'
#' @param object a \code{matched.set} object
#' @param ... Optional additional arguments to be passed to the \code{summary} function
#' @param verbose Logical value specifying whether or not a longer, more verbose summary should be calculated and returned. Default is
#' \code{TRUE}.
#'
#' @return list object with either 5 or 1 element(s), depending on whether or not \code{verbose} is set to \code{TRUE} or not.
#' \item{overview}{A \code{data.frame} object containing information about the treated units (unit id, time of treatment), and the number of matched control units with weights zero and above.}
#' \item{set.size.summary}{a \code{summary} object summarizing the minimum, maximum, and IQR of matched set sizes}
#' \item{number.of.treated.units}{The number of unit, time pairs that are considered to be "treated" units}
#' \item{num.units.empty.set}{The number of units treated at a particular time that were not able to be matched to any control units}
#' \item{lag}{The size of the lag window used for matching on treatment history. This affects which treated and control units are matched.}
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' summary(PM.results$att)
#'
#'
#'
#' @method summary matched.set
#' @export
summary.matched.set <- function(object, ..., verbose = TRUE)
{
  set <- object
  Lengthcol <- sapply(set, length)

  ts <- as.integer(sub(".*\\.", "", names(set)))
  ids <- as.integer(sub("\\..*", "", names(set)))

  df <- data.frame(i = ids, t = ts, matched.set.size = Lengthcol)
  colnames(df)[1:2] <- c(attr(set, "id.var"), attr(set, "t.var"))
  rownames(df) <- NULL

  if(verbose)
  {
    summary.result <- list()
    summary.result$overview <- df
    summary.result$set.size.summary <- summary(Lengthcol, ...)
    summary.result$number.of.treated.units <- length(set)
    summary.result$num.units.empty.set <- sum(Lengthcol == 0)
    summary.result$lag <- attr(set, "lag")
    return(summary.result)
  }
  else
  {
    return(df)
  }
}

#' Plot the distribution of the sizes of matched sets.
#'
#'
#' A plot method for creating a histogram of the distribution of the sizes of matched sets.
#' This method accepts all standard optional \code{hist} arguments via the \code{...} argument.
#' By default, empty matched sets (treated units that could not be
#' matched with any control units) are noted as a vertical bar at x = 0 and not included in the
#' regular histogram. See the \code{include.empty.sets} argument for more information about this.
#'
#' @param x a \code{matched.set} object
#' @param ... optional arguments to be passed to \code{hist}
#' @param border default is NA. This is the same argument as the standard argument for \code{hist}
#' @param col default is "grey". This is the same argument as the standard argument for \code{hist}
#' @param ylab default is "Frequency of Size". This is the same argument as the standard argument for \code{hist}
#' @param xlab default is "Matched Set Size". This is the same argument as the standard argument for \code{hist}
#' @param lwd default is NULL. This is the same argument as the standard argument for \code{hist}
#' @param main default is "Distribution of Matched Set Sizes". This is the same argument as the standard argument for \code{hist}
#' @param freq default is TRUE. See \code{freq} argument in \code{hist} function for more.
#' @param include.empty.sets logical value indicating whether or not empty sets should be included in the histogram. default is FALSE. If FALSE, then empty sets will be noted as a separate vertical bar at x = 0. If TRUE, empty sets will be included as normal sets.
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' plot(PM.results$att)
#' plot(PM.results$att, include.empty.sets = TRUE)
#'
#' @method plot matched.set
#' @export
plot.matched.set <- function(x, ..., border = NA, col = "grey", ylab = "Frequency of Size",
                             xlab ="Matched Set Size" , lwd = NULL,
                             main = "Distribution of Matched Set Sizes",
                             freq = TRUE, include.empty.sets = FALSE)
{
    set <- x
    lvec <- sapply(set, length)

    if(include.empty.sets)
    {
      graphics::hist(x = lvec, freq = freq, border = border, col = col, ylab = ylab, xlab = xlab, main = main, ...)
    }
    else
    {
      lvec.nonempty <- lvec[lvec > 0]

      if(sum(lvec == 0) > 0)
      {
        num.empties <- as.character(sum(lvec == 0))
        graphics::hist(x = lvec.nonempty, freq = freq, border = border, col = col, ylab = ylab,
                       xlab = xlab, main = main, ...)
        graphics::lines(x = c(0,0),
              y = c(0, num.empties),
              lwd = 4,
              col = "#ffc6c4")
      }
      else
      {
        graphics::hist(x = lvec.nonempty, freq = freq, border = border, col = col, ylab = ylab,
                       xlab = xlab, main = main, ...)
      }
    }

}

#' Print \code{matched.set} objects.
#'
#' @param x a \code{matched.set} object
#' @param verbose logical indicating whether or not output should be printed in expanded/raw list form.
#' The verbose form is not recommended unless the data set is small. Default is FALSE
#' @param ... additional arguments to be passed to \code{print}
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' print(PM.results$att)
#'
#'
#'
#' @method print matched.set
#' @export
print.matched.set <- function(x, ..., verbose = FALSE)
{
  set <- x
  if(verbose)
  {
    class(set) <- "list"
    print(set, ...)
  }

  else {
    print(summary(set, verbose = F), ...)
  }
}

#' @export
`[.matched.set` <- function(x, i, j = NULL, drop = NULL)
{

  if(!is.null(j)) stop("matched.set object is a list.")
  class(x) <- "list"
  temp <- x[i]
  attr(temp, "lag") <- attr(x, "lag")
  attr(temp, "refinement.method") <- attr(x, "refinement.method")
  attr(temp, "t.var") <- attr(x, "t.var")
  attr(temp, "id.var" ) <- attr(x, "id.var" )
  attr(temp, "treatment.var") <- attr(x, "treatment.var")
  attr(temp, "distances") <- attr(x, "distances")
  attr(temp, "max.match.size") <- attr(x, "max.match.size")
  attr(temp, "covs.formula") <- attr(x, "covs.formula")
  attr(temp, "match.missing") <- attr(x, "match.missing")
  class(temp) <- "matched.set"

  return(temp)
}

#helper function for get_covariate_balance()
build_balance_mats <- function(idx, ordered_expanded_data, msets)
{

  subset.per.matchedset <- function(sub.idx, set)
  {

    wts <- attr(set, "weights")[which(set == ordered_expanded_data[sub.idx[1:(length(sub.idx) - 1)], attr(msets, "id.var")])]
    return(cbind(ordered_expanded_data[sub.idx,], data.frame("weights" = c(wts, Inf))))
  }
  unnest <- function(mset.idx, mset)
  {
    # print(mset.idx)
    # browser()
    lapply(mset.idx, subset.per.matchedset, set = mset)
  }
  result <- mapply(FUN = unnest, mset.idx = idx, mset = msets, SIMPLIFY = FALSE)
  return(result)
}
