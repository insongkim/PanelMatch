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
#' @author Adam Rauh <amrauh@umich..edu>, In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, and Kosuke Imai <imai@harvard.edu>
#' @export
matched_set <- function(matchedsets, id, t, L, t.var, id.var, treatment.var)
{
  if(length(matchedsets) != length(id) || 
     length(matchedsets) != length(t) || 
     length(id) != length(t))
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
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' dem.sub.panel <- PanelData(dem.sub, 'wbcode2', 'year', 'dem', 'y')
#' # create subset of data for simplicity
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "ps.match",
#'                          data = dem.sub.panel, match.missing = TRUE,
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
#'
#' @examples
#'
#' @method plot matched.set
#' @export
plot.matched.set <- function(x, ..., panel.data, type = "weights", 
                             include.missing = TRUE,
                             low.color = "blue", 
                             mid.color = "white", 
                             high.color = "red", 
                             missing.color = "grey50")
{
  
  
  
  attr(panel.data, "unit.id") -> unit.id
  attr(panel.data, "time.id") -> time.id
  
  rownames(panel.data) <- paste0(panel.data[, unit.id], ".", panel.data[, time.id])
  
  if (type == "weights") {
    x.s <- weights(x)
  } else if (type == "distances") 
  {
    x.s <- distances(x)
  } else {
    stop("type is misspecified")
  }
  
  
  meta.dt <- combine_named_vectors(x.s)
  meta.dt[meta.dt == 0] <- NA
  unique.units <- unique(panel.data[, unit.id]) 
  if (include.missing)
  {
    missing.units <- as.character(setdiff(unique.units, colnames(meta.dt)))
    meta.dt[missing.units] <- NA_real_
  }
  
  
  treated.ts <- as.integer(sub(".*\\.", "", rownames(meta.dt)))
  treated.ids <- as.integer(sub("\\..*", "", rownames(meta.dt)))
  
  meta.dt <- meta.dt[order(treated.ts), ]
  
  dt <- data.table::as.data.table(meta.dt, keep.rownames = "Row")
  
  df_long <- data.table::melt(dt, id.vars = "Row", 
                              variable.name = "Units", 
                              value.name = "Weight")
  
  p <- ggplot2::ggplot(df_long, aes(x = Units, y = Row, fill = Weight)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradientn(colors = c(low.color, mid.color, high.color), 
                         na.value = missing.color) +
    ggplot2::theme_minimal() + labs(y = "Treated Observations")
  return(p)
}


combine_named_vectors <- function(list_of_vectors) {
  # Get all unique names from the vectors
  all_names <- unique(unlist(lapply(list_of_vectors, names)))
  
  # Create a list to hold data frames
  list_of_data_frames <- lapply(list_of_vectors, function(vec) {
    # Create a named vector with all names initialized to NA
    filled_vec <- setNames(rep(NA, length(all_names)), all_names)
    
    # Replace the NA values with the actual values from the vector
    filled_vec[names(vec)] <- vec
    
    # Convert the named vector to a data frame with one row
    as.data.frame(t(filled_vec))
  })
  
  # Combine all data frames into a single data frame
  combined_df <- do.call(rbind, list_of_data_frames)
  
  # Set the row names to the names of the list elements
  rownames(combined_df) <- names(list_of_vectors)
  
  return(combined_df)
}

#' Print \code{matched.set} objects.
#'
#' @param x a \code{matched.set} object
#' @param verbose logical indicating whether or not output should be printed in expanded/raw list form.
#' The verbose form is not recommended unless the data set is small. Default is FALSE
#' @param ... Not used. additional arguments to be passed to \code{print}
#'
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' # create subset of data for simplicity
#' dem.sub.panel <- PanelData(dem.sub, 'wbcode2', 'year', 'dem', 'y')
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "ps.match",
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
    print(set)
  }

  else {
    print(summary(set, verbose = FALSE))
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



handle_pm_qoi <- function(pm.in, qoi.in)
{
  if (!is.null(qoi.in) && identical(qoi.in, "ate"))
  {
    stop("Please specify qoi = NULL, att, art, or atc.")
  }
  if (is.null(qoi.in))
  {
    qoi.actual <- attr(pm.in, "qoi")
  } else {
    qoi.actual <- qoi.in
  }
  msets <- pm.in[[qoi.actual]]
  return(msets)
}

#' @export
weights <- function(object) {
  UseMethod("weights")
}

#' Extract the weights of matched.set objects
#'
#' @param object matched.set object, extracted using the \code{get.PanelMatch()} method
#'
#' @return list of named vectors. Each list element corresponds to a particular treated observation and contains the matched control units, along with their weights.
#' @export
#'
#' @examples
weights.matched.set <- function(object)
{
  weight.list <- lapply(object, function(x) attr(x, "weights"))
  names(weight.list) <- names(object)
  return(weight.list)
}

#' @export
distances <- function(object)
{
  UseMethod("distances")
}

#' Extract distances of matched sets
#'
#' @param object a matched.set object
#'
#' @return A named list of named vectors. Each element corresponds to a matched set and will be a named vector, where the names of each element will identify a matched control unit and its distance from the treated observation within a particular matched set. These correspond to the "distances" attribute, which are calculated and included when the \code{verbose} option is set to TRUE in \code{PanelMatch}.
#' @export
distances.matched.set <- function(object)
{
  distance.list <- lapply(object, function(x) attr(x, "distances"))
  if(all(sapply(distance.list, is.null)))
  {
    stop("distances attribute is missing. Please specify verbose = TRUE and/or specify one of the following methods for refinement: mahalanobis, ps.match, CBPS.match")
  } else {
    names(distance.list) <- names(object)
    return(distance.list)
  }
}