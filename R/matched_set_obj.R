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
  temp <-unlist(strsplit(names(set), split = ".", fixed = TRUE))
  ids <- temp[c(T,F)]
  ts <- temp[c(F,T)]
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
    
    lapply(mset.idx, subset.per.matchedset, set = mset)
  }
  result <- mapply(FUN = unnest, mset.idx = idx, mset = msets, SIMPLIFY = FALSE)
  return(result)
}

#' Calculate covariate balance
#' 
#' 
#' Calculate covariate balance for user specified covariates across matched sets. Balance is assessed by taking the average 
#' of the difference between the values of the specified covariates for the treated unit(s) and the weighted average of 
#' the control units across all matched sets. Results are standardized and are expressed in standard deviations. 
#' Balance is calculated for each period in the specified lag window. 
#' @param matched.sets A \code{matched.set} object
#' @param data The time series cross sectional data set (as a \code{data.frame} object) used to produce the \code{matched.set} object. This data set should be identical to the one passed to \code{PanelMatch} and \code{PanelEstimate} to ensure consistent results.
#' @param verbose logical. When TRUE, the function will return more information about the calculations/results. When FALSE, a more compact version of the results/calculations are returned.
#' @param plot logical. When TRUE, a plot showing the covariate balance calculation results will be shown. When FALSE, no plot is made, but the results of the calculations are returned. default is FALSE
#' @param covariates a character vector, specifying the names of the covariates for which the user is interested in calculating balance. 
#' @param reference.line logical indicating whether or not a horizontal line should be present on the plot at y = 0. Default is TRUE.
#' @param legend logical indicating whether or not a legend identifying the variables should be included on the plot. Default is TRUE.
#' @param ylab Label for y axis. Default is "SD". This is the same as the ylab argument to \code{plot}.
#' @param use.equal.weights logical. If set to TRUE, then equal weights will be assigned to control units, rather than using whatever calculated weights have been assigned. This is helpful for assessing the improvement in covariate balance as a result of refining the matched sets.
#' @param ... Additional graphical parameters to be passed to the \code{plot} function in base R.
#' @examples 
#' #add some additional data to data set for demonstration purposes
#' dem$rdata <- runif(runif(nrow(dem))) 
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, time.id = "year", unit.id = "wbcode2", treatment = "dem",
#'                     outcome.var ="y", refinement.method = "mahalanobis", 
#'                     data = dem, match.missing = TRUE,
#'                     covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                     size.match = 5, qoi = "att")
#' get_covariate_balance(pm.obj$att, dem, covariates = c("tradewb", "rdata"), 
#'                          ylim = c(-2,2))
#' get_covariate_balance(pm.obj$att, dem, covariates = c("tradewb", "rdata"), 
#'                          plot = TRUE, ylim = c(-2,2))
#'
#' @export
get_covariate_balance <- function(matched.sets, data,  covariates, use.equal.weights = FALSE,
                                  verbose = TRUE, plot = FALSE, 
                                  reference.line = TRUE, legend = TRUE, ylab = "SD",...)
{
  if(is.null(covariates))
  {
    stop("please specify the covariates for which you would like to check the balance")
  }
  if(!all(covariates %in% colnames(data)))
  {
    stop("Some of the specified covariates are not columns in the data set.")
  }
  if(!any(class(matched.sets) %in% "matched.set")) stop("Please pass a matched.set object")
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  treatment <- attr(matched.sets, "treatment.var")
  if(class(data[, unit.id]) == "factor") stop("please convert unit id column to character, integer, or numeric")
  if(class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data), formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table(testmat), id = unit.id, variable = time.id, value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table(d), data.table(data), all.x = TRUE, by = c(unit.id, time.id))
    data <- as.data.frame(data)
    #data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  ordered.data[, paste0(unit.id, ".int")] <- as.integer(as.factor(data[, unit.id]))
  
  if(class(ordered.data[, unit.id]) == "character") {
    unit.index.map <- data.frame(original.id = make.names(as.character(unique(ordered.data[, unit.id]))), 
                                 new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(ordered.data[, unit.id]) == "integer") {
    unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))), 
                                 new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(ordered.data[, unit.id]) == "numeric") {
    if(all(unique(ordered.data[, unit.id]) == as.integer(unique(ordered.data[, unit.id])))) #actually integers
    {
      unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))),
                                   new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
    }
    else
    {
      stop("Unit ID data appears to be a non-integer numeric. Please convert.")
    }
  }
  else {
    stop("Unit ID Data is not integer, numeric, or character.")
  }
  
  og.unit.id <- unit.id
  #og.time.id <- time.id
  unit.id <- paste0(unit.id, ".int")
  matched.sets <- matched.sets[sapply(matched.sets, length) > 0]
  matched.sets <- encode_index(matched.sets, unit.index.map, unit.id)
  #they will either all have or not have weights, so we can check the first matched set to see if we need to add equal weighting
  #i dont think that its possible for sets to not have any weights now, but don't think it hurts to keep this in
  if(is.null(attr(matched.sets[[1]], "weights")) | use.equal.weights)
  {
    for(i in 1:length(matched.sets))
    {
      attr(matched.sets[[i]], "weights") <- rep(1/length(matched.sets[[i]]), length(matched.sets[[i]]))
      names(attr(matched.sets[[i]], "weights")) <- matched.sets[[i]]
    }
  }
  treated.ts <- as.integer(unlist(strsplit(names(matched.sets), split = "[.]"))[c(F,T)])
  treated.ids <- as.integer(unlist(strsplit(names(matched.sets), split = "[.]"))[c(T,F)])
  tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
  
  idxlist <- get_yearly_dmats(matrix(nrow = 0, ncol = 0), treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
                                                                       ordered.data[, time.id]), matched_sets = matched.sets, lag)
  balance_mats <- build_balance_mats(ordered_expanded_data = ordered.data, idx =  idxlist, msets = matched.sets)
  unlistedmats <- unlist(balance_mats, recursive = F)
  
  plotpoints <- list()
  for(k in 1:(lag+1))
  {	
    var.points <- list()
    for(i in 1:length(covariates))
    {
      variable <- covariates[i]
      
      sd.val <- sd(sapply(unlistedmats[seq(from = 1, to = (length(matched.sets) * (lag + 1)), by = lag + 1)], 
                          function(x){x[nrow(x), variable]}), na.rm = T)
      if(isTRUE(all.equal(sd.val, 0)))
      {
        sd.val <- NA #make everything fail in a standardized way
      }
      
      tprd <- unlistedmats[seq(from = k, to = (length(matched.sets) * (lag + 1)), by = lag + 1)] #should represent the same relative period across all matched sets. each df is a matched set
      
      get_mean_difs <- function(x, variable) #x is an individual data frame
      {
        return( x[nrow(x), variable] - sum(x[1:(nrow(x) -1),"weights"] * x[1:(nrow(x) -1), variable], na.rm = T ))
      }
      diffs <- sapply(tprd, get_mean_difs, variable = variable)
      
      var.points[[i]] <- mean(diffs / sd.val, na.rm = T)
    }
    names(var.points) <- covariates
    plotpoints[[k]] <- var.points
    
  }
  names(plotpoints) <- paste0("t_", lag:0)
  pointmatrix <- apply((as.matrix(do.call(rbind, plotpoints))), 2, function(x){(as.numeric(x))})
  rownames(pointmatrix) <- names(plotpoints)
  
  remove.vars.idx <- apply(apply(pointmatrix, 2, is.nan), 2, any)
  if(sum(remove.vars.idx) > 0)
  {
    removed.vars <- names(which(apply(apply(pointmatrix, 2, is.nan), 2, any)))
    pointmatrix <- pointmatrix[, !remove.vars.idx]  
    warning(paste0("Some variables were removed due to low variation: ", removed.vars))
  }
  
  
  
  
  
  if(!plot) return(pointmatrix)
  
  if(plot)
  {
    graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), lty =1, ylab = ylab, xaxt = "n", ...)
    graphics::axis(side = 1, labels = paste0("t-", (nrow(pointmatrix) - 1):0), at = 1:nrow(pointmatrix))
    if(legend) legend("topleft", legend = colnames(pointmatrix), col=1:ncol(pointmatrix), lty = 1)
    if(reference.line) graphics::abline(h = 0, lty = "dashed")
  }
}

#' balance_scatter
#'
#' Visualizing the standardized mean differences for covariates via a scatter plot. 
#'
#' \code{balance_scatter} visualizes the standardized mean differences for each covariate.
#' Although users can use the scatter plot in a variety of ways, it is recommended that 
#' the x-axis refers to balance for covariates before refinement, and y-axis
#' refers to balance after refinement. Users can utilize parameters powered by \code{plot}
#' in base R to further customize the figure. 
#' @param non_refined_set a \code{matched.set} object produced by setting `refinement.method` to "none" in `PanelMatch` 
#' @param refined_list a list of one or two \code{matched.set} objects
#' @param xlim xlim of the scatter plot. This is the same as the \code{xlim} argument in \code{plot}
#' @param ylim ylim of the scatter plot. This is the same as the \code{ylim} argument in \code{plot}
#' @param main title of the scatter plot. This is the same as the \code{main} argument in \code{plot}
#' @param x.axis.label x axis label
#' @param y.axis.label y axis label
#' @param pchs one or two pch indicators for the symbols on the scatter plot. See \code{plot} for more information
#' @param covariates variables for which balance is displayed
#' @param data the same time series cross sectional data set used to create the matched sets.
#' @param ... optional arguments to be passed to \code{plot} 
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <adamrauh@mit.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples
#' # get a matched set without refinement
#' sets0 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "none",
#'                     data = dem, match.missing = FALSE,
#'                     covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#' 
#' # get a matched set with refinement using CBPS.match, setting the 
#' # size of matched set to 5
#' sets1 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "mahalanobis",
#'                     data = dem, match.missing = FALSE,
#'                     covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#'                     size.match = 5, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#' 
#' # get another matched set with refinement using CBPS.weight
#' sets2 <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                     treatment = "dem", refinement.method = "ps.weight",
#'                     data = dem, match.missing = FALSE,
#'                     covs.formula = ~ I(lag(y, 1:4)) + I(lag(tradewb, 1:4)),
#'                     size.match = 10, qoi = "att",
#'                     outcome.var = "y",
#'                     lead = 0:4, forbid.treatment.reversal = FALSE)
#' 
#' 
#' # use the function to produce the scatter plot
#' balance_scatter(non_refined_set = sets0$att,
#'               refined_list = list(sets1$att, sets2$att),
#'               data = dem,
#'               covariates = c("y", "tradewb"))
#' # add legend 
#' legend(x = 0, y = 0.8,  
#' legend = c("mahalanobis",
#'            "PS weighting"),
#' y.intersp = 0.65,
#' x.intersp = 0.3,
#' xjust = 0,
#' pch = c(1, 3), pt.cex = 1,
#' bty = "n", ncol = 1, cex = 1, bg = "white")
#'
#' 
#'
#' @export
balance_scatter <- function(non_refined_set, refined_list,
                              xlim = c(0, .8),
                              ylim = c(0, .8),
                              main = "Standardized Mean Difference of Covariates",
                              pchs = c(2,3),
                              covariates, data, 
                              x.axis.label = "Before refinement",
                              y.axis.label = "After refinement",
                              ...) {
  # first, get balance for non-refined set
  non_refined_balance <- get_covariate_balance(matched.sets = non_refined_set, 
                                               data = data,
                                               covariates = covariates)
  
  # second, get balance for refined sets
  refined_balance <- list()
  for (i in 1:length(refined_list)) {
    refined_balance[[i]] <-
      get_covariate_balance(matched.sets = refined_list[[i]], 
                            data = data,
                            covariates = covariates)
  }
  
  # extract values for x-axis from the non-refined sets
  benchmark <- non_refined_balance
  benchmark <- as.vector(benchmark[1:(nrow(benchmark)-1),]) # delete balance results after t-1
  
  # extract values for y-axis from refined sets and delete balance results after t-1
  compared <- sapply(refined_balance, function(x) x <- x[1:(nrow(x)-1),])
  
  graphics::plot(abs(as.numeric(benchmark)), 
       abs(as.numeric(compared[,1])), pch = 1,
       xlab = x.axis.label,
       ylab = y.axis.label,
       xlim = xlim,
       ylim = ylim,
       main = main,
       font.main = 1, ...)
  # logical statement for the length of the refined_balance
  if (length(refined_balance) > 1) {
    for (j in 2:length(refined_balance)) {
      graphics::points(abs(as.numeric(benchmark)), 
             abs(as.numeric(compared[,j])),
             pch = pchs[j])  
    }
  } 
  
  graphics::abline(h = 0, lty = "dashed")
  graphics::abline(0, 1, lty = 2, col = "red")
  
  
} 



decode_index <- function(mset, unit.index, original.unit.id)#, original.time.id)
{
  decode.control.units <- function(in_, unit.index)
  {
    news <- unit.index$original.id[match(in_, unit.index$new.id)]
    if(!is.null(attr(in_, "distances")))
    {
      attr(news, "distances") <- attr(in_, "distances")
      names(attr(news, "distances")) <- news
    }
    if(!is.null(attr(in_, "weights")))
    {
      attr(news, "weights") <- attr(in_, "weights")
      names(attr(news, "weights")) <- news
    }
    return(news)
  }
  new.mset <- sapply(mset, decode.control.units, unit.index = unit.index, simplify = F)
  decode.treated.units <- function(ts, ids, unit.index)#, time.index)
  {
    n.ids <- unit.index$original.id[match(ids, unit.index$new.id)]
    
    return(paste0(n.ids, ".", ts))
  }
  treated.ts <- as.numeric(unlist(strsplit(names(mset), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(mset), split = "[.]"))[c(T,F)])
  attributes(new.mset) <- attributes(mset)
  names(new.mset) <- decode.treated.units(treated.ts, treated.ids, unit.index)
  attr(new.mset, "id.var") <- original.unit.id
  
  return(new.mset)
}


encode_index <- function(mset, unit.index, new.unit.id)
{
  encode.control.units <- function(in_, unit.index)
  {
    news <- unit.index$new.id[match(in_, unit.index$original.id)]
    if(!is.null(attr(in_, "distances")))
    {
      attr(news, "distances") <- attr(in_, "distances")
      names(attr(news, "distances")) <- news
    }
    if(!is.null(attr(in_, "weights")))
    {
      attr(news, "weights") <- attr(in_, "weights")
      names(attr(news, "weights")) <- news
    }
    return(news)
  }
  new.mset <- sapply(mset, encode.control.units, unit.index = unit.index, simplify = F)
  encode.treated.units <- function(ts, ids, unit.index)#, time.index)
  {
    n.ids <- unit.index$new.id[match(ids, unit.index$original.id)]
    
    return(paste0(n.ids, ".", ts))
  }
  treated.ts <- (unlist(strsplit(names(mset), split = "[.]"))[c(F,T)])
  treated.ids <- (unlist(strsplit(names(mset), split = "[.]"))[c(T,F)])
  attributes(new.mset) <- attributes(mset)
  names(new.mset) <- encode.treated.units(treated.ts, treated.ids, unit.index)
  attr(new.mset, "id.var") <- new.unit.id
  
  return(new.mset)
}
