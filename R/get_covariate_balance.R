#' Calculate covariate balance measures for refined and unrefined matched sets
#'
#'
#' Calculate covariate balance for user specified covariates across matched sets. Balance is assessed by taking the average
#' of the difference between the values of the specified covariates for the treated unit(s) and the weighted average of
#' the control units across all matched sets. Results are standardized and are expressed in standard deviations.
#' Balance is calculated for each period in the specified lag window.
#' @param ... one or more PanelMatch objects
#' @param panel.data \code{PanelData} object
#' @param covariates a character vector, specifying the names of the covariates for which the user is interested in calculating balance.
#' @param include.unrefined logical. Indicates whether or not covariate balance measures for unrefined matched sets should be included. If TRUE, the function will return covariate balance results for the PanelMatch configurations provided, as well as a set of balance results that assume all matched controls have equal weight (i.e., the matched sets are unrefined). These results are included in addition to whatever PanelMatch configurations are specified to the function. Note that if you provide a PanelMatch object where no refinement is applied (that is, where \code{refinement.method = "none"}) and set this option to TRUE, then both sets of covariate balance results will be identical. If FALSE, then only balance calculations for the provided PanelMatch specifications are performed and returned. 
#' @returns A list of matrices, or a list of lists (if the QOI is ATE). The matrices contain the calculated covariate balance levels for each specified covariate for each period. Each element in the list (whether that be a matrix or a sublist) corresponds to a \code{PanelMatch} configuration specified to the function. Results are returned in the order they were provided. Unrefined results are stored as a parallel list object in an attribute called "unrefined.balance.results". 
#' 
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' # create subset of data for simplicity
#' #add some additional data to data set for demonstration purposes
#' dem.sub$rdata <- runif(runif(nrow(dem.sub)))
#' dem.sub.panel <- PanelData(dem.sub, "wbcode2", "year", "dem", "y")
#' PM.results <- PanelMatch(panel.data = dem.sub.panel, lag = 4, 
#'                          refinement.method = "ps.match", 
#'                          match.missing = TRUE, 
#'                          covs.formula = ~ tradewb + rdata,
#'                          size.match = 5, qoi = "att",
#'                          lead = 0:4, 
#'                          forbid.treatment.reversal = FALSE)
#' get_covariate_balance(PM.results, panel.data = dem.sub.panel, covariates = c("tradewb", "rdata"))
#'
#' @export
get_covariate_balance <- function(..., 
                                  panel.data, 
                                  covariates,
                                  include.unrefined = TRUE)
{
  if (!inherits(panel.data, "PanelData")) stop("Please provide a PanelData object.")
  
  if(is.null(covariates))
  {
    stop("Please specify covariates")
  }
  if(!all(covariates %in% colnames(panel.data)))
  {
    stop("Some of the specified covariates are not columns in the data set.")
  }
  
  pm.objs <- list(...)
  
  get_qois <- function(pm.obj) {
  
    if (inherits(pm.obj, "PanelMatch")) {
      if (attr(pm.obj, "qoi") %in% c("att", "art", "atc"))
      {
        qoi <- attr(pm.obj, "qoi")
      } else if (attr(pm.obj, "qoi") == "ate") {
        qoi <- c("att", "atc")
      } else {
        stop("Error extracting QOI from provided PanelMatch object")
      }
    } else {
      stop("invalid object: not a PanelMatch object")
    }  
    return(qoi)
  }
  
  qoi.sets <- lapply(pm.objs,
                     get_qois)  
  
  
  balance.results <- list()
  unrefined.balance.results <- list()
  for (i in 1:length(pm.objs)) { # for each PM object
    sub.list <- list()
    unrefined.sub.list <- list()
    for (q in qoi.sets[[i]]) { # get the qoi to extract the matched set object. unless qoi = ate, this will just be either att, art, atc. In case of ate, we need to loop over both att and atc
      pm.obj <- pm.objs[[i]]
      matched.set <- pm.obj[[q]]
      sub.list[[q]] <- get_set_covariate_balance(matched.set, 
                                                 panel.data,
                                                 covariates,
                                                 use.equal.weights = FALSE)
      balance.results[[i]] <- sub.list
      if (include.unrefined)
      {
        unrefined.sub.list[[q]] <- get_set_covariate_balance(matched.set, 
                                                   panel.data,
                                                   covariates,
                                                   use.equal.weights = TRUE)
        unrefined.balance.results[[i]] <- unrefined.sub.list
      }
    }
  }
  
  class(balance.results) <- c("PanelBalance", "list")
  attr(balance.results, "treatment") <- attr(panel.data, 'treatment')
  if(length(unrefined.balance.results) == 0)  unrefined.balance.results <- NULL;
  if (!is.null(unrefined.balance.results))
  {
    class(unrefined.balance.results) <- c("PanelBalance", "list")
    attr(unrefined.balance.results, "treatment") <- attr(balance.results, "treatment")
  }
  attr(balance.results, "unrefined.balance.results") <- unrefined.balance.results
  attr(balance.results, "covariates") <- covariates
  return(balance.results)
}




get_set_covariate_balance <- function(matched.sets, 
                                  panel.data, 
                                  covariates,
                                  use.equal.weights = FALSE,
                                  plot = FALSE,
                                  reference.line = TRUE,
                                  legend = TRUE, 
                                  ylab = "SD",
                                  include.treatment.period = TRUE,
                                  legend.position = "topleft")
{
  
  
  if (!inherits(matched.sets, "matched.set")) 
  {
    stop("Error extracting matched.set object")
  }
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  treatment <- attr(matched.sets, "treatment.var")

  
  matched.sets <- matched.sets[sapply(matched.sets, length) > 0]
  
  
  othercols <- colnames(panel.data)[!colnames(panel.data) %in% c(time.id, 
                                                                     unit.id, 
                                                                     treatment)]
  #subset to keep only needed covariates
  othercols <- othercols[othercols %in% covariates]
  
  ordered.data <- panel.data[, c(unit.id, time.id, treatment, othercols), 
                               drop = FALSE] #reorder columns 
  ordered.data <- ordered.data[, unique(c(unit.id, time.id, 
                                          treatment, covariates)), drop = FALSE]
  
  if(is.null(attr(matched.sets[[1]], "weights")) || 
     use.equal.weights)
  {
    for(i in 1:length(matched.sets))
    {
      attr(matched.sets[[i]], "weights") <- rep(1/length(matched.sets[[i]]), 
                                                length(matched.sets[[i]]))
      names(attr(matched.sets[[i]], "weights")) <- matched.sets[[i]]
    }
  }
  treated.ts <- as.integer(sub(".*\\.", "", names(matched.sets)))
  treated.ids <- as.integer(sub("\\..*", "", names(matched.sets)))
  tlist <- expand_treated_ts(lag, treated.ts)
  
  
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, tlist, 
                              matched.sets, lag)
  
  
  balance_mats <- build_balance_mats(ordered_expanded_data = ordered.data, 
                                     idx =  idxlist, msets = matched.sets)
  unlistedmats <- unlist(balance_mats, recursive = F)
  plotpoints <- list()
  for(k in 1:(lag+1))
  {
    var.points <- list()
    for(i in 1:length(covariates))
    {
      variable <- covariates[i]
      
      sd.val <- sd(sapply(unlistedmats[seq(from = k, 
                                           to = (length(matched.sets) * 
                                                   (lag + 1)), 
                                           by = lag + 1)],
                          function(x){x[nrow(x), variable]}), na.rm = T)
      if(isTRUE(all.equal(sd.val, 0)))
      {
        sd.val <- NA #make everything fail in a standardized way
      }
      
      tprd <- unlistedmats[seq(from = k, 
                               to = (length(matched.sets) * (lag + 1)), 
                               by = lag + 1)] 
      #should represent the same relative period across all matched sets 
      # each df is a matched set
      
      get_mean_difs <- function(x, variable) #x is an individual data frame
      {
        return( x[nrow(x), variable] - 
                  sum(x[1:(nrow(x) -1),"weights"] * x[1:(nrow(x) - 1), 
                                                      variable], na.rm = T ))
      }
      diffs <- sapply(tprd, get_mean_difs, variable = variable)
      
      var.points[[i]] <- mean(diffs / sd.val, na.rm = T)
    }
    names(var.points) <- covariates
    plotpoints[[k]] <- var.points
    
  }
  
  names(plotpoints) <- paste0("t_", lag:0)
  pointmatrix <- apply((as.matrix(do.call(rbind, plotpoints))), 
                       2,
                       function(x){(as.numeric(x))}, 
                       simplify = TRUE)
  rownames(pointmatrix) <- names(plotpoints)
  
  remove.vars.idx <- apply(apply(pointmatrix, 2, is.nan), 2, any)
  if(sum(remove.vars.idx) > 0)
  {
    removed.vars <- names(which(apply(apply(pointmatrix, 
                                            2, 
                                            is.nan), 
                                      2, 
                                      any)))
    pointmatrix <- pointmatrix[, !remove.vars.idx]
    warning(paste0("Some variables were removed due to low variation, inadequate data needed for calculation: ", removed.vars))
  }
  
  return(pointmatrix)
  
  
}

build_balance_mats <- function(idx, ordered_expanded_data, msets)
{
  subset.per.matchedset <- function(sub.idx, set)
  {
    
    wts <- attr(set, "weights")[which(set == ordered_expanded_data[sub.idx[1:(length(sub.idx) - 1)], 
                                                                   attr(msets, "id.var")])]
    return(cbind(ordered_expanded_data[sub.idx,], 
                 data.frame("weights" = c(wts, Inf))))
  }
  unnest <- function(mset.idx, mset)
  {
    lapply(mset.idx, subset.per.matchedset, set = mset)
  }
  result <- mapply(FUN = unnest, mset.idx = idx, 
                   mset = msets, SIMPLIFY = FALSE)
  return(result)
}