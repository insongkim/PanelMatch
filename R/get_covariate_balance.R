#' Calculate covariate balance
#'
#'
#' Calculate covariate balance for user specified covariates across matched sets. Balance is assessed by taking the average
#' of the difference between the values of the specified covariates for the treated unit(s) and the weighted average of
#' the control units across all matched sets. Results are standardized and are expressed in standard deviations.
#' Balance is calculated for each period in the specified lag window.
#' @param matched.sets A \code{matched.set} object
#' @param data The time series cross sectional data set (as a \code{data.frame} object) used to produce the \code{matched.set} object. This data set should be identical to the one passed to \code{PanelMatch()} and \code{PanelEstimate()} to ensure consistent results.
#' @param plot logical. When TRUE, a plot showing the covariate balance calculation results will be shown. When FALSE, no plot is made, but the results of the calculations are returned. default is FALSE
#' @param covariates a character vector, specifying the names of the covariates for which the user is interested in calculating balance.
#' @param reference.line logical indicating whether or not a horizontal line should be present on the plot at y = 0. Default is TRUE.
#' @param legend logical indicating whether or not a legend identifying the variables should be included on the plot. Default is TRUE.
#' @param ylab Label for y axis. Default is "SD". This is the same as the ylab argument to \code{plot()}.
#' @param use.equal.weights logical. If set to TRUE, then equal weights will be assigned to control units, rather than using whatever calculated weights have been assigned. This is helpful for assessing the improvement in covariate balance as a result of refining the matched sets.
#' @param include.treatment.period logical. Default is TRUE. When TRUE, covariate balance measures for the period during which treatment occurs is included. These calculations are not included when FALSE. Users may wish to leave this period off in some circumstances. For instance, one would expect covariate balance to be poor during this period when treatment is continuous and a lagged outcome is included in the refinement formula.
#' @param legend.position position of legend. See documentation for graphics::legend. Default is "topleft"
#' @param ... Additional graphical parameters to be passed to the \code{plot} function in base R.
#' @examples
#' dem.sub <- dem[dem[, "wbcode2"] <= 100, ]
#' # create subset of data for simplicity
#' #add some additional data to data set for demonstration purposes
#' dem.sub$rdata <- runif(runif(nrow(dem.sub)))
#' pm.obj <- PanelMatch(lead = 0:3, lag = 4, time.id = "year", unit.id = "wbcode2", treatment = "dem",
#'                     outcome.var ="y", refinement.method = "ps.match",
#'                     data = dem.sub, match.missing = TRUE,
#'                     covs.formula = ~ tradewb + rdata + I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                     size.match = 5, qoi = "att")
#' get_covariate_balance(pm.obj$att, dem.sub, covariates = c("tradewb", "rdata"),
#'                          ylim = c(-2,2))
#'
#' @export
get_covariate_balance <- function(matched.sets, 
                                  data, 
                                  covariates,
                                  use.equal.weights = FALSE,
                                  plot = FALSE,
                                  reference.line = TRUE,
                                  legend = TRUE, 
                                  ylab = "SD",
                                  include.treatment.period = TRUE,
                                  legend.position = "topleft",
                                  ...)
{
  
  if(is.null(covariates))
  {
    stop("Please specify covariates")
  }
  if(!all(covariates %in% colnames(data)))
  {
    stop("Some of the specified covariates are not columns in the data set.")
  }
  
  if (!inherits(matched.sets, "matched.set")) 
  {
    stop("Please pass a matched.set object")
  }
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  treatment <- attr(matched.sets, "treatment.var")
  
  
  if (!inherits(data[, unit.id], "integer") && 
      !inherits(data[, unit.id], "numeric"))
  {
    stop("please convert unit id column to integer or numeric")
  }
  
  
  
  
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data), 
                                 formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table::data.table(testmat), 
                          id = unit.id, 
                          variable = time.id, 
                          value = treatment,
                          variable.factor = FALSE, 
                          value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table::data.table(d), 
                  data.table::data.table(data), 
                  all.x = TRUE, 
                  by = c(unit.id, time.id))
    data <- as.data.frame(data)
    
  }
  
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  ordered.data <- check_time_data(ordered.data, time.id)
  
  
  matched.sets <- matched.sets[sapply(matched.sets, length) > 0]
  
  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, 
                                                                     unit.id, 
                                                                     treatment)]
  #subset to keep only needed covariates
  othercols <- othercols[othercols %in% covariates]
  
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols), 
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
  
  if (!include.treatment.period)
  {
    pointmatrix <- pointmatrix[-nrow(pointmatrix), ,drop = FALSE]
    stop.val <- 1
    start.val <- nrow(pointmatrix)
  } else {
    stop.val <- 0
    start.val <- nrow(pointmatrix) - 1
  }
  
  
  if (!plot) return(pointmatrix)
  
  if (plot)
  {
    treated.included <- treatment %in% colnames(pointmatrix)
    
    if (treated.included)
    {
      treated.data <- pointmatrix[,which(colnames(pointmatrix) == treatment)] # treated data
      pointmatrix <- pointmatrix[,-which(colnames(pointmatrix) == treatment)] #all non-treatment variable data
      graphics::matplot(pointmatrix, pch = 19,
                        type = "b", 
                        col = 1:ncol(pointmatrix), 
                        lty = 1, ylab = ylab, xaxt = "n", ...)
      graphics::lines(x = 1:nrow(pointmatrix), 
                      y = as.numeric(treated.data), 
                      type = "b",
                      lty = 2, lwd = 3)
      graphics::axis(side = 1, labels = paste0("t-", start.val:stop.val), 
                     at = 1:nrow(pointmatrix), ...)  
    } else
    {
      
      graphics::matplot(pointmatrix, type = "b",
                        pch = 19,
                        col = 1:ncol(pointmatrix), 
                        lty = 1, ylab = ylab, xaxt = "n", ...)
      graphics::axis(side = 1, labels = paste0("t-", start.val:stop.val), 
                     at = 1:nrow(pointmatrix), ...)  
    }
    
    if (legend) {
      if (treated.included)
      {
        legend(legend.position, 
               legend = c(colnames(pointmatrix), treatment), 
               col = c(1:ncol(pointmatrix), "black"), 
               lty = c(rep(1, ncol(pointmatrix)), 2))  
      } else {
        legend(legend.position, 
               legend = colnames(pointmatrix), 
               col = 1:ncol(pointmatrix), lty = 1)
      }
      
    }
    if(reference.line) graphics::abline(h = 0, lty = "dashed")
  }
}

#helper function 
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
