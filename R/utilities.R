check_time_data <- function(data, time.id)
{
  #if(class(data[, time.id]) != "integer") stop("time data is not integer")
  if (!inherits(data[, time.id], "integer")) stop("time data is not integer")
  u.times <- unique(data[, time.id])
  increase.by.one <- all(seq(min(u.times), max(u.times), by = 1) %in% u.times)
  if(increase.by.one)
  {
    return(TRUE)
  }
  else
  {
    stop("integer representation of time data has problematic gaps, as it does not increase by one. Perhaps time data for observations is irregular/not uniform across units?")
  }
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
get_covariate_balance <- function(matched.sets, 
                                  data, 
                                  covariates,
                                  use.equal.weights = FALSE,
                                  verbose = TRUE,
                                  plot = FALSE,
                                  reference.line = TRUE,
                                  legend = TRUE, 
                                  ylab = "SD",
                                  ...)
{
  # not ready yet, so provide defaults 
  calculate.network.proportion.balance = FALSE
  calculate.network.count.balance = FALSE
  adjacency.matrix = NULL
  neighborhood.degree = NULL
  continuous.treatment = FALSE
  
  if(is.null(covariates))
  {
    stop("please specify the covariates for which you would like to check the balance")
  }
  if(!all(covariates %in% colnames(data)))
  {
    stop("Some of the specified covariates are not columns in the data set.")
  }
  
  #if(!any(class(matched.sets) %in% "matched.set")) stop("Please pass a matched.set object")
  if (!inherits(matched.sets, "matched.set")) stop("Please pass a matched.set object")
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  treatment <- attr(matched.sets, "treatment.var")
  
  #if (!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
  if (!inherits(data[, unit.id], "integer") && !inherits(data[, unit.id], "numeric")) stop("please convert unit id column to integer or numeric")
  
  #if (class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  if (!inherits(data[, time.id], "integer")) stop("please convert time id to consecutive integers")
  
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
  if(calculate.network.proportion.balance || calculate.network.count.balance)
  {
    if(is.null(adjacency.matrix))
    {
      stop("Please provide adjacency matrix")
    }
    ordered.data <- calculate_neighbor_treatment(data = ordered.data, 
                                                 edge.matrix = adjacency.matrix, 
                                                 n.degree = neighborhood.degree, 
                                                 unit.id = unit.id, time.id = time.id, 
                                                 treatment.variable = treatment)
    
    if(!is.null(calculate.network.proportion.balance))
    {
      covariates <- c(covariates,
                      make.names(paste0('neighborhood_t_prop', '.', 1:neighborhood.degree)))
    } 
    if(!is.null(calculate.network.count.balance))
    {
      covariates <- c(covariates,
                      make.names(paste0('neighborhood_t_count', '.', 1:neighborhood.degree)))
    }
  }
  
  
  matched.sets <- matched.sets[sapply(matched.sets, length) > 0]
  # matched.sets <- encode_index(matched.sets, unit.index.map, unit.id)
  
  othercols <- colnames(ordered.data)[!colnames(ordered.data) %in% c(time.id, unit.id, treatment)]
  #subset to keep only needed covariates
  othercols <- othercols[othercols %in% covariates]
  
  ordered.data <- ordered.data[, c(unit.id, time.id, treatment, othercols), drop = FALSE] #reorder columns 
  ordered.data <- ordered.data[, unique(c(unit.id, time.id, treatment, covariates)), drop = FALSE]
  # Don't think we need to keep any other columns than the ones here...
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
  treated.ts <- as.integer(sub(".*\\.", "", names(matched.sets)))
  treated.ids <- as.integer(sub("\\..*", "", names(matched.sets)))
  tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
  
  
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, tlist, 
                              matched_sets = matched.sets, lag)
  
  
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
                                           to = (length(matched.sets) * (lag + 1)), 
                                           by = lag + 1)],
                          function(x){x[nrow(x), variable]}), na.rm = T)
      if(isTRUE(all.equal(sd.val, 0)))
      {
        sd.val <- NA #make everything fail in a standardized way
      }
      
      tprd <- unlistedmats[seq(from = k, to = (length(matched.sets) * (lag + 1)), 
                               by = lag + 1)] #should represent the same relative period across all matched sets. each df is a matched set
      
      get_mean_difs <- function(x, variable) #x is an individual data frame
      {
        return( x[nrow(x), variable] - sum(x[1:(nrow(x) -1),"weights"] * x[1:(nrow(x) - 1), 
                                                                           variable], na.rm = T ))
      }
      diffs <- sapply(tprd, get_mean_difs, variable = variable)
      
      var.points[[i]] <- mean(diffs / sd.val, na.rm = T)
    }
    names(var.points) <- covariates
    plotpoints[[k]] <- var.points
    
  }
  
  names(plotpoints) <- paste0("t_", lag:0)
  pointmatrix <- apply((as.matrix(do.call(rbind, plotpoints))), 2, function(x){(as.numeric(x))}, simplify = TRUE)
  rownames(pointmatrix) <- names(plotpoints)
  
  remove.vars.idx <- apply(apply(pointmatrix, 2, is.nan), 2, any)
  if(sum(remove.vars.idx) > 0)
  {
    removed.vars <- names(which(apply(apply(pointmatrix, 2, is.nan), 2, any)))
    pointmatrix <- pointmatrix[, !remove.vars.idx]
    warning(paste0("Some variables were removed due to low variation, inadequate data needed for calculation: ", removed.vars))
  }
  
  
  
  # we can remove time of treatment, since we expect a change
  
  pointmatrix <- pointmatrix[-nrow(pointmatrix), ,drop = FALSE]
  
  if (!plot) return(pointmatrix)
  
  if (plot)
  {
    # when binary, no real reason you'd want to calculate this
    treated.included <- treatment %in% colnames(pointmatrix)
    if (!continuous.treatment)
    {
      
      
      if (treated.included)
      {
        treated.data <- pointmatrix[,which(colnames(pointmatrix) == treatment)] # treated data
        pointmatrix <- pointmatrix[,-which(colnames(pointmatrix) == treatment)] #all non-treatment variable data
        graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), 
                          lty = 1, ylab = ylab, xaxt = "n", ...)
        graphics::lines(x = 1:nrow(pointmatrix), 
                        y = as.numeric(treated.data), 
                        type = "l",
                        lty = 2, lwd = 3)
        graphics::axis(side = 1, labels = paste0("t-", (nrow(pointmatrix)):1), 
                       at = 1:nrow(pointmatrix))  
      } else
      {
        graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), 
                          lty = 1, ylab = ylab, xaxt = "n", ...)
        graphics::axis(side = 1, labels = paste0("t-", (nrow(pointmatrix)):1), 
                       at = 1:nrow(pointmatrix))  
      }
      
    } else
    {
      if (treated.included)
      {
        treated.data <- pointmatrix[,which(colnames(pointmatrix) == treatment)] # treated data
        pointmatrix <- pointmatrix[,-which(colnames(pointmatrix) == treatment)] #all non-treatment variable data
        graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), 
                          lty = 1, ylab = ylab, xaxt = "n", ...)
        graphics::lines(x = 1:(nrow(pointmatrix) - 1), 
                        y = as.numeric(treated.data)[1:(nrow(pointmatrix) - 1)], type = "l",
                        lty = 2, lwd = 3)
        graphics::axis(side = 1, labels = paste0("t-", (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix))    
      } else {
        graphics::matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), 
                          lty = 1, ylab = ylab, xaxt = "n", ...)
        graphics::axis(side = 1, labels = paste0("t-", (nrow(pointmatrix)):1), at = 1:nrow(pointmatrix)) 
      }
      
    }
    
    if (legend) {
      if (treated.included)
      {
        legend("topleft", legend = c(colnames(pointmatrix), treatment), 
               col = c(1:ncol(pointmatrix), "black"), lty = c(rep(1, ncol(pointmatrix)), 2))  
      } else {
        legend("topleft", legend = colnames(pointmatrix), 
               col = 1:ncol(pointmatrix), lty = 1)
      }
      
    }
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
#' @param matched_set_list a list of one or more \code{matched.set} objects
#' @param xlim xlim of the scatter plot. This is the same as the \code{xlim} argument in \code{plot}
#' @param ylim ylim of the scatter plot. This is the same as the \code{ylim} argument in \code{plot}
#' @param main title of the scatter plot. This is the same as the \code{main} argument in \code{plot}
#' @param x.axis.label x axis label
#' @param y.axis.label y axis label
#' @param pchs one or more pch indicators for the symbols on the scatter plot. You should specify a phc symbol for each matched.set you specify in matched_set_list. See \code{plot} for more information
#' @param covariates variables for which balance is displayed
#' @param data the same time series cross sectional data set used to create the matched sets.
#' @param ... optional arguments to be passed to \code{plot}
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <amrauh@umich.edu>, and Kosuke Imai <imai@harvard.edu>
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
#'               matched_set_list = list(sets1$att, sets2$att),
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
balance_scatter <- function(matched_set_list,
                            xlim = c(0, .8),
                            ylim = c(0, .8),
                            main = "Standardized Mean Difference of Covariates",
                            pchs = c(2,3),
                            covariates, data,
                            x.axis.label = "Before refinement",
                            y.axis.label = "After refinement",
                            ...) {
  if (length(matched_set_list) < 1) stop("Please provide at least one matched.set object")
  # first, get balance for non-refined set
  # use first matched set and generate unrefined balance using get_covariate_balance
  non_refined_balance <- get_covariate_balance(matched.sets = matched_set_list[[1]],
                                               data = data,
                                               covariates = covariates,
                                               use.equal.weights = TRUE)
  
  # second, get balance for refined sets
  refined_balance <- list()
  for (i in 1:length(matched_set_list)) {
    refined_balance[[i]] <-
      get_covariate_balance(matched.sets = matched_set_list[[i]],
                            data = data,
                            covariates = covariates)
  }
  
  # extract values for x-axis from the non-refined sets
  benchmark <- as.vector(non_refined_balance)
  
  
  compared <- sapply(refined_balance, function(x) x <- x[1:(nrow(x)),])
  
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



calculate_set_effects <- function(pm.obj, data.in, lead)
{
  if (identical(attr(pm.obj, "qoi"), "att"))
  {
    msets <- pm.obj[["att"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
    
  } else if (identical(attr(pm.obj, "qoi"), "atc"))
  {
    msets <- pm.obj[["atc"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
  } else if (identical(attr(pm.obj, "qoi"), "ate"))
  {
    msets <- pm.obj[["att"]]
    msets.atc <- pm.obj[["atc"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
    
  } else if (identical(attr(pm.obj, "qoi"), "art"))
  {
    msets <- pm.obj[["art"]]
    id.var <- attributes(msets)$id.var
    t.var <- attributes(msets)$t.var
  } else {
    stop("invalid qoi")
  }
  
  rownames(data.in) <- paste0(data.in[, id.var], ".", data.in[, t.var])
  
  get_ind_effects <- function(mset, data_in, 
                              lead.val,
                              mset.name,
                              outcome, 
                              use.abs.value = FALSE,
                              is.atc = FALSE)
  {
    
    if ( identical(length(mset), 0L)) return(NA)
    t.val <- as.numeric(sub(".*\\.", "", mset.name))
    id.val <- as.numeric(sub("\\..*", "", mset.name))
    
    
    past.lookups <- paste0(mset, ".", (t.val - 1))
    future.lookups <- paste0(mset, ".", (t.val + lead.val))
    
    t.past.lookup <- paste0(id.val, ".", (t.val - 1))
    t.future.lookup <- paste0(id.val, ".", (t.val + lead.val))
    
    control.diffs <- data_in[future.lookups, outcome] - data_in[past.lookups, outcome]
    
    treat.diff <- data_in[t.future.lookup, outcome] - data_in[t.past.lookup, outcome]
    
    if (is.atc)
    {
      ind.effects <- sum(attr(mset, "weights") * control.diffs) - treat.diff
    } else {
      ind.effects <- treat.diff - sum(attr(mset, "weights") * control.diffs)
    }
    
    
    denom <- attr(mset, "treatment.change")
    if (use.abs.value) denom <- abs(denom)
    if (is.atc) denom <- 1
    ind.effects <- ind.effects / denom
    return(ind.effects)
    
  }
  
  
  if ( identical(attributes(pm.obj)[["qoi"]], "att"))
  { #using simplify = TRUE because we should always expect a vector, so nothing unexpected should happen
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]]),
                      SIMPLIFY = TRUE)
    
    return(effects)
  } else if (identical(attributes(pm.obj)[["qoi"]], "atc"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]],
                                      is.atc = TRUE),
                      SIMPLIFY = TRUE)
    
    return(effects)
  } else if (identical(attr(pm.obj, "qoi"), "ate"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]]),
                      SIMPLIFY = TRUE)
    
    
    effects.atc <- mapply(FUN = get_ind_effects,
                          mset = msets.atc,
                          mset.name = names(msets.atc),
                          MoreArgs = list(lead.val = lead,
                                          data_in = data.in,
                                          outcome = attributes(pm.obj)[["outcome.var"]]),
                          SIMPLIFY = TRUE)
    
    
    return(list(att = effects,
                atc = effects.atc))
    
  } else if (identical(attr(pm.obj, "qoi"), "art"))
  {
    effects <- mapply(FUN = get_ind_effects,
                      mset = msets,
                      mset.name = names(msets),
                      MoreArgs = list(lead.val = lead,
                                      data_in = data.in,
                                      outcome = attributes(pm.obj)[["outcome.var"]],
                                      use.abs.value = TRUE),
                      SIMPLIFY = TRUE)
    return(effects)
  } else
  {
    stop("invalid qoi")
  }
  
  
  
}



#' get_set_treatment_effects
#'
#' Calculates the treatment effect size at the matched set level
#'
#'
#' Calculate the size of treatment effects for each matched set.
#' @param pm.obj an object of class \code{PanelMatch}
#' @param data data.frame with the original data
#' @param lead integer (or integer vector) indicating the time period(s) in the future for which the treatment effect size will be calculated. Calculations will be made for the period t + lead, where t is the time of treatment. If more than one lead value is provided, then calculations will be performed for each value.
#'
#' @examples
#' 
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
#'                          placebo.test = FALSE)
#' set.effects <- get_set_treatment_effects(pm.obj = PM.results, data = dem, lead = 0)
#'
#'
#' @export

get_set_treatment_effects <- function(pm.obj, data, lead)
{
  return(lapply(lead, calculate_set_effects, pm.obj = pm.obj, data.in = data))

}
#' 
#' 
#' placebo_test
#'
#' Calculates results for a placebo test
#'
#'
#' Calculate the results of a placebo test, looking at the change in outcome at time = t-1, compared to other pre-treatment periods in the lag window.
#' @param pm.obj an object of class \code{PanelMatch}
#' @param data data.frame with the original data
#' @param lag.in integer indicating earliest the time period(s) in the future for which the placebo test change in outcome will be calculated. Calculations will be made over the period t - max(lag) to t-2, where t is the time of treatment. The results are similar to those returned by PanelEstimate, except t-1 is used as the period of comparison, rather than the lead window.
#' @param number.iterations integer specifying the number of bootstrap iterations
#' @param confidence.level confidence level for the calculated standard error intervals
#' @param plot logical indicating whether or not a plot should be generated, or just return the raw data from the calculations
#' @param ... extra arguments to be passed to plot
#'
#' @return list with 2 elements: "estimates", which contains the point estimates for the test and "bootstrapped.estimates", containing the bootstrapped point estimates for the test for each specified lag window period.
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2",
#'                          treatment = "dem", refinement.method = "mahalanobis",
#'                          data = dem, match.missing = TRUE,
#'                          covs.formula = ~ I(lag(tradewb, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE,
#'                          placebo.test = TRUE)
#' placebo_test(PM.results, data = dem, number.iterations = 100, plot = FALSE)
#' 
#' 
#' @export
#'
#'
placebo_test <- function(pm.obj,
                        data,
                        lag.in = NULL,
                        number.iterations = 1000,
                        confidence.level = .95,
                        plot = FALSE,
                        ...)
{
  df.adjustment <- FALSE
  qoi.in <- attr(pm.obj, "qoi")
  if (is.null(lag.in))
  {
    if (identical(qoi.in, "ate"))
    {
      #just pick one
      matchedsets <- pm.obj[["att"]]
      lag.in <- attr(matchedsets, "lag")
    } else {
      matchedsets <- pm.obj[[qoi.in]]
      lag.in <- attr(matchedsets, "lag")
    }
  }
  
  if (identical(qoi.in, "ate"))
  {
    #just pick one for boundary check
    matchedsets <- pm.obj[["att"]]
  }
  #matchedsets <- pm.obj[[qoi.in]]
  if (lag.in == 1) stop("placebo test cannot be conducted for lag = 1")
  if (length(lag.in) >1) stop("lag.in should be a single integer")
  if (lag.in > attr(matchedsets, "lag")) stop("provided lag.in value exceeds lag parameter from matching stage. Please specify a valid lag.in value, such that lag.in < lag")

  lag.in <- lag.in:2


  placebo.results.raw <- panel_estimate(sets = pm.obj,
                                        data = data,
                                        number.iterations = number.iterations,
                                        df.adjustment = df.adjustment,
                                        placebo.test = TRUE,
                                        placebo.lead = lag.in)

  if (plot)
  {
    plot(placebo.results.raw, ...)
  } else {
    colnames(placebo.results.raw$bootstrapped.estimates) <- names(placebo.results.raw$estimates)
    ret.results <- list(estimates = placebo.results.raw$estimates,
         bootstrapped.estimates = placebo.results.raw$bootstrapped.estimates)
    return(ret.results)
  }
}
