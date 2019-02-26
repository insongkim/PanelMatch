## anything related to matched set objects, helper functions, etc. should go in this file
#' @export
matched_set <- function(matchedsets, id, t, L, t.var, id.var, treated.var)
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
  attr(matchedsets, "treated.var") <- treated.var
  return(matchedsets)
}

#' @export
summary.matched.set <- function(set, verbose = T)
{
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
    summary.result$set.size.summary <- summary(Lengthcol)
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
#' @export
plot.matched.set <- function(set, border = NA, col = "grey", ylab = "Frequency of Size", 
                             xlab ="Matched Set Size" , lwd = NULL,
                             main = "Distribution of matched set sizes",
                             ...)
{
    lvec <- sapply(set, length)
    
    hist(x = lvec, freq = TRUE, border = border, col = col, ylab = ylab, xlab = xlab, main = main, ...)
    if(sum(lvec == 0) > 0)
    {
      if(is.null(lwd))
      {
        lwd = 4
      }
      lines(x = c(0,0), y = c(0, length(rep(0, sum(lvec == 0) )) ), col = "red", lwd = lwd)  
    }
}
#' @export
extract.set <- function(set, id, t)
{
  strid <- paste0(id, ".", t)
  subset <- set[names(set) == strid]
  if(length(subset) != 1) stop('t,id pair invalid')
  class(subset) <- "matched.set"
  return(subset)
}

#' @export
print.matched.set <- function(set, verbose = F)
{
  if(verbose) 
  {
    class(set) <- "list"
    print(set)
  }
  
  else {
    print(summary(set, verbose = F))
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
  attr(temp, "treated.var") <- attr(x, "treated.var")
  attr(temp, "distances") <- attr(x, "distances")
  attr(temp, "max.match.size") <- attr(x, "max.match.size")
  attr(temp, "covs.formula") <- attr(x, "covs.formula")
  attr(temp, "match.missing") <- attr(x, "match.missing")
  class(temp) <- "matched.set"
  
  return(temp)
}

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


#' Calculate covariate balance for specified covariates across matched sets. Balance is assessed by taking the difference between 
#' the values of the user specified covariates in the treated unit and the weighted average of that across all matched sets. Furthermore, results are standardized and are expressed in standard deviations. 
#' @param matched.sets A matched.set object
#' @param data The data set used to produce the matched.set object. Please make sure this data set is identical to the one passed to PanelMatch/PanelEstimate to ensure consistent results
#' @param verbose When TRUE, the function will return more information about the calculations/results. When FALSE, a more compact version of the results/calculations are returned.
#' @param plot When TRUE, a plot showing the covariate balance calculation results will be shown. When FALSE, no plot is made, but the results of the calculations are still returned. 
#' @param covariates a character vector, specifying the names of the covariates for which the user is interested in calculating balance. 
#' @param reference.line 
#' @param legend
#' @param ... 
#' @export
get_covariate_balance <- function(matched.sets, data,  covariates, verbose = T, plot = F, reference.line = TRUE, legend = TRUE, ylab = "SD",...)
{
  if(is.null(covariates))
  {
    stop("please specify the covariates for which you would like to check the balance")
  }
  unit.id <- attr(matched.sets, "id.var")
  time.id <- attr(matched.sets, "t.var")
  lag <- attr(matched.sets, "lag")
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  ordered.data <- data[order(data[,unit.id], data[,time.id]), ]
  ordered.data[, paste0(unit.id, ".int")] <- as.integer(as.factor(data[, unit.id]))
  
  if(class(ordered.data[, unit.id]) == "character") {
    unit.index.map <- data.frame(original.id = make.names(as.character(unique(ordered.data[, unit.id]))), new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(ordered.data[, unit.id]) == "integer") {
    unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))), new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(ordered.data[, unit.id]) == "numeric") {
    if(all(unique(ordered.data[, unit.id]) == as.integer(unique(ordered.data[, unit.id])))) #actually integers
    {
      unit.index.map <- data.frame(original.id = (as.character(unique(ordered.data[, unit.id]))), new.id = unique(ordered.data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
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
  matched.sets <- encode_index(matched.sets, unit.index.map, unit.id)
  #they will either all have or not have weights, so we can check the first matched set to see if we need to add equal weighting
  if(is.null(attr(matched.sets[[1]], "weights")))
  {
    for(i in 1:length(matched.sets))
    {
      attr(matched.sets[[i]], "weights") <- rep(1/length(matched.sets[[i]]), length(matched.sets[[i]]))
      names(matched.sets[[i]], "weights") <- matched.sets[[i]]
    }
  }
  treated.ts <- as.integer(unlist(strsplit(names(matched.sets), split = "[.]"))[c(F,T)])
  treated.ids <- as.integer(unlist(strsplit(names(matched.sets), split = "[.]"))[c(T,F)])
  tlist <- expand.treated.ts(lag, treated.ts = treated.ts)
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, tlist, paste0(ordered.data[,unit.id], ".", 
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
      
      sd.val <- sd(sapply(unlistedmats[seq(from = 1, to = (length(matched.sets) * (lag + 1)), by = lag + 1)], function(x){x[nrow(x), variable]}), na.rm = T)
      
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
  if(!plot) return(pointmatrix)
  
  if(plot)
  {
    matplot(pointmatrix, type = "l", col = 1:ncol(pointmatrix), lty =1, ylab = ylab, ...)
    if(legend) legend("topleft", legend = colnames(pointmatrix), col=1:ncol(pointmatrix), lty = 1)
    if(reference.line) abline(h = 0, lty = "dashed")
  }
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
    #n.ts <- time.index$original.time.id[match(ts, time.index$new.time.id)]
    return(paste0(n.ids, ".", ts))
  }
  treated.ts <- as.numeric(unlist(strsplit(names(mset), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(mset), split = "[.]"))[c(T,F)])
  attributes(new.mset) <- attributes(mset)
  names(new.mset) <- decode.treated.units(treated.ts, treated.ids, unit.index)
  attr(new.mset, "id.var") <- original.unit.id
  #attr(new.mset, "t.var") <- original.time.id
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
    #n.ts <- time.index$new.time.id[match(ts, time.index$original.time.id)]
    return(paste0(n.ids, ".", ts))
  }
  treated.ts <- (unlist(strsplit(names(mset), split = "[.]"))[c(F,T)])
  treated.ids <- (unlist(strsplit(names(mset), split = "[.]"))[c(T,F)])
  attributes(new.mset) <- attributes(mset)
  names(new.mset) <- encode.treated.units(treated.ts, treated.ids, unit.index)
  attr(new.mset, "id.var") <- new.unit.id
  #attr(new.mset, "t.var") <- new.time.id
  return(new.mset)
}



