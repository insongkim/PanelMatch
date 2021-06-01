handle_network_covariates <- function(data, edge.matrix, n.degree,
                                      unit.id, time.id, covariates)
{
  
  raw.cov.data <- lapply(covariates, handle_network_covariate, data = data,
                         edge.matrix = edge.matrix, n.degree = n.degree, 
                         unit.id = unit.id, time.id = time.id)
  col.names <- lapply(lapply(covariates, function(x) paste0(x, ".", 1:n.degree)), function(x) paste0(x, ".", c("avg.", "change.")))
  ## stitch back into the main data
  
  col.names <- unlist(col.names)
  data[, col.names] <- unlist(unlist(raw.cov.data, 
                                     recursive = FALSE), 
                              recursive = FALSE)
  
  return(data)
}


# should be more generalized version of the calculate_neighbor_treatment function that currently exists
handle_network_covariate <- function(data, edge.matrix, n.degree, 
                                     unit.id, time.id, covariate)
{
  
  edge.matrix <- edge.matrix[ edge.matrix[,3] == 1, ] #probably want to change how this is done
  g1 <- igraph::graph_from_data_frame(edge.matrix, directed = FALSE)
  ref.names <- paste0(data[, unit.id], ".", data[, time.id])
  treatment.vector <- data[, covariate]
  names(treatment.vector) <- ref.names
  
  
  t1 <- data[, covariate]
  t1 <- c(NA, t1[-nrow(data)])
  t1[which(!duplicated(data[, unit.id]))] <- NA
  names(t1) <- ref.names
  t1 <- treatment.vector - t1
  
  adjusted.neighborhood <- function(graph.in, degree)
  {
    straight.neighborhoods <- igraph::ego(graph.in, order = degree)
    
    ret.list <- straight.neighborhoods
    ret.list <- lapply(straight.neighborhoods, function(x){return(igraph::as_ids(x)[-1])})
    names(ret.list) <- igraph::as_ids(igraph::V(graph.in))
    return(ret.list)
  }
  
  neighborhood.lookup <- lapply(1:n.degree, adjusted.neighborhood, graph.in = g1)
  
  get.neighborhood.treatment.per.time <- function(treatment.lookup, neighborhood.vector, 
                                                  time, return.average = TRUE,
                                                  diff.vector)
  {
    
    lookups <- paste0(neighborhood.vector, '.', time)
    
    
    if (return.average) 
    {
      r.vector <- treatment.lookup[lookups]
      prop.treatment <- mean(r.vector)  
      return(prop.treatment)
    }
    else {
      
      sum.changes <- sum(diff.vector[lookups])
      return(sum.changes)
    }
    
  }
  
  get.treatment.prop.per.row <- function(t.id.pair, degree, neighborhood.lookup,
                                         treatment.vector, return.average = TRUE, 
                                         diff.vector)
  {
    
    id <- as.numeric(unlist(strsplit(t.id.pair, split = "[.]"))[c(T,F)])
    t <- as.numeric(unlist(strsplit(t.id.pair, split = "[.]"))[c(F,T)])
    neighbor.res <- neighborhood.lookup[[degree]][[as.character(id)]]
    if(is.null(neighbor.res) | length(neighbor.res) == 0)
    {
      return.proportion <- NA  
    }
    else
    {
      return.proportion <- get.neighborhood.treatment.per.time(treatment.vector, 
                                                               neighbor.res, t, return.average,
                                                               diff.vector)  
    }
    
    return(return.proportion)
  }
  
  avg.data <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                           neighborhood.lookup = neighborhood.lookup, 
                                                           treatment.vector = treatment.vector, degree = x,
                                                           return.average = TRUE,
                                                           diff.vector = t1)})
  
  
  change.data <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                              neighborhood.lookup = neighborhood.lookup, 
                                                              treatment.vector = treatment.vector, degree = x,
                                                              return.average = FALSE,
                                                              diff.vector = t1)})
  
  return(list(average = avg.data, change = change.data))
  
  
}

calculate_neighbor_treatment <- function(data, edge.matrix, n.degree, 
                                         unit.id, time.id, treatment.variable)
{
  
  edge.matrix <- edge.matrix[ edge.matrix[,3] == 1, ] #probably want to change how this is done
  g1 <- igraph::graph_from_data_frame(edge.matrix, directed = FALSE)
  ref.names <- paste0(data[, unit.id], ".", data[, time.id])
  treatment.vector <- data[, treatment.variable]
  names(treatment.vector) <- ref.names
  
  #df$Lagged_Variable <- c(NA, df$Rate[-nrow(df)])
  #df$Lagged_Variable[which(!duplicated(df$Factor))] <- NA
  
  t1 <- data[, treatment.variable]
  t1 <- c(NA, t1[-nrow(data)])
  t1[which(!duplicated(data[, unit.id]))] <- NA
  names(t1) <- ref.names
  t1 <- treatment.vector - t1
  # update this so that it uses <= degree only, rather than separate adjusted ones
  adjusted.neighborhood <- function(graph.in, degree)
  {
    straight.neighborhoods <- igraph::ego(graph.in, order = degree)
    # if(degree > 1)
    # {
    #   t.neighborhood <- igraph::ego(graph.in, order = (degree - 1))
    #   ret.list <- mapply(difference, big = straight.neighborhoods, small = t.neighborhood, SIMPLIFY = FALSE)
    # }
    # else
    # {
    #   ret.list <- straight.neighborhoods
    # }
    ret.list <- straight.neighborhoods
    ret.list <- lapply(straight.neighborhoods, function(x){return(igraph::as_ids(x)[-1])})
    names(ret.list) <- igraph::as_ids(igraph::V(graph.in))
    return(ret.list)
  }

  neighborhood.lookup <- lapply(1:n.degree, adjusted.neighborhood, graph.in = g1)
  # could use previous line to get a count of number of neighbors for each unit -- normalizing the amount of change to be per unit
  get.neighborhood.treatment.per.time <- function(treatment.lookup, neighborhood.vector, 
                                                  time, return.average = TRUE,
                                                  diff.vector)
  {
    ## designed for binary treatment, will need to update the following chunks to accomodate continuous treatment
    lookups <- paste0(neighborhood.vector, '.', time)
    
    #r.vector <- lookups %in% treated.unit.names
    if (return.average) ## feel like these might need to be updated to replace na's with zeroes or something
    {
      r.vector <- treatment.lookup[lookups]
      prop.treatment <- mean(r.vector)  
      return(prop.treatment)
    }
    else {
      
      sum.changes <- sum(diff.vector[lookups])
      #if we are normalizing over neighborhood size, can do that here
      #prop.treatment <- sum(r.vector, na.rm = TRUE)
      return(sum.changes)
    }
    
  }
  
  get.treatment.prop.per.row <- function(t.id.pair, degree, neighborhood.lookup,
                                         treatment.vector, return.average = TRUE, 
                                         diff.vector)
  {
    
    id <- as.numeric(unlist(strsplit(t.id.pair, split = "[.]"))[c(T,F)])
    t <- as.numeric(unlist(strsplit(t.id.pair, split = "[.]"))[c(F,T)])
    neighbor.res <- neighborhood.lookup[[degree]][[as.character(id)]]
    if(is.null(neighbor.res) | length(neighbor.res) == 0)
    {
      return.proportion <- NA  
    }
    else
    {
      return.proportion <- get.neighborhood.treatment.per.time(treatment.vector, 
                                                               neighbor.res, t, return.average,
                                                               diff.vector)  
    }
    
    return(return.proportion)
  }
  
  ll <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                     neighborhood.lookup = neighborhood.lookup, 
                                                     treatment.vector = treatment.vector, degree = x,
                                                     return.average = TRUE,
                                                     diff.vector = t1)})
  
  data[, make.names(paste0('neighborhood_t_prop', '.', 1:n.degree))] = ll
  
  
  ll <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                     neighborhood.lookup = neighborhood.lookup, 
                                                     treatment.vector = treatment.vector, degree = x,
                                                     return.average = FALSE,
                                                     diff.vector = t1)})
  
  data[, make.names(paste0('neighborhood_t_count', '.', 1:n.degree))] = ll
  rownames(data) <- NULL
  return(data)
}

handle_calipers <- function(plain.ordered.data, caliper.formula, 
                            matched.sets, lag.window, 
                            is.continuous.matching = FALSE,
                            control.threshold = NULL)
{
  
  caliper <- function(y, method, caliper.distance, data_type, 
                      calculation_type, lwindow = lag.window)
  {# redefine the function to extract the method types, variables, caliper distances, might be a better way of doing this
    
    y <- deparse(substitute(y))
    return(c(y, method, caliper.distance, max(lwindow), data_type, calculation_type))
  }
  
  attr(caliper.formula, ".Environment") <- environment()
  caliper.metadata <- model.frame(caliper.formula)
  rownames(caliper.metadata) <- NULL
  
  internal.caliper <- function(x, n = 1L, default = NA)
  {
    if (class(x) == "factor")
    {
      x <- as.numeric(as.character(x))
    }
  
    if (n == 0) return(x)
    xlen <- length(x)
    n <- pmin(n, xlen)
    out <- c(rep(default, n), x[seq_len(xlen - n)])
    attributes(out) <- attributes(x)
    
    return(out)
  }
  
  caliper <- function(y, method, caliper.distance, data_type, calculation_type, lwindow = lag.window)
  {
    
    if(is.null(method) || is.null(caliper.distance))
    {
      stop("arguments missing from caliper function")
    }
    
    return(sapply(lwindow, internal.caliper, x = y))
  }
  
  
  apply_formula <- function(x, form)
  {
  
    attr(form, ".Environment") <- environment()
    tdf <- as.data.frame(as.matrix(model.frame(form, x, na.action = NULL)))
    
    return(cbind(x[, c(1, 2, 3)], tdf))
  }
  
  t.data <- do.call(rbind, by(plain.ordered.data, as.factor(plain.ordered.data[, 1]),
                              FUN = apply_formula, form = caliper.formula))
  #may not be necessary?
  
  t.data <- t.data[order(t.data[,1], t.data[,2]), ]
  rownames(t.data) <- NULL
  
  for(i in 1:ncol(caliper.metadata))
  {
    if (is.continuous.matching)
    {
      col.idx <- i * max(lag.window) + 3
      data.index <- (col.idx - (max(lag.window) - 1)):col.idx
    } else
    {
      col.idx <- i * max(lag.window + 1) + 3
      data.index <- (col.idx - max(lag.window)):col.idx
    }
    #col.idx <- i * max(lag.window + 1) + 3
    
    .IS_FACTOR_VAR <- caliper.metadata[5, i] == "categorical" # or numeric
    .USE_SD_UNITS <- caliper.metadata[6, i] == "sd" # or raw
    
    matched.sets <- handle.single.caliper.per.lag(t.data, matched.sets, caliper.metadata[2,i], 
                                                  as.numeric(caliper.metadata[3, i]), 
                                                  c(1:3, data.index),
                                                  #c(1:3,  (col.idx - max(lag.window)):col.idx),
                                                  .IS_FACTOR_VAR, .USE_SD_UNITS, 
                                                  is.continuous.matching, control.threshold)
  }
  
  return(matched.sets)
  
}

handle.single.caliper.per.lag <- function(plain.ordered.data, matched.sets, caliper.method, caliper.distance, 
                                          data.index, is.factor.var, use.sd.units, do.continuous.matching = FALSE,
                                          control.threshold = NULL)
{
  
  time.var <- attr(matched.sets, "t.var")
  id.var <- attr(matched.sets, "id.var")
  treat.var <- attr(matched.sets, "treatment.var")
  lag.in <- attr(matched.sets, 'lag')
  
  old.lag <- lag.in
  lag.in <- 0 
  #use more efficient regex version
  #treated.ts <- as.numeric(sub(".*\\.", "", names(matched.sets)))
  #treated.ids <- as.numeric(sub("\\..*", "", names(matched.sets)))
  
  ordered.data <- plain.ordered.data[, data.index]
  ##TODO: change everything to numeric matrix? 

  msets <- handle_distance_matrices(as.matrix(ordered.data), matched.sets, caliper.distance,
                           caliper.method, is.factor.var, use.sd.units, id.var, 
                           time.var, lag.in, do.continuous.matching, control.threshold)

  
  lag.in <- old.lag
  
  class(msets) <- c("matched.set", "list")
  attr(msets, 'refinement.method') <- NULL
  attr(msets, "lag") <- lag.in
  attr(msets, "t.var") <- time.var
  attr(msets, "id.var" ) <- id.var
  attr(msets, "treatment.var") <- treat.var
  return(msets)
  
}

handle_perlag_caliper_calculations <- function(nested.list, msets, caliper.method, 
                                               caliper.value, is.factor.var, 
                                               use.sd.units, full.data, 
                                               id.var, time.var, treated_it_pairs)
{
  
  if (length(msets) == 0) return(numeric(0))
  do_calcs <- function(time.df, sd.vals__, is_factor_in)
  {
    
    cov.data <- time.df[, 4:ncol(time.df), drop = FALSE]
    
    do.maths  <- function(control.unit.data, treated.unit.data, sd.vals_, is_factor)
    {
      
      if(!is_factor)
      {
        return(abs(treated.unit.data - control.unit.data) / sd.vals_) 
      } else
      { 
        return(abs(treated.unit.data - control.unit.data))  
      }
    }
    
    results <- apply(cov.data[1:((nrow(cov.data)  - 1)), , drop = FALSE] , MARGIN = 1, FUN  =  do.maths, 
          treated.unit.data = cov.data[nrow(cov.data), ], 
          sd.vals_ = sd.vals__, 
          is_factor = is_factor_in)
    
    return(results)
    
  }
  
  handle_set <- function(sub.list, cal.value, cal.method, standard.deviations, is.factor)
  {
    
    # results.temp <- lapply(sub.list, do_calcs, sd.vals__ = standard.deviations, 
    #                        is_factor_in = is.factor)
    tmat <- do_calcs(sub.list[[1]], sd.vals__ = standard.deviations, is_factor_in = is.factor)
    # tmat <- do.call(rbind, results.temp) 
    ## two things commented out for more efficient refinement implemented 
    colnames(tmat) <- NULL
    # tmat <- tmat[, 1:(ncol(tmat)-1), drop = FALSE]
    meets_caliper <- function(col, cal, cal.method, sd.vals, IS_FACTOR)
    {
      # browser()
      # print(cal)
      if (IS_FACTOR)
      {
        if (cal.method == "max")
        {
          if (any(is.na(col)) || any(col != 0)) #everything must match
          {
            return(FALSE)
          } else 
          {
            return(TRUE)
          }  
        }
        if (cal.method == "average")
        {
          m.val <- mean(col == 0)
          
          if (is.na(m.val) || m.val < cal) #looking at proportion of matches, so higher is more similar
          {
            return(FALSE)
          } else 
          {
            return(TRUE)
          }
        }
      } else
      {
        if(cal.method == "max")
        {
          if(any(is.na(col)) || any(col > cal))
          {
            return(FALSE)
          } else 
          {
            return(TRUE)
          }  
        }
        if(cal.method == "average")
        {
          m.val <- mean(col)
          if(is.na(m.val) || m.val > cal)
          {
            return(FALSE)
          } else 
          {
            return(TRUE)
          }
        }
        else {
          stop("please ensure data for caliper calculations is either a factor, numeric, or integer")
        }
      }
      
      
    }
    unit.index <- apply(tmat, 2, FUN= meets_caliper, 
                        cal = cal.value, 
                        cal.method = cal.method, 
                        sd.vals = standard.deviations,
                        IS_FACTOR = is.factor)
    return(unit.index)
    
  }
  
  #get standard deviation values
  
  row.names(full.data) <- paste0(full.data[, id.var], ".", full.data[, time.var])
  
  tdf <- full.data[treated_it_pairs, ]
  # think this needs to get the full set of treated i,t pairs
  
  if(use.sd.units)
  {
    sd.vals <- apply(tdf[, 4:ncol(tdf)], MARGIN = 2, FUN = sd, na.rm = TRUE)
    if ( 0 %in% sd.vals || any(is.na(sd.vals)) )
    {
      stop("not enough variation for using 'units = sd' in continuous.treatment.info. 
           Try using raw units instead")
    }
  } else {
    sd.vals <- 1
  }
  
  indices.msets <- handle_set(nested.list, caliper.value, caliper.method, sd.vals, is.factor.var)
  #msets <- mapply(function(x, y) return(as.numeric(x[y])), x = msets, y = indices.msets, SIMPLIFY = FALSE)
  #msets <- msets[sapply(msets, length) > 0]
  
  msets <- msets[indices.msets]
  return(msets)
}


handle_network_caliper_and_refinement <- function(network.caliper.info = NULL, network.refinement.info = NULL,
                                                  ordered.data, adjacency.matrix, neighborhood.degree, unit.id, time.id,
                                                  treatment, covs.formula, caliper.formula)
{
  
  if (!is.null(network.caliper.info) || !is.null(network.refinement.info))
  { 
    propstring.formula <- paste0('neighborhood_t_prop', '.', neighborhood.degree)
    countstring.formula <- paste0('neighborhood_t_count', '.', neighborhood.degree)
    if (!is.null(network.refinement.info))
    {
      
      if (network.refinement.info[["use.proportion.data"]])
      {
        if (identical(max(network.refinement.info[['proportion.lags']]), 0))
        {
          if (is.null(covs.formula))
          {
            covs.formula <- reformulate(propstring.formula)
            environment(covs.formula) <- globalenv()
          } else {
            covs.formula <- merge_formula(covs.formula, reformulate(propstring.formula))
          }
          
        }
        else if (max(network.refinement.info[["proportion.lags"]]) > 0)
        {
          propstring.formula.lag <- paste0("I(lag(", propstring.formula, ",", 
                                           deparse(network.refinement.info[["proportion.lags"]]), "))")
          if (is.null(covs.formula))
          {
            covs.formula <- reformulate(propstring.formula.lag)
            environment(covs.formula) <- globalenv()
          } else {
            covs.formula <- merge_formula(covs.formula, reformulate(propstring.formula.lag))
          }
          
        }
        else {
          stop("please enter a valid lag number")
        }
      }
      if(network.refinement.info[["use.count.data"]])
      {
        if(identical(max(network.refinement.info[['count.lags']]), 0))
        {
          if(is.null(covs.formula))
          {
            covs.formula <- reformulate(countstring.formula)
            environment(covs.formula) <- globalenv()
          } else {
            covs.formula <- merge_formula(covs.formula, reformulate(countstring.formula))
          }
          
        }
        else if (max(network.refinement.info[["count.lags"]]) > 0)
        {
          countstring.formula.lag <-paste0("I(lag(", countstring.formula, ",", 
                                           deparse(network.refinement.info[["count.lags"]]), "))")
          if(is.null(covs.formula))
          {
            covs.formula <- reformulate(countstring.formula.lag)
            environment(covs.formula) <- globalenv()
          } else {
            covs.formula <- merge_formula(covs.formula, reformulate(countstring.formula.lag))
          }
          
        }
        else {
          stop("please enter a valid lag number")
        }
      }
      
    }
    if(!is.null(network.caliper.info))
    {
      if(network.caliper.info[["use.proportion.data"]])
      {
        
        propstring.formula.caliper <- paste0("I(caliper(", propstring.formula, ",", 
                                             "'", network.caliper.info[["proportion.caliper.method"]], "'", ",",
                                             network.caliper.info[["proportion.caliper.threshold"]], ",",
                                             '"numeric"', "," ,"'", network.caliper.info[["prop.unit.type"]],"'", 
                                             "))")
        if(is.null(caliper.formula))
        {
          caliper.formula <- reformulate(propstring.formula.caliper)
          environment(caliper.formula) <- environment(covs.formula)
        } else
        {
          caliper.formula <- merge_formula(caliper.formula, reformulate(propstring.formula.caliper))  
        }
        
      }
      if(network.caliper.info[["use.count.data"]])
      {
        countstring.formula.caliper <- paste0("I(caliper(", countstring.formula, ",",
                                              "'", network.caliper.info[["count.caliper.method"]], "'", ",",
                                              network.caliper.info[["count.caliper.threshold"]], ",",
                                              '"numeric"', ",", "'", network.caliper.info[["count.unit.type"]], "'",
                                              "))")
        if(is.null(caliper.formula))
        {
          caliper.formula <- reformulate(countstring.formula.caliper)
          environment(caliper.formula) <- environment(covs.formula)
        } else
        {
          caliper.formula <- merge_formula(caliper.formula, reformulate(countstring.formula.caliper))  
        }
        
      }
    }
    
  }
  return(list(covs.formula, caliper.formula))
  
}


