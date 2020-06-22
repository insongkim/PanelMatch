calculate_neighbor_treatment <- function(data, edge.matrix, n.degree, unit.id, time.id, treatment.variable)
{
  g1 <- igraph::graph_from_data_frame(edge.matrix, directed = FALSE)
  ref.names <- paste0(data[, unit.id], ".", data[, time.id])
  treatment.vector <- ordered.data[, treatment.variable]
  names(treatment.vector) <- ref.names
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
    ret.list <- lapply(ret.list, as.integer)
    names(ret.list) <- as.character(igraph::V(graph.in))
    return(ret.list)
  }
  
  neighborhood.lookup <- lapply(1:n.degree, adjusted.neighborhood, graph.in = g1)
  
  get.neighborhood.treatment.per.time <- function(treatment.lookup, neighborhood.vector, time, return.average = TRUE)
  {
    
    lookups <- paste0(neighborhood.vector, '.', time)
    if(return.average) ## feel like these might need to be updated to replace na's with zeroes or something
    {
      prop.treatment <- mean(treatment.lookup[lookups], na.rm = TRUE)  
      return(prop.treatment)
    }
    else {
      prop.treatment <- sum(treatment.lookup[lookups], na.rm = TRUE)
      return(prop.treatment)
    }
    
  }
  
  get.treatment.prop.per.row <- function(t.id.pair, degree, neighborhood.lookup, treatment.vector, return.average = TRUE)
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
      return.proportion <- get.neighborhood.treatment.per.time(treatment.vector, neighbor.res, t, return.average)  
    }
    
    return(return.proportion)
  }
  
  ll <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                     neighborhood.lookup = neighborhood.lookup, 
                                                     treatment.vector = treatment.vector, degree = x)})
  data[, make.names(paste0('neighborhood_t_prop', '.', 1:n.degree))] = ll
  
  
  ll <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                     neighborhood.lookup = neighborhood.lookup, 
                                                     treatment.vector = treatment.vector, degree = x,
                                                     return.average = FALSE)})
  
  data[, make.names(paste0('neighborhood_t_count', '.', 1:n.degree))] = ll
  
  return(data)
}

handle_calipers <- function(plain.ordered.data, caliper.formula, matched.sets, lag.window)
{
  
  caliper <- function(y, method, caliper.distance, data_type, lwindow = lag.window)
  {# redefine the function to extract the method types, variables, caliper distances, might be a better way of doing this
    y <- deparse(substitute(y))
    return(c(y, method, caliper.distance, max(lwindow), data_type))
  }
  
  attr(caliper.formula, ".Environment") <- environment()
  caliper.metadata <- model.frame(caliper.formula)
  rownames(caliper.metadata) <- NULL
  
  internal.caliper <- function (x, n = 1L, default = NA)
  {
    if(class(x) == "factor")
    {
      x <- as.numeric(as.character(x))
    }
    if (n == 0) return(x)
    xlen <- length(x)
    n <- pmin(n, xlen)
    out <- c(rep(default, n), x[seq_len(xlen - n)])
    attributes(out) <- attributes(x)
    out
  }
  
  caliper <- function(y, method, caliper.distance, data_type, lwindow = lag.window)
  {
    
    if(is.null(method) || is.null(caliper.distance))
    {
      stop("arguments missing from caliper function")
    }
    sapply(lwindow, internal.caliper, x = y)
  }
  
  
  apply_formula <- function(x, form)
  {
    
    attr(form, ".Environment") <- environment()
    tdf <- model.frame(form, x, na.action = NULL)
    cbind(x[, c(1, 2, 3)], model.matrix(form, tdf)[, -1])
  }
  

  t.data <- do.call(rbind, by(plain.ordered.data, as.factor(plain.ordered.data[, 1]), FUN = apply_formula, form = caliper.formula))
  #may not be necessary?
  
  t.data <- t.data[order(t.data[,1], t.data[,2]), ]
  rownames(t.data) <- NULL
  

  
  msets <- matched.sets #initialize the matched sets
  for(i in 1:ncol(caliper.metadata))
  {
    col.idx <- i * max(lag.window + 1) + 3
    
    .IS_FACTOR_VAR <- caliper.metadata[5, i] == "categorical"
    msets <- handle.single.caliper.per.lag(t.data, msets, caliper.metadata[2,i], caliper.metadata[3, i], c(1:3,  (col.idx-max(lag.window)):col.idx), .IS_FACTOR_VAR)
  }
  
  return(msets)
  
}




handle.single.caliper.per.lag <- function(plain.ordered.data, matched.sets, caliper.method, caliper.distance, data.index, is.factor.var)
{
  
  time.var <- attr(matched.sets, "t.var")
  id.var <- attr(matched.sets, "id.var")
  treat.var <- attr(matched.sets, "treatment.var")
  lag.in <- attr(matched.sets, 'lag')
  
  old.lag <- lag.in
  lag.in <- 0 
  
  treated.ts <- as.numeric(unlist(strsplit(names(matched.sets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(matched.sets), split = "[.]"))[c(T,F)])
  
  ordered.data <- plain.ordered.data[, data.index]
  
  tlist <- expand.treated.ts(lag.in, treated.ts = treated.ts)
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, tlist, paste0(ordered.data[,id.var], ".",
                                                                       ordered.data[, time.var]), matched_sets = matched.sets, lag.in)
  mahalmats <- build_maha_mats(ordered_expanded_data = as.matrix(ordered.data), idx =  idxlist)
  msets <- handle_perlag_caliper_calculations(mahalmats, matched.sets, 
                                              caliper.value = caliper.distance, caliper.method = caliper.method, is.factor.var)
  lag.in <- old.lag
  
  return(matched.sets)
  
}

handle_perlag_caliper_calculations <- function(nested.list, msets, caliper.method, caliper.value, is.factor.var)
{
  
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
    results.temp <- lapply(sub.list, do_calcs, sd.vals__ = standard.deviations, is_factor_in = is.factor)
    tmat <- do.call(rbind, results.temp)
    colnames(tmat) <- NULL
    # tmat <- tmat[, 1:(ncol(tmat)-1), drop = FALSE]
    meets_caliper <- function(col, cal, cal.method, sd.vals, IS_FACTOR)
    {
      
      if(IS_FACTOR)
      {
        if(cal.method == "max")
        {
          if(any(is.na(col)) || any(col != 0)) #everything must match
          {
            return(FALSE)
          } else 
          {
            return(TRUE)
          }  
        }
        if(cal.method == "average")
        {
          m.val <- mean(col != 0)
          if(is.na(m.val) || m.val < cal) #looking at proportion of matches, so higher is more similar
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
  tdf <- do.call(rbind, lapply(unlist(nested.list, recursive = F),
                               function(x) return(x[nrow(x), ])))
  sd.vals <- apply(tdf[, 4:ncol(tdf)], MARGIN = 2, FUN = sd, na.rm = TRUE)
  
  indices.msets <- lapply(nested.list, handle_set, 
                          cal.value = caliper.value, 
                          cal.method = caliper.method, 
                          standard.deviations = sd.vals, 
                          is.factor = is.factor.var)
  msets <- mapply(function(x, y) return(x[y]), x = msets, y = indices.msets)
  msets <- msets[sapply(msets, length) > 0]
  return(msets)
}

