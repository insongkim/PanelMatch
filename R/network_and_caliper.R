calculate_neighbor_treatment <- function(data, edge.matrix, n.degree, unit.id, time.id, treatment.variable)
{
  g1 <- graph_from_data_frame(edge.matrix, directed = FALSE)
  
  ref.names <- paste0(data[, unit.id], ".", data[, time.id])
  treatment.vector <- ordered.data[, treatment.variable]
  names(treatment.vector) <- ref.names
  # update this so that it uses <= degree only, rather than separate adjusted ones
  adjusted.neighborhood <- function(graph.in, degree)
  {
    straight.neighborhoods <- ego(graph.in, order = degree)
    if(degree > 1)
    {
      t.neighborhood <- ego(graph.in, order = (degree - 1))
      ret.list <- mapply(difference, big = straight.neighborhoods, small = t.neighborhood, SIMPLIFY = FALSE)
    }
    else
    {
      ret.list <- straight.neighborhoods
    }
    ret.list <- lapply(ret.list, as.integer)
    names(ret.list) <- as.character(V(graph.in))
    return(ret.list)
  }
  
  neighborhood.lookup <- lapply(1:n.degree, adjusted.neighborhood, graph.in = g1)
  
  get.neighborhood.treatment.per.time <- function(treatment.lookup, neighborhood.vector, time)
  {
    
    lookups <- paste0(neighborhood.vector, '.', time)
    prop.treatment <- mean(treatment.lookup[lookups], na.rm = TRUE)
    return(prop.treatment)
  }
  
  get.treatment.prop.per.row <- function(t.id.pair, degree, neighborhood.lookup, treatment.vector)
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
      return.proportion <- get.neighborhood.treatment.per.time(treatment.vector, neighbor.res, t)  
    }
    
    return(return.proportion)
  }
  ll <- lapply(1:n.degree, FUN = function(x) {sapply(ref.names, get.treatment.prop.per.row, 
                                                     neighborhood.lookup = neighborhood.lookup, treatment.vector = treatment.vector, degree = x)})
  data[, make.names(paste0('neighborhood_t_prop', '.', 1:n.degree))] = ll
  return(data)
}


handle.single.caliper.per.lag <- function(plain.ordered.data, matched.sets, caliper.variable, caliper.distance, time.var, id.var, treat.var, lag.in)
{
  
  # time.var <- attr(matched.sets, "t.var")
  # id.var <- attr(matched.sets, "id.var")
  # treat.var <- attr(matched.sets, "treatment.var")
  # lag.in <- attr(matched.sets, 'lag')
  treated.ts <- as.numeric(unlist(strsplit(names(matched.sets), split = "[.]"))[c(F,T)])
  treated.ids <- as.numeric(unlist(strsplit(names(matched.sets), split = "[.]"))[c(T,F)])
  
  ordered.data <- plain.ordered.data[, c(id.var, time.var, treat.var, caliper.variable)]
  
  tlist <- expand.treated.ts(lag.in, treated.ts = treated.ts)
  idxlist <- get_yearly_dmats(as.matrix(ordered.data), treated.ids, tlist, paste0(ordered.data[,id.var], ".",
                                                                       ordered.data[, time.var]), matched_sets = matched.sets, lag.in)
  mahalmats <- build_maha_mats(ordered_expanded_data = as.matrix(ordered.data), idx =  idxlist)
  msets <- handle_perlag_caliper_calculations(mahalmats, matched.sets, caliper.value = caliper.distance)
  lag <- old.lag
  
  return(matched.sets)
  
}

handle_perlag_caliper_calculations <- function(nested.list, msets, caliper.value)
{
  do_calcs <- function(time.df)
  {
    cov.data <- time.df[, 4:ncol(time.df), drop = FALSE]

    result = tryCatch({
      dist(x = cov.data, method = "euclidean")
    }, warning = function(w) {

    }, error = function(e) {
        stop("error in caliper calculations")
    }, finally = {

    })
    result <- as.matrix(result)
    
    #TODO: double check on if euclidean distance with abs() will work
    result <- abs(result[nrow(result), ])
    return(result)
    
  }
  
  handle_set <- function(sub.list, cal.value)
  {
    #TODO: have it check the treated units
    results.temp <- lapply(sub.list, do_calcs)
    tmat <- do.call(rbind, results.temp)
    colnames(tmat) <- NULL
    tmat <- tmat[, 1:(ncol(tmat)-1), drop = FALSE]
    meets_caliper <- function(col, cal)
    {
      if(any(is.na(col)) || any(col > cal))
      {
        return(FALSE)
      }
      else {
        return(TRUE)
      }
    }
    unit.index <- apply(tmat, 2, FUN= meets_caliper, cal = cal.value)
    return(unit.index)
    
  }
  
  indices.msets <- lapply(nested.list, handle_set, cal.value = caliper.value)
  msets <- mapply(function(x, y) return(x[y]), x = msets, y = indices.msets)
  return(msets)
}