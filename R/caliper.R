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
  
  caliper <- function(y, method, 
                      caliper.distance,
                      data_type, 
                      calculation_type, 
                      lwindow = lag.window)
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
  
  t.data <- do.call(rbind, by(plain.ordered.data, 
                              as.factor(plain.ordered.data[, 1]),
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
    
    matched.sets <- handle.single.caliper.per.lag(t.data, matched.sets,
                                                  caliper.metadata[2,i], 
                                                  as.numeric(caliper.metadata[3, i]), 
                                                  c(1:3, data.index),
                                                  .IS_FACTOR_VAR, .USE_SD_UNITS, 
                                                  is.continuous.matching, control.threshold)
  }
  
  return(matched.sets)
  
}

handle.single.caliper.per.lag <- function(plain.ordered.data, 
                                          matched.sets,
                                          caliper.method, 
                                          caliper.distance, 
                                          data.index, 
                                          is.factor.var, 
                                          use.sd.units, 
                                          do.continuous.matching = FALSE,
                                          control.threshold = NULL)
{
  
  time.var <- attr(matched.sets, "t.var")
  id.var <- attr(matched.sets, "id.var")
  treat.var <- attr(matched.sets, "treatment.var")
  lag.in <- attr(matched.sets, 'lag')
  
  old.lag <- lag.in
  lag.in <- 0 
  
  ordered.data <- plain.ordered.data[, data.index]
  
  msets <- apply_calipers(as.matrix(ordered.data), 
                                    matched.sets,
                                    caliper.distance,
                                    caliper.method, 
                                    is.factor.var, 
                                    use.sd.units, 
                                    id.var, 
                                    time.var, 
                                    lag.in,
                                    do.continuous.matching, 
                                    control.threshold)
  
  
  lag.in <- old.lag
  #fix attribute copying
  class(msets) <- c("matched.set", "list")
  attr(msets, 'refinement.method') <- NULL
  attr(msets, "lag") <- lag.in
  attr(msets, "t.var") <- time.var
  attr(msets, "id.var" ) <- id.var
  attr(msets, "treatment.var") <- treat.var
  return(msets)
  
}

meets_caliper <- function(col, 
                          cal, 
                          cal.method,
                          sd.vals,
                          IS_FACTOR)
{
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

do.maths  <- function(control.unit.data, 
                      treated.unit.data, 
                      sd.vals_, 
                      is_factor)
{
  
  if(!is_factor)
  {
    return(abs(treated.unit.data - control.unit.data) / sd.vals_) 
  } else
  { 
    return(abs(treated.unit.data - control.unit.data))  
  }
}

do_calcs <- function(time.df, 
                     sd.vals__, 
                     is_factor_in)
{
  
  cov.data <- time.df[, 4:ncol(time.df), drop = FALSE]

  results <- apply(cov.data[1:((nrow(cov.data)  - 1)),
                            , 
                            drop = FALSE],
                   MARGIN = 1, 
                   FUN  =  do.maths, 
                   treated.unit.data = cov.data[nrow(cov.data), ], 
                   sd.vals_ = sd.vals__, 
                   is_factor = is_factor_in)
  
  return(results)
  
}

handle_set <- function(sub.list, 
                       cal.value, 
                       cal.method, 
                       standard.deviations, 
                       is.factor)
{
  
  tmat <- do_calcs(sub.list[[1]], 
                   sd.vals__ = standard.deviations, 
                   is_factor_in = is.factor)
  colnames(tmat) <- NULL
  unit.index <- apply(tmat, 2, FUN= meets_caliper, 
                      cal = cal.value, 
                      cal.method = cal.method, 
                      sd.vals = standard.deviations,
                      IS_FACTOR = is.factor)
  return(unit.index)
  
}


handle_perlag_caliper_calculations <- function(nested.list, msets, caliper.method, 
                                               caliper.value, is.factor.var, 
                                               use.sd.units, full.data, 
                                               id.var, time.var, treated_it_pairs)
{
  
  if (length(msets) == 0) return(numeric(0))
  
  #get standard deviation values
  
  row.names(full.data) <- paste0(full.data[, id.var], 
                                 ".", 
                                 full.data[, time.var])
  
  tdf <- full.data[treated_it_pairs, ]
  
  if(use.sd.units)
  {
    sd.vals <- apply(tdf[, 4:ncol(tdf)],
                     MARGIN = 2, 
                     FUN = sd, 
                     na.rm = TRUE)
    if ( 0 %in% sd.vals || any(is.na(sd.vals)) )
    {
      stop("not enough variation for using 'units = sd' in continuous.treatment.info. 
           Try using raw units instead")
    }
  } else {
    sd.vals <- 1
  }
  
  indices.msets <- handle_set(nested.list, 
                              caliper.value, 
                              caliper.method, 
                              sd.vals, 
                              is.factor.var)
  msets <- msets[indices.msets]
  return(msets)
}



apply_calipers <- function(ordered_expanded_data,
                                     matched.sets, 
                                     calipervalue,
                                     calipermethod, 
                                     isfactor, 
                                     use.sd, 
                                     id.var,
                                     time.var, 
                                     lag.in, 
                                     continuous.matching = FALSE,
                                     control.threshold = NULL)
{
  
  result <- mapply(FUN = unnest, 
                   matched.set = matched.sets, 
                   treated.unit.info = names(matched.sets),
                   MoreArgs = list(lag.in_ = lag.in,
                                   ordered_expanded_data_ = ordered_expanded_data,
                                   is.continuous.matching = continuous.matching,
                                   idvar = id.var,
                                   timevar = time.var,
                                   all.treated.info = names(matched.sets),
                                   control.threshold = control.threshold,
                                   calipervalue = calipervalue,
                                   calipermethod = calipermethod,
                                   isfactor = isfactor,
                                   use.sd = use.sd,
                                   id.var = id.var,
                                   time.var = time.var,
                                   matched.sets = matched.sets),
                   SIMPLIFY = FALSE)
  
  names(result) <- names(matched.sets)
  if (!continuous.matching)
  {
    result <- result[sapply(result, length) > 0]
  }
  
  return(result)
  
  
}

unnest <- function(matched.set, treated.unit.info,
                   lag.in_, ordered_expanded_data_, 
                   is.continuous.matching = FALSE,
                   idvar, timevar, all.treated.info,
                   control.threshold = NULL, 
                   calipervalue,
                   calipermethod,
                   isfactor,
                   use.sd,
                   id.var,
                   time.var,
                   matched.sets)
{
  treated.ts <- as.integer(sub(".*\\.", "", treated.unit.info))
  treated.ids <- as.integer(sub("\\..*", "", treated.unit.info))
  
  all.treated.ts <- as.integer(sub(".*\\.", "", all.treated.info))
  all.treated.ids <- as.integer(sub("\\..*", "", all.treated.info))
  
  tlist <- expand.treated.ts(lag.in_, treated.ts = treated.ts)
  
  idxlist <- get_yearly_dmats(ordered_expanded_data_, treated.ids, tlist,
                              matched_sets = list(matched.set), lag.in_)
  rr <- lapply(unlist(idxlist, recursive = FALSE), function(x) {ordered_expanded_data_[x, ]})
  
  tset <- handle_perlag_caliper_calculations(rr, matched.set, calipermethod, calipervalue,
                                             isfactor, use.sd, 
                                             ordered_expanded_data_, 
                                             id.var, time.var, names(matched.sets)) 
  return(tset)
}