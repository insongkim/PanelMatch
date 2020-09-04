#' findBinaryTreated
#' 
#' \code{findBinaryTreated} is used to identify t,id pairs of units for which a matched set might exist. 
#' More precisely, it finds units for which at time t, the specified treatment has been applied, but at time t - 1, the treatment has not.
#'
#' @param dmat Data frame or matrix containing data used to identify potential treated units. Must be specified in such a way that a combination of time and id variables will correspond to a unique row. Must also contain at least a binary treatment variable column as well. 
#' @param treatedvar Character string that identifies the name of the column in \code{dmat} that provides information about the binary treatment variable
#' @param time.var Character string that identifies the name of the column in \code{dmat} that contains data about the time variable. This data must be integer that increases by one.
#' @param unit.var Character string that identifies the name of the column in \code{dmat} that contains data about the variable used as a unit id. This data must be integer
#' @param hasbeensorted variable that only has internal usage for optimization purposes. There should be no need for a user to toggle this
#' 
#' @return \code{findBinaryTreated} returns a subset of the data in the \code{dmat} data frame, containing only treated units for which a matched set might exist
#'
#' @keywords internal
#' 

findBinaryTreated <- function(dmat, treatedvar, time.var, unit.var, hasbeensorted = FALSE)
{
  dmat <- dmat[, c(unit.var, time.var, treatedvar)]
  #subset the columns to just the data needed for this operation
  
  colidx <- which(colnames(dmat) == treatedvar)
  if (length(colidx) > 1) stop("error in column naming scheme")
  uidx <- which(colnames(dmat) == unit.var)
  if (length(uidx) > 1) stop("error in column naming scheme")
  
  if (hasbeensorted)
  {
    odf <- dmat
  }
  else
  {
    odf <- dmat[order(dmat[,unit.var], dmat[,time.var]), ]  
  }
  classes <- sapply(dmat, class)
  # Verifying that time and unit data are integers -- perhaps later shift this to an automatic conversion process?
  if (classes[time.var] != "integer")
  {
    warning("time variable data provided not integer. Automatic conversion attempted")
    dmat[, time.var] <- as.integer(dmat[, time.var])
  }
  if (classes[unit.var] != "integer")
  {
    stop("unit id variable data provided not integer")
    #dmat[, unit.var] <- as.integer(dmat[, unit.var])
  }
  
  t.history <- odf[,treatedvar]
  t.idxs <- which(t.history == 1)
  
  num.df <- as.matrix(odf)
  if (!is.numeric(num.df)) stop("data in treated, time, or id columns is not numeric")
  
  ind <- get_treated_indices(num.df, t.idxs - 1, colidx - 1, uidx - 1)
  treated.unit.indices <- t.idxs[ind]
  odf <- odf[treated.unit.indices, ]
  rownames(odf) <- NULL
  return(odf)
}


#' get.matchedsets
#' 
#' \code{get.matchedsets} is used to identify matched sets for a given unit with a specified i, t.
#' 
#' @param t integer vector specifying the times of treated units for which matched sets should be found. This vector should be the same length as the following \code{id} parameter -- the entries at corresponding indices in each vector should form the t,id pair of a specified treatment unit. 
#' @param id integer vector specifying the unit ids of treated units for which matched sets should be found. note that both \code{t} and \code{id} can be of length 1
#' @param L An integer value indicating the length of treatment history to be matched
#' @param data data frame containing the data to be used for finding matched sets.
#' @param t.column Character string that identifies the name of the column in \code{data} that contains data about the time variable. Each specified entry in \code{t} should be somewhere in this column in the data. This data must be integer that increases by one.
#' @param id.column Character string that identifies the name of the column in \code{data} that contains data about the unit id variable. Each specified entry in \code{id} should be somewhere in this column in the data. This data must be integer.
#' @param treatedvar Character string that identifies the name of the column in \code{data} that contains data about the binary treatment variable. 
#' @param hasbeensorted variable that only has internal usage for optimization purposes. There should be no need for a user to toggle this
#' @param match.on.missingness TRUE/FALSE indicating whether or not the user wants to "match on missingness." That is, should units with NAs in their treatment history windows be matched with control units that have NA's in corresponding places?
#' @param matching logical indicating whether or not the treatment history should be used for matching. This should almost always be set to TRUE, except for specific situations where the user is interested in particular diagnostic questions.
#' @return \code{get.matchedsets} returns a "matched set" object, which primarily contains a named list of vectors. Each vector is a "matched set" containing the unit ids included in a matched set. The list names will indicate an i,t pair (formatted as "<i variable>.<t variable>") to which the vector/matched set corresponds.
#'
#' 
#' @keywords internal
#' 
get.matchedsets <- function(t, id, data, L, t.column, id.column, treatedvar, 
                            hasbeensorted = FALSE, match.on.missingness = TRUE, 
                            matching = TRUE, continuous = FALSE, 
                            continuous.treatment.info = NULL) 
{
  
  if (length(t) == 0 | length(id) == 0)
  {
    stop("time and/or unit information missing")
  }
  if (!hasbeensorted)
  {
    data <- data[order(data[,id.column], data[,t.column]), ]  
    
  }
  # Verifying that time and unit data are integers -- perhaps later shift this to an automatic conversion process?
  classes <- sapply(data, class)
  if (classes[t.column] != "integer")
  {
    warning("time variable data provided not integer. Automatic conversion attempted")
    data[, t.column] <- as.integer(data[, t.column])
  }
  if (classes[id.column] != "integer")
  {
    stop("unit id variable data provided not integer")
    #data[, id.column] <- as.integer(data[, id.column])
  }
  
  d <- data[, c(id.column, t.column, treatedvar)]
  d <- as.matrix(d)
  if (!is.numeric(d)) stop('data in treated, time, or id columns is not numeric')
  
  if (!continuous)
  {
    
    #CHECK TO MAKE SURE COLUMNS ARE IN ORDER!!!
    #fix factor conversion to be more smooth and allow for different data types.
    
    #dont think these are necessary anymore? since we balance up front...
    compmat <- data.table::dcast(data.table::as.data.table(d), formula = paste0(id.column, "~", t.column), 
                                 value.var = treatedvar) #reshape the data so each row corresponds to a unit, columns specify treatment over time
    # d <- data.table::melt(compmat, id = id.column, variable = t.column, value = treatedvar, 
    #                       variable.factor = FALSE, value.name = treatedvar)
    # #ensuring that we have integers
    # #newdata <- as.integer(as.character(unlist(d[, ..t.column])))
    # newdata <- as.integer(d[, get(t.column)])
    # d[, (t.column) := newdata]
    # #d <-  d[order(d[,..id.column], d[,..t.column]), ] #cast -> melt fills in missing data with NA's but order is not preserved, so second sort necessary
    # d <- d[order(d[, get(id.column)], d[,get(t.column)]), ]
    # d <- data.matrix(d)
    # ## try deleting all of the above up to note above
    if (match.on.missingness)
    {
      d[is.na(d[,treatedvar]), treatedvar] <- -1
      compmat[is.na(compmat)] <- -1  
    }
  
    control.histories <- get_comparison_histories(d, t, id, which(colnames(d) == t.column) - 1 ,
                                                  which(colnames(d) == id.column) - 1, L,
                                                  which(colnames(d) == treatedvar) - 1) #control histories should be a list
    
    if (!matching & !match.on.missingness)
    {
      tidx <- !unlist(lapply(control.histories, function(x)return(any(is.na(x)))))
      control.histories <- control.histories[tidx]
      t2 <- t[!tidx]
      id2 <- id[!tidx]
      t <- t[tidx]
      id <- id[tidx]
      t.map <- match(t, unique(d[, t.column]))
      sets <- non_matching_matcher(control.histories, as.matrix(compmat), t.map, id, L = 1, missing_window = L)
      if (any(!tidx))
      {
        l2 <- replicate(sum(!tidx), numeric())
        named.sets <- matched_set(matchedsets = c(sets, l2), id = c(id, id2), t = c(t, t2), L = L, 
                                  t.var = t.column, id.var = id.column, treatment.var = treatedvar)
      } else
      {
        named.sets <- matched_set(matchedsets = sets, id = id, t = t, L = L, 
                                  t.var = t.column, id.var = id.column, treatment.var = treatedvar)
      }
      
    } 
    else
    {
      t.map <- match(t, unique(d[, t.column]))
      sets <- get_msets_helper(control.histories, as.matrix(compmat), t.map, id, L)
      named.sets <- matched_set(matchedsets = sets, id = id, t = t, L = L, 
                                t.var = t.column, id.var = id.column, treatment.var = treatedvar)  
    }
    
    return(named.sets)
  } else # continuous
  {
    ##convert to caliper formula here...
    
    continuous.treatment.formula <- as.formula(paste0("~ I(caliper(", treatedvar,",", "'", continuous.treatment.info[["method"]], "'", 
                                                      ",", continuous.treatment.info[["matching.threshold"]],
                                                      ",", "'",continuous.treatment.info[["type"]],"'",
                                                      ",", "'", continuous.treatment.info[["units"]],"'" ,"))"))
    #~ I(caliper(cal.data,"max", 3, "numeric", "raw"))
    attr(continuous.treatment.formula, ".Environment") <- environment()
    matched.sets <- vector("list", length(id))
    attr(matched.sets, "")
    names(matched.sets) <- paste0(id, ".", t)
    
    attr(matched.sets, "lag") <- L
    attr(matched.sets, "t.var") <- t.column
    attr(matched.sets, "id.var" ) <- id.column
    attr(matched.sets, "treatment.var") <- treatedvar
    
    continuous.matched.sets <- handle_calipers(plain.ordered.data = d, caliper.formula = continuous.treatment.formula, 
                            matched.sets = matched.sets, lag.window = 0:L, is.continuous.matching = TRUE)
    return(continuous.matched.sets) #should maybe attach extra attributes for continuously matched msets
  }
}

  
#' findContinuousTreated
#' 
#' \code{findContinuousTreated} is used to identify t,id pairs of units for which a matched set might exist. 
#' 
#'
#' @param dmat Sorted data frame or matrix containing data used to identify potential treated units. Must be specified in such a way that a combination of time and id variables will correspond to a unique row. Must also contain at least a continuous treatment variable column as well. 
#' @param treatedvar Character string that identifies the name of the column in \code{dmat} that provides information about the continuous treatment variable
#' @param time.var Character string that identifies the name of the column in \code{dmat} that contains data about the time variable. This data must be integer that increases by one.
#' @param unit.var Character string that identifies the name of the column in \code{dmat} that contains data about the variable used as a unit id. This data must be integer
#' 
#' @return \code{findContinuousTreated} returns a subset of the data in the \code{dmat} data frame, containing only treated units for which a matched set might exist
#'
#' @keywords internal
#' 
findContinuousTreated <- function(dmat, treatedvar, time.var, unit.var,
                                  qoi, continuous.treatment.info)
{
  identifyContinuousIndex <- function(x, treatedvar.in, 
                                      qoi, threshold)
  {
    if (qoi == "att")
    {
      unitIndex <- which(diff(x[, treatedvar]) >= threshold) + 1 ## need to add one since it does not pad with NA values
    } else if (qoi == "atc") 
    {
      unitIndex <- which(diff(x[, treatedvar]) <= threshold) + 1 ## need to add one since it does not pad with NA values
    } else {
      warning("Undefined qoi for continuous matching!")
    }
    
    return(x[unitIndex, ])
  }
  # TODO: fix this
  ##evaluate the function here to extract threshold
  treatment.threshold <- continuous.treatment.info[["treatment.threshold"]] #F(continuous.treatment.formula)
  treatedUnits <- by(dmat, INDICES = dmat[, unit.var], FUN = identifyContinuousIndex, 
     treatedvar.in = treatedvar, qoi = qoi, threshold = treatment.threshold, simplify = FALSE)
  
  treatedDF <- do.call(rbind, treatedUnits)
  if (nrow(treatedDF) == 0) stop("No viable treated units for continuous matching specification")
  rownames(treatedDF) <- NULL
  return(treatedDF)
}
  




