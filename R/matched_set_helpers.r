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

findBinaryTreated <- function(dmat, qoi.in,
                              treatedvar, 
                              time.var, 
                              unit.var, 
                              hasbeensorted = FALSE)
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
  
  if (classes[time.var] != "integer")
  {
    warning("time variable data provided not integer. Automatic conversion attempted")
    dmat[, time.var] <- as.integer(dmat[, time.var])
  }
  if (classes[unit.var] != "integer")
  {
    stop("unit id variable data provided not integer")
  }
  
  if (identical(qoi.in, "atc")) 
  {
    
    t.history <- odf[, treatedvar]
    id.history <- odf[, unit.var]
    c.idxs <- which(t.history == 0)
    t1 <- c.idxs - 1
    
    t1[t1 == 0] <- NA
    
    idx <- t.history[c.idxs] == 0 & 
      t.history[t1] == 0 & 
      (id.history[c.idxs] == id.history[t1])
    idx[is.na(idx)] <- FALSE
    odf <- odf[c.idxs[idx],]
    rownames(odf) <- NULL
  } else {
    t.history <- odf[,treatedvar]
    t.idxs <- which(t.history == 1)
    
    num.df <- as.matrix(odf)
    if (!is.numeric(num.df)) stop("data in treated, time, or id columns is not numeric")
    
    ind <- get_treated_indices(num.df, t.idxs - 1, colidx - 1, uidx - 1)
    treated.unit.indices <- t.idxs[ind]
    odf <- odf[treated.unit.indices, ]
    rownames(odf) <- NULL
  }
  
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
                            matching = TRUE, 
                            qoi.in,
                            restrict.control.period = NULL,
                            continuous.treatment.info = NULL)
{
  continuous <- !is.null(continuous.treatment.info)
  if (length(t) == 0 || length(id) == 0)
  {
    stop("time and/or unit information missing")
  }
  if (!hasbeensorted)
  {
    data <- data[order(data[,id.column], data[,t.column]), ]

  }
  classes <- sapply(data, class)
  if (classes[t.column] != "integer")
  {
    warning("time variable data provided not integer. Automatic conversion attempted")
    data[, t.column] <- as.integer(data[, t.column])
  }
  if (classes[id.column] != "integer")
  {
    stop("unit id variable data provided not integer")
  }

  d <- data[, c(id.column, t.column, treatedvar)]
  d <- as.matrix(d)
  if (!is.numeric(d)) stop('data in treated, time, or id columns is not numeric')
  
  if (!continuous)
  {
    compmat <- data.table::dcast(data.table::as.data.table(d), 
                                 formula = paste0(id.column, "~", t.column),
                                 value.var = treatedvar) #reshape the data so each row corresponds to a unit, columns specify treatment over time
    
    
    if (match.on.missingness)
    {
      d[is.na(d[,treatedvar]), treatedvar] <- -1
      compmat[is.na(compmat)] <- -1
    }
    
    control.histories <- get_comparison_histories(d, t, id, which(colnames(d) == t.column) - 1 ,
                                                  which(colnames(d) == id.column) - 1, L,
                                                  which(colnames(d) == treatedvar) - 1,
                                                  identical(qoi.in, "atc")) #control histories should be a list
    
    if (!is.null(restrict.control.period))
    {
      indx <- enforce_strict_histories(control.histories, 
                                       restrict.control.period)
      control.histories <- control.histories[indx]
      t <- t[indx]
      id <- id[indx]
    }
    
    if (!matching & !match.on.missingness)
    {
      tidx <- !unlist(lapply(control.histories, 
                             function(x)return(any(is.na(x)))))
      control.histories <- control.histories[tidx]
      t2 <- t[!tidx]
      id2 <- id[!tidx]
      t <- t[tidx]
      id <- id[tidx]
      t.map <- match(t, unique(d[, t.column]))
      sets <- non_matching_matcher(control.histories, as.matrix(compmat), t.map, 
                                   id, L = 1, missing_window = L)
      if (any(!tidx))
      {
        l2 <- replicate(sum(!tidx), numeric())
        named.sets <- matched_set(matchedsets = c(sets, l2), 
                                  id = c(id, id2), 
                                  t = c(t, t2), L = L,
                                  t.var = t.column, 
                                  id.var = id.column, 
                                  treatment.var = treatedvar)
      } else
      {
        named.sets <- matched_set(matchedsets = sets, 
                                  id = id, t = t, L = L,
                                  t.var = t.column, id.var = id.column, 
                                  treatment.var = treatedvar)
      }
    }
    else
    {
      t.map <- match(t, unique(d[, t.column]))
      sets <- get_msets_helper(control.histories, 
                               as.matrix(compmat), 
                               t.map, id, L)
      named.sets <- matched_set(matchedsets = sets, 
                                id = id, t = t, L = L,
                                t.var = t.column, 
                                id.var = id.column, 
                                treatment.var = treatedvar)
    }
    attr(named.sets, "restrict.control.period") <- restrict.control.period
    return(named.sets)
  } else #continuous
  {
    
    
    continuous.matched.sets <- continuous_treatment_matching(d,id.column,
                                  t.column,
                                  treatedvar,
                                  id,
                                  t,
                                  L, 
                                  continuous.matching.info = continuous.treatment.info)

    named.sets <- matched_set(matchedsets = continuous.matched.sets, 
                              id.t.pairs = names(continuous.matched.sets),
                              L = L,
                              t.var = t.column, 
                              id.var = id.column, 
                              treatment.var = treatedvar)
    return(named.sets) 
  }
  
   
}


extract.baseline.treatment <- function(matched.sets,
                                       data.in,
                                       id.variable,
                                       time.variable,
                                       treatment.variable)
{
  treated.t <- as.integer(sub(".*\\.", "", names(matched.sets)))
  treated.id <- as.integer(sub("\\..*", "", names(matched.sets)))
  treated.tm1 <- treated.t - 1
  
  rownames(data.in) <- paste0(data.in[,id.variable],
                              ".",
                              data.in[,time.variable])
  
  idx <- paste0(treated.id, ".", treated.tm1)
  
  t.starting <- data.in[idx, treatment.variable]
  names(t.starting) <- NULL
  for (i in 1:length(matched.sets)) {
    attr(matched.sets[[i]], "treatment.baseline") <- t.starting[i]
  }
  
  return(matched.sets)
}

############################################################
#### This function calculates the differences from t-1 to 1 for treated and control units in the treatment variable
#### While functionality is somewhat trivial for current implementation of package, it will be needed for multicategorical treatment version.
############################################################
extract.differences <- function(indexed.data, matched.set, 
                                treatment.variable,
                                qoi)
{
  treated.t <- as.integer(sub(".*\\.", "", names(matched.set)))
  treated.id <- as.integer(sub("\\..*", "", names(matched.set)))
  treated.tm1 <- treated.t - 1
  
  treated.key.t <- names(matched.set)
  treated.key.tm1 <- paste0(treated.id, ".", treated.tm1)
  
  multi.factor <- 1 
  if(length(matched.set[[1]]) > 0)
  {
    control.keys.t <- paste0(matched.set[[1]], ".", treated.t)
    control.keys.tm1 <- paste0(matched.set[[1]], ".", treated.tm1)
    
    keys.t <- c(treated.key.t, control.keys.t)
    keys.tm1 <- c(treated.key.tm1, control.keys.tm1)
    
    differences <- as.numeric(indexed.data[keys.t, treatment.variable] - 
                                indexed.data[keys.tm1, treatment.variable])
    
    attr(matched.set[[1]], "treatment.change") <- differences[1] * multi.factor
     
    attr(matched.set[[1]], "control.change") <- differences[2:length(differences)] * multi.factor
    
  } else {
    differences <- as.numeric(indexed.data[treated.key.t, treatment.variable] - 
                                indexed.data[treated.key.tm1,
                                             treatment.variable])
    attr(matched.set[[1]], "treatment.change") <- differences[1] * multi.factor
  }
  return(matched.set[[1]])
  
}

extractDifferences <- function(indexed.data, 
                               data.keys,
                               matched.sets, 
                               id.variable,
                                treatment.variable, 
                               time.variable)
{
  # use data.table to do grouped lag variable
  
  indexed.data[, paste0(treatment.variable, "_l1")] <- unlist(tapply(indexed.data[, treatment.variable], 
                                                                     indexed.data[, id.variable], 
                                                                     function(x) c(NA, x[-length(x)])))
  
  diff.vec <- indexed.data[, treatment.variable] - 
    indexed.data[, paste0(treatment.variable, "_l1")]

  names(diff.vec) <- data.keys
  
  if (length(matched.sets) == 0)
  {
    stop("no matched sets")
  }
  for (i in 1:length(matched.sets)) {
    
    treated.t <- as.integer(sub(".*\\.", "", names(matched.sets)[i]))
    attr(matched.sets[[i]], "control.change") <- diff.vec[paste0(matched.sets[[i]], 
                                                                 ".", treated.t)]
    attr(matched.sets[[i]], "treatment.change") <- diff.vec[names(matched.sets)[i]]
  }
  return(matched.sets)
}


identifyDirectionalChanges <- function(msets, ordered.data, id.var, time.var,
                                       treatment.var, qoi)
{
  rownames(ordered.data) <- paste0(ordered.data[, id.var], 
                                   ".", 
                                   ordered.data[, time.var])
  msets <- extractDifferences(ordered.data,
                              rownames(ordered.data),
                              msets,
                              id.var,
                              treatment.var,
                              time.var)
  
  # for (i in 1:length(msets)) {
  #     
  #   msets[[i]] <- extract.differences(ordered.data, msets[i], 
  #                                     treatment.var, qoi)
  #     
  # }
  return(msets)
}

# builds a list that contains all times in a lag window that correspond to a particular treated unit. This is structured as a list of vectors. Each vector is lag + 1 units long. The overall list will
# be the same length as the number of matched sets
expand.treated.ts <- function(lag, treated.ts)
{
  helper <- function(treated.t)
  {
    return(seq(from =  (treated.t - lag), 
               to = treated.t, by = 1))
  }
  lapply(treated.ts, helper)
}

#' findContinuousTreated
#'
#' \code{findContinuousTreated} is used to identify t,id pairs of units for which a matched set might exist.
#'
#'
#' @param dmat Sorted data frame or matrix containing data used to identify potential treated units. Must be specified in such a way that a combination of time and id variables will correspond to a unique row. Must also contain at least a continuous treatment variable column as well. This data should be balanced as well. 
#' @param treatedvar Character string that identifies the name of the column in \code{dmat} that provides information about the continuous treatment variable
#' @param time.var Character string that identifies the name of the column in \code{dmat} that contains data about the time variable. This data must be integer that increases by one.
#' @param unit.var Character string that identifies the name of the column in \code{dmat} that contains data about the variable used as a unit id. This data must be integer
#'
#' @return \code{findContinuousTreated} returns a subset of the data in the \code{dmat} data frame, containing only treated units for which a matched set might exist
#'
#' @keywords internal
#'
findContinuousTreated <- function(dmat, treatedvar, time.var, unit.var,
                                  qoi, continuous.treatment.info, 
                                  lag = NULL)
{
  identifyContinuousIndex <- function(x, treatedvar.in,
                                      qoi, threshold)
  {
    if (identical(qoi, "att"))
    {
      unitIndex <- which( diff(x[, treatedvar]) >= threshold) + 1 ## need to add one since it does not pad with NA values
    } else if (identical(qoi, "art"))
    {
      # assuming that for ART, the matching threshold parameter provided is always positive.
      unitIndex <- which( diff(x[, treatedvar]) <= -1 * threshold) + 1 ## need to add one since it does not pad with NA values
    } else {
      warning("Undefined qoi for continuous matching!")
    }
    
    return(x[unitIndex, ])
  }
  
  
  
  treatment.threshold <- continuous.treatment.info[["treatment.threshold"]] #(continuous.treatment.formula)
  treatedUnits <- by(dmat, INDICES = dmat[, unit.var], FUN = identifyContinuousIndex,
                     treatedvar.in = treatedvar, qoi = qoi, threshold = treatment.threshold, simplify = FALSE)
  
  temp.treateds <- do.call(rbind, treatedUnits)
  if (nrow(temp.treateds) == 0) stop("No viable treated units for continuous matching specification")
  
  
  if(!is.null(continuous.treatment.info[["minimum.treatment.value"]]))
  {
    indx <- temp.treateds[, treatedvar] >= continuous.treatment.info[["minimum.treatment.value"]]
    temp.treateds <- temp.treateds[indx,]
  }
  if(!is.null(continuous.treatment.info[["maximum.treatment.value"]]))
  {
    indx <- temp.treateds[, treatedvar] <= continuous.treatment.info[["maximum.treatment.value"]]
    temp.treateds <- temp.treateds[indx,]
  }
  
  if (!is.null(lag))
  {
    min.per.unit <- min(dmat[,time.var])
    indx <- temp.treateds[, time.var] - lag >= min.per.unit #eliminate anything where lag is too large to get a full treatment history.
    temp.treateds <- temp.treateds[indx,]
  }
  
  rownames(temp.treateds) <- NULL
  return(temp.treateds)
}


prepContinuousControlUnits <- function(ordered.data,
                                       idvar, time.var,
                                       all.treated.ids,
                                       all.treated.ts,
                                       treated.ids,
                                       treated.ts,
                                       control.threshold
)
{
  full.controls <- unique(ordered.data[, idvar])
  
  not.valid.ids <- all.treated.ids[all.treated.ts %in% treated.ts] #cant include other treated units from the same time
  matched.set <- full.controls[!full.controls %in% not.valid.ids] # start with everything, then remove any units that are treated
  # during the same period as the current t/id pair under consideration
  rownames(ordered.data) <- paste0(ordered.data[, idvar], ".", ordered.data[, time.var])
  controls.to.check.t <- paste0(matched.set, ".", treated.ts)
  controls.to.check.t1 <- paste0(matched.set, ".", treated.ts - 1)
  # assuming third column is treated column
  diffs <- as.numeric(ordered.data[controls.to.check.t, 3]) - as.numeric(ordered.data[controls.to.check.t1, 3])
  matched.set <- matched.set[abs(diffs) <= control.threshold]
  matched.set <- matched.set[!is.na(matched.set)] #final cleanup, remove the NAs 
  return(matched.set)
}


filterContinuousTreated <- function(msets, e.sets, qoi,
                                    control.threshold = NULL)
{
  
  
  # if (identical(qoi, "att"))
  # {
  #   idx <- sapply(msets, function(x) attr(x, "treatment.change")) >= 0
  #   msets <- msets[idx]
  #   if (length(e.sets) > 0) 
  #   {
  #     idx <- sapply(e.sets, function(x) attr(x, "treatment.change")) >= 0
  #     e.sets <- e.sets[idx]  
  #   }
  #   
  # } else if (identical(qoi, "art"))
  # {
  #   idx <- sapply(msets, function(x) attr(x, "treatment.change")) <= 0
  #   msets <- msets[idx]
  #   if (length(e.sets) > 0) 
  #   {
  #     idx <- sapply(e.sets, function(x) attr(x, "treatment.change")) <= 0
  #     e.sets <- e.sets[idx]  
  #   }
  # } else 
  # {
  #   stop("Invalid QOI for continuous treatment!")
  # }  
  
  ## remove the invalid controls
  if (length(msets) > 0)
  {
    for (i in 1:length(msets)) {
      idx.vector <- abs(attr(msets[[i]], "control.change")) <= control.threshold
      
      clean_vector <- function(vec) {
        sapply(vec, function(x) {
          if (is.na(x)) {
            return(FALSE)
          } else {
            return(x)
          }
        })
      }
      idx.vector.cleaned <- clean_vector(idx.vector)
      control.changes <- attr(msets[[i]], "control.change")[idx.vector.cleaned]
      treatment.baseline <- attr(msets[[i]], "treatment.baseline") 
      treatment.change <- attr(msets[[i]], "treatment.change") 
      msets[[i]] <- msets[[i]][idx.vector.cleaned]
      attr(msets[[i]], "control.change") <- control.changes
      attr(msets[[i]], "treatment.baseline") <- treatment.baseline
      attr(msets[[i]], "treatment.change") <- treatment.change
      
    }  
  }
  
  # if (length(e.sets) > 0)
  # {
  #   for (i in 1:length(e.sets)) {
  #     e.sets[[i]][abs(attr(e.sets[[i]], "control.change")) <= control.threshold]
  #   }
  # }
  e.sets <- msets[sapply(msets, length) == 0]
  msets <- msets[sapply(msets, length) > 0 ]
  
  return(list("sets" = msets, 
              "empty.sets" = e.sets))
}


continuous_treatment_matching <- function(expanded_data,
                                          unit.id.variable,
                                          time.variable,
                                          treatment.variable,
                                          treated_ids,
                                          treated_ts,
                                          lag,
                                          continuous.matching.info)
{
  
  dt <- expanded_data[, c(unit.id.variable,
                          time.variable,
                          treatment.variable)]
  
  rownames(dt) <- paste0(expanded_data[, unit.id.variable],
                         '.',
                         expanded_data[, time.variable])
  ft <- function(x)
  {
    as.matrix(dist(x[,treatment.variable,
                     drop = FALSE], 
                   method = "manhattan", 
                   diag = TRUE, 
                   upper = TRUE))
  }
  
  dl <- by(dt, 
           INDICES = as.factor(dt[, time.variable]), 
           FUN = ft, 
           simplify = FALSE)
  
  if (length(dl) > 0) {
    # we should verify that rownames/colnames are identical in each
    treated.ids <- as.integer(sub("\\..*", "", 
                                  rownames(dl[[1]])))
    id.row.map <- data.frame(treated.id = treated.ids,
                             row.col.idx = 1:length(treated.ids))
    
  }
  
  year.dists <- data.frame(t.var = names(dl),
                           idx = 1:length(dl))
  
  yt <- data.frame(treated.id = treated_ids,
                   treated.t = as.character(treated_ts), # since we use this as grouping variable
                   L = lag)
  
  temp1 <- merge(yt,
                 year.dists,
                 by.x = "treated.t",
                 by.y = "t.var",
                 all.x = TRUE)
  temp2 <- merge(temp1,
                 id.row.map,
                 by = "treated.id",
                 all.x= TRUE)
  # re sort these all later
  fl <- function(xrow)
  {
    #browser()
    idx.row <- as.numeric(xrow['row.col.idx'])
    dist.mat.idcs <- (as.numeric(xrow['idx']) - as.numeric(xrow['L'])):(as.numeric(xrow['idx']) - 1)
    
    get_units <- function(i, j, caliper = continuous.matching.info[["matching.threshold"]])
    {
      #browser()
      ds <- dl[[i]]
      id.t <- colnames(ds)[which(ds[j, ] <= caliper)]
      treated.ids <- as.integer(sub("\\..*", "", id.t))
      return(treated.ids)
    }
    
    res <- lapply(dist.mat.idcs,
           FUN = get_units, 
           j = idx.row)
    
    attr(res, "treated.obs.time") <- xrow["treated.t"]
    return(res)
    
  }

  sdf <-apply(temp2,
        MARGIN = 1,
        FUN = fl,
        simplify = FALSE)

  # Function to find the intersection of numeric vectors in a sublist
  find_intersection <- function(sublist) {
    Reduce(intersect, sublist)
  }
  
  # Apply the function to each sublist in the outer list
  intersections <- lapply(sdf, 
                          find_intersection)

  filter_mset <- function(reduced.mset,
                          treated.t,
                          treated.observations)
    
  {
    to.check <- paste0(reduced.mset, '.', treated.t) 
    return(reduced.mset[!to.check %in% treated.observations])
  }
  treated.obs <- paste0(treated_ids, '.', treated_ts)
  matched.sets <- mapply(filter_mset, 
         reduced.mset = intersections,
         treated.t = temp2$treated.t,
         MoreArgs = list(treated.observations = treated.obs))
  names(matched.sets) <- paste0(temp2[, "treated.id"], 
                                '.', 
                                temp2[, "treated.t"])
  return(matched.sets)
  
}