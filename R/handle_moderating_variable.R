#' handle_moderating_variable
#' 
#' handles moderating variable calculations: In practice, this just involves slicing the data up according to the moderator, calling PanelEstimate() and putting everything back together This function creates the sets of objects on which PanelEstimate() will be called. It identifies the set of valid values the moderating variable can take on.
#' @param ordered.data data.frame
#' @param att.sets matched.set object for the ATT or ART
#' @param atc.sets matched.set object for the ATC
#' @param PM.object PanelMatch object
#' @param moderator string specifying the name of the moderating variable
#' @param unit.id string specifying the unit id variable
#' @param time.id string specifying the time id variable
#' @param qoi.in string specifying the QOI
#'
#' @return Character vector of valid moderating variable values
#' @keywords internal
handle_moderating_variable <- function(ordered.data, att.sets, 
                                       atc.sets, PM.object,
                                       moderator, unit.id, time.id, qoi.in)
{
  
  .reconstruct_pm_objects <- function(att.set = NULL, 
                                      atc.set = NULL, 
                                      PM.object_)
  {
    t.pm.object <- list()
    if(!is.null(att.set))
    {
      
      t.pm.object[[qoi.in]] <- att.set
    }
    if(!is.null(atc.set))
    {
      t.pm.object[["atc"]] <- atc.set
    }
    attrib <- names(attributes(PM.object_))[names(attributes(PM.object_)) != "names"]
    for(tatt in attrib)
    {
      attr(t.pm.object, tatt) <- attr(PM.object_, tatt)
    }
    return(t.pm.object)
  }
  
  ref.names <- paste0(ordered.data[, unit.id], ".", ordered.data[, time.id])
  moderator.vector <- ordered.data[, moderator]
  names(moderator.vector) <- ref.names
  moderated.sets.att <- list()
  moderated.sets.atc <- list()
  subset.list <- list()
  moderating.values <- unique(ordered.data[, moderator])
  valid.moderating.values <- c()
  for(val in as.vector(na.omit(moderating.values)))
  {
    #make sure we handle empty set situation
    if(!is.null(att.sets))
    {
      indx.set <- moderator.vector[names(att.sets)] == val
      t.set <- att.sets[indx.set]
      t.sum <- summary(t.set)
      if(length(t.set) > 0 && 
         (t.sum$num.units.empty.set < nrow(t.sum$overview)) )
      {
        moderated.sets.att[[make.names(val)]] <- t.set
        valid.moderating.values <- append(valid.moderating.values, val)
      }
      else
      {
        moderated.sets.att[[make.names(val)]] <- NULL
      }
    }
    if(!is.null(atc.sets))
    {
      indx.set <- moderator.vector[names(atc.sets)] == val
      t.set <- atc.sets[indx.set]
      t.sum <- summary(t.set)
      if(length(t.set) > 0 && 
         (t.sum$num.units.empty.set < nrow(t.sum$overview)) )
      {
        moderated.sets.atc[[make.names(val)]] <- t.set
        valid.moderating.values <- append(valid.moderating.values, val)
      }
      else
      {
        moderated.sets.atc[[make.names(val)]] <- NULL
      }
      
    }
  }
  if(!is.null(att.sets) & !is.null(atc.sets))
  {
    if(!identical(names(moderated.sets.att), names(moderated.sets.atc)))
    {
      stop("Data insufficient for calculating ATE with moderating variables")
    }
  }
  
  if(is.null(atc.sets))
  {
    #do att
    ret.obj <-lapply(moderated.sets.att,
                     FUN = .reconstruct_pm_objects,
                     PM.object_ = PM.object, atc.set = NULL)
  }
  else if(is.null(att.sets))
  {
    #do atc
    ret.obj <-lapply(moderated.sets.atc,
                     FUN = .reconstruct_pm_objects,
                     PM.object_ = PM.object, att.set = NULL)
  }
  else
  {
    
    ret.obj <- mapply(FUN = .reconstruct_pm_objects, SIMPLIFY = FALSE,
                      att.set = moderated.sets.att, atc.set = moderated.sets.atc,
                      MoreArgs = list(PM.object_ = PM.object))
  }
  
  names(ret.obj) <- as.character(valid.moderating.values)
  return(ret.obj)
}