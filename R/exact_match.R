do_exact_matching <- function(sets, balanced.panel.data, exact.match.vars)
{
  L <- attr(sets, "lag")

  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])
  make.years <- function(t, repnum)
  {
    q <- rep(c((t - 0):t), repnum)
    return(q)
  }
  lsets <- sapply(sets, length)
  ts <- mapply(FUN = make.years, t = ts, repnum = lsets, SIMPLIFY = F)
  iddata <- lapply(sets, as.numeric)
  names(iddata) <- NULL
  create.keys <- function(t, id)
  {
    return(paste0(id,'.',t))
  }
  control.data <- (mapply(create.keys, t = ts, id = iddata, SIMPLIFY = F))
  treatment.data <- names(sets)
  rowkeys <- paste0(balanced.panel.data[,1], '.', balanced.panel.data[, 2]) #think we can assume we have unit + time columns first
  
  colidx <- which(colnames(balanced.panel.data) %in% exact.match.vars)
  bpd <- as.matrix(balanced.panel.data[, c(1,2, colidx)])
  expanded.lists <- do_exact_matching_refinement(bpd, L, rowkeys, control.data,
                               treatment.data, (3:ncol(bpd) - 1))
  expanded.list <- unlist(expanded.lists, recursive = F)
  condensed.list <- list()
  for (i in 1:length(sets)) {
    tlist <- expanded.list[seq(from = i, to = length(expanded.list), by = length(sets))]
    matrix.set <- do.call(rbind, tlist)
    condensed.list[[i]] <- apply(matrix.set, MARGIN = 2, all) #gives us consolidated t/f index for each matched set
  }
  for (i in 1:length(sets)) {
    sets[[i]] <- sets[[i]][condensed.list[[i]]]
  }
  
  return(sets)
}
