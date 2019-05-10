exact.match.pcs <- function(sets)
{
  browser()
  L <- attr(sets, "lag")
  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])
  make.years <- function(t, repnum)
  {
    q <- rep(c((t - L):t), repnum)
    return(q)
  }
  lsets <- sapply(sets, length)
  ts <- mapply(FUN = make.years, t = ts, repnum = lsets, SIMPLIFY = F)
  iddata <- lapply(sets, function(x){sapply(x, function(y){rep(y, L + 1)})})
  iddata <- lapply(iddata, as.numeric)
  names(iddata) <- NULL
  create.keys <- function(t, id)
  {
    return(paste0(id,'.',t))
  }
  return(mapply(create.keys, t = ts, id = iddata))
}

exact.match.pts <- function(sets, lead.in = 0)
{
  include <- sapply(sets, length) > 0
  num.empty <- sum(!include)
  tids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])[include]
  tids <- rep(tids, rep(2, length(tids)))
  ts <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(F,T)])[include]
  make.years <- function(t, lead, repnum)
  {
    q <- rep( c((t - 1), (t + lead)), repnum)
    return(q)
  }
  ts <- unlist(mapply(FUN = make.years, t = ts, lead = lead.in, repnum = 1, SIMPLIFY = F))
  data.frame(t = ts, id = tids)
}