PS_DID <- function (L, F, M = 1, time.id = "year", unit.id = "ccode", treatment, 
                    covariate, dependent, data) 
{
  L <- L
  F <- F
  M <- M
  data <- data[c(unit.id, time.id,  treatment, dependent, covariate)]
  dlist <- lapply(1:L, 
                  function (i) slide(data = data, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                     NewVar = paste("dependent_l", i, sep="")))
  data <- Reduce(function(x, y) {merge(x, y)}, dlist)
  varnames <- c(time.id, unit.id, treatment, dependent, c(covariate, names(data[, (4 + length(covariate) + 1):length(data)])))
  d2 <- na.omit(data[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x)))
  dmatrix <- as.matrix(d2)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L, 
                                                                       F, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  smallerlist <- lapply(smallerlist, dframelist.rb_dup)
  smallerlist <- Filter(function(x) nrow(x) > 2 * (L + F + 
                                                     1), smallerlist)
  even_smaller1 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
                                          x$V1 == unique(x$V1)[1], ]$V3 == 0, smallerlist)
  Fs <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 == max(unique(x$V1)), ]
    return(x)
  })
  even_smaller1 <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 %in% sort(unique(x$V1))[1:(L+1)], ]
    return(x)
  })
  
  even_smaller1 <- lapply(even_smaller1, function (x) {
    colnames(x) <- varnames
    return(x)}
  )
  Fs <- lapply(Fs, function (x) {
    colnames(x) <- varnames
    return(x)}
  )
  
  pooled <- rbindlist(even_smaller1)
  
  fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
              family = binomial(link = "logit"), data = pooled)
  pooled$ps <- fit0$fitted.values
  
  aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
  newlist <- lapply(even_smaller1, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
  newlist <- Map(rbind.fill, newlist, Fs)
  # superlist <- mapply(function(x, y) rbind.fill(x,y), newlist, Fs)
  
  all.diffs_PS <- lapply(newlist, PS_m_did, L = L, F, M = M)
  ATT <- mean(unlist(all.diffs_PS), na.rm = T)
  
  dmatrix[, 3] <- ifelse(dmatrix[, 3] == 1, 0, 1)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L = L, 
                                                                       F, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  smallerlist <- lapply(smallerlist, dframelist.rb_dup)
  smallerlist <- Filter(function(x) nrow(x) > 2 * (L + F + 
                                                     1), smallerlist)
  even_smaller2 <- Filter(function(x) x[x$V2 == unique(x$V2)[2] & 
                                          x$V1 == unique(x$V1)[1], ]$V3 == 0, smallerlist)
  
  Fs <- lapply(even_smaller2, function(x) {
    x <- x[x$V1 == max(unique(x$V1)), ]
    return(x)
  })
  even_smaller2 <- lapply(even_smaller2, function(x) {
    x <- x[x$V1 != max(unique(x$V1)), ]
    return(x)
  })
  
  even_smaller2 <- lapply(even_smaller2, function (x) {
    colnames(x) <- varnames
    return(x)}
  )
  Fs <- lapply(Fs, function (x) {
    colnames(x) <- varnames
    return(x)}
  )
  
  pooled <- rbindlist(even_smaller2)
  
  fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
              family = binomial(link = "logit"), data = pooled)
  pooled$ps <- fit0$fitted.values
  
  aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
  newlist <- lapply(even_smaller2, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
  newlist <- Map(rbind.fill, newlist, Fs)
  # superlist <- mapply(function(x, y) rbind.fill(x,y), newlist, Fs)
  
  all.diffs_PS2 <- lapply(newlist, PS_m_did, L = L, F, M = M)
  
  ATC <- mean(unlist(all.diffs_PS2), na.rm = T)
  
  ATC <- -(ATC)
  DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(na.omit(all.diffs_PS)) + 
                                        ATC * length(na.omit(all.diffs_PS2)))/(length(na.omit(all.diffs_PS)) + 
                                                                                 length(na.omit(all.diffs_PS2))), ATT)
  return(DID_ATE)
}
