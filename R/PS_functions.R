PS_m_did <- function(x, L, FORWARD, M = M) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  colnames(x)[1:5] <- c("V2", "V1", "ps", "V4", "V5")
  x <- x[!duplicated(x[c("V2", "V1")]),]
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x[x$V4 == 1 & x$V1 == (max(x$V1)-FORWARD), ]$V2
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-((length(unique(x$V1)) - FORWARD):length(unique(x$V1)))]
  matched_set <- x[x$V1 %in% timeid, ]
  
  PS_distance <- abs(tapply(matched_set$ps, matched_set$V2, mean) - mean(matched_set$ps[which(matched_set$V2 == treated.id)]))
  
  first.diff <- x[x$V2 == treated.id & x$V1 == max(x$V1), ]$V5 - 
    x[x$V2 == treated.id & x$V1 == sort(unique(x$V1))[L], ]$V5
  
  if (M < length(testid)) {
    matchid <- as.numeric(names(sort(PS_distance))[2:(M+1)])
  } else {
    matchid <- testid[testid != treated.id]
  }

  
  second.diff <- mean(x[x$V2 %in% matchid & 
                          x$V1 == sort(unique(x$V1))[L + 1 + FORWARD], ]$V5 -
                        x[x$V2 %in% matchid & 
                            x$V1 == sort(unique(x$V1))[L], ]$V5) 
  
  return(first.diff - second.diff)
}

PS_m_weight <- function(x, L, FORWARD, M = M, d2 = d2,
                        unit.id = unit.id, time.id = time.id) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  colnames(x)[1:5] <- c("V2", "V1", "ps", "V4", "V5")
  x <- x[!duplicated(x[c("V2", "V1")]),]
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x[x$V4 == 1 & x$V1 == (max(x$V1)-FORWARD), ]$V2
  testid <- unique(x$V2)
  timeid_later <- unique(x$V1)
  timeid <- unique(x$V1)[-((length(unique(x$V1)) - FORWARD):length(unique(x$V1)))]
  matched_set <- x[x$V1 %in% timeid, ]
  
  PS_distance <- abs(tapply(matched_set$ps, matched_set$V2, mean) - mean(matched_set$ps[which(matched_set$V2 == treated.id)]))
  
  if (M < length(testid) - 1) {
    matchid <- as.numeric(names(sort(PS_distance))[2:(M+1)])
    weights <- as.data.frame(rbind(cbind(1/M, matchid), cbind(w.weight = 1, treated.id)))
  } else {
    matchid <- as.numeric(names(sort(PS_distance)))[-1]
    weights <- as.data.frame(rbind(cbind(1/(length(testid)-1), matchid), cbind(w.weight = 1, treated.id)))
  }
  
  
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  colnames(weights)[1] <- "w.weight"
  merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
  merged <- merged[order(merged$V2, merged$V1), ]
  merged$w.weight <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == treated.id, 1, 
                            ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 == treated.id, 1, 
                                   ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[testid != treated.id], merged$w.weight, 
                                          ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 %in% testid[testid != treated.id], -merged$w.weight, 0) 
                                   )))
  
  # if (MSMD == 0.0) {
  #   merged$w.weight <- 0
  # }
  new.W <- d2[c(unit.id, time.id)] # create a new data.frame as large as the dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
  total2$dit <- 0
  total2$dit[which(total2$V2 == testid[2] & total2$V1 == max(timeid_later))] <- 1
  return (list("w.weight" = total2$w.weight, "dit" = total2$dit)) # return the weight variable
} 

PS_m_weights <- function (L, FORWARD, M = 1, time.id = "year", unit.id = "ccode",
                          treatment, 
                          covariate, dependent, data, qoi = "att") 
{
  
  data <- data[c(unit.id, time.id,  treatment, dependent, covariate)]
  dlist <- lapply(1:L, 
                  function (i) slide(data = data, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                     NewVar = paste("dependent_l", i, sep="")))
  data <- Reduce(function(x, y) {merge(x, y)}, dlist)
  # to include ldvs in varnames
  varnames <- c(time.id, unit.id, treatment, dependent, c(covariate, colnames(data)[(4 + length(covariate) + 1):length(data)]))
  d2 <- na.omit(data[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x)))
  d2 <- d2[order(d2$unit, d2$time), ]
  dmatrix <- as.matrix(d2)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L, 
                                                                       F = FORWARD, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  
  # take the last time period from each subset:
  Fs <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 == max(unique(x$V1)), ]
    return(x)
  })
  
  # to only include L and the first treatment period
  even_smaller1 <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 %in% sort(unique(x$V1))[1:(L+1)], ]
    return(x)
  })
  
  # add varnames
  even_smaller1 <- lapply(even_smaller1, function (x) {
    colnames(x) <- varnames
    return(x)
  })
  
  # add varnames to the last-time-period-subset
  Fs <- lapply(Fs, function (x) {
    colnames(x) <- varnames
    return(x)
  })
  
  pooled <- rbindlist(even_smaller1) # get a dataset for propensity score generation
  
  # get propensity scores
  fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
              family = binomial(link = "logit"), data = pooled)
  pooled$ps <- fit0$fitted.values
  
  # aggregate to delete duplicates
  aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
  
  newlist <- lapply(even_smaller1, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
  newlist <- Map(rbind.fill, newlist, Fs)
  # superlist <- mapply(function(x, y) rbind.fill(x,y), newlist, Fs)
  
  # calibrating weights
  weights_and_dits <- lapply(newlist, PS_m_weight, L = L, FORWARD = FORWARD, M = M,
                             d2 = d2, unit.id = unit.id, time.id = time.id)
  all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
  all.dits <- (lapply(weights_and_dits, function (x) x$dit))
  d2$weights_att <- Reduce("+", all.weights)
  d2$dit_att <- Reduce("+", all.dits)
  if (qoi == "att") {
    return(d2)
  } else {
    
    dmatrix[, 3] <- ifelse(dmatrix[, 3] == 1, 0, 1)
    smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L = L, 
                                                                         F = FORWARD, dmatrix)), delete.NULLs)
    smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
    even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
    
    Fs <- lapply(even_smaller2, function(x) {
      x <- x[x$V1 == max(unique(x$V1)), ]
      return(x)
    })
    
    even_smaller2 <- lapply(even_smaller2, function(x) {
      x <- x[x$V1 %in% sort(unique(x$V1))[1:(L+1)], ]
      return(x)
    })
    
    even_smaller2 <- lapply(even_smaller2, function (x) {
      colnames(x) <- varnames
      return(x)
    })
    
    Fs <- lapply(Fs, function (x) {
      colnames(x) <- varnames
      return(x)
    })
    
    pooled <- rbindlist(even_smaller2)
    
    fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
                family = binomial(link = "logit"), data = pooled)
    pooled$ps <- fit0$fitted.values
    
    aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
    newlist <- lapply(even_smaller2, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
    newlist <- Map(rbind.fill, newlist, Fs)
    # superlist <- mapply(function(x, y) rbind.fill(x,y), newlist, Fs)
    
    # calibrating weights
    weights_and_dits <- lapply(newlist, PS_m_weight, L = L, FORWARD = FORWARD, M = M,
                               d2 = d2, unit.id = unit.id, time.id = time.id)
    all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    d2$weights_atc <- Reduce("+", all.weights)
    d2$dit_atc <- Reduce("+", all.dits)
    
    return(d2)
  }
}



PS_plot <- function(x, L, FORWARD, M = M,
                    show.covariate = FALSE,
                    post.treatment = FALSE) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  colnames(x)[1:5] <- c("V2", "V1", "ps", "V4", "V5")
  x <- x[!duplicated(x[c("V2", "V1")]),]
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x[x$V4 == 1 & x$V1 == (max(x$V1)-FORWARD), ]$V2
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-((length(unique(x$V1)) - FORWARD):length(unique(x$V1)))]
  matched_set <- x[x$V1 %in% timeid, ]
  
  PS_distance <- abs(tapply(matched_set$ps, matched_set$V2, mean) - mean(matched_set$ps[which(matched_set$V2 == treated.id)]))
  
  if (M < length(testid)) {
    matchid <- as.numeric(names(sort(PS_distance))[2:(M+1)])
  } else {
    matchid <- testid[testid != treated.id]
  }
  
  if (post.treatment == FALSE) {
    if (show.covariate == FALSE) {
      return(list("gap" = as.vector(tapply(matched_set$V5[which(matched_set$V2 %in% matchid)],
                                           matched_set[matched_set$V2 %in% matchid, ]$V1, mean)) -
                    x[x$V2 == treated.id & x$V1 %in% timeid, ]$V5, 
                  "unit.id" = paste(treated.id, unique(x$V1)[(L+FORWARD+1)], sep = ",")))
    } else {
      return(list("gap" = as.vector(tapply(matched_set$x[which(matched_set$V2 %in% matchid)],
                                           matched_set[matched_set$V2 %in% matchid, ]$V1, mean)) -
                    x[x$V2 == treated.id & x$V1 %in% timeid, ]$x, 
                  "unit.id" = paste(treated.id, unique(x$V1)[(L+FORWARD+1)], sep = ",")))
    }
  } else {
    if (show.covariate == FALSE) {
      return(list("gap" = as.vector(tapply(x$V5[which(x$V2 %in% matchid)], 
                                           x$V1[which(x$V2 %in% matchid)], mean)) -
                    x$V5[which(x$V2 == treated.id)],
                  "unit.id" = paste(treated.id, unique(x$V1)[(L+FORWARD+1)], sep = ",")))
    } else {
      return(list("gap" = as.vector(tapply(x$x[which(x$V2 %in% matchid)], 
                                           x$V1[which(x$V2 %in% matchid)], mean)) -
                    x$x[which(x$V2 == treated.id)],
                  "unit.id" = paste(treated.id, unique(x$V1)[(L+FORWARD+1)], sep = ",")))
    }
  }
  
}



