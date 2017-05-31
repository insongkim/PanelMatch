cscwdid <- function(x, L, FORWARD) {
  L <- L; FORWARD <- FORWARD
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  dataprep.out <- dataprep(foo = x, 
                           dependent = "V5",
                           unit.variable = "V2",
                           # unit.names.variable = "unit.name",
                           time.variable = "V1",
                           treatment.identifier = testid[2], # the regionno of the treated unit
                           controls.identifier = testid[-2],
                           time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                           time.predictors.prior = min(timeid):max(timeid-F-1),
                           predictors = "V4")
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
  # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
  # then making it a data.frame so as to merge it in the next step
  weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  x <- merge(x, weights, by = "V2") # merge it with the data.frame (
  first.diff <- x[x$V2 == testid[2], ]$V5[L+1+F] - x[x$V2 == testid[2], ]$V5[L]
  second.diff.1 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$V5*
                         x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$w.weight)
  second.diff.2 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$V5*
                         x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$w.weight)
  second.diff <- second.diff.1 - second.diff.2                    
  return(first.diff - second.diff)
}

cscwdid_tmp <- function(x, L, FORWARD) {
  L <- L; FORWARD <- FORWARD
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  if (nrow(x) > 2*(L+F+1)) {
    dataprep.out <- dataprep(foo = x, 
                             dependent = "V5",
                             unit.variable = "V2",
                             # unit.names.variable = "unit.name",
                             time.variable = "V1",
                             treatment.identifier = testid[2], # the regionno of the treated unit
                             controls.identifier = testid[-2],
                             time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                             time.predictors.prior = min(timeid):max(timeid-F-1),
                             predictors = "V4")
    synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
    # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
    # then making it a data.frame so as to merge it in the next step
    weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
    colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
    x <- merge(x, weights, by = "V2") # merge it with the data.frame (
    first.diff <- x[x$V2 == testid[2], ]$V5[L+1+F] - x[x$V2 == testid[2], ]$V5[L]
    second.diff.1 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$V5*
                           x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$w.weight)
    second.diff.2 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$V5*
                           x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$w.weight)
    second.diff <- second.diff.1 - second.diff.2                    
    return(first.diff - second.diff)
  } else {
    first.diff <- x[x$V2 == testid[2], ]$V5[L+1+F] - x[x$V2 == testid[2], ]$V5[L]
    second.diff <- x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$V5 -
      x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$V5
    return(first.diff - second.diff)
  }
}



### writing a function to apply synthetic control and generate a weight vector of the same length as the dataset
callSynth <- function (x, unit.id, time.id, d2) {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  dataprep.out <- dataprep(foo = x, 
                           dependent = "V5",
                           unit.variable = "V2",
                           # unit.names.variable = "unit.name",
                           time.variable = "V1",
                           treatment.identifier = testid[2], # the regionno of the treated unit
                           controls.identifier = testid[-2],
                           time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                           time.predictors.prior = min(timeid):max(timeid-F-1),
                           predictors = "V4")
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
  # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
  # then making it a data.frame so as to merge it in the next step
  weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
  merged$w.weight <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                            ifelse(merged$V1 == max(timeid) - F - 1 & merged$V2 == testid[2],1, 
                                   ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                          ifelse(merged$V1 == max(timeid) - F - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                   )))
  new.W <- d2[c(unit.id, time.id)] # create a new data.frame as large as the dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
  total2$dit <- 0
  total2$dit[which(total2$V2 == testid[2] & total2$V1 == max(timeid))] <- 1
  return (list("w.weight" = total2$w.weight, "dit" = total2$dit)) # return the weight variable
} 

callSynth_tmp <- function (x, unit.id, time.id, L, FORWARD, d2) {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  if (nrow(x) > 2*(L+FORWARD+1)) {
    dataprep.out <- dataprep(foo = x, 
                             dependent = "V5",
                             unit.variable = "V2",
                             # unit.names.variable = "unit.name",
                             time.variable = "V1",
                             treatment.identifier = testid[2], # the regionno of the treated unit
                             controls.identifier = testid[-2],
                             time.optimize.ssr = min(timeid):max(timeid-FORWARD-1), # the pre-treatment preiod
                             time.predictors.prior = min(timeid):max(timeid-FORWARD-1),
                             predictors = "V4")
    synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
    # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
    # then making it a data.frame so as to merge it in the next step
    weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
    colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
    merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
    merged$w.weight <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid) - FORWARD - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid) - FORWARD - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    new.W <- d2[c(unit.id, time.id)] # create a new data.frame as large as the dataset
    names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
    total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
    total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
    total2$dit <- 0
    total2$dit[which(total2$V2 == testid[2] & total2$V1 == max(timeid))] <- 1
    return (list("w.weight" = total2$w.weight, "dit" = total2$dit)) # return the weight variable
  } else {
    merged <- x
    merged$w.weight <- 1
    merged$w.weight <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid) - FORWARD - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid) - FORWARD - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    new.W <- d2[c(unit.id, time.id)] # create a new data.frame as large as the dataset
    names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
    total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
    total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
    total2$dit <- 0
    total2$dit[which(total2$V2 == testid[2] & total2$V1 == max(timeid))] <- 1
    return (list("w.weight" = total2$w.weight, "dit" = total2$dit)) # return the weight variable
  }
}



# diagnostics
cscwplot <- function(x, L, FORWARD, Main = "ATT") {
  L <- L; FORWARD <- FORWARD
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  if (nrow(x) > 2*(L+FORWARD+1)) {
    dataprep.out <- dataprep(foo = x, 
                             dependent = "V5",
                             unit.variable = "V2",
                             # unit.names.variable = "unit.name",
                             time.variable = "V1",
                             treatment.identifier = testid[2], # the regionno of the treated unit
                             controls.identifier = testid[-2],
                             time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                             time.predictors.prior = min(timeid):max(timeid-F-1),
                             predictors = "V4")
    synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
    # plot the pre-treatment gaps
    return(list("gap" = x$V5[which(x$V2 == testid[2] & x$V1 %in% timeid[-length(timeid)])] - dataprep.out$Y0plot %*% synth.out$solution.w, 
                "unit.id" = paste(testid[2], timeid[L + FORWARD + 1], sep = ",")))
  } else {
    return(NULL)
  }
}

cscMSPE <- function(x, L, FORWARD) {
  L <- L; FORWARD <- FORWARD
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  if (nrow(x) > 2*(L+FORWARD+1)) {
    dataprep.out <- dataprep(foo = x, 
                             dependent = "V5",
                             unit.variable = "V2",
                             # unit.names.variable = "unit.name",
                             time.variable = "V1",
                             treatment.identifier = testid[2], # the regionno of the treated unit
                             controls.identifier = testid[-2],
                             time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                             time.predictors.prior = min(timeid):max(timeid-F-1),
                             predictors = "V4")
    synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
    # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
    # then making it a data.frame so as to merge it in the next step
    return(synth.out$loss.v)} else {
      return(NULL)
    }
}

