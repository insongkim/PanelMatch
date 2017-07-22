Panel_vit <- function(x, lag, lead, scheme, M) {
  if (scheme == "Synth") {
    return(synth_vit(x, lag = lag, lead = lead))
  } else if(scheme == "Maha"){
    return(Maha_vit(x, lag, lead, M = M))
  } else if(scheme == "Pscore") {
    return(PS_vit(x, lag, lead, M = M))
  } else {
    return("WRONG")
  }
}



synth_vit <- function(x, lag, lead) {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  covariate.names <- colnames(x)[5:length(x)]
  if (nrow(x) > 2*(lag + lead + 1)) {
    if (is.na(colnames(x)[5:length(x)][1])) {
      V2 = x$V2; V1 = x$V1
      control_data <- reshape2::dcast(x[V2 != testid[2], ], V2 ~ V1)
      
      # treat_data <- x %>% 
      #   filter(V2 == testid[2]) %>% 
      #   select_("V4") %>% 
      #   as.data.frame()
      
      treat_data <- as.data.frame(x$V4[which(x$V2 == testid[2])])
      
      synth_out <- synth_constReg_weight(
        Y_t = as.vector(treat_data[,1]), 
        Y_c = as.matrix(control_data[,-1]), 
        T0 = (lag)
      )
      
      weights <- as.data.frame(rbind(cbind(synth_out$weight, testid[-2]), cbind(w.weight = 1, testid[2])))
      colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
      merged <- merge(x, weights, by = "V2") # merge it with the data.frame
      merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                                ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 == testid[2],1, 
                                       ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                              ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                       )))
      return(merged)
    } else {
      dataprep.out <- dataprep(foo = x, 
                               dependent = "V4",
                               unit.variable = "V2",
                               # unit.names.variable = "unit.name",
                               time.variable = "V1",
                               treatment.identifier = testid[2], # the regionno of the treated unit
                               controls.identifier = testid[-2],
                               time.optimize.ssr = min(timeid):max(timeid-lead-1), # the pre-treatment preiod
                               time.predictors.prior = min(timeid):max(timeid-lead-1),
                               predictors = covariate.names)
      synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
      # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
      # then making it a data.frame so as to merge it in the next step
      weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
      colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
      merged <- merge(x, weights, by = "V2") # merge it with the data.frame
      merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                                ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 == testid[2],1, 
                                       ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                              ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                       )))
      return(merged)
    }
  } else {
    merged <- x
    merged$w.weight <- 1
    merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid) - lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    return(merged)
  }
}

Maha_vit <- function(x, lag, lead, M = 3) {
  testid <- unique(x$V2)
  timeid_later <- unique(x$V1)
  timeid <- timeid_later[1:lag]
  
  if (nrow(x) > 2*(lag + lead + 1)) {
    matched_set <- x[x$V1 %in% timeid, ]
    MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    MSMD <- append(Reduce("+", MSMDlist)/length(timeid), 0, after = 1)
    
    # first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + lead] - x[x$V2 == testid[2], ]$V4[L]
    
    if (M < length(testid)-1) {
      matchid <- order(MSMD)[2:(M+1)]
      weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, testid[2])))
    } else {
      matchid <- order(MSMD)
      weights <- as.data.frame(rbind(cbind(1/(length(testid)-1), testid[matchid[-1]]), cbind(w.weight = 1, testid[2])))
      
    } 
    
    # second.diff <- mean(x[x$V2 %in% testid[matchid] & 
    #                         x$V1 == unique(x$V1)[L + 1 + lead], ]$V4 -
    #                       x[x$V2 %in% testid[matchid] & 
    #                           x$V1 == unique(x$V1)[L], ]$V4, na.rm = T)
    # 
    # 
    
    # weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, testid[2])))
    
    colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
    merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
    merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    return(merged)
  } else {
    merged <- x
    merged$w.weight <- 1
    merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    return(merged)
  }
}


PS_vit <- function(x, lag, lead, M = M) {

  colnames(x)[1:5] <- c("V2", "V1", "ps", "V4", "V5")
  x <- x[!duplicated(x[c("V2", "V1")]),]
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x[x$V4 == 1 & x$V1 == (max(x$V1)-lead), ]$V2
  testid <- unique(x$V2)
  timeid_later <- unique(x$V1)
  timeid <- unique(x$V1)[-((length(unique(x$V1)) - lead):length(unique(x$V1)))]
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
  colnames(merged)[c(4,5)] <- c("V3", "V4")
  merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == treated.id, 1, 
                            ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 == treated.id, 1, 
                                   ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[testid != treated.id], merged$w.weight, 
                                          ifelse(merged$V1 == max(timeid_later) - lead - 1 & merged$V2 %in% testid[testid != treated.id], -merged$w.weight, 0) 
                                   )))
  return(merged) # return the weight variable
} 


PanelWit <- function(data, unit.id, time.id, matched_set,
                     lag, lead) {
  # set testid and timeid
  treated.time <- max(matched_set$V1)-lead
  treated.id <- matched_set[matched_set$V3 == 1 & 
                              matched_set$V1 == treated.time, ]$V2
  
  new.W <- data[c(unit.id, time.id)] # create a new data.frame as large as the *cleaned* dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, matched_set, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$wit <- ifelse(is.na(total2$wit), 0, total2$wit) # turn NAs into zero
  total2$dit <- 0
  # IMPORTANT: CHANGE HERE IF THINGS GO WRONG FOR lag >0
  total2$dit[which(total2$V2 == treated.id & total2$V1 == treated.time)] <- 1
  ######################################################
  return (list("wit" = total2$wit, "dit" = total2$dit)) 
}



PanelDiDResult <- function(x, lag, lead){
  testid <- unique(x$V2)
  location <- match(x$V2[x$V3 == 1 & x$V1 == (max(x$V1)-lead)], testid)
  
  first.diff <- x[x$V2 == testid[location], ]$V4[lag+1+lead] - x[x$V2 == testid[location], ]$V4[lag]
  second.diff.1 <- sum(x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag+1+lead], ]$V4*
                         x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag+1+lead], ]$w.weight)
  second.diff.2 <- sum(x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag], ]$V4*
                         x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag], ]$w.weight)
  second.diff <- second.diff.1 - second.diff.2  
  return(first.diff - second.diff)
  
}


# newlist$`Matched sets for ATC`[[1]]
# PanelDiDResult(x = merged, 
#                lag = 4, lead = 2)








