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
      return(merged)
    }
  } else {
    merged <- x
    merged$w.weight <- 1
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
    return(merged)
  } else {
    merged <- x
    merged$w.weight <- 1
    return(merged)
  }
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






