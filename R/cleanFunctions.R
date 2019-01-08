gaps_plot_tmp <- function(x, lag, lead, data, dependent,
                          qoi, refinement = TRUE, 
                          covariate_names,
                          method,
                          treated_set,
                          covariate = NULL) {
  treated_set <- treated_set
  x <- as.data.frame(x)
  # colnames(data) <- c("time.id", "unit.id",
  #                     "treatment", "dependent", covariate_names)
  # if (method == "Maha"|method == "Synth") {
  #   x <- as.data.frame(x)
  #   colnames(x)[5:length(x)] <- c(covariate_names,
  #                                 "dependent",
  #                                 "w.weight")
  #   colnames(treated_set)[5:(length(treated_set)-3)] <- 
  #     covariate_names
  # }
  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
  if (is.null(covariate)) {
    treated_set$V5 <- treated_set$V4
    if (refinement == TRUE) {
      gap <- x$V4[x$V2 == treated.id] - 
        tapply(x$V4[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
               x$V1[x$V2 != treated.id], sum)
    } else {
      gap <- x$V4[x$V2 == treated.id] - 
        tapply(x$V4[x$V2 != treated.id], 
               x$V1[x$V2 != treated.id], mean)
    }
    
  } else {
    x$V5 <- as.numeric(x[,c(covariate)])
    treated_set$V5 <- as.numeric(treated_set[,c(covariate)])
    if (refinement == TRUE){
      gap <- x$V5[x$V2 == treated.id] - 
        tapply(x$V5[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
               x$V1[x$V2 != treated.id], sum)
    } else {
      gap <- x$V5[x$V2 == treated.id] - 
        tapply(x$V5[x$V2 != treated.id], 
               x$V1[x$V2 != treated.id], mean)
    }
    
  }
  # if (is.null(covariate) == FALSE) {
  #   data$dependent <- data[covariate][,1]
  # }
  # divide the gap by the sd of the outcome variable among all treated units 
  # in that treatment time period.
  if (qoi == "att") {
    # overall <- rep(NA, (lag+1+lead))
    overall <- tapply(treated_set$V5, treated_set$big_L, sd, na.rm = T)
    # sub.data <- data[which(data$time.id <= max(x$V1) & 
    #                          data$time.id >= (max(x$V1)-lead-lag)),]
    # index.l <- 
    #   as.numeric(rownames(sub.data[which(sub.data$time.id == (max(x$V1)-lead) & sub.data$treatment == 1), ]))
    # for (i in 1:(lag + 1 + lead)) {
    #   overall[i] <- sd(sub.data$dependent[rownames(sub.data) %in% (index.l-lag-1 + i)])
    #   overall[i] <- ifelse(overall[i] == 0, NA, overall[i]) # prevent inf
    # }
    
  } else {
    overall <- tapply(treated_set$V5, treated_set$big_L, sd, na.rm = T)
    # sub.data <- data[which(data$time.id <= max(x$V1) & 
    #                          data$time.id >= (max(x$V1)-lead-lag)),]
    # index.l <- 
    #   as.numeric(rownames(sub.data[which(sub.data$time.id == (max(x$V1)-lead) & 
    #                                                 sub.data$treatment == 0), ]))
    # for (i in 1:(lag + 1 + lead)) {
    #   overall[i] <- sd(sub.data$dependent[rownames(sub.data) %in% (index.l-lag-1 + i)])
    #   overall[i] <- ifelse(overall[i] == 0, NA, overall[i]) # prevent inf
    # }
  }
  
  return(list("gap" = gap/overall,
              "unit" = paste(treated.id, unique(x$V1)[lag + 1], sep = ",")))
}


## Caliper
gaps_caliper <- function(x, lag, covariate = NULL, data,
                         balance_type,
                         qoi) {
  
  colnames(data) <- c("time.id", "unit.id",
                      "treatment", "dependent", covariate)
  
  
  treated.id <- x$V2[which(x$V3 == 1 & x$V1 == (min(x$V1)+lag))] # check this
  x <- x[x$V1 < (min(x$V1)+lag), ]
  
  if (is.null(covariate)) {
    gap <- x$V4[x$V2 == treated.id] - 
      tapply(x$V4[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
             x$V1[x$V2 != treated.id], sum)
  } else {
    gap <- x$V5[x$V2 == treated.id] - 
      tapply(x$V5[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
             x$V1[x$V2 != treated.id], sum)
  }
  
  if (is.null(covariate) == FALSE) {
    data$dependent <- data[c(covariate)]
  }
  if (qoi == "att") {
    # for all the time periods that matched with pre-treatment L and 
    # and are treated, gimme the standard deviation
    overall <- tapply(data$dependent[data$time.id %in% unique(x$V1) & data$treatment == 1],
                      data$time.id[data$time.id %in% unique(x$V1) & data$treatment == 1], sd)
  } else {
    overall <- tapply(data$dependent[data$time.id %in% unique(x$V1) & data$treatment == 0],
                      data$time.id[data$time.id %in% unique(x$V1) & data$treatment == 0], sd)
  }
  if (balance_type == "parallel") {
    return(abs(mean(gap/overall)))
  } else {
    return(abs(sd(gap/overall)))
  }
  
  
}


## WANT TO DELETE THESE AFTER PANELESTIMATE2 IS FINISHED
### updating (adding) wits
apply_leads <- function(data, unit.id, time.id, matched_set,
                        lag, lead, estimator = c("did", "matching"),
                        inference = c("bootstrap", "wfe")) {
  unit.id = unit.id; time.id = time.id;
  #inference = inference; estimator = estimator
  # set testid and timeid
  treated.time <- min(matched_set$V1) + lag
  treated.id <- matched_set[matched_set$V3 == 1 & 
                              matched_set$V1 == treated.time, ]$V2
  
  testid <- unique(matched_set$V2)
  
  if (estimator == "did" & inference == "bootstrap") {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
                              ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, -1, 
                                     ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 
                                            ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 0) 
                                     )))
    
  } else if (estimator == "did" & inference == "wfe") {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
                              ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, 1, 
                                     ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 
                                            ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 0) 
                                     )))
  } else if (estimator == "matching" & inference == "bootstrap") {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
           ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 != treated.id, -matched_set$w.weight, 0))
                 
  } else if (estimator == "matching" & inference == "wfe")  {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
                              ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 != treated.id, matched_set$w.weight, 0))
  }

  
  new.W <- data[c(unit.id, time.id)] # create a new data.frame as large as the *cleaned* dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, matched_set, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$wit <- ifelse(is.na(total2$wit), 0, total2$wit) # turn NAs into zero
  total2$dit <- 0
  # IMPORTANT: CHANGE HERE IF THINGS GO WRONG FOR lead >0
  total2$dit[which(total2$V2 == treated.id & total2$V1 == treated.time)] <- 1
  return (list("wit" = total2$wit, "dit" = total2$dit)) 
}

# this will give us wit and dit at the dataset level for each matched set
lapply_leads <- function(x, leads, data, unit.id = unit.id, time.id = time.id, lag,
                         estimator = estimator, inference = inference){
  return(lapply(leads, apply_leads, unit.id = unit.id, time.id = time.id,
         matched_set = x, lag = lag, estimator = estimator, inference = inference,
         data = data))
}

extract_objects <- function(x, objective = c("wit", "dit")) {
  if (objective == "wit") {
    lapply(x, function(x) x$wit)
  } else {
    lapply(x, function(x) x$dit)
  }
  
}

# to eliminate "[[1]]"
each_lead <- function(x, lead) {
  lapply(x, function(a) a[[lead]])
}



