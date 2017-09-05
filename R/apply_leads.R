### updating (adding) wits
apply_leads <- function(data, unit.id, time.id, matched_set,
                        lag, lead) {
  unit.id = unit.id; time.id = time.id
  # set testid and timeid
  treated.time <- min(matched_set$V1) + lag
  treated.id <- matched_set[matched_set$V3 == 1 & 
                              matched_set$V1 == treated.time, ]$V2
  
  testid <- unique(matched_set$V2)
  
  
  matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
                              ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, -1, 
                                     ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 
                                            ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 0) 
                                     )))
  
  
  new.W <- data[c(unit.id, time.id)] # create a new data.frame as large as the *cleaned* dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, matched_set, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$wit <- ifelse(is.na(total2$wit), 0, total2$wit) # turn NAs into zero
  total2$dit <- 0
  # IMPORTANT: CHANGE HERE IF THINGS GO WRONG FOR lead >0
  total2$dit[which(total2$V2 == treated.id & total2$V1 == treated.time)] <- 1
  return (list("wit" = total2$wit, "dit" = total2$dit)) 
}

lapply_leads <- function(x, leads, data, unit.id = unit.id, time.id = time.id, lag){
  return(lapply(leads, apply_leads, unit.id = unit.id, time.id = time.id,
         matched_set = x, lag = lag,
         data = data))
}

extract_objects <- function(x, objective = c("wit", "dit")) {
  if (objective == "wit") {
    lapply(x, function(x) x$wit)
  } else {
    lapply(x, function(x) x$dit)
  }
  
}

each_lead <- function(x, lead) {
  lapply(x, function(a) a[[lead]])
}
