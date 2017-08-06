MSMD_result <- function(x, L, FORWARD, M = 3) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[-length(unique(x$V1))]
  
  if (nrow(x) > 2*(L+FORWARD+1)) {
    matched_set <- x[x$V1 %in% timeid, ]
    MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    MSMD <- append(Reduce("+", MSMDlist)/length(timeid), 0, after = 1)
    
    first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + FORWARD] - x[x$V2 == testid[2], ]$V4[L]
    
    if (M < length(testid)) {
      matchid <- order(MSMD)[2:(M+1)]
    } else {
      matchid <- order(MSMD)[-1]
    }
    
    second.diff <- mean(x[x$V2 %in% testid[matchid] & 
                            x$V1 == unique(x$V1)[L + 1 + FORWARD], ]$V4 -
                          x[x$V2 %in% testid[matchid] & 
                              x$V1 == unique(x$V1)[L], ]$V4, na.rm = T) 
    return(first.diff - second.diff)
  } else {
    first.diff <- x[x$V2 == testid[2], ]$V4[L+1+FORWARD] - x[x$V2 == testid[2], ]$V4[L]
    second.diff <- x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+FORWARD], ]$V4 -
      x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$V4
    return(first.diff - second.diff)
  }
  
}

MSMD_each <- function(time_id, matched_set, testid) {
  sub_sub <- matched_set[matched_set$V1 == time_id, ]
  cov_matrix <- cov(as.matrix(sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)]))
  if(all(eigen(cov_matrix)$values > .00001)) {
    if (length(sub_sub) == 4) {
      return((sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)] - 
                as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)])) * cov_matrix^(-1) *
               (sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)] - 
                  as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)]))
      )
    } else {
      return(mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)], 
                         center = as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)]),
                         cov = cov_matrix))
    }
   
  } else {
    if (length(sub_sub) == 4) {
      return((sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)] - 
                as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)])) * cov_matrix^(-1) *
               (sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)] - 
                  as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)])))
    } else {
      return(tryCatch(mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)], 
                         center = as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)]),
                         cov = diag((length(sub_sub) - 3)) * cov_matrix), 
             error = function(e) mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[-2], 4:length(sub_sub)], 
                                             center = as.numeric(sub_sub[sub_sub$V2 == testid[2], 4:length(sub_sub)]),
                                             cov = diag((length(sub_sub) - 3)))))
    }
   
  }
}




MSMD_weight <- function(x, L, FORWARD, M = 3, d2 = d2, unit.id = unit.id, time.id = time.id) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  
  testid <- unique(x$V2)
  timeid_later <- unique(x$V1)
  timeid <- timeid_later[-length(timeid_later)]
  
  if (nrow(x) > 2*(L+FORWARD+1)) {
    matched_set <- x[x$V1 %in% timeid, ]
    MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    MSMD <- append(Reduce("+", MSMDlist)/length(timeid), 0, after = 1)
    
    # first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + FORWARD] - x[x$V2 == testid[2], ]$V4[L]
    
    if (M < length(testid)-1) {
      matchid <- order(MSMD)[2:(M+1)]
      weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, testid[2])))
    } else {
      matchid <- order(MSMD)
      weights <- as.data.frame(rbind(cbind(1/(length(testid)-1), testid[matchid[-1]]), cbind(w.weight = 1, testid[2])))
      
    } 
    
    # second.diff <- mean(x[x$V2 %in% testid[matchid] & 
    #                         x$V1 == unique(x$V1)[L + 1 + FORWARD], ]$V4 -
    #                       x[x$V2 %in% testid[matchid] & 
    #                           x$V1 == unique(x$V1)[L], ]$V4, na.rm = T)
    # 
    # 
    
    # weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, testid[2])))
    
    colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
    merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
    merged$w.weight <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
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
  } else {
    merged <- x
    merged$w.weight <- 1
    merged$w.weight <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
                              ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 == testid[2],1, 
                                     ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
                                            ifelse(merged$V1 == max(timeid_later) - FORWARD - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                     )))
    new.W <- d2[c(unit.id, time.id)] # create a new data.frame as large as the dataset
    names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
    total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
    total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
    total2$dit <- 0
    total2$dit[which(total2$V2 == testid[2] & total2$V1 == max(timeid_later))] <- 1
    return (list("w.weight" = total2$w.weight, "dit" = total2$dit)) # return the weight variable
  }
  
}

MSMD_DID_weights <- function(L, FORWARD, time.id = "year", qoi = "ate",
                             unit.id = "ccode", M = 3,
                             treatment, covariate, dependent, d) {
  
  d <- d # set dataset
  
  varnames <- c(time.id, unit.id, treatment, dependent, covariate)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
  
  ### finding matches using the cpp function ###
  # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
  ### finding matches using the cpp function ###
  
  ### cleaning the output from cpp ###
  # delete both higher level and lower level null entries
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = L, F = FORWARD, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  # # subset out any dataframe that have 2 or fewer than 2 units
  # even_smaller1 <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), even_smaller1)
  # 
  # # only focus on ATT
  # even_smaller1 <- Filter(function (x) x[x$V2 == unique(x$V2)[2] & x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
  # 
  # d.sum1 <- length(even_smaller1) 
  
  # calibrating weights
  weights_and_dits <- lapply(even_smaller1, MSMD_weight, L = L, FORWARD = FORWARD, M = M,
                             d2 = d2, unit.id = unit.id, time.id = time.id)
  all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
  all.dits <- (lapply(weights_and_dits, function (x) x$dit))
  d2$weights_att <- Reduce("+", all.weights)
  d2$dit_att <- Reduce("+", all.dits)
  if (qoi == "att") {
    return(d2)
  } else {
    ### from 1 to zero ###
    dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
    
    ### finding matches using the cpp function ###
    # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
    ### finding matches using the cpp function ###
    
    ### cleaning the output from cpp ###
    # delete both higher level and lower level null entries
    smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = L, F = FORWARD, dmatrix)), delete.NULLs) 
    # further cleaning
    smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
    # use function dframelist.rb_dup to turn every list element into a data.frame
    even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
    # subset out any dataframe that have 2 or fewer than 2 units
    # even_smaller2 <- Filter(function (x) nrow(x) > 2*(L+FORWARD+1), even_smaller2)
    # 
    # # only focus on ATC
    # even_smaller2 <- Filter(function (x) x[x$V2 == unique(x$V2)[2] & x$V1 == unique(x$V1)[L], ]$V3 == 0, smallerlist)
    # 
    # d.sum2 <- length(even_smaller2) 
    
    # calibrating weights
    weights_and_dits <- lapply(even_smaller2, MSMD_weight, L = L, FORWARD = FORWARD, M = M,
                               d2 = d2, unit.id = unit.id, time.id = time.id)
    all.weights <- (lapply(weights_and_dits, function (x) x$w.weight))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    d2$weights_atc <- Reduce("+", all.weights)
    d2$dit_atc <- Reduce("+", all.dits)
    
    return(d2)
  }
  
}

MSMD_plot <- function(x, L, FORWARD, M = 3, unit.id = unit.id, time.id = time.id,
                      post.treatment = FALSE, show.covariate = FALSE) {
  L <- L
  FORWARD <- FORWARD
  M <- M
  
  testid <- unique(x$V2)
  timeid_later <- unique(x$V1)
  timeid <- timeid_later[1:L]
  
  if (nrow(x) > 2*(L+FORWARD+1)) {
    matched_set <- x[x$V1 %in% timeid, ]
    MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    MSMD <- append(Reduce("+", MSMDlist)/length(timeid), 0, after = 1)
    
    # first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + FORWARD] - x[x$V2 == testid[2], ]$V4[L]
    
    if (M < length(testid)-1) {
      matchid <- order(MSMD)[2:(M+1)]
    } else {
      matchid <- order(MSMD)
    } 
    
    if(post.treatment == FALSE)
    {
      if (show.covariate == FALSE) {
        # taking average with M
        return(list("gaps" = (as.vector(tapply(matched_set[matched_set$V2 %in% testid[matchid], ]$V4,
                                              matched_set[matched_set$V2 %in% testid[matchid], ]$V1,
                                              mean)) - x$V4[which(x$V2 == testid[2] &
                                                                    x$V1 %in% timeid_later[1:L])])/sd(x$V4[which(x$V2 == testid[2] & x$V1 %in% timeid)]),
                    "unit.id" = paste(testid[2], timeid_later[L + FORWARD + 1], sep = ",")))
        
      } else {
        return(list("gaps" = (as.vector(tapply(matched_set[matched_set$V2 %in% testid[matchid], ]$V5,
                                              matched_set[matched_set$V2 %in% testid[matchid], ]$V1,
                                              mean)) - x$V5[which(x$V2 == testid[2] &
                                                                    x$V1 %in% timeid_later[1:L])])/sd(x$V5[which(x$V2 == testid[2] & x$V1 %in% timeid)]),
                    "unit.id" = paste(testid[2], timeid_later[L + FORWARD + 1], sep = ",")))
        
        
      }
      
    } else {
      if (show.covariate == FALSE) {
        return(list("gaps" = (as.vector(tapply(x[x$V2 %in% testid[matchid], ]$V4,
                                              x[x$V2 %in% testid[matchid], ]$V1,
                                              mean)) - x$V4[which(x$V2 == testid[2])])/sd(x$V4[which(x$V2 == testid[2] & x$V1 %in% timeid)]),
                    "unit.id" = paste(testid[2], timeid_later[L + FORWARD + 1], sep = ",")))
      } else {
        return(list("gaps" = (as.vector(tapply(x[x$V2 %in% testid[matchid], ]$V5,
                                              x[x$V2 %in% testid[matchid], ]$V1,
                                              mean)) - x$V5[which(x$V2 == testid[2])])/sd(x$V5[which(x$V2 == testid[2] & x$V1 %in% timeid)]),
                    "unit.id" = paste(testid[2], timeid_later[L + FORWARD + 1], sep = ",")))
        
      }
    }
    
    
  } else {
    return(NULL)
  }
  
}

