Panel_vit <- function(x, lag, max.lead, method, M, weighting) {
  if (method == "Synth"|method == "SynthPscore"|method == "SynthCBPS") {
    return(synth_vit(x, lag = lag, max.lead = max.lead, method = method))
  } else if(method == "Maha"){
    return(Maha_vit(x, lag, max.lead, M = M))
  } else if(method == "Pscore"|method == "CBPS") {
    return(PS_vit(x, lag, max.lead, M = M, weighting = weighting))
  } else {
    return("WRONG")
  }
}

# get_treated_set <- function(x) {
#   treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
#   return(x[x$V2 == treated.id,])
#   
# }

# synth_constReg_weight is borrowed from the package "synthR" written 
# by Soichiro Yamauchi 
synth_constReg_weight <- function (Y_t, Y_c, T0, init = NULL, method = "solnp") 
{
  obj_func <- function(par, Y_t, Y_c) {
    W_c <- par
    value <- sum((Y_t - apply(Y_c, 2, function(x) x %*% W_c))^2)
    value <- value/1000
    return(value)
  }
  const_func <- function(par, Y_t, Y_c) {
    sum(par)
  }
  N_c <- nrow(Y_c)
  T_total <- ncol(Y_c)
  T_post <- T_total - T0
  Y_c_pre <- Y_c[, 1:T0]
  Y_t_pre <- Y_t[1:T0]
  convergence <- 1
  while (convergence != 0) {
    if (is.null(init)) {
      W_c_init <- rdirichlet(n = 1, alpha = rep(1, 
                                                          N_c))
    }
    else {
      W_c_init <- init
    }
    if (method == "solnp") {
      opt_out <- solnp(pars = W_c_init, fun = obj_func, 
                               Y_t = Y_t_pre, Y_c = Y_c_pre, eqfun = const_func, 
                               eqB = 1, LB = rep(0, N_c), UB = rep(1, N_c) 
      )
      MSE <- opt_out$value
      weight_est <- opt_out$par
      convergence <- ifelse(is.null(init), opt_out$convergence, 
                            0)
      if (convergence == 0) {
        cat("-----------------\n")
        cat("Converged!\n")
        cat("-----------------\n")
      }
      else {
        cat("--------------------------\n")
        cat("Did not convergence...\n")
        cat("Trying with different intial values...\n")
        cat("--------------------------\n")
      }
    }
    else if (method == "nlopt") {
      opt_out <- nloptr(x0 = W_c_init, eval_f = obj_func, 
                                lb = rep(0, N_c), ub = rep(1, N_c), eval_g_eq = const_func, 
                                opts = list(algorithm = "NLOPT_LD_AUGLAG_EQ"), 
                                Y_t = Y_t_pre, Y_c = Y_c_pre)
      MSE <- opt_out$objective
      weight_est <- opt_out$solution
      convergence <- 0
    }
  }
  Y_synth <- t(Y_c) %*% weight_est
  return(list(weight = weight_est, value = MSE, Y_synth = Y_synth, 
              init = W_c_init))
}

synth_vit <- function(x, lag, max.lead, method = "Synth") {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  if(method == "Synth") {
    covariate.names <- colnames(x)[5:length(x)]
  } else {
    covariate.names <- "ps"
  }
  
  if (nrow(x) > 2*(lag + max.lead + 1)) {
    if (is.na(colnames(x)[5:length(x)][1])) {
      V2 = x$V2; V1 = x$V1
      control_data <- dcast(x[V2 != testid[2], ], V2 ~ V1)
      
      # treat_data <- x %>% 
      #   filter(V2 == testid[2]) %>% 
      #   select_("V4") %>% 
      #   as.data.frame()
      
      treat_data <- as.data.frame(x$V4[which(x$V2 == testid[2])])
      
      invisible(capture.output(synth_out <- synth_constReg_weight(
        Y_t = as.vector(treat_data[,1]), 
        Y_c = as.matrix(control_data[,-1]), 
        T0 = (lag),
        init = rep(1/(length(testid)-1), length(testid)-1)
      )))
      
      weights <- as.data.frame(rbind(cbind(synth_out$weight, testid[-2]), cbind(w.weight = 1, testid[2])))
      colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
      merged <- merge(x, weights, by = "V2") # merge it with the data.frame
      # if (max.lead > 0) {
      #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
      #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],-1, 
      #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], -merged$w.weight, 
      #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], merged$w.weight, 0) 
      #                               )))
      # } else {
      #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
      #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],1, 
      #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
      #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
      #                               )))
      # }
    
      return(merged)
    } else {
      invisible(capture.output(dataprep.out <- dataprep(foo = x, 
                               dependent = "V4",
                               unit.variable = "V2",
                               # unit.names.variable = "unit.name",
                               time.variable = "V1",
                               treatment.identifier = testid[2], # the regionno of the treated unit
                               controls.identifier = testid[-2],
                               time.optimize.ssr = min(timeid):max(timeid-max.lead-1), # the pre-treatment preiod
                               time.predictors.prior = min(timeid):max(timeid-max.lead-1),
                               predictors = covariate.names)))
      invisible(capture.output(synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS"))) # calibrate the weights
      # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
      # then making it a data.frame so as to merge it in the next step
      weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
      colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
      merged <- merge(x, weights, by = "V2") # merge it with the data.frame
      # if (max.lead > 0) {
      #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
      #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],-1, 
      #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], -merged$w.weight, 
      #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], merged$w.weight, 0) 
      #                               )))
      # } else {
      #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
      #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],1, 
      #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
      #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
      #                               )))
      # }

      return(merged)
    }
  } else {
    merged <- x
    merged$w.weight <- 1
    # if (max.lead > 0){
    #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],-1, 
    #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], -merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], merged$w.weight, 0) 
    #                               )))
    # } else {
    #   merged$wit <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 == testid[2],1, 
    #                               ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid) - max.lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
    #                               )))
    # }
  
    return(merged)
  }
}

Maha_vit <- function(x, lag, max.lead, M = 3) {
  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-max.lead), ]$V2
  x <- na.omit(x)
  x <- x[x$V2 %in% as.numeric(names(which(table(x$V2) == max.lead + 1 + lag))),]
  x <- rbind(x[x$V2 == treated.id, ], x[x$V2 != treated.id, ])
  testid <- unique(x$V2)
  timeid <- unique(x$V1)[1:(lag+1)]
  
  if(treated.id %in% testid == FALSE){
    return(NULL)
  }
  
  if (nrow(x) > 2*(lag + max.lead + 1)) {
    matched_set <- x[x$V1 %in% timeid, ]
    # matched_set <- as.data.frame(apply(matched_set, 2, function(x){
    #   miss_idx    <- is.na(x)
    #   mean_x      <- mean(x, na.rm = TRUE)
    #   x[miss_idx] <- mean_x
    #   return(x)
    # }))
    
    MSMDlist <- suppressWarnings(MSMD_each(timeid[1], matched_set = matched_set, testid = testid, 
                          treated.id = treated.id))
    for (i in 2:length(unique(timeid))) {
      MSMDlist <- MSMDlist + suppressWarnings(MSMD_each(timeid[i], matched_set = matched_set, testid = testid, 
                                       treated.id = treated.id))
    }
    if (length(x) <= 5) {
      MSMDlist <- as.numeric(unlist(as.list(MSMDlist)))
    }
    # MSMDlist <- lapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    # MSMDlist <- sapply(unique(timeid), MSMD_each, matched_set = matched_set, testid = testid)
    MSMD <- append(MSMDlist, 0, after = match(treated.id, testid)-1)
   # MSMD <- append(Reduce("+", lapply(MSMDlist, as.data.frame))/length(timeid), 0, after = 1)
    # MSMD <- append(Reduce("+", lapply(MSMDlist, function(x) as.matrix(x)))/length(timeid), 0, after = 1)
    
    # first.diff <- x[x$V2 == testid[2], ]$V4[L + 1 + max.lead] - x[x$V2 == testid[2], ]$V4[L]
    
    if (M < length(testid)-1) {
      matchid <- order(MSMD)[2:(M+1)]
      weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, treated.id)))
    } else {
      matchid <- order(MSMD)
      weights <- as.data.frame(rbind(cbind(1/(length(testid)-1), testid[matchid[-1]]), cbind(w.weight = 1, treated.id)))
      
    } 
    
    # second.diff <- mean(x[x$V2 %in% testid[matchid] & 
    #                         x$V1 == unique(x$V1)[L + 1 + max.lead], ]$V4 -
    #                       x[x$V2 %in% testid[matchid] & 
    #                           x$V1 == unique(x$V1)[L], ]$V4, na.rm = T)
    # 
    # 
    
    # weights <- as.data.frame(rbind(cbind(1/M, testid[matchid]), cbind(w.weight = 1, testid[2])))
    
    colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
    merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
    # if (max.lead > 0) {
    #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == testid[2],-1, 
    #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], -merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[-2], merged$w.weight, 0) 
    #                               )))
    # } else {
    #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == testid[2],1, 
    #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
    #                               )))
    # }
   
    return(merged)
  } else {
    merged <- x
    merged$w.weight <- 1
    # if (max.lead > 0) {
    #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == testid[2],-1, 
    #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], -merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[-2], merged$w.weight, 0) 
    #                               )))
    # } else {
    #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == testid[2], 1, 
    #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == testid[2],1, 
    #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[-2], merged$w.weight, 
    #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
    #                               )))
    # }
    
    return(merged)
  }
}

MSMD_each <- function(time_id, matched_set, testid, treated.id = treated.id) {
  
  sub_sub <- matched_set[matched_set$V1 == time_id, ]
  # sub_sub[, colSums(is.na(sub_sub)) != 0]
  cov_matrix <- cov(as.matrix(sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)]))
  #cov_matrix <- cov(as.matrix(sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)]),
  #              use = "pairwise.complete.obs") # to avoid NAs causing problems from
  # lagged outcomes
  # cov_matrix <- cov_matrix[rowSums(is.na(cov_matrix))!=dim(cov_matrix)[1],
  #            colSums(is.na(cov_matrix))!=dim(cov_matrix)[2]]
  # 
  if(all(eigen(cov_matrix)$values > .00001)) {
    if (length(sub_sub) == 5) {
      return((sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)] - 
                as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)])) * cov_matrix^(-1) *
               (sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)] - 
                  as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)]))
      )
    } else {
      return(mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)], 
                         center = as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)]),
                         cov = cov_matrix))
    }
    
  } else {
    if (length(sub_sub) == 5) {
      return((sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)] - 
                as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)])) * cov_matrix^(-1) *
               (sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)] - 
                  as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)])))
    } else {
      return(tryCatch(mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)], 
                                  center = as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)]),
                                  cov = diag((length(sub_sub) - 4)) * cov_matrix), 
                      error = function(e) mahalanobis(x = sub_sub[sub_sub$V2 %in% testid[testid != treated.id], 5:length(sub_sub)], 
                                                      center = as.numeric(sub_sub[sub_sub$V2 == treated.id, 5:length(sub_sub)]),
                                                      cov = diag((length(sub_sub) - 4)))))
    }
    
  }
}


PS_vit <- function(x, lag, max.lead, M = M, weighting = FALSE) {
  
#  colnames(x)[1:4] <- c("V1", "V2", "V3", "V4")
  # x <- x[!duplicated(x[c("V2", "V1")]),]
  x <- x[order(x$V2, x$V1), ]
  treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-max.lead), ]$V2
  testid <- unique(x$V2)
  treated.set <- x[x$V1 == (unique(x$V1)[lag+1]), ]
  # timeid_later <- unique(x$V1)
  # timeid <- unique(x$V1)[-((length(unique(x$V1)) - max.lead):length(unique(x$V1)))]
  # matched_set <- x[x$V1 %in% timeid, ]
  
  if (weighting == TRUE) {
    m_ps <- treated.set$ps[which(treated.set$V2!=treated.id)]
    
    weights <- as.data.frame(rbind(cbind((m_ps/(1-m_ps))/sum(m_ps/(1-m_ps)), 
                                         testid[!testid == treated.id]), cbind(w.weight = 1, treated.id)))
  } else {
    PS_distance <- abs(treated.set$ps[which(treated.set$V2 == treated.id)] - 
                         treated.set$ps[which(treated.set$V2 != treated.id)])
    names(PS_distance) <- testid[testid != treated.id]
    if (M < length(testid) - 1) {
      matchid <- as.numeric(names(sort(PS_distance[!names(PS_distance) == treated.id]))[1:M])
      weights <- as.data.frame(rbind(cbind(1/M, matchid), cbind(w.weight = 1, treated.id)))
    } else {
      matchid <- as.numeric(names(sort(PS_distance[!names(PS_distance) == treated.id])))
      weights <- as.data.frame(rbind(cbind(1/(length(testid)-1), matchid), cbind(w.weight = 1, treated.id)))
    }
    
  }
  
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  colnames(weights)[1] <- "w.weight"
  merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
  merged <- merged[order(merged$V2, merged$V1), ]
 # colnames(merged)[c(4,5)] <- c("V3", "V4")
  # if (max.lead > 0) {
  #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == treated.id, 1, 
  #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == treated.id, -1, 
  #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[testid != treated.id], -merged$w.weight, 
  #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[testid != treated.id], merged$w.weight, 0) 
  #                               )))
  # } else {
  #   merged$wit <- ifelse(merged$V1 == max(timeid_later) & merged$V2 == treated.id, 1, 
  #                        ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 == treated.id, 1, 
  #                               ifelse(merged$V1 == max(timeid_later) & merged$V2 %in% testid[testid != treated.id], merged$w.weight, 
  #                                      ifelse(merged$V1 == max(timeid_later) - max.lead - 1 & merged$V2 %in% testid[testid != treated.id], -merged$w.weight, 0) 
  #                               )))
  # }
  # 
  return(merged) # return the weight variable
} 


# PanelWit <- function(data, unit.id, time.id, matched_set,
#                      lag, lead) {
#   # set testid and timeid
#   treated.time <- max(matched_set$V1)-lead
#   treated.id <- matched_set[matched_set$V3 == 1 & 
#                               matched_set$V1 == treated.time, ]$V2
#   
#   new.W <- data[c(unit.id, time.id)] # create a new data.frame as large as the *cleaned* dataset
#   names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
#   total2 <- merge(new.W, matched_set, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
#   total2$wit <- ifelse(is.na(total2$wit), 0, total2$wit) # turn NAs into zero
#   total2$dit <- 0
#   # IMPORTANT: CHANGE HERE IF THINGS GO WRONG FOR lead >0
#   total2$dit[which(total2$V2 == treated.id & total2$V1 == treated.time)] <- 1
#   ######################################################
#   return (list("wit" = total2$wit, "dit" = total2$dit)) 
# }

take_out <- function(matched_set, lag, lead) {
  matched_set <- matched_set[matched_set$V1 <= (min(matched_set$V1) + lag + max(lead)),]
}

# PanelWit2 <- function(data, unit.id, time.id, matched_set,
#                      lag, lead) {
#   # set testid and timeid
#   treated.time <- min(matched_set$V1) + lag
#   treated.id <- matched_set[matched_set$V3 == 1 & 
#                               matched_set$V1 == treated.time, ]$V2
#   
#   testid <- unique(matched_set$V2)
#   
#   if (lead > 0) {
#     matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
#                          ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, -1, 
#                                 ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 
#                                        ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 0) 
#                                 )))
#   } else {
#     matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
#                          ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, 1, 
#                                 ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 
#                                        ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 0) 
#                                 )))
#   }
#   
#   new.W <- data[c(unit.id, time.id)] # create a new data.frame as large as the *cleaned* dataset
#   names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
#   total2 <- merge(new.W, matched_set, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
#   total2$wit <- ifelse(is.na(total2$wit), 0, total2$wit) # turn NAs into zero
#   total2$dit <- 0
#   # IMPORTANT: CHANGE HERE IF THINGS GO WRONG FOR lead >0
#   total2$dit[which(total2$V2 == treated.id & total2$V1 == treated.time)] <- 1
#   ######################################################
#   return (list("wit" = total2$wit, "dit" = total2$dit)) 
# }
# 




# PanelDiDResult <- function(x, lag, lead){
#   testid <- unique(x$V2)
#   location <- match(x$V2[x$V3 == 1 & x$V1 == (max(x$V1)-lead)], testid)
#   
#   first.diff <- x[x$V2 == testid[location], ]$V4[lag+1+lead] - x[x$V2 == testid[location], ]$V4[lag]
#   second.diff.1 <- sum(x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag+1+lead], ]$V4*
#                          x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag+1+lead], ]$w.weight)
#   second.diff.2 <- sum(x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag], ]$V4*
#                          x[x$V2 %in% testid[-location] & x$V1 == unique(x$V1)[lag], ]$w.weight)
#   second.diff <- second.diff.1 - second.diff.2  
#   return(first.diff - second.diff)
#   
# }


# newlist$`Matched sets for ATC`[[1]]
# PanelDiDResult(x = merged, 
#                lag = 4, lead = 2)


gaps_plot_tmp <- function(x, lag, lead, data, dependent,
                          qoi, refinement = TRUE, 
                          covariate_names,
                          method,
                          treated_set,
                          covariate = NULL) {
  treated_set <- treated_set
  # colnames(data) <- c("time.id", "unit.id",
  #                     "treatment", "dependent", covariate_names)
  if (method == "Maha"|method == "Synth") {
    x <- as.data.frame(x)
    colnames(x)[5:length(x)] <- c(covariate_names,
                                  "dependent",
                                  "w.weight")
    colnames(treated_set)[5:(length(treated_set)-3)] <- 
      covariate_names
  }
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


# parallel_trends <- function(x, lag, lead, adjustment) {
#   
#   treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
#   treated.outcome <- x$V4[which(x$V2 == treated.id)]
#   if (adjustment == FALSE) {
#     control.outcome <- tapply(x$V4[which(x$V2 != treated.id)],
#                               x$V1[which(x$V2 != treated.id)], mean)
#   } else {
#     control.outcome <- tapply(x$V4[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
#                               x$V1[x$V2 != treated.id], sum)
#   }
#   
#   
#   return(list("treated.outcome" = treated.outcome,
#               "control.outcome" = control.outcome))
# }
# 




# gaps_plot <- function(x, lag, lead, covariate = NULL) {
#   treated.id <- x[x$V3 == 1 & x$V1 == (max(x$V1)-lead), ]$V2 # check this
#   if (is.null(covariate)) {
#     return(list("gap" = x$V4[x$V2 == treated.id] - 
#                   tapply(x$V4[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
#                          x$V1[x$V2 != treated.id], sum),
#                 "unit" = paste(treated.id, unique(x$V1)[lag + 1], sep = ",")))
#   } else {
#     return(list("gap" = x$V5[x$V2 == treated.id] - 
#                   tapply(x$V5[x$V2 != treated.id] * x$w.weight[x$V2 != treated.id], 
#                          x$V1[x$V2 != treated.id], sum),
#                 "unit" = paste(treated.id, unique(x$V1)[lag + 1], sep = ",")))
#   }
#   
# }


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


equality_four <- function(x, y, z){
  return(sum(x*y)/sum(z))
}



delete.NULLs <- function(x.list) {
  x.list[unlist(lapply(x.list, length) != 0)]
}
dframelist.rb_dup <- function(x) {
  test <- data.table::rbindlist(lapply(x, as.data.frame))
  return(test[!duplicated(test), ])
}

grab_index <- function(x) {
  return(list("unit.index" = unique(x$V2), 
              "time.index" = unique(x$V1)))
}
# indexes <- lapply(matched_sets$`Matched sets for ATT`, grab_index)

gen_sets <- function(x, data, unit.id, time.id) {
  return(data[data[unit.id][[1]] %in% x$unit.index &
                data[time.id][[1]] %in% x$time.index, ])
}

# new.matched_sets <- lapply(indexes, gen_sets, data = Data.obs, unit.id = "unit",
 #                          time.id = "time")


# 
# merge.formula <- function(form1, form2, ...){
#   
#   # get character strings of the names for the responses 
#   # (i.e. left hand sides, lhs)
#   lhs1 <- deparse(form1[[2]])
#   lhs2 <- deparse(form2[[2]])
#   if(lhs1 != lhs2) stop('both formulas must have the same response')
#   
#   # get character strings of the right hand sides
#   rhs1 <- strsplit(deparse(form1[[3]]), " \\+ ")[[1]]
#   rhs2 <- strsplit(deparse(form2[[3]]), " \\+ ")[[1]]
#   
#   # create the merged rhs and lhs in character string form
#   rhs <- c(rhs1, rhs2)
#   lhs <- lhs1
#   
#   # put the two sides together with the amazing 
#   # reformulate function
#   out <- reformulate(rhs, lhs)
#   
#   # set the environment of the formula (i.e. where should
#   # R look for variables when data aren't specified?)
#   environment(out) <- parent.frame()
#   
#   return(out)
# }

# this is how you get the addition operator working for formulas
Ops.formula <- function(e1, e2){
  FUN <- .Generic
  if(FUN == '+'){
    out <- merge(e1, e2)
    environment(out) <- parent.frame()
    return(out)
  }
  else stop('can not yet subtract formula objects')
}



### updating (adding) wits
apply_leads <- function(data, unit.id, time.id, matched_set,
                        lag, lead, estimator = c("did", "matching")) {
  unit.id = unit.id; time.id = time.id
  # set testid and timeid
  treated.time <- min(matched_set$V1) + lag
  treated.id <- matched_set[matched_set$V3 == 1 & 
                              matched_set$V1 == treated.time, ]$V2
  
  testid <- unique(matched_set$V2)
  
  if (estimator == "did") {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
                              ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 == treated.id, -1, 
                                     ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 %in% testid[testid != treated.id], -matched_set$w.weight, 
                                            ifelse(matched_set$V1 == (treated.time+lead) - lead - 1 & matched_set$V2 %in% testid[testid != treated.id], matched_set$w.weight, 0) 
                                     )))
    
  } else {
    matched_set$wit <- ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 == treated.id, 1, 
           ifelse(matched_set$V1 == (treated.time+lead) & matched_set$V2 != treated.id, -matched_set$w.weight, 0))
                 
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
                         estimator = estimator){
  return(lapply(leads, apply_leads, unit.id = unit.id, time.id = time.id,
         matched_set = x, lag = lag, estimator = estimator,
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
