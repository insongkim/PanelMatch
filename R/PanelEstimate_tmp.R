PanelEstimate_tmp <- function(lead, 
                          covariate.only = FALSE,
                          inference = c("wfe", "bootstrap"), 
                          ITER = 1000, matched_sets = NULL,
                          plot = FALSE) {
  
  # stop if lead > max.lead
  if (lead > matched_sets$max.lead) 
    stop(paste("The number of leads you choose 
               has exceeded the maximum number you set
               when finding matched sets, which is", matched_sets$max.lead))
  
  lag = matched_sets$lag
  data <- matched_sets$data
  dependent = matched_sets$dependent
  treatment = matched_sets$treatment
  covariate.only = matched_sets$covariate.only
  unit.id = matched_sets$unit.id
  
  if (is.null(matched_sets$`Matched sets for ATT`) == FALSE) {
    matched_sets$`Matched sets for ATT` <- lapply(matched_sets$`Matched sets for ATT`, 
                                                  take_out, lag = lag, lead = lead)
  }
  
  if (is.null(matched_sets$`Matched sets for ATC`) == FALSE) {
    matched_sets$`Matched sets for ATC` <- lapply(matched_sets$`Matched sets for ATC`, 
                                                  take_out, lag = lag, lead = lead)
  }
  

  
  if (matched_sets$qoi == "att") {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATT`,
                               PanelWit2, 
                               unit.id = matched_sets$unit.id, 
                               time.id = matched_sets$time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    # cat("doing it")
    data$weights_att <- Reduce("+", all.weights)
    data$dit_att <- Reduce("+", all.dits)
  } else if (matched_sets$qoi == "atc") {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATC`,
                               PanelWit2, 
                               unit.id = matched_sets$unit.id, 
                               time.id = matched_sets$time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_atc <- Reduce("+", all.weights)
    data$dit_atc <- Reduce("+", all.dits)
  } else if (matched_sets$qoi == "ate") {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATT`,
                               PanelWit2, 
                               unit.id = matched_sets$unit.id, 
                               time.id = matched_sets$time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_att <- Reduce("+", all.weights)
    data$dit_att <- Reduce("+", all.dits)
    
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATC`,
                               PanelWit2, 
                               unit.id = matched_sets$unit.id, 
                               time.id = matched_sets$time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_atc <- Reduce("+", all.weights)
    data$dit_atc <- Reduce("+", all.dits)
  }
  
  
  # ATT
  if (matched_sets$qoi == "att") {
    if (inference == "wfe"){
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = matched_sets$unit.id,
                      time.index = matched_sets$time.id, method = "unit", 
                      qoi = "att", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
      
    } else if (inference == "bootstrap"){
      all.diffs.weighted <- sapply(matched_sets$`Matched sets for ATT`, 
                                   PanelDiDResult, lag = lag, lead = lead)
      
      coefs <- rep(NA, ITER) 
      dit.atts <- rep(NA, ITER)
      wit.atts <- list()
      
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        if (lead > 0) {
          att.new <- sum(d.sub1$weights_att*d.sub1$dv)/sum(d.sub1$dit_att)
        } else {
          att.new <- sum(d.sub1$weights_att*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_att)
        }
        
        coefs[k] <- att.new
        dit.atts[k] <- sum(d.sub1$dit_att)
        wit.atts[[k]] <- d.sub1$weights_att
      }
      
      return(list("o.coef" = mean(all.diffs.weighted, na.rm = T),
                  "boots" = coefs))
    }
    
    # ATC
  } else if (matched_sets$qoi == "atc"){
    if (inference == "wfe") {
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = matched_sets$unit.id,
                      time.index = matched_sets$time.id, method = "unit", 
                      qoi = "atc", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
    } else if (inference == "bootstrap") {
      all.diffs.weighted <- sapply(matched_sets$`Matched sets for ATC`, 
                                   PanelDiDResult, lag = lag, lead = lead)
      
      coefs <- rep(NA, ITER) 
      dit.atcs <- rep(NA, ITER)
      wit.atcs <- list()
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        if (lead > 0) {
          atc.new <- -sum(d.sub1$weights_atc*d.sub1$dv)/sum(d.sub1$dit_atc)
        } else {
          atc.new <- sum(d.sub1$weights_atc*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_atc)
        }
        coefs[k] <- atc.new
        dit.atcs[k] <- sum(d.sub1$dit_atc)
        wit.atcs[[k]] <- d.sub1$weights_atc
      }
      
      return(list("o.coef" = -mean(all.diffs.weighted, na.rm = T),
                  "boots" = coefs))
    }
    
  } else if (matched_sets$qoi == "ate") {
    if (inference == "wfe"){
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = matched_sets$unit.id,
                      time.index = matched_sets$time.id, method = "unit", 
                      qoi = "ate", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
    } else if (inference == "bootstrap"){
      all.diffs.weighted1 <- sapply(matched_sets$`Matched sets for ATT`, 
                                    PanelDiDResult, lag = lag, lead = lead)
      ATT <- mean(all.diffs.weighted1, na.rm = T)
      all.diffs.weighted2 <- sapply(matched_sets$`Matched sets for ATC`, 
                                    PanelDiDResult, lag = lag, lead = lead)
      ATC <- -mean(all.diffs.weighted2, na.rm = T)
      DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(all.diffs.weighted1) + 
                                            ATC * length(all.diffs.weighted2))/(length(all.diffs.weighted1) + 
                                                                                  length(all.diffs.weighted2)), ATT)
      coefs <- rep(NA, ITER) 
      dit.atts <- rep(NA, ITER)
      dit.atcs <- rep(NA, ITER)
      wit.atts <- list()
      wit.atcs <- list()
      
      for (k in 1:ITER) {
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        
        if (lead > 0) {
          att.new <- sum(d.sub1$weights_att*d.sub1$dv)/sum(d.sub1$dit_att)
        } else {
          att.new <- sum(d.sub1$weights_att*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_att)
        }
        
        d.sub1$control <- ifelse(d.sub1$treatment == 1, 0, 1)
        if (lead > 0) {
          atc.new <- -sum(d.sub1$weights_atc*d.sub1$dv)/sum(d.sub1$dit_atc)
        } else {
          atc.new <- sum(d.sub1$weights_atc*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_atc)
        }
        
        coefs[k] <- (att.new*sum(d.sub1$dit_att) + atc.new*sum(d.sub1$dit_atc))/(sum(d.sub1$dit_att) + sum(d.sub1$dit_atc))
        dit.atts[k] <- sum(d.sub1$dit_att)
        dit.atcs[k] <- sum(d.sub1$dit_atc)
        wit.atts[[k]] <- d.sub1$weights_att
        wit.atcs[[k]] <- d.sub1$weights_atc
      }
      return(list("o.coef" = DID_ATE, "boots" = coefs))
      
    }
    
    
  }
}