PanelEstimate <- function(lag, lead, scheme, 
                          data, treatment, qoi, M,
                          time.id, unit.id, covariate.only = FALSE,
                          inference = c("wfe", "bootstrap"),
                          covariate, dependent, matched_sets = NULL) {
  if(is.null(matched_sets)) {
    matched_sets <- PanelMatch(lag = lag, lead = lead, data = data, 
                          treatment = treatment,
                          dependent = dependent, 
                          covariate = covariate, qoi = qoi, 
                          scheme = scheme,
                          M = M, covariate.only = covariate.only,
                          time.id = time.id, unit.id = unit.id)
  } 
  
  data <- matched_sets$data
  if (length(matched_sets) == 3 &
      "Matched sets for ATT" %in% names(matched_sets)) {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATT`,
                               PanelWit, 
                               unit.id = unit.id, 
                               time.id = time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    # cat("doing it")
    data$weights_att <- Reduce("+", all.weights)
    data$dit_att <- Reduce("+", all.dits)
  } else if (length(matched_sets) == 3 &
             "Matched sets for ATC" %in% names(matched_sets)) {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATC`,
                               PanelWit, 
                               unit.id = unit.id, 
                               time.id = time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_atc <- Reduce("+", all.weights)
    data$dit_atc <- Reduce("+", all.dits)
  } else if (length(matched_sets) == 5) {
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATT`,
                               PanelWit, 
                               unit.id = unit.id, 
                               time.id = time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_att <- Reduce("+", all.weights)
    data$dit_att <- Reduce("+", all.dits)
    
    weights_and_dits <- lapply(matched_sets$`Matched sets for ATC`,
                               PanelWit, 
                               unit.id = unit.id, 
                               time.id = time.id,
                               data = matched_sets$data, 
                               lag = lag, lead = lead)
    
    all.weights <- (lapply(weights_and_dits, function (x) x$wit))
    all.dits <- (lapply(weights_and_dits, function (x) x$dit))
    data$weights_atc <- Reduce("+", all.weights)
    data$dit_atc <- Reduce("+", all.dits)
  }
  
  
  # ATT
  if (length(matched_sets) == 3 & 
      "Matched sets for ATT" %in% names(matched_sets)) {
    if (inference == "wfe"){
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = unit.id,
                      time.index = time.id, method = "unit", 
                      qoi = "att", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      return(fit)
    } else if (inference == "bootstrap"){
      all.diffs.weighted <- sapply(matched_sets$`Matched sets for ATT`, 
                                   PanelDiDResult, lag = lag, lead = lead)
      
      return(list("ATT" = mean(all.diffs.weighted, na.rm = T)))
    }
   
  # ATC
  } else if (length(matched_sets) == 3 & 
             "Matched sets for ATC" %in% names(matched_sets)){
    if (inference == "wfe") {
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = unit.id,
                      time.index = time.id, method = "unit", 
                      qoi = "atc", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      return(fit)
    } else if (inference == "bootstrap") {
      all.diffs.weighted <- sapply(matched_sets$`Matched sets for ATC`, 
                                   PanelDiDResult, lag = lag, lead = lead)
      
      return(list("ATC" = mean(all.diffs.weighted, na.rm = T)))
    }
   
  } else if (length(matched_sets) == 5) {
    if (inference == "wfe"){
      fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                      treat = treatment, unit.index = unit.id,
                      time.index = time.id, method = "unit", 
                      qoi = "ate", estimator = "did", 
                      hetero.se = TRUE, 
                      auto.se = TRUE, White = TRUE,  
                      data = data)
      return(fit)
    } else if (inference == "bootstrap"){
      all.diffs.weighted1 <- sapply(matched_sets$`Matched sets for ATT`, 
                                    PanelDiDResult, lag = lag, lead = lead)
      ATT <- mean(all.diffs.weighted1, na.rm = T)
      all.diffs.weighted2 <- sapply(matched_sets$`Matched sets for ATC`, 
                                    PanelDiDResult, lag = lag, lead = lead)
      ATC <- mean(all.diffs.weighted2, na.rm = T)
      ATC <- -(ATC)
      DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(all.diffs.weighted1) + 
                                            ATC * length(all.diffs.weighted2))/(length(all.diffs.weighted1) + 
                                                                                  length(all.diffs.weighted2)), ATT)
      return(list("ATE" = DID_ATE))
    }
  
   
  }
}