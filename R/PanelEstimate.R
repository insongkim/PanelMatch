PanelEstimate <- function(lag, lead, scheme, 
                          data, treatment, qoi, M,
                          time.id, unit.id,
                          covariate, dependent, matched_set = NULL) {
  if(is.null(matched_set)) {
    matched_set <- PanelMatch(lag = lag, lead = lead, data = data, 
                          treatment = treatment,
                          dependent = dependent, 
                          covariate = covariate, qoi = qoi, 
                          scheme = scheme,
                          M = M,
                          time.id = time.id, unit.id = unit.id)
  } 
  # ATT
  if (length(matched_set) == 2 & 
      "Matched sets for ATT" %in% names(matched_set)) {
    all.diffs.weighted <- sapply(matched_set$`Matched sets for ATT`, 
                                 PanelDiDResult, lag = lag, lead = lead)
    
    return(list("ATT" = mean(all.diffs.weighted, na.rm = T)))
  # ATC
  } else if (length(matched_set) == 2 & 
             "Matched sets for ATC" %in% names(matched_set)){
    all.diffs.weighted <- sapply(matched_set$`Matched sets for ATC`, 
                                 PanelDiDResult, lag = lag, lead = lead)
    return(list("ATC" = mean(all.diffs.weighted, na.rm = T)))
  } else if (length(matched_set) == 4) {
    all.diffs.weighted1 <- sapply(matched_set$`Matched sets for ATT`, 
                                 PanelDiDResult, lag = lag, lead = lead)
    ATT <- mean(all.diffs.weighted1, na.rm = T)
    all.diffs.weighted2 <- sapply(matched_set$`Matched sets for ATC`, 
                                  PanelDiDResult, lag = lag, lead = lead)
    ATC <- mean(all.diffs.weighted2, na.rm = T)
    ATC <- -(ATC)
    DID_ATE <- ifelse(length(ATC) > 0, (ATT * length(all.diffs.weighted1) + 
                                          ATC * length(all.diffs.weighted2))/(length(all.diffs.weighted1) + 
                                                                            length(all.diffs.weighted2)), ATT)
    return(list("ATE" = DID_ATE))
  }
}