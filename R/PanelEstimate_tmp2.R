#' Illustration of PanelEstimate_tmp2
#' 


PanelEstimate_tmp2 <- function(lead, 
                               covariate.only = FALSE,
                               inference = c("wfe", "bootstrap"), 
                               ITER = 1000, matched_sets = NULL,
                               plot = FALSE,
                               estimator = "did",
                               df.adjustment = FALSE,
                               qoi = NULL,
                               CI = .90) {
  
  # stop if lead > max.lead
  if (max(lead) > matched_sets$max.lead) 
    stop(paste("The number of leads you choose 
               has exceeded the maximum number you set
               when finding matched sets, which is", matched_sets$max.lead))
  
  lag = matched_sets$lag
  data <- matched_sets$data
  dependent = matched_sets$dependent
  treatment = matched_sets$treatment
  covariate.only = matched_sets$covariate.only
  unit.id = matched_sets$unit.id
  time.id = matched_sets$time.id
  method = matched_sets$method
  if (is.null(qoi)) {
    qoi = matched_sets$qoi
  } else {
    qoi = qoi
  }
  
  
  
  if (is.null(matched_sets$`ATT_matches`) == FALSE) {
    matched_sets$`ATT_matches` <- lapply(matched_sets$`ATT_matches`, 
                                         take_out, lag = lag, lead = max(lead))
    if (length(matched_sets$ATT_matches) == 0) {
      qoi <- "atc"
    }
  }
  
  if (is.null(matched_sets$`ATC_matches`) == FALSE) {
    matched_sets$`ATC_matches` <- lapply(matched_sets$`ATC_matches`, 
                                         take_out, lag = lag, lead = max(lead))
    if (length(matched_sets$ATC_matches) == 0) {
      qoi <- "att"
    }
  }
  
  
  if (qoi == "att") {
    newlist <- lapply(matched_sets$`ATT_matches`, lapply_leads, unit.id = unit.id, 
                      time.id = time.id, lag = lag, estimator = estimator,
                      data = data, leads = lead)
    
    W_it_by_lead <- lapply(newlist, extract_objects, objective = "wit")
    dits <- lapply(newlist, extract_objects, objective = "dit")
    dits <- each_lead(dits, lead = 1)
    lead2 <- 1:length(lead)
    W_it_by_lead <- lapply(lead2, each_lead, x = W_it_by_lead)
    W_it_by_lead <- lapply(W_it_by_lead, function(x) Reduce("+", x))
    
    data[, (length(data) + 1):(length(data) + length(W_it_by_lead))] <- unlist(W_it_by_lead)
    colnames(data)[match(tail(colnames(data), n = length(W_it_by_lead)), colnames(data))] <- sapply(lead, function(x) paste0("Wit_att", x))
    data[, length(data) +1] <- Reduce("+", dits)
    colnames(data)[length(data)] <- "dits_att"
    data$`Wit_att-1` <- 0
  } else if (qoi == "atc") {
    newlist <- lapply(matched_sets$`ATC_matches`, lapply_leads, unit.id = unit.id, 
                      time.id = time.id, lag = lag, estimator = estimator,
                      data = data, leads = lead)
    
    W_it_by_lead <- lapply(newlist, extract_objects, objective = "wit")
    dits <- lapply(newlist, extract_objects, objective = "dit")
    dits <- each_lead(dits, lead = 1)
    lead2 <- 1:length(lead)
    W_it_by_lead <- lapply(lead2, each_lead, x = W_it_by_lead)
    W_it_by_lead <- lapply(W_it_by_lead, function(x) Reduce("+", x))
    
    data[, (length(data) + 1):(length(data) + length(W_it_by_lead))] <- unlist(W_it_by_lead)
    colnames(data)[match(tail(colnames(data), n = length(W_it_by_lead)), colnames(data))] <- sapply(lead, function(x) paste0("Wit_atc", x))
    data[, length(data) +1] <- Reduce("+", dits)
    colnames(data)[length(data)] <- "dits_atc"
    data$`Wit_atc-1` <- 0
  } else if (qoi == "ate") {
    # ATT
    newlist <- lapply(matched_sets$`ATT_matches`, lapply_leads, unit.id = unit.id, 
                      time.id = time.id, lag = lag, estimator = estimator,
                      data = data, leads = lead)
    
    W_it_by_lead <- lapply(newlist, extract_objects, objective = "wit")
    dits <- lapply(newlist, extract_objects, objective = "dit")
    dits <- each_lead(dits, lead = 1)
    lead2 <- 1:length(lead)
    W_it_by_lead <- lapply(lead2, each_lead, x = W_it_by_lead)
    W_it_by_lead <- lapply(W_it_by_lead, function(x) Reduce("+", x))
    
    data[, (length(data) + 1):(length(data) + length(W_it_by_lead))] <- unlist(W_it_by_lead)
    colnames(data)[match(tail(colnames(data), n = length(W_it_by_lead)), colnames(data))] <- sapply(lead, function(x) paste0("Wit_att", x))
    data[, length(data) +1] <- Reduce("+", dits)
    colnames(data)[length(data)] <- "dits_att"
    # ATC
    newlist <- lapply(matched_sets$`ATC_matches`, lapply_leads, unit.id = unit.id, 
                      time.id = time.id, lag = lag, estimator = estimator, 
                      data = data, leads = lead)
    
    W_it_by_lead <- lapply(newlist, extract_objects, objective = "wit")
    dits <- lapply(newlist, extract_objects, objective = "dit")
    dits <- each_lead(dits, lead = 1)
    lead2 <- 1:length(lead)
    W_it_by_lead <- lapply(lead2, each_lead, x = W_it_by_lead)
    W_it_by_lead <- lapply(W_it_by_lead, function(x) Reduce("+", x))
    
    data[, (length(data) + 1):(length(data) + length(W_it_by_lead))] <- unlist(W_it_by_lead)
    colnames(data)[match(tail(colnames(data), n = length(W_it_by_lead)), colnames(data))] <- sapply(lead, function(x) paste0("Wit_atc", x))
    data[, length(data) +1] <- Reduce("+", dits)
    colnames(data)[length(data)] <- "dits_atc"
    data$`Wit_att-1` <- 0
    data$`Wit_atc-1` <- 0
  }
  
  
  # ATT
  if (qoi == "att") {
    if (inference == "wfe"){
      if(length(lead) != 1 || lead != 0)
        stop("The wfe option can only take lead = 0")
      data$Wit_att0 <- ifelse(data$dits_att == 1, -1, data$Wit_att0)
      data$Wit_att0 <- -(data$Wit_att0)
      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "unit", 
                        qoi = "att", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "time", 
                        qoi = "att", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      }
      
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
      
    } else if (inference == "bootstrap"){
      o.coefs <- sapply(data[, sapply(lead, function(x) paste0("Wit_att", x)), drop = FALSE],
                        equality_four,
                        y = data[c(dependent)][,1],
                        z = data$dits_att)
      
      if (length(lead[lead<0]) > 1) {
        names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                         length(o.coefs)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                         length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
      
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(W_it_by_lead))
      # dit.atts <- rep(NA, ITER)
      # wit.atts <- list()
      
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        att_new <-  sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1$dv,
                           z = d.sub1$dits_att)
        # if (lead > 0) {
        #   att.new <- sum(d.sub1$weights_att*d.sub1$dv)/sum(d.sub1$dit_att)
        # } else {
        #   att.new <- sum(d.sub1$weights_att*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_att)
        # }
        # 
        coefs[k,] <- att_new
        # dit.atts[k] <- sum(d.sub1$dit_att)
        # wit.atts[[k]] <- d.sub1$weights_att
      }
      # changed return to class
      # return(list("o.coef" = mean(all.diffs.weighted, na.rm = T),
      #             "boots" = coefs))
      z <- list("o.coef" = o.coefs,
                "boots" = coefs, "ITER" = ITER,
                "method" = method, "lag" = lag,
                "lead" = lead, "CI" = CI, "qoi" = qoi)
      class(z) <- "PanelEstimate_tmp2"
      z
    }
    
    # ATC
  } else if (qoi == "atc"){
    if (inference == "wfe") {
      if(length(lead) != 1 || lead != 0)
        stop("The wfe option can only take lead = 0")
      data$Wit_atc0 <- ifelse(data$dits_atc == 1, -1, data$Wit_atc0)
      data$Wit_atc0 <- -(data$Wit_atc0)
      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "unit", 
                        qoi = "atc", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "time", 
                        qoi = "atc", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      }
      
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
    } else if (inference == "bootstrap") {
      o.coefs <-  -sapply(data[, sapply(lead, function(x) paste0("Wit_atc", x)), drop = FALSE],
                          equality_four,
                          y = data[c(dependent)][,1],
                          z = data$dits_atc)
      
      if (length(lead[lead<0]) > 1) {
        names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                         length(o.coefs)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                         length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(W_it_by_lead))
      # dit.atts <- rep(NA, ITER)
      # wit.atts <- list()
      
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        atc_new <- -sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_atc", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1$dv,
                           z = d.sub1$dits_atc)
        # if (lead > 0) {
        #   att.new <- sum(d.sub1$weights_att*d.sub1$dv)/sum(d.sub1$dit_att)
        # } else {
        #   att.new <- sum(d.sub1$weights_att*(2*d.sub1$treatment-1)*d.sub1$dv)/sum(d.sub1$dit_att)
        # }
        # 
        coefs[k,] <- atc_new
        # dit.atts[k] <- sum(d.sub1$dit_att)
        # wit.atts[[k]] <- d.sub1$weights_att
      }
      # changed return to class
      # return(list("o.coef" = mean(all.diffs.weighted, na.rm = T),
      #             "boots" = coefs))
      z <- list("o.coef" = o.coefs,
                "boots" = coefs, "ITER" = ITER,
                "method" = method, "lag" = lag,
                "lead" = lead, "CI" = CI, "qoi" = qoi)
      class(z) <- "PanelEstimate_tmp2"
      z
      # return(list("o.coef" = -mean(all.diffs.weighted, na.rm = T),
      #             "boots" = coefs))
    }
    
  } else if (qoi == "ate") {
    if (inference == "wfe"){
      if(length(lead) != 1 || lead != 0)
        stop("The wfe option can only take lead = 0")
      data$Wit_att0 <- ifelse(data$dits_att == 1, -1, data$Wit_att0)
      data$Wit_att0 <- -(data$Wit_att0)
      
      data$Wit_atc0 <- ifelse(data$dits_atc == 1, -1, data$Wit_atc0)
      data$Wit_atc0 <- -(data$Wit_atc0)

      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "unit", 
                        qoi = "ate", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = matched_sets$unit.id,
                        time.index = matched_sets$time.id, method = "time", 
                        qoi = "ate", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
      }
      
      if (plot == TRUE) {
        fit$matched_sets <- matched_sets
        return(fit)
      } else {
        return(fit)
      }
    } else if (inference == "bootstrap"){
      o.coefs_att <-  sapply(data[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                  drop = FALSE],
                             equality_four,
                             y = data[c(dependent)][,1],
                             z = data$dits_att)
      
      o.coefs_atc <-  -sapply(data[, sapply(lead, function(x) paste0("Wit_atc", x)), 
                                   drop = FALSE],
                              equality_four,
                              y = data[c(dependent)][,1],
                              z = data$dits_atc)
      
      o.coefs_ate <- (o.coefs_att*sum(data$dits_att) + o.coefs_atc*sum(data$dits_atc))/
        (sum(data$dits_att) + sum(data$dits_atc))
      
      if (length(lead[lead<0]) > 1) {
        names(o.coefs_ate)[(length(o.coefs_ate)-max(lead[lead>=0])):
                         length(o.coefs_ate)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs_ate)[(length(o.coefs_ate)-length(lead) + 1):
                         length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else {
        names(o.coefs_ate) <- sapply(lead, function(x) paste0("t+", x))
      }
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(W_it_by_lead))
      
      # dit.atts <- rep(NA, ITER)
      # dit.atcs <- rep(NA, ITER)
      # wit.atts <- list()
      # wit.atcs <- list()
      
      for (k in 1:ITER) {
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        # create bootstap sample with sapply
        df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        d.sub1 <- data[unlist(df.bs),]
        colnames(d.sub1)[3:4] <- c("treatment", "dv")
        
        att_new <-sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                drop = FALSE],
                         equality_four,
                         y = d.sub1$dv,
                         z = d.sub1$dits_att)
        
        atc_new <- -sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_atc", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1$dv,
                           z = d.sub1$dits_atc)
        
        coefs[k,] <- (att_new*sum(d.sub1$dits_att) + atc_new*sum(d.sub1$dits_atc))/
          (sum(d.sub1$dits_att) + sum(d.sub1$dits_atc))
        
        
        # dit.atts[k] <- sum(d.sub1$dit_att)
        # dit.atcs[k] <- sum(d.sub1$dit_atc)
        # wit.atts[[k]] <- d.sub1$weights_att
        # wit.atcs[[k]] <- d.sub1$weights_atc
      }
      # return(list("o.coef" = DID_ATE, "boots" = coefs))
      z <- list("o.coef" = o.coefs_ate,
                "boots" = coefs, "ITER" = ITER,
                "method" = method, "lag" = lag,
                "lead" = lead, "CI" = CI, "qoi" = qoi)
      class(z) <- "PanelEstimate_tmp2"
      z
      
    }
    
    
  }
}


summary.PanelEstimate_tmp2 <- function(object) {
  if(object$method == "Maha"){
    cat("Weighted Difference-in-Differences with Mahalanobis Distance\n")
  } else if (object$method == "Pscore") {
    cat("Weighted Difference-in-Differences with Propensity Score\n")
  } else if (object$method == "Synth") {
    cat("Weighted Difference-in-Differences with Synthetic Control\n")
  } else if (object$method == "CBPS") {
    cat("Weighted Difference-in-Differences with Covariate Balancing Propensity Score\n")
  }
  cat("Matches created with", object$lag, "lags\n")
  cat("\nStandard errors computed with", object$ITER, "Weighted bootstrap samples\n")
  # cat("\nTotal effect in", object$lead, "periods after the treatment:")
  qoi <- ifelse(object$qoi == "ate", "Average Treatment Effect (ATE)", 
                ifelse(object$qoi == "att", "Average Treatment Effect on the Treated (ATT)", 
                       "Average Treatment Effect on the Control (ATC)"))
  cat("\nEstimate of", qoi, "by Period:")
  df <- rbind(t(as.data.frame(object$o.coef)), # point estimate
              colSds(object$boots, na.rm = T), # bootstrap se
              
              # Efron & Tibshirani 1993 p170 - 171
              t(colQuantiles(object$boots,
                           probs = c((1-object$CI)/2, object$CI+(1-object$CI)/2), 
                           na.rm = T, drop = FALSE)), # percentile CI
              # Efron & Tibshirani 1993 p138
              2*object$o.coef - colMeans(object$boots, na.rm = T), # bc point estimate
              
              t(colQuantiles(2*matrix(nrow = object$ITER, ncol = length(object$o.coef), 
                                      object$o.coef, byrow = TRUE) - object$boots,
                             probs = c((1-object$CI)/2, object$CI+(1-object$CI)/2), 
                             na.rm = T, drop = FALSE)) # bc percentile CI
              )
  rownames(df) <- c("Point Estimate(s)", "Standard Error(s)", 
                    paste("Lower Limit of", object$CI*100, "% Regular Confidence Interval"),
                    paste("Upper Limit of", object$CI*100, "% Regular Confidence Interval"),
                    "Bias-corrected Estimate(s)", 
                    paste("Lower Limit of", object$CI*100, "% Bias-corrected Confidence Interval"),
                    paste("Upper Limit of", object$CI*100, "% Bias-corrected Confidence Interval"))
  print(kable(df))
  # cat("Bias:", object$o.coef - mean(object$boots, na.rm = T), "\n")
  # cat("Standard Error:", 
  #     sd(object$boots), "\n")
}
