#' @export
PanelEstimate2 <- function(lead, #probably want to swap the order of these around to be more intuitive
                          inference = c("wfe", "bootstrap"),
                          ITER = 1000, matched_sets = NULL,
                          estimator = "did",
                          df.adjustment = FALSE, qoi = NULL,
                          CI = .95,
                          sets,
                          data,
                          outcome.variable) {

  if (inference == "wfe" & length(lead) > 1) 
    stop("When inference method is wfe, please only supply 1 lead at a time. 
         For example, please call this function with `lead` = 1 and then call it with `lead` = 2,
         rather than supplying `lead`` = 1:2")
  lag <- attr(sets, "lag")
  dependent = outcome.variable
  treatment <- attr(sets, "treated.var")
  unit.id <- attr(sets, "id.var")
  time.id <- attr(sets, "t.var")
  method <- attr(sets, "refinement.method")
  restricted <- attr(sets, "restricted") # this doesnt exist yet, not sure what it means.
  #check the leads, update sets, and weights
  #browser()
  if(is.null(restricted)){restricted <- FALSE}
  
  sets <- prep_for_leads(sets, data, max(lead), time.id, unit.id, outcome.variable)
  sets <- sets[sapply(sets, length) > 0]
  treated.unit.ids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])
  # WHAT TO DO WHEN EMPTY SETS ARE PRODUCED AS A RESULT OF THIS?? REMOVE THEM I ASSUME?
  
  if (is.null(qoi)) {
    qoi = matched_sets$qoi
  } else {
    qoi = qoi
  }
  
  
  # DONT KNOW WHAT THESE ARE DOING
  if (qoi == "att") {
    
    data[, paste0("Wit_att", lead)] <- do.call(cbind, lapply(lead, FUN = getWits, data = data, matched_sets = sets))
    data$dit_att <- getDits(matched_sets = sets, data = data)
    colnames(data)[length(data)] <- "dits_att"
    data$`Wit_att-1` <- 0
    ##NOTE THE COMMENT/ASSUMPTION
    data[, dependent][is.na(data[, dependent])] <- 0 #replace the NAs with zeroes. I think this is ok because the dits should always be zero for these, so the value is irrelevant. this just makes the implementation a little bit easier 
    
  } else if (qoi == "atc") {
    newlist <- lapply(matched_sets$`ATC_matches`, lapply_leads, unit.id = unit.id, 
                      time.id = time.id, lag = lag, estimator = estimator,
                      inference = inference,
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
                      inference = inference,
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
                      inference = inference,
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
  #browser()
  if (qoi == "att") {
    if (inference == "wfe"){

      
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]


      #browser()
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
      
      #browser()
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
      coefs <- matrix(NA, nrow = ITER, ncol = length(lead))
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        while(all(!units %in% treated.unit.ids)) #while none of the units are treated units, resample
        {
          units <- sample(clusters, size = length(clusters), replace=T)
        }
        # create bootstap sample with sapply
        d.sub1 <- data[ data[,unit.id] %in% units, ]
        #df.bs <- lapply(units, function(x) which(data[,unit.id]==x))
        #d.sub1 <- data[unlist(df.bs),]
        #colnames(d.sub1)[3:4] <- c("treatment", "dv")
        att_new <-  sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1[,outcome.variable],
                           z = d.sub1$dits_att)
        #if(is.na(att_new)) browser()
        coefs[k,] <- att_new
      }
      # changed return to class
      z <- list("o.coef" = o.coefs,
                "boots" = coefs, "ITER" = ITER,
                "method" = method, "lag" = lag,
                "lead" = lead, "CI" = CI, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PE"
      z
    }
    
    # ATC
  } else if (qoi == "atc"){
    if (inference == "wfe") {

      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
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
       
        coefs[k,] <- atc_new
        
        
      }
      z <- list("o.coef" = o.coefs,
                "boots" = coefs, "ITER" = ITER,
                "lead" = lead, "CI" = CI, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PE"
      z

    }
    
  } else if (qoi == "ate") {
    if (inference == "wfe"){
      
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]
      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
      
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
        
      
      }
      # return(list("o.coef" = DID_ATE, "boots" = coefs))
      z <- list("o.coef" = o.coefs_ate,
                "boots" = coefs, "ITER" = ITER,
                "lead" = lead, "CI" = CI, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PE"
      return(z)
    }
    
    
  }
}

#' Get summaries of PanelEstimate objects
#'
#' \code{summary.PanelEstimate()} takes an object returned by
#' \code{PanelEstimate}, and returns a summary table of point
#' estimates and the confidence intervales.
#'
#' @usage \method{summary}{PanelEstimate}(object, ...)
#' @param object A PanelEstimate object
#' @param ... Further arguments to be passed to \code{summary.PanelEstimate()}.
#' 
#' @export
summary.PE <- function(object, verbose = TRUE, bias.corrected = FALSE, ...) {
  
  if(verbose)
  {
      if(attr(object$matched.sets, "refinement.method") == "Maha")
        {
          cat("Weighted Difference-in-Differences with Mahalanobis Distance\n")
        } 
      else if (attr(object$matched.sets, "refinement.method") == "Pscore") 
        {
          cat("Weighted Difference-in-Differences with Propensity Score\n")
        } 
      else if (attr(object$matched.sets, "refinement.method") == "Synth") 
        {
          cat("Weighted Difference-in-Differences with Synthetic Control\n")
        } 
      else if (attr(object$matched.sets, "refinement.method") == "CBPS") 
        {
          cat("Weighted Difference-in-Differences with Covariate Balancing Propensity Score\n")
        }
    cat("Matches created with", attr(object$matched.sets, "lag"), "lags\n")
    cat("\nStandard errors computed with", object$ITER, "Weighted bootstrap samples\n")
    # cat("\nTotal effect in", object$lead, "periods after the treatment:")
    qoi <- ifelse(object$qoi == "ate", "Average Treatment Effect (ATE)", 
                ifelse(object$qoi == "att", "Average Treatment Effect on the Treated (ATT)", 
                       "Average Treatment Effect on the Control (ATC)"))
    cat("\nEstimate of", qoi, "by Period:\n")
  }
  if(bias.corrected)
  {
    df <- rbind(t(as.data.frame(object$o.coef)), # point estimate
              
              apply(object$boots, 2, sd, na.rm = T), # bootstrap se
              
              # Efron & Tibshirani 1993 p170 - 171
              apply(object$boots, 2, quantile, probs = c( (1-object$CI)/2, object$CI+(1-object$CI)/2 ), na.rm = T), # percentile CI
              # Efron & Tibshirani 1993 p138
              2*object$o.coef - colMeans(object$boots, na.rm = T), # bc point estimate
              
              apply( (2*matrix(nrow = object$ITER, ncol = length(object$o.coef), object$o.coef, byrow = TRUE) - object$boots), 2, quantile, 
                             probs = c((1-object$CI)/2, object$CI+(1-object$CI)/2), 
                             na.rm = T) ) # bc percentile CI)
  rownames(df) <- c("estimate", "std.error", 
                    paste0((1-object$CI)/2 * 100, "%"),
                    paste0( (object$CI+(1-object$CI)/2) * 100, "%"),
                    "estimate(bias corrected)", 
                    paste0((1-object$CI)/2 * 100, "%", "(bias corrected)"),
                    paste0((object$CI+(1-object$CI)/2) * 100, "%", "(bias corrected)"))
  }
  else
  {
    df <- rbind(t(as.data.frame(object$o.coef)), # point estimate
              
              apply(object$boots, 2, sd, na.rm = T), # bootstrap se
              
              # Efron & Tibshirani 1993 p170 - 171
              apply(object$boots, 2, quantile, probs = c( (1-object$CI)/2, object$CI+(1-object$CI)/2 ), na.rm = T))
    rownames(df) <- c("estimate", "std.error", 
                    paste0((1-object$CI)/2 * 100, "%"),
                    paste0( (object$CI+(1-object$CI)/2) * 100, "%"))
  }
  
  tdf <- t(df)
  if(!verbose) return(t(df))
  return(list("summary" = tdf, "lag" = attr(object$matched.sets, "lag"), "iterations" = object$ITER, "qoi" = object$qoi) ) 
  # cat("Bias:", object$o.coef - mean(object$boots, na.rm = T), "\n")
  # cat("Standard Error:", 
  #     sd(object$boots), "\n")
}


#' @export
plot.PE <- function(pe.object, ylab = "Estimated Effect of Treatment", xlab = "Time", main = "Estimated Effects of Treatment Over Time",...)
{
  plot.data <- summary(pe.object, verbose = F, bias.corrected = F)
  plot(x = 1:5,y = plot.data[, 1], ylim = c(min(plot.data[,3]) - .1, max(plot.data[,4]) + .1), pch = 16, 
       xaxt = "n", ylab = ylab, xlab = xlab, main = main, ...)
  axis(side = 1, at = 1:nrow(plot.data), labels = rownames(plot.data))
  arrows(1:5, plot.data[,3], 1:5, plot.data[,4], length=0.05, angle=90, code=3)
  abline(h = 0, lty = "dashed")
}