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
  #method = inference
  method <- attr(sets, "refinement.method")
  restricted <- attr(sets, "restricted") # this doesnt exist yet, not sure what it means.
  
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  #do we need to order the data?
  data <- data[order(data[,unit.id], data[,time.id]), ]
  
  if(is.null(restricted)){restricted <- FALSE}

  if (qoi == "att" | qoi == "ate") 
  {
    sets <- prep_for_leads(sets, data, max(lead), time.id, unit.id, outcome.variable)
    sets <- sets[sapply(sets, length) > 0]
    treated.unit.ids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])
    
    data[, paste0("Wit_att", lead)] <- do.call(cbind, lapply(lead, FUN = getWits, data = data, matched_sets = sets, estimation.method = inference))
    data$dit_att <- getDits(matched_sets = sets, data = data)
    colnames(data)[length(data)] <- "dits_att"
    data$`Wit_att-1` <- 0
    
  } 
  if (qoi == "atc" | qoi == "ate") 
  {
    #first we need to "flip" the treatment variable, then re-run panelmatch with this as the new treatment variable
    data$atc_variable <- ifelse(data[, treatment] == 1, 0, 1)
    sets2 <- PanelMatch2(lag = lag, time.id = time.id, unit.id = unit.id, treatment = "atc_variable", outcome = outcome.variable,
                                     refinement.method = attr(sets, "refinement.method"),
                                     size.match = attr(sets, "max.match.size"),
                                     data = data,
                                     match.missing = attr(sets, "match.missing"),
                                     covs.formula = attr(sets, "covs.formula"),
                                     verbose = FALSE)
    sets2 <- prep_for_leads(sets2, data, max(lead), time.id, unit.id, outcome.variable)
    sets2 <- sets2[sapply(sets2, length) > 0]
    treated.unit.ids2 <- as.numeric(unlist(strsplit(names(sets2), split = "[.]"))[c(T,F)])
    
    data[, paste0("Wit_atc", lead)] <- do.call(cbind, lapply(lead, FUN = getWits, data = data, matched_sets = sets2, estimation.method = inference))
    data$dit_atc <- getDits(matched_sets = sets2, data = data)
    colnames(data)[length(data)] <- "dits_atc"
    data$`Wit_atc-1` <- 0
    
  } 
  #NOTE THE COMMENT/ASSUMPTION
  if(inference == "bootstrap")
  {
    data[, dependent][is.na(data[, dependent])] <- 0 #replace the NAs with zeroes. I think this is ok because the dits should always be zero for these, so the value is irrelevant. this just makes the implementation a little bit easier   
  }
  
  
  # ATT
  if (qoi == "att") {
    if (inference == "wfe"){

      
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]

      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "unit", 
                        qoi = "att", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = sets, "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "time", 
                        qoi = "att", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = sets, "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
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
        att_new <-  sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1[,outcome.variable],
                           z = d.sub1$dits_att)
        coefs[k,] <- att_new
      }
      # changed return to class
      z <- list("coefficients" = o.coefs,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER,
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PanelEstimate"
      return(z)
    }
    
    # ATC
  } else if (qoi == "atc"){
    if (inference == "wfe") {

      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "unit", 
                        qoi = "atc", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = sets2, "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "time", 
                        qoi = "atc", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = sets2, "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
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
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(lead))
      
      for (k in 1:ITER) {  
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        while(all(!units %in% treated.unit.ids2)) #while none of the units are treated units, resample
        {
          units <- sample(clusters, size = length(clusters), replace=T)
        }
        d.sub1 <- data[ data[,unit.id] %in% units, ]
        atc_new <- -sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_atc", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1[,outcome.variable],
                           z = d.sub1$dits_atc)
      
        coefs[k,] <- atc_new
        
        
      }
      z <- list("coefficients" = o.coefs,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = sets2)
      class(z) <- "PanelEstimate"
      return(z)

    }
    
  } else if (qoi == "ate") {
    if (inference == "wfe"){
      
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]
      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
      
      if (estimator == "did") {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "unit", 
                        qoi = "ate", estimator = "did", 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = list(sets,sets2), "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
      } else {
        fit <- PanelWFE(formula = as.formula(paste(dependent, "~", treatment)), 
                        treat = treatment, unit.index = unit.id,
                        time.index = time.id, method = "time", 
                        qoi = "ate", estimator = NULL, 
                        df.adjustment = df.adjustment,
                        hetero.se = TRUE, 
                        auto.se = TRUE, White = TRUE,  
                        data = data)
        coef <- as.numeric(fit$coefficients)
        names(coef) <- paste0("t+", lead)
        z <- list("coefficients" = coef,
                  "method" = method, "lag" = lag,
                  "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                  "standard.error" = fit$se,"matched.sets" = list(sets,sets2), "df" = fit$df)
        class(z) <- "PanelEstimate"
        return(z)
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
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(lead))
      
      
      for (k in 1:ITER) {
        # make new data
        clusters <- unique(data[, unit.id])
        units <- sample(clusters, size = length(clusters), replace=T)
        while(all(!units %in% treated.unit.ids) | all(!units %in% treated.unit.ids2)) #while none of the units are treated units (att and atc), resample
        {
          units <- sample(clusters, size = length(clusters), replace=T)
        }
        # create bootstap sample with sapply
        d.sub1 <- data[ data[,unit.id] %in% units, ]
        
        att_new <-sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_att", x)), 
                                drop = FALSE],
                         equality_four,
                         y = d.sub1[,outcome.variable],
                         z = d.sub1$dits_att)
        
        atc_new <- -sapply(d.sub1[, sapply(lead, function(x) paste0("Wit_atc", x)), 
                                  drop = FALSE],
                           equality_four,
                           y = d.sub1[,outcome.variable],
                           z = d.sub1$dits_atc)
        coefs[k,] <- (att_new*sum(d.sub1$dits_att) + atc_new*sum(d.sub1$dits_atc))/
          (sum(d.sub1$dits_att) + sum(d.sub1$dits_atc))
        
      
      }
      # return(list("o.coef" = DID_ATE, "boots" = coefs))
      z <- list("coefficients" = o.coefs_ate,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = list(sets, sets2))
      class(z) <- "PanelEstimate"
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
summary.PanelEstimate <- function(object, verbose = TRUE, bias.corrected = FALSE, ...) {
  
  if(verbose)
  {
    if(object$qoi == "ate")
    {
      refinement.method <- attr(object$matched.set[[1]], "refinement.method")
      lag <- attr(object$matched.set[[1]], "lag")
    }
    else
    {
      refinement.method <- attr(object$matched.set, "refinement.method")
      lag <- attr(object$matched.set, "lag")
    }
    if(refinement.method == "mahalanobis")
        {
          cat("Weighted Difference-in-Differences with Mahalanobis Distance\n")
        } 
    else if (refinement.method == "ps.weight" | refinement.method == "ps.match") 
        {
          cat("Weighted Difference-in-Differences with Propensity Score\n")
    } else(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match") 
        {
          cat("Weighted Difference-in-Differences with Covariate Balancing Propensity Score\n")
        }
    cat("Matches created with", attr(object$matched.sets, "lag"), "lags\n")
    if(!is.null(object$bootstrap.iterations))
    {
      cat("\nStandard errors computed with", object$bootstrap.iterations, "Weighted bootstrap samples\n")  
    }
    
    # cat("\nTotal effect in", object$lead, "periods after the treatment:")
    qoi <- ifelse(object$qoi == "ate", "Average Treatment Effect (ATE)", 
                ifelse(object$qoi == "att", "Average Treatment Effect on the Treated (ATT)", 
                       "Average Treatment Effect on the Control (ATC)"))
    cat("\nEstimate of", qoi, "by Period:\n")
  }
  if(bias.corrected)
  {
    if(is.null(object$bootstrap.iterations)) stop("bias corrected estimates only available for bootstrap method currently")
    df <- rbind(t(as.data.frame(object$coefficients)), # point estimate
              
              apply(object$bootstrapped.coefficients, 2, sd, na.rm = T), # bootstrap se
              
              # Efron & Tibshirani 1993 p170 - 171
              apply(object$bootstrapped.coefficients, 2, quantile, probs = c( (1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2 ), na.rm = T), # percentile CI
              # Efron & Tibshirani 1993 p138
              2*object$coefficients - colMeans(object$bootstrapped.coefficients, na.rm = T), # bc point estimate
              
              apply( (2*matrix(nrow = object$bootstrap.iterations, ncol = length(object$coefficients), object$coefficients, byrow = TRUE) - object$bootstrapped.coefficients), 2, quantile, 
                             probs = c((1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2), 
                             na.rm = T) ) # bc percentile CI)
          rownames(df) <- c("estimate", "std.error", 
                    paste0((1-object$confidence.level)/2 * 100, "%"),
                    paste0( (object$confidence.level+(1-object$confidence.level)/2) * 100, "%"),
                    "estimate(bias corrected)", 
                    paste0((1-object$confidence.level)/2 * 100, "%", "(bias corrected)"),
                    paste0((object$confidence.level+(1-object$confidence.level)/2) * 100, "%", "(bias corrected)"))
  }
  else
  {
    if(!is.null(object$bootstrap.iteration))
    {
      df <- rbind(t(as.data.frame(object$coefficients)), # point estimate
                  
                  apply(object$bootstrapped.coefficients, 2, sd, na.rm = T), # bootstrap se
                  
                  # Efron & Tibshirani 1993 p170 - 171
                  apply(object$bootstrapped.coefficients, 2, quantile, probs = c( (1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2 ), na.rm = T))
      rownames(df) <- c("estimate", "std.error", 
                        paste0((1-object$confidence.level)/2 * 100, "%"),
                        paste0( (object$confidence.level+(1-object$confidence.level)/2) * 100, "%"))
      tdf <- t(df)
      if(!verbose) return(t(df))
      return(list("summary" = tdf, "lag" = lag, "iterations" = object$bootstrap.iterations, "qoi" = object$qoi) )   
    }
    else #wfe method
    {
      critical.vals <- rep(qt((1 - object$confidence.level) / 2, object$df), 2) * c(1, -1)
      quants <- object$coefficients + critical.vals * object$standard.error
      df <- as.data.frame(c(object$coefficients, # point estimate
                  object$standard.error,
                  quants))
      rownames(df) <- c("estimate", "std.error", 
                        paste0((1-object$confidence.level)/2 * 100, "%"),
                        paste0( (object$confidence.level+(1-object$confidence.level)/2) * 100, "%"))
      tdf <- t(df)
      rownames(tdf) <- names(object$coefficients)
      if(!verbose) return(t(df))
      return(list("summary" = tdf, "lag" = lag, "qoi" = object$qoi) )
    }
    
  }
  
   
  
}

#' plot the point estimates and standard errors from a PanelEstimate calculation. The only mandatory argument is an object of the PanelEstimate class
#' Use standard arguments to the \code{plot} function to modify the plot as needed.
#' @param pe.object a PanelEstimate object
#' @export
plot.PanelEstimate <- function(pe.object, ylab = "Estimated Effect of Treatment", xlab = "Time", main = "Estimated Effects of Treatment Over Time",...)
{
  
  plot.data <- summary(pe.object, verbose = F, bias.corrected = F)
  # plot(x = 1:5,y = plot.data[, 1], ylim = c(min(plot.data[,3]) - .1, max(plot.data[,4]) + .1), pch = 16, 
  #      xaxt = "n", ylab = ylab, xlab = xlab, main = main, ...)
  plot(x = 1:(nrow(plot.data)),y = plot.data[, 1], pch = 16, 
       xaxt = "n", ylab = ylab, xlab = xlab, main = main, ...)
  axis(side = 1, at = 1:nrow(plot.data), labels = rownames(plot.data))
  arrows(1:(nrow(plot.data)), plot.data[,3], 1:(nrow(plot.data)), plot.data[,4], length=0.05, angle=90, code=3)
  abline(h = 0, lty = "dashed")
}


