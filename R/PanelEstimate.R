#' PanelEstimate
#'
#' \code{PanelEstimate} estimates a causal quantity of interest, including the average treatment effect for treated or control units (att and atc, respectively), or average treatment effect (ate), as specified in \code{PanelMatch}
#' This is done by estimating the counterfactual outcomes for each treated unit using
#' matched sets. Users will specify matched sets that were obtained by the
#' \code{PanelMatch} function and obtain point estimates via weighted fixed effects regressions or via
#' weighted average computation with weighted bootstrap standard errors.
#' 
#' @param inference One of ``wfe'' (weighted fixed effects) or
#' ``bootstrap'' methods for standard error calculation. The default
#' is \code{bootstrap}.
#' @param ITER An integer value indicating the number of bootstrap
#' iteration. The default is 1000.
#' @param sets A list of class `PanelMatch' attained by
#' \code{PanelMatch}.
#' @param df.adjustment A logical value indicating whether a
#' degree-of-freedom adjustment should be performed for standard error
#' calculation. The default is \code{FALSE}.
#' @param CI A numerical value specifying the range of interval
#' estimates for statistical inference. The default is .95.
#' @param data The same time series cross sectional data set provided to the PanelMatch function to produce the \code{sets}
#' @return \code{PanelEstimate} returns a list of class
#' `PanelEstimate' containing the following components:
#' \item{coefficients}{the point estimates of the quantity of interest}
#' \item{bootstrapped.coefficients}{the bootstrapped coefficients, if applicable}
#' \item{bootstrap.iterations}{the number of iterations, if applicable}
#' \item{method}{refinement method used to create the matched sets from which the estimates were calculated}
#' \item{lag}{See PanelMatch argument \code{lag} for more information.}
#' \item{lead}{The lead window sequence for which PanelEstimate is producing point estimates and standard errors.}
#' \item{confidence.level}{the confidence interval range}
#' \item{qoi}{the quantity of interest}
#' \item{matched.sets}{the refined matched sets used to produce the estimations}
#' \item{df}{if \code{inference} is "wfe", the degrees of freedom}
#' \item{standard.error}{the standard error of the point estimates}
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <adamrauh@mit.edu>, and Kosuke Imai <kimai@Princeton.edu>
#'
#' @examples \dontrun{
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = T, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = TRUE)
#' PE.results <- PanelEstimate(inference = "bootstrap", sets = PM.results, data = dem)
#' 
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = T, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0, forbid.treatment.reversal = TRUE)
#' PE.results <- PanelEstimate(inference = "wfe", sets = PM.results, data = dem)
#' }
#' @export
PanelEstimate <- function(inference = c("wfe", "bootstrap"),
                          ITER = 1000,
                          df.adjustment = FALSE,
                          CI = .95,
                          sets,
                          data) 
{
  if(class(inference) == "list" & class(ITER) == "list" & class(df.adjustment) == "list" & class(CI) == "list" & class(sets) == "list")
  {
    if(length(unique(length(inference), length(ITER), length(df.adjustment), length(CI), length(sets))) == 1)
    {
      res = mapply(FUN = panel_estimate, inference = inference, ITER = ITER, 
                   df.adjustment = df.adjustment, CI= CI, sets = sets, 
                   MoreArgs = list(data = data),
                   SIMPLIFY = FALSE)
    }
    else {
      stop("arguments are not provided in equal length lists")
    }
  }
  else 
  {
    res = panel_estimate(inference = inference, ITER = ITER, df.adjustment = df.adjustment, CI = CI, sets = sets, data = data)
  }
  return(res)
}




panel_estimate <- function(inference = c("wfe", "bootstrap"),
                           ITER = 1000,
                           df.adjustment = FALSE,
                           CI = .95,
                           sets,
                           data)
{
  
  lead <- attr(sets, "lead")
  outcome.variable <- attr(sets, "outcome.var")
  if(!attr(sets, "forbid.treatment.reversal") & inference == "wfe")
  {
    stop("WFE cannot be used unless treatment stability is guaranteed. Please re-run PanelMatch with forbid.treatment.reversal = TRUE")
  }
  if(class(sets) != "PanelMatch") stop("sets is not a PanelMatch object")
  qoi <- attr(sets, "qoi")
  if(qoi == "ate")
  {
    temp.sets <- sets
    sets <- sets[["att"]] #just picking one of the two because they should be the same
  }
  else
  {
    sets <- sets[[qoi]]  
  }
  sets <- sets[sapply(sets, length) > 0]
  if (inference == "wfe" & length(lead) > 1) 
    stop("When inference method is wfe, please only supply 1 lead at a time. 
         For example, please call this function with `lead` = 1 and then call it with `lead` = 2,
         rather than supplying `lead`` = 1:2. Please re-run PanelMatch in accordance with this.")
  
  lag <- attr(sets, "lag")
  dependent = outcome.variable
  treatment <- attr(sets, "treated.var")
  unit.id <- attr(sets, "id.var")
  time.id <- attr(sets, "t.var")
  #method = inference
  method <- attr(sets, "refinement.method")
  
  forbid.treatment.reversal <- attr(sets, "forbid.treatment.reversal") # this doesnt exist yet, not sure what it means.
  #add in checks about forbid.treatment.reversal and wfe, etc. 
  
  if(!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  
  if(class(data[, unit.id]) == "factor") stop("please convert unit id column to character, integer, or numeric")
  if(class(data[, time.id]) != "integer") stop("please convert time id to consecutive integers")
  
  if(any(table(data[, unit.id]) != max(table(data[, unit.id]))))
  {
    testmat <- data.table::dcast(data.table::as.data.table(data), formula = paste0(unit.id, "~", time.id),
                                 value.var = treatment)
    d <- data.table::melt(data.table(testmat), id = unit.id, variable = time.id, value = treatment,
                          variable.factor = FALSE, value.name = treatment)
    d <- data.frame(d)[,c(1,2)]
    class(d[, 2]) <- "integer"
    data <- merge(data.table(d), data.table(data), all.x = TRUE, by = c(unit.id, time.id))
    data <- as.data.frame(data)
    #data <- make.pbalanced(data, balance.type = "fill", index = c(unit.id, time.id))
  }
  check_time_data(data, time.id)
  
  data <- data[order(data[,unit.id], data[,time.id]), ]
  if(any(is.na(data[, unit.id]))) stop("Cannot have NA unit ids")
  data[, paste0(unit.id, ".int")] <- as.integer(as.factor(data[, unit.id]))
  
  if(class(data[, unit.id]) == "character") {
    unit.index.map <- data.frame(original.id = make.names(as.character(unique(data[, unit.id]))), new.id = unique(data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(data[, unit.id]) == "integer") {
    unit.index.map <- data.frame(original.id = (as.character(unique(data[, unit.id]))), new.id = unique(data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
  }
  else if(class(data[, unit.id]) == "numeric") {
    if(all(unique(data[, unit.id]) == as.integer(unique(data[, unit.id])))) #actually integers
    {
      unit.index.map <- data.frame(original.id = (as.character(unique(data[, unit.id]))), new.id = unique(data[, paste0(unit.id, ".int")]), stringsAsFactors = F)
    }
    else
    {
      stop("Unit ID data appears to be a non-integer numeric. Please convert.")
    }
  }
  else {
    stop("Unit ID Data is not integer, numeric, or character.")
  }
  og.unit.id <- unit.id
  unit.id <- paste0(unit.id, ".int")
  othercols <- colnames(data)[!colnames(data) %in% c(time.id, unit.id, treatment)]
  data <- data[, c(unit.id, time.id, treatment, othercols)] #reorder columns 
  
  if(qoi == "att")
  {
    sets <- encode_index(sets, unit.index.map, unit.id)
  }
  if(qoi == "atc")
  {
    sets2 <- encode_index(sets, unit.index.map, unit.id)
  }
  if(qoi == "ate")
  {
    sets <- encode_index(temp.sets$att, unit.index.map, unit.id)
    sets2 <- encode_index(temp.sets$atc, unit.index.map, unit.id)
  }
  
  if (qoi == "att" | qoi == "ate") 
  {
    treated.unit.ids <- as.numeric(unlist(strsplit(names(sets), split = "[.]"))[c(T,F)])
    data[, paste0("Wit_att", lead)] <- do.call(cbind, lapply(lead, FUN = getWits, data = data, matched_sets = sets, estimation.method = inference))
    data$dit_att <- getDits(matched_sets = sets, data = data)
    colnames(data)[length(data)] <- "dits_att"
    data$`Wit_att-1` <- 0
    
  } 
  if (qoi == "atc" | qoi == "ate") 
  {
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
  
  
  if (qoi == "att") 
  {
    if (inference == "wfe")
    {
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]
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
      sets <- decode_index(sets, unit.index.map, og.unit.id)
      z <- list("coefficients" = coef,
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                "standard.error" = fit$se,"matched.sets" = sets, "df" = fit$df)
      class(z) <- "PanelEstimate"
      return(z)
      
    } else if (inference == "bootstrap")
    {
      o.coefs <- sapply(data[, sapply(lead, function(x) paste0("Wit_att", x)), drop = FALSE],
                        equality_four,
                        y = data[c(dependent)][,1],
                        z = data$dits_att)
      
      if (length(lead[lead<0]) > 1) 
      {
        names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                         length(o.coefs)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                         length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else 
      {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
      coefs <- matrix(NA, nrow = ITER, ncol = length(lead))
      
      for (k in 1:ITER) 
      {  
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
      sets <- decode_index(sets, unit.index.map, og.unit.id)
      # changed return to class
      z <- list("coefficients" = o.coefs,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PanelEstimate"
      return(z)
    }
  } else if (qoi == "atc")
  {
    if (inference == "wfe") 
    {
      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
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
      sets2 <- decode_index(sets2, unit.index.map, og.unit.id)
      z <- list("coefficients" = coef,
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                "standard.error" = fit$se,"matched.sets" = sets2, "df" = fit$df)
      class(z) <- "PanelEstimate"
      return(z)
      
    } else if (inference == "bootstrap") 
    {
      o.coefs <-  -sapply(data[, sapply(lead, function(x) paste0("Wit_atc", x)), drop = FALSE],
                          equality_four,
                          y = data[c(dependent)][,1],
                          z = data$dits_atc)
      
      if (length(lead[lead<0]) > 1) 
      {
        names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                         length(o.coefs)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                         length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else 
      {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
      
      
      coefs <- matrix(NA, nrow = ITER, ncol = length(lead))
      
      for (k in 1:ITER) 
      {  
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
      sets2 <- decode_index(sets2, unit.index.map, og.unit.id)
      z <- list("coefficients" = o.coefs,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = sets2)
      class(z) <- "PanelEstimate"
      return(z)
      
    }
  } else if (qoi == "ate") 
  {
    if (inference == "wfe")
    {
      
      data$Wit_att0 <- data[c(paste0("Wit_att", lead))][,1]
      data$Wit_atc0 <- data[c(paste0("Wit_atc", lead))][,1]
      
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
      sets <- decode_index(sets, unit.index.map, og.unit.id)
      sets2 <- decode_index(sets2, unit.index.map, og.unit.id)
      z <- list("coefficients" = coef,
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, 
                "standard.error" = fit$se,"matched.sets" = list(sets,sets2), "df" = fit$df)
      class(z) <- "PanelEstimate"
      return(z)
      
    } else if (inference == "bootstrap")
    {
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
      
      if (length(lead[lead<0]) > 1) 
      {
        names(o.coefs_ate)[(length(o.coefs_ate)-max(lead[lead>=0])):
                             length(o.coefs_ate)] <- sapply(lead[lead>=0], function(x) paste0("t+", x))
        names(o.coefs_ate)[(length(o.coefs_ate)-length(lead) + 1):
                             length(lead[lead<0])] <- sapply(lead[lead<0], function(x) paste0("t", x))
        
      } else 
      {
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
      sets <- decode_index(sets, unit.index.map, og.unit.id)
      sets2 <- decode_index(sets2, unit.index.map, og.unit.id)
      z <- list("coefficients" = o.coefs_ate,
                "bootstrapped.coefficients" = coefs, "bootstrap.iterations" = ITER, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lead" = lead, "confidence.level" = CI, "qoi" = qoi, "matched.sets" = list(sets, sets2))
      class(z) <- "PanelEstimate"
      return(z)
    }
  }
  
  
}


#' Get summaries of PanelEstimate objects
#' \code{summary.PanelEstimate()} takes an object returned by
#' \code{PanelEstimate}, and returns a summary table of point
#' estimates and the confidence intervales.
#'
#' @param object A PanelEstimate object
#' @param verbose logical indicating whether or not output should be output in an expanded form.
#' @param bias.corrected logical indicating whether or not bias corrected estimates should be provided.
#' @param ... Further arguments to be passed to \code{summary.PanelEstimate()}.
#'
#' @method summary PanelEstimate
#' @export
summary.PanelEstimate <- function(object,..., verbose = TRUE, bias.corrected = FALSE) {
  
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
    
    if(object$qoi == "att")
    {
      qoi <- "Average Treatment Effect on the Treated (ATT)"
    }
    if(object$qoi == "atc")
    {
      qoi <- "Average Treatment Effect on the Control (ATC)"
    }
    if(object$qoi == "ate")
    {
      qoi <- "Average Treatment Effect (ATE)"
    }
    
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

#' @param x a PanelEstimate object
#' @param ylab default is "Estimated Effect of Treatment. This is the same argument as the standard argument for \code{plot}
#' @param xlab default is "Time". This is the same argument as the standard argument for \code{plot}
#' @param main default is "Estimated Effects of Treatment Over Time". This is the same argument as the standard argument for \code{plot}
#' @param ylim default is NULL. This is the same argument as the standard argument for \code{plot}
#' @param ... Additional optional arguments to \code{plot}.
#' 
#' @export
plot.PanelEstimate <- function(x, ylab = "Estimated Effect of Treatment", 
                               xlab = "Time", main = "Estimated Effects of Treatment Over Time", ylim = NULL, ...)
{
  
  pe.object <- x
  plot.data <- summary(pe.object, verbose = F, bias.corrected = F)
  if(is.null(ylim))
  {
    ylim <- c(min(plot.data[, 3]) - abs(mean(plot.data[, 3])), max(plot.data[, 4]) + abs(mean(max(plot.data[, 4]))))
  }
  graphics::plot(x = 1:(nrow(plot.data)),y = plot.data[, 1], pch = 16, cex = 1.5,
       xaxt = "n", ylab = ylab, xlab = xlab, main = main, ylim = ylim, ...)
  graphics::axis(side = 1, at = 1:nrow(plot.data), labels = rownames(plot.data))
  graphics::segments(1:(nrow(plot.data)), plot.data[,3], 1:(nrow(plot.data)), plot.data[,4])#, length=0.05, angle=90, code=3)
  graphics::abline(h = 0, lty = "dashed")
}


