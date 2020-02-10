#' PanelEstimate
#'
#' \code{PanelEstimate} estimates a causal quantity of interest, including the average treatment effect for 
#' treated or control units (att and atc, respectively), or average treatment effect (ate), as specified in \code{PanelMatch}.
#' This is done by estimating the counterfactual outcomes for each treated unit using
#' matched sets. Users will provide matched sets that were obtained by the
#' \code{PanelMatch} function and obtain point estimates via a
#' weighted average computation with weighted bootstrap standard errors. Point estimates and standard errors will be 
#' produced for each period in the lead window specified by the \code{lead} argument from \code{PanelMatch}. 
#' Users may run multiple estimations by providing lists of each argument to the function. 
#' However, in this format, every argument must be explicitly specified in each configuration 
#' and must adhere to the same data types/structures outlined below. See the included code examples for more about 
#' how this functionality works.
#' 
#' @param number.iterations An integer value indicating the number of bootstrap
#' iterations. The default is 1000.
#' @param sets A \code{PanelMatch} object attained via the
#' \code{PanelMatch} function.
#' @param df.adjustment A logical value indicating whether or not a
#' degree-of-freedom adjustment should be performed for the standard error
#' calculation. The default is \code{FALSE}.
#' @param confidence.level A numerical value specifying the confidence level and range of interval
#' estimates for statistical inference. The default is .95.
#' @param data The same time series cross sectional data set provided to the PanelMatch function used to produce 
#' the matched sets
#' @return \code{PanelEstimate} returns a list of class
#' `PanelEstimate' containing the following components:
#' \item{estimates}{the point estimates of the quantity of interest for the lead periods specified}
#' \item{bootstrapped.estimates}{the bootstrapped point estimate values}
#' \item{bootstrap.iterations}{the number of iterations used in bootstrapping}
#' \item{method}{refinement method used to create the matched sets from which the estimates were calculated}
#' \item{lag}{See PanelMatch argument \code{lag} for more information.}
#' \item{lead}{The lead window sequence for which \code{PanelEstimate} is producing point estimates and standard errors.}
#' \item{confidence.level}{the confidence level}
#' \item{qoi}{the quantity of interest}
#' \item{matched.sets}{the refined matched sets used to produce the estimations}
#' \item{standard.error}{the standard error(s) of the point estimates}
#' 
#' @references Imai, Kosuke, In Song Kim, and Erik Wang (2018)
#' @author In Song Kim <insong@mit.edu>, Erik Wang
#' <haixiao@Princeton.edu>, Adam Rauh <adamrauh@mit.edu>, and Kosuke Imai <imai@harvard.edu>
#'
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)), 
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = TRUE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem)
#'
#' 
#'
#' @export
PanelEstimate <- function(sets,
                          number.iterations = 1000,
                          df.adjustment = FALSE,
                          confidence.level = .95,
                          data) 
{
  inference <- "bootstrap"
  if(inference == "wfe") stop("wfe is no longer supported. Please specify inference = 'bootstrap'")
  if(class(number.iterations) == "list" & class(df.adjustment) == "list" & class(confidence.level) == "list" & class(sets) == "list")
  {
    if(length(unique(length(inference), length(number.iterations), length(df.adjustment), length(confidence.level), length(sets))) == 1)
    {
      res = mapply(FUN = panel_estimate, number.iterations = number.iterations, 
                   df.adjustment = df.adjustment, confidence.level= confidence.level, sets = sets, 
                   MoreArgs = list(data = data, inference = inference),
                   SIMPLIFY = FALSE)
    }
    else {
      stop("arguments are not provided in equal length lists")
    }
  }
  else 
  {
    res = panel_estimate(inference = inference, number.iterations = number.iterations, df.adjustment = df.adjustment, confidence.level = confidence.level, sets = sets, data = data)
  }
  return(res)
}




panel_estimate <- function(inference = "bootstrap",
                           number.iterations = 1000,
                           df.adjustment = FALSE,
                           confidence.level = .95,
                           sets,
                           data)
{
  
  lead <- attr(sets, "lead")
  outcome.variable <- attr(sets, "outcome.var")
  
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
  
  
  lag <- attr(sets, "lag")
  dependent = outcome.variable
  treatment <- attr(sets, "treatment.var")
  unit.id <- attr(sets, "id.var")
  time.id <- attr(sets, "t.var")
  #method = inference
  method <- attr(sets, "refinement.method")
  
  forbid.treatment.reversal <- attr(sets, "forbid.treatment.reversal") # this doesnt exist yet, not sure what it means.
  #add in checks about forbid.treatment.reversal and wfe, etc. 
  
  if(!"data.frame" %in% class(data)) stop("please convert data to data.frame class")
  
  if(!class(data[, unit.id]) %in% c("integer", "numeric")) stop("please convert unit id column to integer or numeric")
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
    
    for(j in lead)
    {
      dense.wits <- getWits(lead = j, data = data, matched_sets = sets, estimation.method = inference)
      data = merge(x = data, y = dense.wits, all.x = TRUE, by.x = colnames(data)[1:2], by.y = c("id", "t"))
      colnames(data)[length(data)] <- paste0("Wit_att", j)
      data[is.na(data[, length(data)]), length(data)] <- 0 #replace NAs with zeroes
    }
    
    data$dit_att <- getDits(matched_sets = sets, data = data)
    colnames(data)[length(data)] <- "dits_att"
    data$`Wit_att-1` <- 0
    
  } 
  if (qoi == "atc" | qoi == "ate") 
  {
    treated.unit.ids2 <- as.numeric(unlist(strsplit(names(sets2), split = "[.]"))[c(T,F)])
    
    for(j in lead)
    {
      dense.wits <- getWits(lead = j, data = data, matched_sets = sets2, estimation.method = inference)
      data = merge(x = data, y = dense.wits, all.x = TRUE, by.x = colnames(data)[1:2], by.y = c("id", "t"))
      colnames(data)[length(data)] <- paste0("Wit_atc", j)
      data[is.na(data[, length(data)]), length(data)] <- 0 #replace NAs with zeroes
    }
    
    
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
    if (inference == "bootstrap")
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
      coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
      
      for (k in 1:number.iterations) 
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
      z <- list("estimates" = o.coefs,
                "bootstrapped.estimates" = coefs, "bootstrap.iterations" = number.iterations, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "method" = method, "lag" = lag,
                "lead" = lead, "confidence.level" = confidence.level, "qoi" = qoi, "matched.sets" = sets)
      class(z) <- "PanelEstimate"
      return(z)
    }
  } else if (qoi == "atc")
  {
    if (inference == "bootstrap") 
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
      
      
      coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
      
      for (k in 1:number.iterations) 
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
      z <- list("estimates" = o.coefs,
                "bootstrapped.estimates" = coefs, "bootstrap.iterations" = number.iterations, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lead" = lead, "confidence.level" = confidence.level, "qoi" = qoi, "matched.sets" = sets2)
      class(z) <- "PanelEstimate"
      return(z)
      
    }
  } else if (qoi == "ate") 
  {
    if (inference == "bootstrap")
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
      
      coefs <- matrix(NA, nrow = number.iterations, ncol = length(lead))
      
      
      for (k in 1:number.iterations) {
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
      z <- list("estimates" = o.coefs_ate,
                "bootstrapped.estimates" = coefs, "bootstrap.iterations" = number.iterations, "standard.error" = apply(coefs, 2, sd, na.rm = T),
                "lead" = lead, "confidence.level" = confidence.level, "qoi" = qoi, "matched.sets" = list(sets, sets2))
      class(z) <- "PanelEstimate"
      return(z)
    }
  }
  
  
}


#' Get summaries of PanelEstimate objects/calculations
#' 
#' 
#' \code{summary.PanelEstimate} takes an object returned by
#' \code{PanelEstimate}, and returns a summary table of point
#' estimates and confidence intervals
#'
#' @param object A PanelEstimate object
#' @param verbose logical indicating whether or not output should be printed in an expanded form. Default is TRUE
#' @param bias.corrected logical indicating whether or not bias corrected estimates should be provided. Default is FALSE
#' @param ... optional additional arguments. Currently, no additional arguments are supported. 
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem)
#' summary(PE.results)
#' 
#'
#' 
#' @method summary PanelEstimate
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
    if (refinement.method == "ps.weight" | refinement.method == "ps.match") 
        {
          cat("Weighted Difference-in-Differences with Propensity Score\n")
        } 
    if(refinement.method == "CBPS.weight" | refinement.method == "CBPS.match") 
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
    df <- rbind(t(as.data.frame(object$estimates)), # point estimate
              
              apply(object$bootstrapped.estimates, 2, sd, na.rm = T), # bootstrap se
              
              # Efron & Tibshirani 1993 p170 - 171
              apply(object$bootstrapped.estimates, 2, quantile, probs = c( (1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2 ), na.rm = T), # percentile confidence.level
              # Efron & Tibshirani 1993 p138
              2*object$estimates - colMeans(object$bootstrapped.estimates, na.rm = T), # bc point estimate
              
              apply( (2*matrix(nrow = object$bootstrap.iterations, ncol = length(object$estimates), object$estimates, byrow = TRUE) - object$bootstrapped.estimates), 2, quantile, 
                             probs = c((1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2), 
                             na.rm = T) ) # bc percentile confidence.level)
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
      df <- rbind(t(as.data.frame(object$estimates)), # point estimate
                  
                  apply(object$bootstrapped.estimates, 2, sd, na.rm = T), # bootstrap se
                  
                  # Efron & Tibshirani 1993 p170 - 171
                  apply(object$bootstrapped.estimates, 2, quantile, probs = c( (1-object$confidence.level)/2, object$confidence.level+(1-object$confidence.level)/2 ), na.rm = T))
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
      quants <- object$estimates + critical.vals * object$standard.error
      df <- as.data.frame(c(object$estimates, # point estimate
                  object$standard.error,
                  quants))
      rownames(df) <- c("estimate", "std.error", 
                        paste0((1-object$confidence.level)/2 * 100, "%"),
                        paste0( (object$confidence.level+(1-object$confidence.level)/2) * 100, "%"))
      tdf <- t(df)
      rownames(tdf) <- names(object$estimates)
      if(!verbose) return(t(df))
      return(list("summary" = tdf, "lag" = lag, "qoi" = object$qoi) )
    }
    
  }
  
   
  
}


#' Plot point estimates and standard errors from a PanelEstimate calculation.
#' 
#' 
#' The \code{plot.PanelEstimate} method takes an object returned by the \code{PanelEstimate} function and plots the calculated 
#' point estimates and standard errors over the specified \code{lead} time period. 
#' The only mandatory argument is an object of the \code{PanelEstimate} class.
#'
#' @param x a \code{PanelEstimate} object
#' @param ylab default is "Estimated Effect of Treatment. This is the same argument as the standard argument for \code{plot}
#' @param xlab default is "Time". This is the same argument as the standard argument for \code{plot}
#' @param main default is "Estimated Effects of Treatment Over Time". This is the same argument as the standard argument for \code{plot}
#' @param ylim default is NULL. This is the same argument as the standard argument for \code{plot}
#' @param ... Additional optional arguments to be passed to \code{plot}.
#' @examples
#' PM.results <- PanelMatch(lag = 4, time.id = "year", unit.id = "wbcode2", 
#'                          treatment = "dem", refinement.method = "mahalanobis", 
#'                          data = dem, match.missing = TRUE, 
#'                          covs.formula = ~ I(lag(tradewb, 1:4)) + I(lag(y, 1:4)),
#'                          size.match = 5, qoi = "att",
#'                          outcome.var = "y", lead = 0:4, forbid.treatment.reversal = FALSE)
#' PE.results <- PanelEstimate(sets = PM.results, data = dem)
#' plot(PE.results)
#'
#'
#' @method plot PanelEstimate
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
  graphics::segments(1:(nrow(plot.data)), plot.data[,3], 1:(nrow(plot.data)), plot.data[,4])
  graphics::abline(h = 0, lty = "dashed")
}


