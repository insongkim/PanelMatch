#' calculate_point_estimates
#' Helper function that calculates the point estimates for the specified QOI
#'
#' @param qoi.in string specifying the QOI
#' @param data.in data.frame providing the processed/parsed data to be used for calculations
#' @param lead see PanelMatch() documentation
#' @param outcome.variable string specifying the outcome variable
#' @param pooled Logical. See PanelEstimate() documentation.
#'
#' @return A named vector of point estimates
#' @keywords internal
calculate_point_estimates <- function(qoi.in, data.in, lead,
                                      outcome.variable,
                                      pooled = FALSE)
{
  if ( identical(qoi.in, "att") ||
       identical(qoi.in, "atc") ||
       identical(qoi.in, "art"))
  {
    
    col.idx <- sapply(lead, function(x) paste0("Wit_", qoi.in, x))
    x.in <- data.in[, col.idx, drop = FALSE]
    y.in <- data.in[c(outcome.variable)][,1]
    z.in <- data.in[, paste0("dits_", qoi.in)]
    
    o.coefs <- sapply(x.in, equality_four,
                      y = y.in,
                      z = z.in)
    
    #do coefficient flip for atc
    if (identical(qoi.in, "atc")) o.coefs <- -o.coefs
    
    
    
    if (all(lead >= 0)) 
    {
      if (length(lead[lead < 0]) > 1)
      {
        names(o.coefs)[(length(o.coefs) - 
                          max(lead[lead >= 0])):length(o.coefs)] <-
          sapply(lead[lead >= 0], function(x) paste0("t+", x))
        names(o.coefs)[(length(o.coefs) - 
                          length(lead) + 1):length(lead[lead < 0])] <-
          sapply(lead[lead < 0], function(x) paste0("t", x))
        
      } else
      {
        names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
      }
    }
    
    if (pooled)
    {
      o.coefs <- mean(o.coefs, na.rm = TRUE)
      names(o.coefs) <- NULL
    }
  } else if (identical(qoi.in, "ate")) {
    o.coefs_att <-  sapply(data.in[, sapply(lead, function(x) paste0("Wit_att", x)),
                                   drop = FALSE],
                           equality_four,
                           y = data.in[c(outcome.variable)][,1],
                           z = data.in$dits_att)
    
    o.coefs_atc <-  -sapply(data.in[, sapply(lead, function(x) paste0("Wit_atc",x)),
                                    drop = FALSE],
                            equality_four,
                            y = data.in[c(outcome.variable)][,1],
                            z = data.in$dits_atc)
    
    o.coefs <- (o.coefs_att*sum(data.in$dits_att) + 
                  o.coefs_atc*sum(data.in$dits_atc))/
      (sum(data.in$dits_att) + sum(data.in$dits_atc))
    
    
    
    if (length(lead[lead<0]) > 1)
    {
      names(o.coefs)[(length(o.coefs)-max(lead[lead>=0])):
                       length(o.coefs)] <- sapply(lead[lead>=0], 
                                                  function(x) paste0("t+", 
                                                                     x))
      names(o.coefs)[(length(o.coefs)-length(lead) + 1):
                       length(lead[lead<0])] <- sapply(lead[lead<0], 
                                                       function(x) paste0("t", 
                                                                          x))
      
    } else
    {
      names(o.coefs) <- sapply(lead, function(x) paste0("t+", x))
    }
    
    if (pooled)
    {
      o.coefs <- mean(o.coefs, na.rm = TRUE)
      names(o.coefs) <- NULL
    }
  }
  return(o.coefs)
}


#' equality_four
#' Small helper function implementing estimation function from Imai, Kim, and Wang (2023)
#' @return Returns numeric vector of results.
#' @keywords internal
equality_four <- function(x, y, z){
  
  return(sum(x*y)/sum(z))
}