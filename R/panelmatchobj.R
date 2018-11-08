panelmatch <- function(method = NULL, treatment, qoi, dependent, covariate, unit.id, time.id, M, 
                           covariate.only, lag, max.lead, data, ATT_matches, ATC_matches, restricted, covariate_names)
{
    
}

print.panelmatch <- function(x, ...)
{
  if(!hasArg(verbose) | (hasArg(verbose) & verbose == FALSE))
  {
    print(x[names(x) != "ATT_matches" & names(x) != "ATC_matches" & names(x) != "data"])
    s <- paste0("Number of ATT matched sets: ", length(x$ATT_matches))
    print(s)
    s <- paste0("Number of ATC matches sets: ", length(x$ATC_matches))
    print(s)  
  }
  else
  {
    print(x)
  }
  
    
}