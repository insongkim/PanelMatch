# The calculation of point estimates and standard errors first 
# requires the calculation of a variety of different weights, parameters, and 
# indicator variables.
# This function prepares the data within PanelEstimate() such that the estimates
# can be calculated easily.
prepare_data <- function(data.in, lead, sets.att = NULL,
                        sets.atc = NULL,
                        qoi.in, dependent.variable)
{
  if ( identical(qoi.in, "att") || 
       identical(qoi.in, "art") || 
       identical(qoi.in, "ate"))
  {
    if (!identical(qoi.in, "ate")) qoi.t <- qoi.in
    if (identical(qoi.in, "ate")) qoi.t <- "att"
    for (j in lead)
    {
      
      dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.att)
      data.in = merge(x = data.in, y = dense.wits, all.x = TRUE,
                   by.x = colnames(data.in)[1:2], by.y = c("id", "t"))
      colnames(data.in)[length(data.in)] <- paste0("Wit_", qoi.t, j)
      data.in[is.na(data.in[, length(data.in)]), 
              length(data.in)] <- 0 #replace NAs with zeroes
    }

    data.in[, paste0("dit_", qoi.t)] <- getDits(matched_sets = sets.att,
                                                data = data.in)
    colnames(data.in)[length(data.in)] <- paste0("dits_", qoi.t)
    data.in[, paste0("Wit_", qoi.t, "-1")] <- 0

  }
  if (identical(qoi.in, "atc") || identical(qoi.in, "ate"))
  {

    for (j in lead)
    {
      
      dense.wits <- getWits(lead = j, data = data.in, matched_sets = sets.atc)
      data.in = merge(x = data.in, y = dense.wits, all.x = TRUE,
                   by.x = colnames(data.in)[1:2], by.y = c("id", "t"))
      colnames(data.in)[length(data.in)] <- paste0("Wit_atc", j)
      data.in[is.na(data.in[, length(data.in)]), 
              length(data.in)] <- 0 #replace NAs with zeroes
    }

    data.in$dit_atc <- getDits(matched_sets = sets.atc, data = data.in)
    colnames(data.in)[length(data.in)] <- "dits_atc"
    data.in$`Wit_atc-1` <- 0

  }
  data.in[, dependent.variable][is.na(data.in[, dependent.variable])] <- 0 
  #replace the NAs with zeroes.
  #This is ok because the dits should always be zero for corresponding entries
  # so the value is irrelevant.
  return(data.in)

}


