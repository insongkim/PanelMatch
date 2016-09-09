####### Here begins coding up the estimator for "Synthetic Control for Multiple
####### Treated Units", the first approach in the TSCS paper updated on Aug 4

findMatched2 <- function(unit.id, time.id, treatment, covariate, unit.name, dependent, data) {
  # a list of lists to contain matched matrix for each treated observation
  master_list <- list()
  # a function to delete NULLs from a list
  delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
  }
  
  # load the dataframe that the user specifies
  data <- na.omit(data[c(unit.id, time.id, treatment, covariate, unit.name, dependent)]) # omit NAs 
  colnames(data) <- c("unit.id", "time.id", "treatment", "covariate", "unit.name", "dependent") # rename variables to match with the object names in the loop below
  data <- data[order(data$unit.id, data$time.id), ] # order by unit and time. This is important as the loop below works with 
  # this order
  
  print("finding a matched set of control units for each treated observation...")
  
  # the loop
  for (i in 1:nrow(data)) {
    if (data$treatment[i] == 1) { # if the observation gets treated 
      # record the time and unit of the treated observation
      treated.time <- data$time.id[i]
      treated.unit <- data$unit.id[i]
      # creating treated data and untreated data based on what's given in the math
      untreated.data <- data[data$time == treated.time & data$unit.id != treated.unit, ]
      treated.data <- data[data$time == treated.time & data$unit.id == treated.unit, ]
      id <- unique(untreated.data$unit.id) # extract unit ids of units that are untreated at time t
      result_list <- list() # create a list to store data that will combine both treated and untreated data
      for (k in 1:length(id)) {
        result_list[[k]] <- rbind(untreated.data[untreated.data$unit.id == id[k] & untreated.data$treatment ==0, ], treated.data)
      }
      master_list[[i]] <- result_list[unlist(lapply(result_list, nrow) == 2)] # save the entire result_list to the master_list, indexed by i
      # here I use function "delete.NULLs" to delete all the empty elements
      # in the result_list
    } else {
      next
    }
  }
  
  newlist <- delete.NULLs(master_list)

  print("computing regression weights...")
  
  # create a data.frame of weights that takes directly the unit.id and time.id from "data"
  new.W <- data.frame(data$unit.id, data$time.id)
  names(new.W)[1:2] <- c("unit.id", "time.id")
  
  # initialize weights of zeroes. This variable is a vector of length N*T
  new.W$weights <- 0  
  
  ## loop through 
  for (h in 1:length(newlist)) {
    test <- rbindlist(newlist[[h]]) # rbindlist to turn a list into a dataset
    # since each element in test has the same treated observation, there 
    # are duplicates
    test <- test[!duplicated(test), ] 
    # extract all unique unit.id s
    testid <- unique(test$unit.id)
    
    # calling package "synth" to create weights using synthetic control
    dataprep.out <- dataprep(foo = test, 
                             dependent = "dependent",
                             unit.variable = "unit.id",
                             unit.names.variable = "unit.name",
                             time.variable = "time.id",
                             treatment.identifier = testid[2], # the id of the treated unit
                             controls.identifier = testid[-2], # the ids of all control units
                             time.optimize.ssr = mean(test$time.id), # the pre-treatment preiod
                             time.predictors.prior = mean(test$time.id),
                             predictors = "covariate")
    
    # extract weights: "solution.w" from "synth" output is the weights. 
    # we use cbind and rbind to combine these weights with the treated observation, which receives a weight of 1
    # and a time variable that takes the same value across all observations in this small dataset.
    # the value is going to be the year of the treated observation.
    weights <- cbind(rbind(as.data.frame(cbind(synth(data.prep.obj = dataprep.out, method = "BFGS")$solution.w, testid[-2])), 
                           c(1, testid[2])), 
                     "time" = mean(test$time.id))
    names(weights)[2] <- "unit.id" # change the name to "unit.id" just for later
                                   # merging purposes
    
    # importantly, we merge the smaller dataset "weights" with new.W
    # so that the merged dataset has the nrow of N*T. This will make
    # summing up the weights easier.
    total2 <- merge(new.W, weights, by = "unit.id", all.x = TRUE, all.y = TRUE)
    total2 <- total2[order(total2$unit.id, total2$time.id), ] # sort it by unit and time
    total2$w.weight[is.na(total2$w.weight)] <- 0 # set all NAs to zero in w.weight
    
    # most important step: for all observations with the same time identifier as 
    # the one for the treated observations, add the values in variable w.weights from 
    # the merged dataset (created in the previous step) to their original values in new.W$weights
    # for those with different time identifiers, let their weights remain what they are
    # THIS WAY WE CONSTANTLY UPDATE THE WEIGHTS IN new.W AS WE LOOP THROUGH THE newlist
    new.W$weights <- ifelse(new.W$time.id == mean(test$time.id), new.W$weights + total2$w.weight, new.W$weights)
  }
  
  # add weights to data
  data$weights <- new.W$weights 
  # return data
  return(data) # so that the returned object from this function will be 
               # a ready-to-use dataframe with a variable storing all the weights
}
