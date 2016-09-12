####### Here begins coding up the estimator for "Synthetic Control for Multiple
####### Treated Units", the first approach in the TSCS paper updated on Aug 4

findMatched2_vectorized <- function(unit.id, time.id, treatment, covariate, dependent, unit.name, data) {
  # a list of lists to contain matched matrix for each treated observation
  master_list <- list()
  # a function to delete NULLs from a list
  delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
  }
  
  # load the dataframe that the user specifies
  data <- na.omit(data[c(unit.id, time.id, treatment, covariate, dependent, unit.name)]) # omit NAs 
  
  llist <- list()
  for (i in 1:length(covariate)) {
    llist[i] <- paste("covariate",i, sep = "")
  }
  
  colnames(data) <- c("unit.id", "time.id", "treatment", llist, "dependent", "unit.name") # rename variables to match with the object names in the loop below
  
  data$unit.name <- as.character(data$unit.name)
  data[1:(length(data)-1)] <- lapply(data[1:(length(data)-1)], function(x) as.numeric(as.character(x)))
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
  
test.list <- lapply(lapply(newlist, function(y) rbindlist(y)), function(x) x[!duplicated(x), ])
testid.list <- lapply(test.list, function(x) unique(x$unit.id))


weights.list <- apply(mapply(function(x, y) {
  dataprep(foo = x,
           dependent = "dependent",
           unit.variable = "unit.id",
           unit.names.variable = "unit.name",
           time.variable = "time.id",
           treatment.identifier = y[2], # the regionno of the treated unit
           controls.identifier = y[-2],
           time.optimize.ssr = mean(x$time.id), # the pre-treatment preiod
           time.predictors.prior = mean(x$time.id),
           predictors = llist)
}, test.list, testid.list), 
                      2, function (x) synth(x, method = "BFGS"))

list_of_merged_weights <- mapply(function(x, y, z)
  cbind(rbind(as.data.frame(cbind(x$solution.w, y[-2])), 
              c(1, y[2])), 
        "time" = mean(z$time.id)), weights.list, testid.list, test.list)


temp2 <- lapply(lapply(apply(list_of_merged_weights, 2, function(x) rbind(x)),
                       function(x) as.data.frame(cbind(as.data.frame(x[[1]]),
                                                              as.data.frame(x[[2]]),
                                                              as.data.frame(x[[3]])))),
                setNames, nm = c("weights", "unit.id", "time.id"))

merged <- lapply(temp2, function (x)
  merge(new.W, x, by = "unit.id", all.x = TRUE, all.y = TRUE))

merged2 <- lapply(merged, 
                  function(x) x[order(x$unit.id, x$time.id.x), ]) # sort it by unit and time


merged3 <- lapply(merged2, function(x) {
  x$weights[is.na(x$weights)] <- 0
  return(x)
})

##############

merged4 <- lapply(merged3, function(x){
  x$time.id.y <- mean(x$time.id.y, na.rm = TRUE)
  return(x)
})

merged5 <- lapply(merged4, function(x) {
  x$weights <- ifelse(x$time.id.x == mean(x$time.id.y), x$weights, 0)
  return(x)
})


data$weights <- Reduce("+", lapply(merged5, function(x) return(x$weights)))

  # return data
  return(data) # so that the returned object from this function will be 
  # a ready-to-use dataframe with a variable storing all the weights
}
