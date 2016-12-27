syn_DID <- function(L, F, time.id = "year", 
                    unit.id = "ccode",
                    treatment, covariate, dependent) {

L <- L # set past history
F <- F # set the future

varnames <- c(time.id, unit.id, treatment, covariate, dependent)

# a function to delete null/empty entries in a list
delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(lapply(x.list, length) != 0)]
}

# a function to make every element in a list a data.frame while deleting duplicates
dframelist.rb_dup <- function (x) {
  test <- rbindlist(lapply(x, as.data.frame))
  return(test[!duplicated(test), ])
}

cscwdid <- function(x) {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  dataprep.out <- dataprep(foo = x, 
                           dependent = "V5",
                           unit.variable = "V2",
                           # unit.names.variable = "unit.name",
                           time.variable = "V1",
                           treatment.identifier = testid[2], # the regionno of the treated unit
                           controls.identifier = testid[-2],
                           time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                           time.predictors.prior = min(timeid):max(timeid-F-1),
                           predictors = "V4")
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
  # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
  # then making it a data.frame so as to merge it in the next step
  weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  x <- merge(x, weights, by = "V2") # merge it with the data.frame (
  first.diff <- x[x$V2 == testid[2], ]$V5[L+1+F] - x[x$V2 == testid[2], ]$V5[L]
  second.diff.1 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$V5*
                         x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L+1+F], ]$w.weight)
  second.diff.2 <- sum(x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$V5*
                         x[x$V2 %in% testid[-2] & x$V1 == unique(x$V1)[L], ]$w.weight)
  second.diff <- second.diff.1 - second.diff.2                    
  return(first.diff - second.diff)
}

### writing a function to apply synthetic control and generate a weight vector of the same length as the dataset
callSynth <- function (x, unit.id, time.id) {
  testid <- unique(x$V2)
  timeid <- unique(x$V1)
  dataprep.out <- dataprep(foo = x, 
                           dependent = "V5",
                           unit.variable = "V2",
                           # unit.names.variable = "unit.name",
                           time.variable = "V1",
                           treatment.identifier = testid[2], # the regionno of the treated unit
                           controls.identifier = testid[-2],
                           time.optimize.ssr = min(timeid):max(timeid-F-1), # the pre-treatment preiod
                           time.predictors.prior = min(timeid):max(timeid-F-1),
                           predictors = "V4")
  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS") # calibrate the weights
  # use rbind and cbind to add to the weights a weight vector of 1s for the treated observation,
  # then making it a data.frame so as to merge it in the next step
  weights <- as.data.frame(rbind(cbind(synth.out$solution.w, testid[-2]), cbind(w.weight = 1, testid[2])))
  colnames(weights)[2] <- "V2" # give the critical column the unit.name, so can merge
  merged <- merge(x, weights, by = "V2") # merge it with the data.frame (smaller data.frame as a list element)
  merged$w.weight <- ifelse(merged$V1 == max(timeid) & merged$V2 == testid[2], 1, 
                            ifelse(merged$V1 == max(timeid) - F - 1 & merged$V2 == testid[2],1, 
                                   ifelse(merged$V1 == max(timeid) & merged$V2 %in% testid[-2], merged$w.weight, 
                                          ifelse(merged$V1 == max(timeid) - F - 1 & merged$V2 %in% testid[-2], -merged$w.weight, 0) 
                                   )))
  new.W <- data.frame(unit.id, time.id) # create a new data.frame as large as the dataset
  names(new.W)[1:2] <- c("V2", "V1") # assigning names so as to merge it next
  total2 <- merge(new.W, merged, by = c("V2", "V1"), all= T) # merge, so now we have a full data.frame with weights
  total2$w.weight <- ifelse(is.na(total2$w.weight), 0, total2$w.weight) # turn NAs into zero
  return (total2$w.weight) # return the weight variable
} 

# subsetting the data.frame to include only relevant variables
d2 <- na.omit(d[varnames])
d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric

### from zero to 1 ###
# as.matrix it so that it can work with the cpp function
dmatrix <- as.matrix(d2)
# dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]

### finding matches using the cpp function ###
biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
### finding matches using the cpp function ###

### cleaning the output from cpp ###
# delete both higher level and lower level null entries
smallerlist <- lapply(Filter(function (x) !is.null(x), biglist), delete.NULLs) 
# further cleaning
smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
# use function dframelist.rb_dup to turn every list element into a data.frame
smallerlist <- lapply(smallerlist, dframelist.rb_dup)
# subset out any dataframe that have 2 or fewer than 2 units
smallerlist <- Filter(function (x) nrow(x) > 2*(L+F+1), smallerlist)

# only focus on ATT
even_smaller1 <- Filter(function (x) x[x$V2 == unique(x$V2)[2] & x$V1 == unique(x$V1)[1], ]$V3 == 0, smallerlist)

# brute forcing qoi
all.diffs.weighted <- lapply(even_smaller1, cscwdid)
qoi <- Reduce("+", all.diffs.weighted)/length(all.diffs.weighted)
return(qoi)
}
