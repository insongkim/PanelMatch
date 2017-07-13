syn_DID_MSPE <- function(L, F, time.id = "year", qoi = "ATE",
                             unit.id = "ccode",
                             treatment, covariate, dependent, d) {
  
  L <- L # set past history
  F <- F # set the future
  FORWARD <- F
  d <- d # set dataset
  
  varnames <- c(time.id, unit.id, treatment, covariate, dependent)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
  
  ### finding matches using the cpp function ###
  # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
  ### finding matches using the cpp function ###
  
  ### cleaning the output from cpp ###
  # delete both higher level and lower level null entries
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L, F, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  
  if (qoi == "ATT") {
    return(sapply(even_smaller1, cscMSPE, L = L, FORWARD = FORWARD))
  } else {
    ### from 1 to zero ###
    dmatrix[,3] <- ifelse(dmatrix[,3] == 1, 0, 1)
    
    ### finding matches using the cpp function ###
    # biglist <- findDDmatched2(L = L, F, dmatrix) # first argument is matched treatment history
    ### finding matches using the cpp function ###
    
    ### cleaning the output from cpp ###
    # delete both higher level and lower level null entries
    smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L = L, F, dmatrix)), delete.NULLs) 
    # further cleaning
    smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
    # use function dframelist.rb_dup to turn every list element into a data.frame
    even_smaller2 <- lapply(smallerlist, dframelist.rb_dup)
    
    return(list("MSPE for ATT" = sapply(even_smaller1, cscMSPE, L = L, FORWARD = FORWARD),
                "MSPE for ATC" = sapply(even_smaller2, cscMSPE, L = L, FORWARD = FORWARD)))
    
  }
  
}

syn_DID_gapsplot <- function(L, FORWARD, time.id = "year", 
                             xlab = "Time periods",
                             ylab = "Differences between the treated units and their synthetic control units",
                             unit.id = "ccode", 
                             legend.position = "right",
                             treatment, covariate, dependent, d,
                             show.covariate = FALSE,
                             post.treatment = FALSE, vline = TRUE,
                             colour = "red", linetype = "longdash") {
  
  
  varnames <- c(time.id, unit.id, treatment, covariate, dependent)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
  
  ### finding matches using the cpp function ###
  # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
  ### finding matches using the cpp function ###
  
  ### cleaning the output from cpp ###
  # delete both higher level and lower level null entries
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L, FORWARD, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  
  
  # execute cscwplot to obtain w.weight and unit.id
  plot.materials <- delete.NULLs(lapply(even_smaller1, cscwplot, L = L, FORWARD = FORWARD,
                                        show.covariate = show.covariate,
                                        post.treatment = post.treatment,
                                        covariate = covariate))
  
  if (post.treatment == FALSE) {
    df <- data.frame(x=rep(-L:-1, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L))
  } else {
    df <- data.frame(x=rep(-L:FORWARD, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L+1+FORWARD))
  }
  
  
  t <- tapply(df$val, df$x, quantile)
  
  # aggregate(val ~ x, data = df, FUN = quantile)
  y2 <- apply(do.call(rbind,t),2, FUN = unlist)[,1]
  for (i in 2:5) {
    y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  }
  
  if (post.treatment == FALSE) {
    df2 <- data.frame(x2=rep(-L:-1,5), 
                      y2=y2,
                      g2=gl(5,(L), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  } else {
    df2 <- data.frame(x2=rep(-L:FORWARD,5), 
                      y2=y2,
                      g2=gl(5,(L+1+FORWARD), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  }
  
  if (post.treatment == FALSE) {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:-1,
                       labels=paste("t", -L:-1, sep = "")) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) 
  } else {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:FORWARD,
                       labels=c(paste("t", -L:-1, sep = ""), paste("t", 0, sep = ""), 
                                paste("t+", 1:FORWARD, sep = ""))
      ) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) +
      if(vline == TRUE) {
        geom_vline(xintercept = L+1, linetype = linetype, 
                   colour = colour)
      } else {
        geom_vline()
      }
    
  }
  
}

MSMD_gapsplot <- function(L, FORWARD, time.id = "year", 
                          xlab = "Time periods",
                          ylab = "Differences between the treated units and their synthetic control units",
                          unit.id = "ccode", M = 3,
                          legend.position = "right",
                          treatment, covariate, dependent, d,
                          show.covariate = FALSE,
                          post.treatment = FALSE, vline = TRUE,
                          colour = "red", linetype = "longdash") {
  
  
  varnames <- c(time.id, unit.id, treatment, covariate, dependent)
  
  # subsetting the data.frame to include only relevant variables
  d2 <- na.omit(d[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x))) # as.numeric
  
  ### from zero to 1 ###
  # as.matrix it so that it can work with the cpp function
  dmatrix <- as.matrix(d2)
  # dmatrix2 <- dmatrix[dmatrix[, 3] == 1, ]
  
  ### finding matches using the cpp function ###
  # biglist <- findDDmatched2(L, F, dmatrix) # first argument is matched treatment history
  ### finding matches using the cpp function ###
  
  ### cleaning the output from cpp ###
  # delete both higher level and lower level null entries
  smallerlist <- lapply(Filter(function (x) !is.null(x), findDDmatched2(L, FORWARD, dmatrix)), delete.NULLs) 
  # further cleaning
  smallerlist <- Filter(function (x) length(x) > 0, smallerlist)
  # use function dframelist.rb_dup to turn every list element into a data.frame
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  # execute cscwplot to obtain w.weight and unit.id
  plot.materials <- delete.NULLs(lapply(even_smaller1, MSMD_plot,
                                        L = L, FORWARD = FORWARD,
                                        M = M,
                                        show.covariate = show.covariate,
                                        post.treatment = post.treatment))
  
  if (post.treatment == FALSE) {
    df <- data.frame(x=rep(-L:-1, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L))
  } else {
    df <- data.frame(x=rep(-L:FORWARD, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L+1+FORWARD))
  }
  
  
  t <- tapply(df$val, df$x, quantile)
  
  # aggregate(val ~ x, data = df, FUN = quantile)
  y2 <- apply(do.call(rbind,t),2,  FUN = unlist)[,1]
  for (i in 2:5) {
    y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  }
  
  if (post.treatment == FALSE) {
    df2 <- data.frame(x2=rep(-L:-1,5), 
                      y2=y2,
                      g2=gl(5,(L), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  } else {
    df2 <- data.frame(x2=rep(-L:FORWARD,5), 
                      y2=y2,
                      g2=gl(5,(L+1+FORWARD), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  }
  
  if (post.treatment == FALSE) {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:-1,
                       labels=paste("t", -L:-1, sep = "")) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) 
  } else {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:FORWARD,
                       labels=c(paste("t", -L:-1, sep = ""), paste("t", 0, sep = ""), 
                                paste("t+", 1:FORWARD, sep = ""))
      ) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) +
      if(vline == TRUE) {
        geom_vline(xintercept = L+1, linetype = linetype, 
                   colour = colour)
      } else {
        geom_vline()
      }
    
  }
  
  
  
}

PS_gapsplot <- function(L, FORWARD, time.id = "year", 
                        xlab = "Time periods",
                        ylab = "Differences between the treated units and their synthetic control units",
                        unit.id = "ccode", M = 3,
                        legend.position = "right",
                        treatment, covariate, dependent, data,
                        show.covariate = FALSE,
                        post.treatment = FALSE, vline = TRUE,
                        colour = "red", linetype = "longdash") {
  
  
  data <- data[c(unit.id, time.id,  treatment, dependent, covariate)]
  dlist <- lapply(1:L, 
                  function (i) slide(data = data, Var = dependent, GroupVar = unit.id, TimeVar = time.id, slideBy = -(i),
                                     NewVar = paste("dependent_l", i, sep="")))
  data <- Reduce(function(x, y) {merge(x, y)}, dlist)
  
  # to include ldvs in varnames
  varnames <- c(time.id, unit.id, treatment, dependent, c(covariate, colnames(data)[(4 + length(covariate) + 1):length(data)]))
  
  d2 <- na.omit(data[varnames])
  d2[1:(length(d2))] <- lapply(d2[1:(length(d2))], function(x) as.numeric(as.character(x)))
  d2 <- d2[order(d2[,2], d2[,1]), ]
  dmatrix <- as.matrix(d2)
  smallerlist <- lapply(Filter(function(x) !is.null(x), findDDmatched2(L, 
                                                                       F = FORWARD, dmatrix)), delete.NULLs)
  smallerlist <- Filter(function(x) length(x) > 0, smallerlist)
  even_smaller1 <- lapply(smallerlist, dframelist.rb_dup)
  
  # take the forward periods from each subset:
  # IMPORTANT
  Fs <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 %in% unique(x$V1)[(L+2):(L+1+FORWARD)], ]
    return(x)
  })
  
  # to only include L and the first treatment period
  even_smaller1 <- lapply(even_smaller1, function(x) {
    x <- x[x$V1 %in% sort(unique(x$V1))[1:(L+1)], ]
    return(x)
  })
  
  # add varnames
  even_smaller1 <- lapply(even_smaller1, function (x) {
    colnames(x) <- varnames
    return(x)
  })
  
  # add varnames to the forward-periods-subset
  Fs <- lapply(Fs, function (x) {
    colnames(x) <- varnames
    return(x)
  })
  
  pooled <- rbindlist(even_smaller1) # get a dataset for propensity score generation
  
  # get propensity scores
  fit0 <- glm(reformulate(response = treatment, termlabels = c(covariate, names(pooled[, (4 + length(covariate) + 1):length(pooled)]))), 
              family = binomial(link = "logit"), data = pooled)
  pooled$ps <- fit0$fitted.values
  
  # aggregate to delete duplicates
  aggregated <- aggregate(reformulate(response = "ps", termlabels = c(unit.id, time.id)), pooled, mean)
  
  newlist <- lapply(even_smaller1, function(x) merge(aggregated, x, by = c(unit.id, time.id)))
  newlist <- Map(rbind.fill, newlist, Fs)
  # superlist <- mapply(function(x, y) rbind.fill(x,y), newlist, Fs)
  
  
  plot.materials <- delete.NULLs(lapply(newlist, PS_plot,
                                        L = L, FORWARD = FORWARD,
                                        M = M,
                                        show.covariate = show.covariate,
                                        post.treatment = post.treatment))
  
  if (post.treatment == FALSE) {
    df <- data.frame(x=rep(-L:-1, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L))
  } else {
    df <- data.frame(x=rep(-L:FORWARD, length(plot.materials)), 
                     val=unlist(lapply(plot.materials, function(x) x$gap)), 
                     variable=rep(paste0("unit", 
                                         unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                  each=L+1+FORWARD))
  }
  
  
  t <- tapply(df$val, df$x, quantile)
  
  # aggregate(val ~ x, data = df, FUN = quantile)
  y2 <- apply(do.call(rbind,t),2,  FUN = unlist)[,1]
  for (i in 2:5) {
    y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  }
  
  if (post.treatment == FALSE) {
    df2 <- data.frame(x2=rep(-L:-1,5), 
                      y2=y2,
                      g2=gl(5,(L), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  } else {
    df2 <- data.frame(x2=rep(-L:FORWARD,5), 
                      y2=y2,
                      g2=gl(5,(L+1+FORWARD), 
                            labels=c("Min", "1st Qu.", "Median",
                                     "3rd Qu.", "Max.")))
  }
  
  if (post.treatment == FALSE) {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:-1,
                       labels=paste("t", -L:-1, sep = "")) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) 
  } else {
    ggplot() + 
      geom_line(aes(x=factor(x), y=val, group = variable,
                    colour = variable), df) +
      labs(x = xlab, y = ylab) + 
      scale_x_discrete(breaks = -L:FORWARD,
                       labels=c(paste("t", -L:-1, sep = ""), paste("t", 0, sep = ""), 
                                paste("t+", 1:FORWARD, sep = ""))
      ) +
      geom_line(aes(x=factor(x2), y=y2, group = g2), 
                colour = "black",
                df2) + theme(legend.position=legend.position) +
      if(vline == TRUE) {
        geom_vline(xintercept = L+1, linetype = linetype, 
                   colour = colour)
      } else {
        geom_vline()
      }
    
  }
  
  
  
}


