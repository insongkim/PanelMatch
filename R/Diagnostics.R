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
                             xlab = "Treatment periods",
                             ylab = "Gaps between real and synthetic",
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
