syn_DID_gapsplot <- function(L, F, time.id = "year", 
                             xlab = "Pre-treatment periods",
                             ylab = "Gaps between real and synthetic",
                             unit.id = "ccode", 
                             legend.position = "right",
                             treatment, covariate, dependent, d) {
  
  FORWARD <- F
  
  
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
  # execute cscwplot
  plot.materials <- delete.NULLs(lapply(even_smaller1, cscwplot, L = L, FORWARD = FORWARD))

  df <- data.frame(x=rep(-L:-1, length(plot.materials)), 
                   val=unlist(lapply(plot.materials, function(x) x$gap)), 
                   variable=rep(paste0("unit", 
                                       unlist(lapply(plot.materials, function(x) x$unit.id))), 
                                each=L))

  t <- tapply(df$val, df$x, quantile)
  
 # aggregate(val ~ x, data = df, FUN = quantile)
  y2 <- apply(do.call(rbind,t),2,  FUN = unlist)[,1]
  for (i in 2:5) {
    y2 <- c(y2, apply(do.call(rbind,t),2,  FUN = unlist)[,i])
  }
  df2 <- data.frame(x2=rep(-L:-1,5), 
                    y2=y2,
                    g2=gl(5,(L), 
                    labels=c("Min", "1st Qu.", "Median",
                             "3rd Qu.", "Max.")))

  ggplot() + 
    geom_line(aes(x=factor(x), y=val, group = variable,
                  colour = variable), df) +
    labs(x = xlab, y = ylab) + 
    scale_x_discrete(breaks = -L:-1,
                     labels=paste("t", -L:-1, sep = "")) +
    geom_line(aes(x=factor(x2), y=y2, group = g2), 
              colour = "black",
              df2) + theme(legend.position=legend.position) 
  
 
}

# 
# 
# a <- tapply(df$val, df$x, quantile)
# df2 <- data.frame(x2=rep(1:L,5), y2=unlist(a))
# 
# d1 <- data.frame(x1=rep(1:10,3), y1=rnorm(10*3), g1=gl(3,10,labels=letters[1:3]))
# d2 <- data.frame(x2=rep(1:10,3), y2=rnorm(10*3), g2=gl(3,10, labels=letters[4:6]))
# 
# ggplot() + 
#   geom_line(aes(x1, y1, colour=g1), d1) +  
#   geom_line(aes(x2, y2, colour=g2), d2)


