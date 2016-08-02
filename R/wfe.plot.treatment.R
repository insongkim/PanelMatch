### write the function
wfe.plot.treatment <- function(unit, time, treatment, graph.label = NULL, ylab = NULL) {
  # unit, time and treatment are variables taken directly from the dataset
  # graph.label is simply a character string 
  # ylab is another character string
  unit.names <- as.character(unique(unit)) # getting names of units to be "pasted" into the 
  # X-axis
  times <- unique(time)      # unique time values
  n.unit <- length(unit.names) # number of N
  n.times <- length(times) # length of T
  d <- na.omit(as.data.frame(cbind(unit, time, treatment))) # this deals with NAs
  d$time <- as.numeric(as.character(d$time)) # convert factor into numeric
  d$treatment <- as.numeric(as.character(d$treatment)) # convert factor into numeric
  
  par(mar = c(4, 4, .8, 0))
  plot(1, 1, type = "n", axes = FALSE, xlim = c(-10, n.unit  -10), ylim = c(.5, n.times + .5),
       xlab = "", ylab = ylab, main = graph.label)
  axis(1, at = -10:(n.unit -10), labels = FALSE)
  if (n.unit <= 200){
    mtext(text = unit.names, side = 1, at = -9.5:(n.unit-10.5), las = 2, line = 1, 
          cex = sqrt(sqrt(174/n.unit))*.3)
  }
  begin.time <- min(times)
  interval <- max(times) - begin.time
  ticks <- c(begin.time, ceiling(begin.time + interval/5), ceiling(begin.time + 2*interval/5), 
             ceiling(begin.time + 3*interval/5), ceiling(begin.time + 4*interval/5), max(times))
  axis(2, at = ticks - (begin.time-1), labels = ticks, line = 1)
  polygon.color = c("lightgrey", "dimgrey")
  
  
  for (i in 1:n.unit) {
    for (t in 1:n.times) {
      if (length(d$treatment[d$unit == unit.names[i] &
                             d$time == times[t]]) > 0) {
        temp <- d$treatment[d$unit == unit.names[i] &
                              d$time == times[t]]
        if (is.na(temp)) next
        polygon(y = c(.5+t-1, .5+t, .5+t, .5+t-1), 
                x = c(.5+i-11.5, .5+i-11.5, .5+i-10.5, .5+i-10.5), density = NA, border = "white",
                col = polygon.color[temp + 1])
      } else {
        next
      }
    }
  }
  
}
