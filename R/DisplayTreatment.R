
DisplayTreatment <- function(unit.id, time.id, treatment, data, 
                        color.of.treated = "red",
                        color.of.untreated = "blue", 
                        title = "Treatment Distribution \n Across Units and Time",
                        legend.position= "none", xlab = "time", ylab = "unit",
                        x.size = 10, y.size = 5)

{
  # load the dataframe that the user specifies
  data <- na.omit(data[c(unit.id, time.id, treatment)])
  
  # rename variables to match with the object names in the loop below
  colnames(data) <- c("unit.id", "time.id", "treatment")
  
  # make unit.id a character: this is useful when the unit the user
  # passes to the function is numeric (e.g. dyad id)
  # data$unit.id <- as.character(data$unit.id)
  
  ## Sorting units by treatment intensity 
 
  data$trintens <- tapply(data$treatment, data$unit.id, mean, na.rm = T)[as.character(data$unit.id)]
  data <- data[order(data$trintens), ]
  data$old.index <- data$unit.id
  data$unit.id <- match(data$unit.id, unique(data$unit.id) ) 
  data$unit.id <- factor(data$unit.id, levels=unique(as.character(data$unit.id)))
  
      # use ggplot2
      p <- ggplot(data, aes(unit.id, time.id)) + geom_tile(aes(fill = treatment),
                                                           colour = "white") +
        scale_fill_gradient(low = color.of.untreated,
                            high = color.of.treated, guide = "legend", 
                            breaks = c(0,1), labels = c("not treated", "treated")) +
        theme_bw() +
        labs(list(title = title, x = ylab, y = xlab, fill = "")) +
        theme(axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.position = legend.position,
              panel.background = element_blank(), 
              axis.text.x = element_text(angle=90, size = x.size, vjust=0.5),
              axis.text.y = element_text(size = y.size),
              plot.title = element_text(hjust = 0.5)) +
        scale_x_discrete(expand = c(0, 0), labels = unique(as.character(data$old.index))) +
        scale_y_continuous(limits = c(min(data$time.id)-1,max(data$time.id) + 1), expand = c(0,0)) +
        coord_flip()
        return(p) # return the plot
}
