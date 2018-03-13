#' Illustration of DisplayTreatment
#' 
#' \code{DisplayTreatment} visualizes the treatment distribution across unit and time in a panel dataset
#' 
#' @author Erik Wang, \email{haixiaow@@princeton.edu}
#' 
#' 
#' 
#' 
#' @param unit.id A numeric vector of unit identifiers
#' @param time.id A numeric vector of time identifiers
#' @param treatment Name of the treatment variable in character class
#' @param data The data frame in data.frame class
#' @param color.of.treated Color of the treated observations in character. Default is red.
#' @param color.of.untreated Color of the untreated observations in character. Defalt is blue.
#' @param title Title of the plot in character
#' @param legend.position Position of the legend with the same choice set as in ggplot2
#' @param xlab Character label of the x-axis
#' @param ylab Character label of the y-axis
#' @param x.size Numeric size of the text for xlab. Default is 10.
#' @param x.size Numeric size of the text for ylab. Default is 5. 
#' @param letend.labels Character vector of length two describing the labels of the legend to be shown in the plot
#' @examples 
#' set.seed(123)
#' alphai <- rnorm(5000, mean = 0, sd = 1)
#' gammat <- rnorm(100, mean = 0, sd = 1)
#' 
#' sim_data <- function (N = 100, Time = 20, lag.one = 4, lag.two = 6,
#'                       lead = 0, M = 1,
#'                       rho_1 = .6, rho_t_1 = .6, rho_tt_1 = .4, 
#'                       rho_x = .4, rho_x2 = 0, lagTreOutc = .4, 
#'                       beta = 1, beta_x = .4, beta_x2 = 0, ephi = 0,
#'                       rho_2 = .3, hetereo = T,phi = .5, interc = -4,
#'                       ITER = 500) {
#'   y <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   x <- matrix(rnorm(N*(Time + 50), 0.5,1), ncol=N)
#'   x2 <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   x1.t0 <- rnorm(1, 0.5, 1)
#'   x2.t0 <- rnorm(1, 0.5, 1)
#'   for (i in 1:N) {
#'     x2[1, i] <- phi * x2.t0 + rnorm(1, mean = 0.5, sd = 1)
#'     for (t in 2:(Time + 50)){
#'       x2[t, i] <- phi*x2[t-1, i] + rnorm(1, mean = 0.5, sd = 1)
#'     }
#'   }
#'   treat <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   y.lagged <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   treat.lagged <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   
#'   eps <- matrix(NA, ncol = N, nrow = (Time + 50))
#'   
#'   for (i in 1:N) {
#'     #' 100 means t0
#'     y.lagged[1, i] <- 0
#'     treat.lagged[1,i] <- 0
#'     treat.error <- interc
#'     prob <- exp(rho_t_1*y.lagged[1,i] + alphai[i] + rho_tt_1*treat.lagged[1,i] + rho_x*x[1,i] + rho_x2*x2[1,i] + gammat[1] + treat.error)/
#'       (1+exp(rho_t_1*y.lagged[1,i] + alphai[i] + rho_tt_1*treat.lagged[1,i] + rho_x*x[1,i] + rho_x2*x2[1,i] + gammat[1] + treat.error))
#'     treat[1,i] <- rbinom(1,1, prob)
#'     #' error for the first period
#'     if(hetereo == T){
#'       eps[1, i] <- rnorm(1, 0, sd = runif(1, 0, 0.5 + treat[1,i]))
#'     } else {
#'       eps[1, i] <- rnorm(1, 0, 1)
#'     }
#'     
#'     
#'     y[1,i] <- rho_1*y.lagged[1,i] + lagTreOutc * treat.lagged[1,i] + 
#'       beta*treat[1,i] + beta_x*x[1,i] + beta_x2*x2[1,i] +
#'       alphai[i] + gammat[1] + eps[1, i]
#'     
#'     for (t in 2:(Time + 50)) {
#'       prob <- exp(rho_t_1*y[t-1,i] + alphai[i] + rho_tt_1*treat[t-1,i] + rho_x*x[t,i] + rho_x2*x2[t,i] +gammat[t] + treat.error)/
#'         (1+exp(rho_t_1*y[t-1,i] + alphai[i] + rho_tt_1*treat[t-1,i] + rho_x*x[t,i] + rho_x2*x2[t,i] + gammat[t] + treat.error))
#'       treat[t,i] <- rbinom(1,1, prob)
#'       treat.lagged[t,i] <- treat[t-1,i]
#'       
#'       if(hetereo == T) {
#'         eps[t, i] <- ephi*eps[t-1, i] + rnorm(1, 0, sd = runif(1, 0, .5 + treat[t, i]))
#'       } else {
#'         eps[t, i] <- ephi*eps[t-1, i] + rnorm(n = 1, mean = 0, sd = 1)
#'       }
#'       
#'       #' truth:
#'       y[t, i] <- rho_1*y[t-1, i] + beta*treat[t,i] + lagTreOutc*treat[t-1,i] + beta_x*x[t,i] + beta_x2*x2[t,i] + alphai[i] + gammat[t] + eps[t, i] #' the current period
#'       
#'       y.lagged[t,i] <- y[t-1,i]
#'     }
#'     
#'   }
#'   y.vec <- c(y)
#'   y.lagged.vec <- c(y.lagged)
#'   treat.vec <- c(treat)
#'   treat.lagged.vec <- c(treat.lagged)
#'   x.vec <- c(x)
#'   #'y2.vec <- c(y2)
#'   #'#' generate unit and (Time + 50) index
#'   unit.index <- rep(1:N, each = (Time + 50))
#'   time.index <- rep(1:(Time + 50), N)
#'   #' Data.str <- as.data.frame(cbind(y.vec, treat.vec, unit.index, x1.vec, x2.vec))
#'   #' colnames(Data.str) <- c("y", "tr", "strata.id", "x1", "x2")
#'   Data.obs <- as.data.frame(cbind(time.index, unit.index, y.vec, y.lagged.vec,
#'                                   treat.vec, treat.lagged.vec, x.vec))
#'   colnames(Data.obs) <- c("time", "unit", "y", "y_l1", "treat", "treat_l1", "x")
#'   Data.obs <- Data.obs[which(Data.obs$time > 50), ]
#'   return(Data.obs)
#' }
#' 
#' library(ggplot2)
#' out <- sim_data(N = 100, T = 20, 
#'                 rho_1 = 0.6, rho_t_1 = 0.6, 
#'                 ephi = 0,  hetereo = F)
#' DisplayTreatment(unit.id = "unit",
#'                  time.id = "time",
#'                  treatment = "treat",
#'                  data = out)
#'                  
#' @export
DisplayTreatment <- function(unit.id, time.id, treatment, data, 
                        color.of.treated = "red",
                        color.of.untreated = "blue", 
                        title = "Treatment Distribution \n Across Units and Time",
                        legend.position= "none", xlab = "time", ylab = "unit",
                        x.size = 10, y.size = 5,
                        legend.labels = c("not treated", "treated"))

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
                            breaks = c(0,1), labels = legend.labels) +
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
