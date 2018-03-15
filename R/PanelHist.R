#' Illustration of PanelHist
#' 
#' \code{PanelHist} visualizes the distribution of the number of controls in each matched set
#' 
#' @author Erik Wang, \email{haixiaow@@princeton.edu}
#' 

#' @param matched_sets Name of the object returned from \code{\link{PanelMatch.R}}
#' @param qoi Quantity of interest, must be either "att" or "atc"
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
#'     # 100 means t0
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
#'       # truth:
#'       y[t, i] <- rho_1*y[t-1, i] + beta*treat[t,i] + lagTreOutc*treat[t-1,i] + beta_x*x[t,i] + beta_x2*x2[t,i] + alphai[i] + gammat[t] + eps[t, i] # the current period
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
#'   #y2.vec <- c(y2)
#'   ## generate unit and (Time + 50) index
#'   unit.index <- rep(1:N, each = (Time + 50))
#'   time.index <- rep(1:(Time + 50), N)
#'   # Data.str <- as.data.frame(cbind(y.vec, treat.vec, unit.index, x1.vec, x2.vec))
#'   # colnames(Data.str) <- c("y", "tr", "strata.id", "x1", "x2")
#'   Data.obs <- as.data.frame(cbind(time.index, unit.index, y.vec, y.lagged.vec,
#'                                   treat.vec, treat.lagged.vec, x.vec))
#'   colnames(Data.obs) <- c("time", "unit", "y", "y_l1", "treat", "treat_l1", "x")
#'   Data.obs <- Data.obs[which(Data.obs$time > 50), ]
#'   return(Data.obs)
#' }
#' Matches_Synth <- PanelMatch(lag = 4, max.lead = 0, time.id = "time",
#'                                             unit.id = "unit",
#'                                             treatment = "treat",
#'                                             formula = y~treat + x,
#'                                             method = "Synth",
#'                                             qoi = "ate", 
#'                                             data = Data.obs)
#'                                  
#' PanelMatch::PanelHist(Matches_Synth, qoi = "att",
#'                       main = "Distribution of the number of control units",
#'                       xlab = "")
#' 
#' PanelHist(Matches_Synth, qoi = "atc",
#'                       main = "Distribution of the number of control units",
#'                       xlab = "")
#' @export
PanelHist <- function(matched_sets, qoi = "att",...) {
  if (qoi == "att") {
    hist(unlist(matched_sets$NC_ATT), ...)
  } else if (qoi == "atc") {
    hist(unlist(matched_sets$NC_ATC), ...)
  } else {
    stop("please select either att or atc for qoi")
  }
  
}

