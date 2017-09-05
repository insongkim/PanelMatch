options(digits=6)
# clean slate
# rm(list=ls())
# library(devtools)
# install_github("insongkim/wfe", ref = "new_branch")
# load packages
pkg <- c("foreign","ggplot2")
lapply(pkg, require, character.only = TRUE)

alphai <- rnorm(n =1000, mean = 10, sd = 6)
gammat <- rnorm(n = 20, mean = 10, sd = 6)
sim_wfe2_data <- function (N = 100, Time = 20, lag.one = 4, lag.two = 6,
                           lead = 0,
                           rho_1 = .6, rho_t_1 = .6, rho_tt_1 = .2, 
                           rho_x = .4, rho_x2 = 0, lagTreOutc = .2, 
                           beta = 1, beta_x = .2, beta_x2 = 0, 
                           phi = .3, rho_t_2 = .3, ephi = .5,
                           rho_2 = .3, M = 1, hetereo = T,
                           x_fe = 0, frac = 1.3,
                           ITER = 500) {
  y <- matrix(NA, ncol = N, nrow = Time)
  eps <- matrix(NA, ncol = N, nrow = Time)
  treat <- matrix(NA, ncol = N, nrow = Time)
  y.lagged <- matrix(NA, ncol = N, nrow = Time)
  treat.lagged <- matrix(NA, ncol = N, nrow = Time)
  
  
  x <- matrix(rep(NA, N*Time), ncol = N)
  for (i in 1:N) {
    x[1, i] <- rnorm(1, 3, 1) + x_fe*gammat[1] + x_fe*alphai[i]
    for(t in 2:Time){
      x[t, i] <- phi*x[t-1, i] + rnorm(1, 3, 1) + x_fe*gammat[t] + x_fe*alphai[i]
    }
    
  }
  
  x2 <- matrix(rep(NA, N*Time), ncol = N)
  for (i in 1:N) {
    x2[1, i] <- rnorm(1, 0.5, 1) + gammat[1] + alphai[i]
    for(t in 2:Time){
      x2[t, i] <- phi*x2[t-1, i] + rnorm(1, 0.5, 1) + gammat[t] + alphai[i]
    }
    
  }
  
  for (i in 1:N) {
    y.lagged[1, i] <- rnorm(1) + alphai[i] + gammat[1]
    treat.lagged[1,i] <- rbinom(1,1,exp(rho_x*x[1,i] + rho_x2*x2[1,i] + alphai[i] + gammat[1])/(1+exp(rho_x*x[1,i] + rho_x2*x2[1,i] + alphai[i] + gammat[1])))
    # think about commenting out treat.error
    # treat.error <- -4 #rnorm(1,-4)
    prob <- exp(rho_t_1*y.lagged[1,i] + alphai[i] + rho_tt_1*treat.lagged[1,i] + rho_x*x[1,i] + rho_x2*x2[1,i] + gammat[1])/
      (1+exp(rho_t_1*y.lagged[1,i] + alphai[i] + rho_tt_1*treat.lagged[1,i] + rho_x*x[1,i] + rho_x2*x2[1,i] + gammat[1]))
    treat[1,i] <- rbinom(1,1, prob/frac)
    eps[1, i] <- rnorm(1, 0, 1)
    y[1,i] <- rho_1*y.lagged[1,i] + alphai[i] + gammat[1] + 
      beta*treat[1,i] + lagTreOutc*treat.lagged[1,i] + beta_x*x[1,i] + beta_x2*x2[1,i] +
      eps[1,i]
    
    for (t in 2:Time) {
      # treat.error <- -4 #rnorm(1,-4)
      prob <- exp(rho_t_1*y[t-1,i] + alphai[i] + rho_tt_1*treat[t-1,i] + rho_x*x[t,i] + rho_x2*x2[t,i] +gammat[t])/
        (1+exp(rho_t_1*y[t-1,i] + alphai[i] + rho_tt_1*treat[t-1,i] + rho_x*x[t,i] + rho_x2*x2[t,i] + gammat[t]))
      treat[t,i] <- rbinom(1,1, prob/frac)
      treat.lagged[t,i] <- treat[t-1,i]
      
      if(hetereo == T) {
        eps[t, i] <- ephi*eps[t-1, i] + rnorm(n = 1, mean = 0, sd = runif(1, 0.5, 1.5))
      } else {
        eps[t, i] <- ephi*eps[t-1, i] + rnorm(n = 1, mean = 0, sd = 1)
      }
      
      # truth:
      y[t, i] <- rho_1*y[t-1, i] + beta*treat[t,i] + lagTreOutc*treat[t-1,i] + beta_x*x[t,i] + beta_x2*x2[t,i] + 
        alphai[i] + gammat[t] + eps[t, i] # the current period
      
      y.lagged[t,i] <- y[t-1,i]
    }
    
  }
  y.vec <- c(y)
  y.lagged.vec <- c(y.lagged)
  treat.vec <- c(treat)
  treat.lagged.vec <- c(treat.lagged)
  x.vec <- c(x)
  #y2.vec <- c(y2)
  ## generate unit and time index
  unit.index <- rep(1:N, each = Time)
  time.index <- rep(1:Time, N)
  # Data.str <- as.data.frame(cbind(y.vec, treat.vec, unit.index, x1.vec, x2.vec))
  # colnames(Data.str) <- c("y", "tr", "strata.id", "x1", "x2")
  Data.obs <- as.data.frame(cbind(time.index, unit.index, y.vec, y.lagged.vec,
                                  treat.vec, treat.lagged.vec, x.vec))
  colnames(Data.obs) <- c("time", "unit", "y", "y_l1", "treat", "treat_l1", "x")
  
  return(Data.obs)
}

Data.obs <- sim_wfe2_data(N = 100, T = 20)

# Matches_CBPSL4_2 <- PanelMatch(lag = 4, max.lead = 5, time.id = "time", 
#                                unit.id = "unit",
#                                treatment = "treat", 
#                                formula = y~ treat,
#                                method = "CBPS", weighting = TRUE,
#                                qoi = "ate",  M = 3,
#                                data = Data.obs)
# 
# 
# result.ate <- PanelEstimate_tmp2(lead = 0:5, CI = .95, 
#                                  matched_sets = Matches_CBPSL4_2,
#                                  inference = "bootstrap",
#                                  ITER = 100)
# 
# PanelEffectPlot(result.ate)