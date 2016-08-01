

pwfe <- function (formula, treat = "treat.name", outcome, data, pscore = NULL,
                  unit.index, time.index = NULL, method = "unit", within.unit = TRUE,
                  qoi = c("ate", "att"), estimator = NULL, C.it = NULL,
                  White = TRUE, White.alpha = 0.05,
                  hetero.se = TRUE, auto.se = TRUE, unbiased.se = FALSE,
                  verbose = TRUE) {


  pwfe.call <- match.call()
  tn.row <- nrow(data) # total number of rows in data



### Warnings
  ## Warning for missing unit & time index
  if (missing(unit.index))
    stop("'unit.index' or index for strata should be provided")

  if (is.null(time.index) & method == "time")
    stop("'time.index' should be provided")

  ## Warning for methods
  if(method=="time" && !is.null(estimator) && estimator == "fd")
    stop("First Difference is not compatible with 'time' method")

  if(is.null(time.index) && !is.null(estimator) && estimator  == "fd")
    stop("First Difference cannot calculate when 'time.index' is missing")

  if(is.null(time.index) && !is.null(estimator) && estimator  == "did")
    stop("Difference-in-Differences is not compatible with pwfe")

  ## Warning for C.it
  if (!is.null(C.it)){
    Cit <- data[,C.it]
    if (!is.numeric(Cit) && length(Cit)!= tn.row)
      stop("'C.it' must be a numeric vector with length equal to number of observations")
    if ( sum(Cit < 0) > 0 )
      stop("'C.it' must be a non-negative numeric vector")
  }

  ## Waring for pscore (if propensity score is provided by researcher)
  p.score <- data[,pscore]
  if (!is.null(pscore) && !is.numeric(p.score) && length(p.score)!= tn.row)
    stop("'pscore' must be a numeric vector with length equal to number of observations")
  if (!is.null(pscore)) {
    if (!((0 < p.score ) && (p.score < 1)))
      stop("'pscore' must be a bounded away from zero and one")
  }
  if (!is.null(pscore) && !missing(formula))
    stop("'formula' should not be provided when pscore is specified")




  ## cat("warnings done:\n")
  ## C.it
  ## Default for ATE
  if (is.null(C.it)){
    data$C.it <- as.integer(rep(1, nrow(data)))
  }

  ## White.alpha
  if (is.null(White.alpha)){
    White.alpha <- 0.05
  } else {
    White.alpha <- White.alpha
  }

  ## Warning for binary treatment

  ## treatment variable
  data$TR <- as.numeric(data[,treat])

  if (length(unique(data$TR)) !=2 ){
    stop("'treat' must be a binary vector")
  }
  if (sum(unique(data$TR)) !=1) {
    stop("'treat' must be either 0 or 1 where 1 indicates treated")
  }


  
### Unit and Time index
  ## Creating time index for strata fixed effect analysis


  ## unit index
  
  numeric.u.index <- as.numeric(as.factor(data[,unit.index]))
  numeric.u.index[is.na(numeric.u.index)] <- 0
  ## handling missing unit index
  uniq.u <- unique(na.omit(numeric.u.index))
  uniq.u <- sort(uniq.u[!(uniq.u %in% 0)])
  J.u <- length(uniq.u)

  
  data$u.index <- Index(numeric.u.index, uniq.u, J.u, tn.row)

  ## time index
  if (is.null(time.index)) {
    data$t.index <- GenTime(data$u.index, tn.row, length(uniq.u))
    numeric.t.index <- as.numeric(as.factor(data$t.index))
  } else {
    numeric.t.index <- as.numeric(as.factor(data[,time.index]))
    numeric.t.index[is.na(numeric.t.index)] <- 0
    ## handling missing time index
    uniq.t <- unique(na.omit(numeric.t.index))
    uniq.t <- sort(uniq.t[!(uniq.t %in% 0)])
    ## needs to sort for unbalnced panel, See Index()
    J.t <- length(uniq.t)
    data$t.index <- Index(numeric.t.index, uniq.t, J.t, tn.row)
  }
  uniq.t <- unique(data$t.index)

  

  ## unique unit number and sorting data
  if (method == "unit"){
    unit.number <- length(uniq.u)
  } else if (method == "time"){
    unit.number <- length(uniq.t)
  } else {
    stop("method should be either unit or time")
  }

  if (verbose)
    cat(" \nNumber of unique", method, "is", unit.number, "\n")


  ## order data 
  data <- data[order(data$u.index, data$t.index),]

  ## saving unit index for each unit
  name.unit <- unique(data[,unit.index])
  number.unit <- unique(numeric.u.index)

  ## saving time index for each time
  name.time <- unique(data[,time.index])
  number.time <- unique(numeric.t.index)

  
  
  ## quantity of interest: qoi

  if (missing(qoi)){
    causal <- "ate"
    ate.n <- 1
  }

  ate.n <-  att.n <- 0
  if (qoi == "ate"){
    causal <- "ATE (Average Treatment Effect)"
    ate.n <- 1
  }
  if (qoi == "att"){
    causal <- "ATT (Average Treatment Effect for the Treated)"
    att.n <- 1
  }

  if (is.null(estimator)) {
    est <- "NULL"
  }
  if (!is.null(estimator) && estimator == "fd"){
    est <- "FD (First-Difference)"
  }
  if (!is.null(estimator) && estimator == "did"){
    est <- "DID (Difference-in-Differences)"
  }


### propensity score estimation

  if (!is.null(pscore)) {
    data$y.star <- Transform(data[,outcome], data[,treat], data[,pscore])
  }
  if (is.null(pscore)) {

    ## formula for propensity score estimation
    f.left <- as.character(treat)
    f.right <- unlist(strsplit(as.character(formula), "~"))[2]
    c.formula <- paste(f.left,"~ -1 + ", f.right)
    
    
    ps.formula <- as.formula(c.formula) # formula for propensity score estimation


    data$index <- 1:nrow(data) # indexing for saving y.star
    
    if (verbose){
      cat("\nPropensity Score estimation started based on the following model:\n", c.formula, "\n")
    }
    flush.console()
    ## transform Y_{it}
    if (method=="time" & within.unit == TRUE){
      ## regression for propensity score estimation is run within unit
      Data.ps <- c()
      for (j in 1:length(uniq.t)){
        ## print(j)
        if (verbose){
          cat(".")
          flush.console()
        }

        pool <- data[data$t.index == uniq.t[j], ] # subset data by time j
        n.pool <- nrow(pool) # gives total number of units for the year j
        n.treat <- sum(pool[,treat]) # number of treated units in a given year t
        n.cont <- sum(1-pool[,treat]) # number of controlled units in a given year t
        ## estimate propensity score within each time t
        fit.ps <- bayesglm(ps.formula, family=binomial(link="logit"),
                           data = pool)
        pool$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
        pool.u.index <- unique(pool$u.index)
        pool$ystar <- Transform(pool[,outcome], pool[,treat], pool$fitted.ps)
        Data.ps <- rbind(Data.ps, pool)
      }
      data$y.star[Data.ps$index] <- Data.ps$ystar

    }

    if (method=="unit" & within.unit == TRUE) {
      ## regression for propensity score estimation is run within unit
      Data.ps <- c()
      for (i in 1:length(uniq.u)){

        if (verbose){
          if (i %% 100 == 0)
            cat(".")
          flush.console()
        }

        pool <- data[data$u.index == uniq.u[i], ] # subset data by unit i
        n.pool <- nrow(pool) # gives total number of times for the unit i
        n.treat <- sum(pool[,treat]) # number of treated times in a given unit i
        n.cont <- sum(1-pool[,treat]) # number of controlled times in a given unit i
        ## estimate propensity score within each time t
        fit.ps <- bayesglm(ps.formula, family=binomial(link="logit"), data
                           = pool)
        pool$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
        pool.t.index <- unique(pool$t.index)
        pool$ystar <- Transform(pool[,outcome], pool[,treat], pool$fitted.ps)
        Data.ps <- rbind(Data.ps, pool)
      }
      data$y.star[Data.ps$index] <- Data.ps$ystar

    }

    

    ## transform Y_{it}
    if (within.unit == FALSE){
      ## regression for propensity score estimation is on entire data
      n.treat <- sum(data[,treat]) # number of treated units 
      n.cont <- sum(1-data[,treat]) # number of controlled units 
      ## estimate propensity score within each time t
      fit.ps <- bayesglm(ps.formula, family=binomial(link="logit"),
                         data = data)
      data$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
      data$y.star <- Transform(data[,outcome], data[,treat], data$fitted.ps)
    }

  }
  
  if (verbose)
    cat("\nPropensity Score estimation done \n")
  flush.console()
  


  if (verbose)
    cat("\nWeight calculation started ")
  flush.console()
  ## Unit fixed effects models
  if ( (method=="unit" & qoi=="ate" & is.null(estimator) ) | (method=="unit" & qoi=="att" & is.null(estimator)) ) {
    W <- GenWeightsUnit(data$u.index, data$t.index, data$TR, data$C.it, tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, length(uniq.u)*length(uniq.t), verbose)
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)
    data$W.it <- VectorizeC(as.matrix(W), data$t.index, data$u.index, tn.row)
  }

  ## Time fixed effects models
  if ( (method=="time" & qoi=="ate" & is.null(estimator) ) | (method=="time" & qoi=="att" & is.null(estimator)) ) {
    W <- GenWeightsTime(data$t.index, data$u.index, data$TR, data$C.it, tn.row, length(uniq.t), length(uniq.u), ate.n, att.n, length(uniq.t)*length(uniq.u), verbose)
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)
    data$W.it <- VectorizeC(as.matrix(W), data$t.index, data$u.index, tn.row)
  }


  ## Within Unit First Difference
  if(( (method=="unit") && (qoi == "ate") && (!is.null(estimator) && estimator == "fd")) | ((method == "unit") && (qoi =="att") && (!is.null(estimator) && estimator == "fd"))) {
    W <- GenWeightsFD(data$u.index, data$t.index, data$TR, data$C.it, tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, verbose)
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)
    data$W.it <- VectorizeC(as.matrix(W), data$t.index, data$u.index, tn.row)
  }
  if (verbose)
    cat(" Weight calculation done \n")
  ## print(Sys.time())
  flush.console()
  ##cat("weights:", print(W), "\n")


  ## e <- environment()
  ## save(file = "temp.RData", list = ls(), env = e)


  Data.final <- as.data.frame(cbind(data$y.star, data$TR, data$W.it, data$u.index, data$t.index))
  colnames(Data.final) <- c("y.star", "TR", "W.it", "u.index", "t.index")



  ## Creating the matrix for final analysis
  Data.dm <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(Data.final)))
  Data.wdm <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(Data.final)))


  ## column names
  ## print(colnames(Data.final))
  colnames(Data.dm) <- colnames(Data.wdm) <- c("dep.var", "treat")

  ## de-meaning for fast calculation
  for (k in 1:2){
    if (method == "unit"){
      w.wdemean <- WWDemean(Data.final[,k], Data.final$W.it, Data.final$u.index, unit.number, tn.row)
      demean <- Demean(Data.final[,k], Data.final$u.index, unit.number, tn.row)
    }
    if (method == "time"){
      w.wdemean <- WWDemean(Data.final[,k], Data.final$W.it, Data.final$t.index, unit.number, tn.row)
      demean <- Demean(Data.final[,k], Data.final$u.index, unit.number, tn.row)
    }
    Data.wdm[,k] <- as.vector(w.wdemean)
    Data.dm[,k] <- as.vector(demean)
  }

  Data.dm$unit <- Data.final$u.index
  Data.dm$time <- Data.final$t.index
  Data.wdm$unit <- Data.final$u.index
  Data.wdm$time <- Data.final$t.index

  ## Print de(weighted)-meaned data

  ## cat("#####", "\n","De-meaned Data:","\n")
  ## print(Data.dm)
  ## cat("#####", "\n","sqrt(W) x (weighted demeaned):","\n")
  ## print(Data.wdm)

  ## final regression on weighted demeaned data
  fit.final <- lm(dep.var ~ -1 + treat, data = Data.wdm)
  fit.ols <- lm(dep.var ~ -1 + treat, data = Data.dm)

  ## print(summary(fit.final))

  ## set up data frame, with support for standard and modified responses
  
  mf.final <- model.frame(formula="dep.var ~ -1 + treat", data=Data.wdm)
  X <- model.matrix(attr(mf.final, "terms"), data=mf.final)
  Y <- model.response(mf.final, "numeric")
  p <- ncol(X)

  coef.wls <- fit.final$coef
  coef.ols <- fit.ols$coef
  d.f <- fit.final$df - unit.number

  sigma2 <- sum(resid(fit.final)^2)/d.f
  vcov.wls <- vcov(fit.final)*((fit.final$df)/(fit.final$df - unit.number))
  var.cov <- vcov.wls
  vcov.ols <- vcov(fit.ols)

  ## e <- environment()
  ## save(file = "temp.RData", list = ls(), env = e)

  
  ## residual <-  resid(fit.final)*1/sqrt(Data.final$W.it)
  residual <- c(data[,outcome] -  data[,treat]%*%t(coef.wls))
  resid.ols <- resid(fit.ols)

  ## save residuals
  Data.wdm$u.tilde <-  sqrt(Data.final$W.it)*resid(fit.final)
  Data.dm$u.hat <- u.hat <- resid(fit.ols)



### Robust Standard Errors

  ## contructing de(weighted)-meaned X matrix
  X.tilde <- as.matrix(Data.wdm[,2]) # NT x 1 (where p is number of covariates)
  X.hat <- as.matrix(Data.dm[,2]) # NT x 1 (where p is number of covariates)

  ## residuals
  u.tilde <- as.matrix(Data.wdm$u.tilde) # NT x 1
  u.hat <- as.matrix(Data.dm$u.hat)  # NT x 1

  ## unit vector (the length should be same as nrow(X.tilde) and length(u.tilde)
  wdm.unit <- as.vector(Data.wdm$unit)
  uniq.u.pw <- length(unique(wdm.unit)) # number of uniq units with positive weights


  ## cat("3: Robust error calculation started\n")
  ## print(Sys.time())

  ginv.XX.tilde <- ginv(crossprod(X.tilde, X.tilde))
  ginv.XX.hat <- ginv(crossprod(X.hat, X.hat))

  diag.ee.tilde <- c(u.tilde^2)
  diag.ee.hat <- c(u.hat^2)

  
  if ((hetero.se == TRUE) & (auto.se == TRUE)) {# Default is Arellano
    std.error <- "Heteroscedastic / Autocorrelation Robust Standard Error"

    ## 1. arbitrary autocorrelation as well as heteroskedasticity (Eq 12)

    Omega.hat.HAC <- OmegaHatHAC(nrow(X.tilde), p, wdm.unit, J.u, X.tilde, u.tilde)
    Omega.hat.HAC <- matrix(Omega.hat.HAC, nrow = p, ncol = p)
    Omega.hat.HAC <- (1/(nrow(X.tilde)))* Omega.hat.HAC # without degree of freedom adjustment 
    ## Omega.hat.HAC <- (1/(nrow(X.tilde) - J.u - p))* Omega.hat.HAC
    
    Omega.hat.fe.HAC <- OmegaHatHAC(nrow(X.hat), p, wdm.unit, J.u, X.hat, u.hat)
    Omega.hat.fe.HAC <- matrix(Omega.hat.fe.HAC, nrow = p, ncol = p)
    Omega.hat.fe.HAC <- (1/(nrow(X.hat))) * Omega.hat.fe.HAC # without degree of freedom adjustment
    ## Omega.hat.fe.HAC <- (1/(nrow(X.hat) - J.u - p)) * Omega.hat.fe.HAC

    Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Omega.hat.HAC %*% (nrow(X.tilde) * ginv.XX.tilde)
    Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (nrow(X.hat) * ginv.XX.hat)

  } else if ( (hetero.se == TRUE) & (auto.se == FALSE) ) {# independence across observations but heteroskedasticity
    std.error <- "Heteroscedastic Robust Standard Error"

    ## 2. independence across observations but heteroskedasticity (Eq 11)
    
    ## Omega.hat.HC <- OmegaHatHC(nrow(X.tilde), p, wdm.unit, J.u, X.tilde, u.tilde)
    ## Omega.hat.HC <- matrix(Omega.hat.HC, nrow = p, ncol = p)
    ## ## same as the following matrix multiplication but slower
    ## ## Omega.hat.he <- t(X.tilde) %*% diag(c(u.tilde)^2, nrow=nrow(X.tilde)) %*% X.tilde
    ## Omega.hat.HC <- (1/(nrow(X.tilde) - J.u - p)) * Omega.hat.HC


    Omega.hat.HC <- (1/(nrow(X.tilde) - J.u - p))*(crossprod((X.tilde*diag.ee.tilde), X.tilde)) 
    
    Omega.hat.fe.HC <- (1/(nrow(X.hat) - J.u - p))*(crossprod((X.hat*diag.ee.hat), X.hat))  
    
    
### Stock-Watson (Econometrica 2008: Eq(6)): Bias-asjusted for balance panel
    if (unbiased.se == TRUE) {
      ## check if panel is balanced
      if (sum(as.numeric(apply(matrix(table(wdm.unit)), 1, mean) != mean(matrix(table(wdm.unit)))))  == 0 ){
        std.error <- "Heteroskedastic Standard Error (Stock-Watson Biased Corrected)"
        
        B.hat <- matrix(0, nrow=dim(X.tilde)[2], ncol=dim(X.tilde)[2])
        B2.hat <- matrix(0, nrow=dim(X.hat)[2], ncol=dim(X.hat)[2])
        for (i in 1:J.u) {
          ## cat("unit", i, "\n")
          X.i <- X.tilde[Data.wdm$unit == i,]
          X2.i <- X.hat[Data.wdm$unit == i,]
          ## print(sum(as.numeric(Data.wdm$unit == i)))
          ## print(X.i)
          if (sum(as.numeric(Data.wdm$unit == i)) > 1) { 
            u.i <- u.tilde[Data.wdm$unit == i]
            u2.i <- u.hat[Data.wdm$unit == i]
            ## print(length(u.i))
            flush.console()
            XX.i <- crossprod(X.i, X.i)
            XX2.i <- crossprod(X2.i, X2.i)
            B.hat <- B.hat + ((1/J.t)* XX.i*(1/(J.t-1))*sum(u.i^2))
            B2.hat <- B2.hat + ((1/J.t)* XX2.i*(1/(J.t-1))*sum(u2.i^2))
          }
        }

        B.hat <- B.hat * (1/J.u)
        B2.hat <- B2.hat * (1/J.u)

        cat("time", J.t, "\n")
        Sigma_HRFE <- ((J.t-1)/(J.t-2))*(Omega.hat.HC - (1/(J.t-1))*B.hat)
        Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Sigma_HRFE %*% (nrow(X.tilde) * ginv.XX.tilde)

        Sigma2_HRFE <- ((J.t-1)/(J.t-2))*(Omega.hat.fe.HC - (1/(J.t-1))*B2.hat)
        Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Sigma2_HRFE %*% (nrow(X.hat) * ginv.XX.hat)
      } else { 
        stop ("unbiased.se == TRUE is allowed only when panel is balanced")
      }
    }
    
    Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Omega.hat.HC %*% (nrow(X.tilde) * ginv.XX.tilde)
    
    Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Omega.hat.fe.HC %*% (nrow(X.hat) * ginv.XX.hat)


  } else if ( (hetero.se == FALSE) & (auto.se == FALSE) ) {# indepdence and homoskedasticity
    std.error <- "Homoskedastic Standard Error"
    
    Psi.hat.wfe <- nrow(X.tilde) * (sigma2 * ginv.XX.tilde)
    
    ## same as the following

    ## Psi.hat.wfe2 <- J.u * sigma2 * solve(XX.tilde)
    ## cat("compare homoskedasticity s.e\n")
    ## print(sqrt(diag(Psi.hat.wfe)))
    ## print(sqrt(diag(Psi.hat.wfe2)))
    
    Psi.hat.fe <- nrow(X.hat) * vcov(fit.ols)

  } else if ( (hetero.se == FALSE) & (auto.se == TRUE) ) {# Kiefer

    stop ("robust standard errors with autocorrelation and homoskedasiticy is not supported")

  }

  var.cov <- Psi.hat.wfe * (1/nrow(X.tilde))
  var.cov.fe <- Psi.hat.fe *(1/nrow(X.hat))



### White (1980) Test: Theorem 4

      if (White == TRUE){
        
        diag.ee <- c(u.hat) * c(u.tilde)
        
        Lambda.hat1 <-  1/((nrow(X.hat)))* (crossprod((X.hat*diag.ee), X.tilde))  
        Lambda.hat2 <-  1/((nrow(X.tilde)))* (crossprod((X.tilde*diag.ee), X.hat))  

        Phi.hat <- Psi.hat.wfe + Psi.hat.fe - (nrow(X.hat)*ginv.XX.hat) %*% Lambda.hat1 %*% (nrow(X.tilde)*ginv.XX.tilde) - (nrow(X.tilde)*ginv.XX.tilde) %*% Lambda.hat2 %*% (nrow(X.hat)*ginv.XX.hat)

        ## White test: null hypothesis is ``no misspecification''

        white.stat <- as.double(Re(nrow(X.hat) * t(coef.ols - coef.wls) %*% ginv(Phi.hat) %*% (coef.ols - coef.wls)))
        test.null <- pchisq(as.numeric(white.stat), df=p, lower.tail=F) < White.alpha
        white.p <- pchisq(as.numeric(white.stat), df=p, lower.tail=F)
        flush.console()

        ## if (verbose) {
        ##   cat("\nWhite calculation done")
        ##   flush.console()
        ## }

      } else {
        white.stat <- "NULL"
        test.null <- "NULL"
        white.p <- "NULL"
      }

  mf <- as.data.frame(cbind(data[,outcome], data[,treat]))
  colnames(mf) <- c(outcome, treat)
  y <- as.data.frame(data[,outcome])
  colnames(y) <- outcome
  x <- as.data.frame(data[,treat])
  colnames(x) <- treat
  

### Saving results
  z <- list(coefficients = coef.wls,
            x = X,
            y = y,
            mf = mf,
            call = pwfe.call,
            vcov = var.cov,
            se = sqrt(diag(var.cov)),
            sigma = sqrt(sigma2),
            df = d.f,
            residuals = residual,
            W = W,
            unit.name = name.unit,
            unit.index = number.unit,
            time.name = name.time,
            time.index = number.time,
            method = method,
            causal = causal,
            est = est,
            std.error = std.error,
            White.pvalue = white.p,
            White.alpha = White.alpha,
            White.stat = white.stat,
            White.test = test.null)
  class(z) <- "pwfe"
  z



}


print.pwfe <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nStd.Err:\n")
  print(x$se)
}



summary.pwfe <- function(object, signif.stars = getOption("show.signif.stars"),...){
  se <- object$se
  sigma <- object$sigma
  df <- object$df
  tval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               Std.Err = se,
               t.value = tval,
               p.value = 2*pt(-abs(tval), df = object$df))
  res <- list(call = object$call,
              coefficients = TAB,
              sigma = object$sigma,
              df = object$df,
              W = object$weights,
              residuals = object$residuals,
              method = object$method,
              causal = object$causal,
              estimator = object$est,
              std.error = object$std.error,
              White.pvalue = object$White.pvalue,
              White.alpha = object$White.alpha,
              White.stat = object$White.stat,
              White.test = object$White.test)
  class(res) <- "summary.pwfe"
  res
}


print.summary.pwfe <- function(x, ...){
  cat("\nMethod:", x$method, "Fixed Effects (Propensity Score)\n")
  cat("\nQuantity of Interest:", x$causal)
  cat("\nEstimator:", x$estimator)
  cat("\nStandard Error:", x$std.error)
  cat("\n")
  cat("\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
  cat("\nResidual standard error:", format(signif(x$sigma,
                                                  4)), "on", x$df, "degrees of freedom")
  cat("\nWhite statistics for functional misspecification:", x$White.stat, "with Pvalue=", x$White.pvalue)
  cat("\nReject the null of NO misspecification:", x$White.test)
  cat("\n")
}
