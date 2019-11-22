PanelWFE <- function (formula, data, treat = "treat.name",
                      unit.index, time.index = NULL, method = "unit",
                      dyad1.index = NULL, dyad2.index = NULL,
                      qoi = "ate", estimator = NULL, C.it = NULL,
                      hetero.se = TRUE, auto.se = TRUE,
                      dyad.se = FALSE,
                      White = TRUE,
                      verbose = FALSE, unbiased.se = FALSE, unweighted = FALSE,
                      store.wdm = FALSE,
                      tol = sqrt(.Machine$double.eps)){
  if (qoi == "att") {
    data$big_W_it <- data$Wit_att0
  } else if (qoi == "atc"){
    data$big_W_it <- data$Wit_atc0
    qoi = "att"
  } else if (qoi == "ate") {
    data$big_W_it <- data$Wit_att0 + data$Wit_atc0
  }
  
  wfe.call <- match.call()
  ## set up data frame, with support for standard and modified responses
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf, "numeric")
  tn.row <- nrow(mf) # total number of rows in data
  
  class(data) <- "data.frame"
  
  ## ## remove missing variables: removing rows with missing values in either y or treat
  remove.indices <- which(!rownames(data) %in% rownames(mf))
  
  
  if (length(remove.indices) > 0){
    data <- data[-remove.indices,]
    
  }
  
  data$y <- y
  
  ## Creating dummies variables for White test in the end
  X <- as.data.frame(x[,-1])
  p <- ncol(X)
  
  ## C.it
  ## Default for ATE
  if (is.null(C.it)){
    data$C.it <- as.integer(rep(1, nrow(data)))
  }
  
  ## --------------------------------------------------------
  ## Unit and Time index
  ## -------------------------------------------------------- 
  ## Creating time index for strata fixed effect analysis
  
  ## storing original unit and time index
  orig.unit.idx <- as.character(data[,unit.index])
  
  ## unit index
  numeric.u.index <- as.numeric(as.factor(data[,unit.index]))
  numeric.u.index[is.na(numeric.u.index)] <- 0
  ## handling missing unit index
  uniq.u <- unique(na.omit(numeric.u.index))
  uniq.u <- sort(uniq.u[!(uniq.u %in% 0)])
  J.u <- length(uniq.u)
  data$u.index <- Index(numeric.u.index, uniq.u, J.u, tn.row)
  
  ## for dyadic data: unit 1 and unit 2 that compose each dyad
  if(dyad.se == TRUE){
    ## Warning 
    if(is.null(dyad1.index) | is.null(dyad2.index)){
      stop("Warning: For dyadic data, two separate unit indices for the members of each dyad -- dyad1.index, dyad2.indx -- should be provided")
    }
    
    data$dyad <- data[, unit.index]
    data$c1 <- data[, dyad1.index]
    data$c2 <- data[, dyad2.index]
    
  }
  
  ## time index
  if (is.null(time.index)) {
    data$t.index  <- GenTime(data$u.index, tn.row, length(uniq.u))
    numeric.t.index <- as.numeric(as.factor(data$t.index))
    
    ## storing original unit and time index
    orig.time.idx <- as.character(data$t.index)
  } else {
    ## storing original unit and time index
    orig.time.idx <- as.character(data[,time.index])
    
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

  
  ## order data by unit index
  
  tmp <- cbind(X, data$u.index, data$t.index) 
  tmp <- tmp[order(tmp[,(p+1)], tmp[,(p+2)]),]
  X <- as.data.frame(tmp[,-((p+1):(p+2))])
  colnames(X) <- colnames(x)[-1]
  
  data <- data[order(data$u.index, data$t.index),]
  y <- data$y[order(data$u.index, data$t.index)]
  
  
  ## saving unit index for each unit
  name.unit <- unique(data[,unit.index])
  number.unit <- unique(numeric.u.index)
  units <- as.data.frame(cbind(number.unit,as.character(name.unit)))
  colnames(units) <- c("unit.index", "unit")
  rownames(units) <- seq(1:nrow(units))
  
  ## saving time index for each time
  
  if (is.null(time.index)) {
    name.time  <- unique(data$t.index)
  } else {
    name.time <- unique(data[,time.index])
  }
  number.time <- unique(numeric.t.index)
  times <- cbind(number.time, name.time)
  times <- as.data.frame(times[order(times[,1]),])
  colnames(times) <- c("time.index", "time")
  
  
  ## new model frame order by u.index
  mf.col <- colnames(mf)
  mf.sorted <- cbind(y,X)
  colnames(mf.sorted) <- mf.col
  
  ## treatment variable
  data$TR <- as.numeric(data[,treat])
  
  if (length(unique(data$TR)) !=2)
    stop("'treat' must be a binary vector: there are more than two values of treatment")
  
  if (sum(unique(data$TR)) !=1)
    stop("'treat' must be a either zero or one where one indicates treatment")
  
  ### --------------------------------------------------------
  ### Quantity of interest: qoi
  ### -------------------------------------------------------- 
  
  if (missing(qoi)){
    qoi <- "ate"
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

  
  ### Weights calculation
    ## Differences-in-difference
  if(( (method=="unit") & (qoi == "ate") & (!is.null(estimator) & estimator == "did")) |
      ( (method == "unit") & (qoi =="att") & (!is.null(estimator) & estimator == "did")) |
      ( (method=="unit") & (qoi == "ate") & (!is.null(estimator) & estimator == "Mdid")) |
      ( (method == "unit") & (qoi =="att") & (!is.null(estimator) & estimator == "Mdid"))       
    ) {
      
      data$W.it <- data$big_W_it
      
      ## creating index for sparse dummy matrix
      
      u <- as.matrix(table(data$u.index))
      
      Udummy.i <- seq(1:sum(u))
      Udummy.j <- c()
      
      for (j in 1:length(uniq.u)) {
        Udummy.j <- c(Udummy.j, rep(j, u[j,1]))
      }
      Udummy  <- Matrix::sparseMatrix(x=1, i=Udummy.i, j=Udummy.j)
      
      t <- as.matrix(table(data$u.index, data$t.index))
      
      ## checking panel structure
      if (verbose)
        if (length(which(t==0)) > 0){
          cat("\nUnbalanced Panel Data\n")
        }
      if ( length(which(t>1)) > 0 ){
        stop ("\nunit-time pair is not unique\n")
      }
      flush.console()
      
      ## this takes time: should be made more efficient
      Tdummy <- array(0,dim=c(0,length(uniq.t)))
      for (j in 1:nrow(t)) {
        Tdummy <- rbind(Tdummy,  Diagonal(x = t[j,], n=length(uniq.t)))
      }
      
      ## when panel is unbalanced there will be rows of zeros
      zero <- which(apply(Tdummy, 1, mean) == 0)
      if(length(zero) > 0) {
        Tdummy <- Tdummy[-zero,]
      }
      
      
      a <- unlist(strsplit(as.character(formula), "~"))
      formula.ni <- as.formula(paste(a[2], "~ -1 + ",  a[3]))
      
      ############################################################
    
      ## subset observations with non-zero weights
      if (White == TRUE){
        nz.index <- seq(1,tn.row)
      } else { # cannot calculate White statistics 
        ## exclude zero-weights observations for efficient calculation
        nz.index <- data$W.it !=0
        tn.row <- length(which(data$W.it !=0))
      }
      
      nz.obs <- sum(as.numeric(data$W.it !=0))
      
      X <- as.matrix(X)[nz.index,]
      Y <- as.matrix(data$y)[nz.index]
      
      
      ## removing zero weights rows
      Udummy <- Udummy[nz.index,]
      Tdummy <- Tdummy[nz.index,]
      
      
      ### removing zero columns for full-rank (after deleting zero weights observations)
      
      ## 1. Unit dummies
      
      u.zero <- try(which(as(apply(Udummy, 2, sum), "sparseVector") == 0), silent=TRUE)
      
      if (class(u.zero) == "try-error") {
        u.zero <- c()
        for (i in 1:ncol(Udummy)) {
          if (sum(Udummy[,i]) == 0){
            temp <- i
            u.zero <- c(u.zero, temp)
          }
        }
        if (verbose)
          cat("\nReached Memory Limit: White == FALSE option is recommended\n")
        flush.console()
        
        if (length(u.zero) > 0) {
          Udummy <- Udummy[,-u.zero]
          gc()
        }
      } else {
        if (length(u.zero) > 0) {
          Udummy <- Udummy[,-u.zero]
          gc()
        }
      }
      
      
      ## 2. Time dummies
      
      if (length(which(as(apply(Tdummy, 2, sum), "sparseVector")==0)) > 0) {
        zero <- which(as(apply(Tdummy, 2, sum), "sparseVector") == 0)
        n.zero <- length(zero)
        Tdummy <- Tdummy[,-zero]
        ## delete last column
        ## last <- ncol(Tdummy) 
        ## Tdummy <- Tdummy[,-last]
        ## adding a column of 1000's for numerical stability of ginv
        Tdummy <- cbind(Tdummy, rep(1000, nrow(Tdummy)))
        gc()
      } else {
        last <- ncol(Tdummy)
        Tdummy <- Tdummy[,-last] # for identification exclude the last year dummy
        ## adding a column of 1's for numerical stability of ginv
        Tdummy <- cbind(Tdummy, rep(1000, nrow(Tdummy)))
        gc()
      }
      
      ## number of columns for Unit/Time dummy matrix
      n.Udummy <- ncol(Udummy)
      n.Tdummy <- ncol(Tdummy)
      
      
      ## combining unit dummy matrix and X matrix
      D <- Matrix(cbind(Udummy, Tdummy))
      ## Note: Tdummy part does not have zero columns, but Udummy part has it
      ## will be addressed thie issue below by n.zero
      
      ## final number of dummies

      
      Transformed <- complex_project(data, D, Tdummy, Udummy, n.Tdummy, n.Udummy, nz.index, tn.row, X, Y, J.u, tol, White)
      browser()
      y.tilde <- Transformed[,1]
      X.tilde <- as.matrix(Transformed[,-1])
      
      ## save weighted demeaned dataframe
      if (store.wdm == TRUE){
        Y.wdm <- y.tilde
        X.wdm <- X.tilde
      } else {
        Y.wdm <- NULL
        X.wdm <- NULL
      }
      
      rm(Transformed)
      gc()
      
      if (ncol(X.tilde) == 1) {
        colnames(X.tilde) <- a[3]
      } else {
        colnames(X.tilde) <- colnames(X)
      }
      
      ginv.XX.tilde <- MASS::ginv(crossprod(X.tilde))
      betaT <- ginv.XX.tilde%*% crossprod(X.tilde, y.tilde)
      if (length(betaT) == 1) {
        colnames(betaT) <- a[3]
      }
      ## print(betaT)
      coef.wls <- matrix(as.double(Re(betaT)))
      rownames(coef.wls) <- colnames(X.tilde)
      #NOTE: the results of coef.wls match the coefficient calculated by the bootstrap implementation
      browser()
      
      #NOTE: here is where the problem begins. e tilde is huge.      
      ## weighted residuals
      e.tilde <- (y.tilde - X.tilde %*% betaT)
      colnames(e.tilde) <- "e.tilde"
      ## true residuals
      w.sqrt <- sapply(data$W.it[nz.index], Im.sqrt)
      
      resid <- try(1/w.sqrt * (y.tilde - X.tilde %*% betaT))
      
      ## in case zero weights observations are not excluded
      if (White == TRUE){
        if (sum(as.numeric(w.sqrt==0)) > 0) { 
          zero.index <- data$W.it ==0
          resid[zero.index] <- 0
        }
      }
      rm(w.sqrt)
      gc()
      

      diag.ee.tilde <- as.vector(e.tilde * e.tilde)
      

      #ginv.XX.hat <- MASS::ginv(crossprod(X.hat, X.hat))
      
      d.f <- length(y.tilde)
      
      
      sigma2 <- as.double(Re(sum(resid^2)/d.f))
      
      ## Remove observations with zero weights
      ## data backup
      data.zero <- data
      zero.ind <- which(data$W.it==0)
      if(length(zero.ind) > 0){
        data.nonzero <- data[-zero.ind, ]
      } else {
        data.nonzero <- data
      }
      n.units <- length(unique(data$u.index))
      n.times <- length(unique(data$t.index))
      
      Mstar <- nrow(data.nonzero)            
      
                    
      if(unweighted == FALSE){
        n.nonzero.units <- length(unique(data.nonzero$u.index))
        n.nonzero.times <- length(unique(data.nonzero$t.index))
      } else {
        n.nonzero.units <- n.units
        n.nonzero.times <- n.times
      }
      
      x.vars <- colnames(x)
      x.vars <- x.vars[-grep("Intercept", x.vars)]
      nK <- length(x.vars) 
      variables <- c("y", x.vars)
      
      
      
      ## #######################################################################            
      ## (Robust) standard errors (GMM asymptotic variance for wfe)
      ## calculating GMM standard errors
      ## #######################################################################
      #### CONTAINS DEGREE OF FREEDOM ADJUSTMENTS 
      if ((hetero.se == TRUE) & (auto.se == TRUE)){
        
        ## 1. arbitrary autocorrelation as well as heteroskedasticity (Eq 12)
        ## -----------------------------------------------------
        ## vcov matrix for WFE 
        ## -----------------------------------------------------
        browser()
        ## degrees of freedom adjustment
        df_wfe2 <- (Mstar/(Mstar-1))*((Mstar-nK)/(Mstar- n.nonzero.units - n.nonzero.times - nK))
        #df_wfe2 <- 1 / (Mstar - n.nonzero.times - n.nonzero.units - nK)

        ##NOTE: XeeX, e.tilde are extremely large
        XeeX <- as.double(comp_OmegaHAC(c(X.tilde), e.tilde, c(X.tilde), e.tilde,
                                        dim(X.tilde)[1], dim(X.tilde)[2], data$u.index, J.u))
        XeeX.wfe <- matrix(XeeX, nrow=ncol(X.tilde), ncol=ncol(X.tilde), byrow=T)
        
        Psi.hat.wfe <- df_wfe2*((ginv.XX.tilde %*% XeeX.wfe %*% ginv.XX.tilde))

        
        
      } 
      
      ## storing standard errors
      vcov.wfe <- Psi.hat.wfe 
      se.did <- as.double(Re(sqrt(diag(vcov.wfe))))
      
      ### Saving results
      
      z <- list(coefficients = coef.wls,
                se = se.did,
                sigma = try(sqrt(sigma2))
                )
      
      class(z) <- "wfedid"
      z
    }
    
  
}



### Matrix multiplication: two sparse matrix : A=R1+I1i, B=R2+I2i

## A%*%B: result is a list where [[1]] is real part [[2]] is imaginary part of complex matrix multiplication
Sparse_compMatrixMultiply <- function(R1,I1,R2,I2) {
  result <- list()
  result[[1]] <- drop0(R1%*%R2) - drop0(I1%*%I2)
  result[[2]] <- drop0(R1%*%I2) + drop0(I1%*%R2)
  result
}

## A %*% t(B): result is a list where [[1]] is real part [[2]] is imaginary part of complex matrix multiplication
Sparse_compMatrix_tcrossprod <- function(R1,I1,R2,I2) {
  ## real part
  result <- list()
  result[[1]] <- drop0(tcrossprod(R1, R2)) - drop0(tcrossprod(I1, I2))
  result[[2]] <- drop0(tcrossprod(R1, I2)) + drop0(tcrossprod(I1, R2))
  result
}

## t(A) %*% B: result is a list where [[1]] is real part [[2]] is imaginary part of complex matrix multiplication
Sparse_compMatrix_crossprod <- function(R1,I1,R2,I2) {
  ## real part
  result <- list()
  result[[1]] <- drop0(crossprod(R1, R2)) - drop0(crossprod(I1, I2))
  result[[2]] <- drop0(crossprod(R1, I2)) + drop0(crossprod(I1, R2))
  result
}

## sqrt of weights (imaginary numbers)
Im.sqrt <- function(weights) {
  if (weights >= 0) # real number
    im.w <- complex(real=sqrt(weights), imaginary =0)
  if (weights < 0) # imaginary number
    im.w <- complex(real=0, imaginary = sqrt(-weights))
  invisible(im.w)
}

## sum of sqrt(W.it) across years for each dyad
general.inv <- function(weight, tol = tol){
  if(abs(weight) < tol) {
    out <- 0
  } else {
    out <- 1/weight
  }
  out
}

complex_project <- function(data, D, Tdummy, Udummy, n.Tdummy, n.Udummy, nz.index, tn.row, X, Y, J.u, tol, White)
{
  ### Projection onto Complex-plane
  
  ## vector of sqrt(W.it)
  w.sqrt <- sapply(data$W.it[nz.index], Im.sqrt)
  
  
  #########################################################################
  
  
  R1 <- Diagonal(x = Re(w.sqrt), n=tn.row)
  I1 <- Diagonal(x = Im(w.sqrt), n=tn.row)
  
  yL <- list()
  yL[[1]] <- Matrix(Y)
  yL[[2]] <- Matrix(0, nrow=nrow(yL[[1]]),  ncol=ncol(yL[[1]]))
  y.starL <- Sparse_compMatrixMultiply(R1,I1, yL[[1]], yL[[2]])
  rm(yL)
  gc()
  
  xL <- list()
  xL[[1]] <- Matrix(X)
  xL[[2]] <- Matrix(0, nrow=nrow(xL[[1]]),  ncol=ncol(xL[[1]]))
  x.starL <- Sparse_compMatrixMultiply(R1,I1, xL[[1]], xL[[2]])
  rm(xL)
  gc()
  
  #########################################################################
  
  ## create D1, and D2 
  R2 <- Matrix(D)
  I2 <- Matrix(0, nrow=nrow(D), ncol=ncol(D))
  
  D.starL <- Sparse_compMatrixMultiply(R1,I1,R2,I2)
  
  ## e <- environment()
  ## save(file = "temp.RData", list = ls(), env = e)
  
  rm(R2, I2)
  gc()
  
  
  D1.starL <- D2.starL <- list()
  
  ## e <- environment()
  ## save(file = "temp.RData", list = ls(), env = e, compress = TRUE)
  
  
  D1.starL[[1]] <- drop0(D.starL[[1]][,1:n.Udummy])
  D1.starL[[2]] <- drop0(D.starL[[2]][,1:n.Udummy])
  D2.starL[[1]] <- drop0(D.starL[[1]][,((n.Udummy+1):(n.Udummy+n.Tdummy))])
  D2.starL[[2]] <- drop0(D.starL[[2]][,((n.Udummy+1):(n.Udummy+n.Tdummy))])
  ## cat("Number of columns for D1", n.Udummy, "\n")
  ## cat("Number of columns for D2", n.Tdummy, "\n")
  
  rm(D.starL)
  gc()
  #########################################################################
  
  
  
  sum.sqrtW <- complex(real=apply(Sparse_compMatrixMultiply(R1,I1,D1.starL[[1]],D1.starL[[2]])[[1]], 2, sum), imaginary=rep(0, length(J.u)))
  rm(R1, I1)
  gc()
  
  
  inv.weight <- sapply(sum.sqrtW, general.inv, tol)
  
  rm(sum.sqrtW)
  gc()
  
  ## e <- environment()
  ## save(file = "temp.RData", list = ls(), env = e)
  
  ginvW <- list()
  ginvW[[1]] <- Diagonal(x = Re(inv.weight))
  ginvW[[2]] <- Diagonal(x = Im(inv.weight))
  
  rm(inv.weight)
  gc()
  
  Dginv <- Sparse_compMatrixMultiply(D1.starL[[1]], D1.starL[[2]], ginvW[[1]], ginvW[[2]])
  rm(ginvW)
  gc()
  
  P1L <-  Sparse_compMatrix_tcrossprod(Dginv[[1]], Dginv[[2]], D1.starL[[1]], D1.starL[[2]])
  
  rm(Dginv, D1.starL)
  gc()
  
  Q1L <- list()
  Q1L[[1]] <- drop0(Diagonal(x=1, n=nrow(P1L[[1]])) - P1L[[1]])
  Q1L[[2]] <- drop0(Diagonal(x=0, n=nrow(P1L[[1]])) - P1L[[2]])
  
  
  
  Q <- try(Sparse_compMatrixMultiply(Q1L[[1]], Q1L[[2]], D2.starL[[1]], D2.starL[[2]]), silent = TRUE)
  if ((class(Q) == "try-error") & (White == TRUE)) {
    stop ("Insufficient memory. White = FALSE option is needed")
  }
  rm(Q1L, D2.starL)
  gc()
  
  
  Q.matrix <- matrix(complex(real=as.matrix(Q[[1]]), imaginary=as.matrix(Q[[2]])), nrow=nrow(Q[[1]]))
  
  QQ.inv <- list()
  QQ.inv[[1]] <- drop0(Matrix(Re(MASS::ginv(crossprod(Q.matrix)))))
  QQ.inv[[2]] <- drop0(Matrix(Im(MASS::ginv(crossprod(Q.matrix)))))
  rm(Q.matrix)
  gc()
  
  PL <- list()
  
  Q.QQinv <- Sparse_compMatrixMultiply(Q[[1]], Q[[2]], QQ.inv[[1]], QQ.inv[[2]]) 
  rm(QQ.inv)
  gc()
  
  
  
  #########################################################################
  ### Fast Projection in R
  
  
  YX.starL <- list()
  YX.starL[[1]] <- cbind(y.starL[[1]], x.starL[[1]])
  YX.starL[[2]] <- cbind(y.starL[[2]], x.starL[[2]])
  rm(y.starL, x.starL)
  gc()
  
  
  P1.YX <- Sparse_compMatrixMultiply(P1L[[1]], P1L[[2]], YX.starL[[1]], YX.starL[[2]])
  rm(P1L)
  gc()
  
  
  Q.YX <- Sparse_compMatrix_crossprod(Q[[1]], Q[[2]], YX.starL[[1]], YX.starL[[2]])
  rm(Q)
  gc()
  
  QQQQ.YX <- Sparse_compMatrixMultiply(Q.QQinv[[1]], Q.QQinv[[2]], Q.YX[[1]], Q.YX[[2]])
  
  
  Transformed <- matrix(complex(real=as.vector(YX.starL[[1]] - P1.YX[[1]] - QQQQ.YX[[1]]), imaginary=as.vector(YX.starL[[2]] - P1.YX[[2]] - QQQQ.YX[[2]])), nrow=tn.row)
  rm(YX.starL, P1.YX, QQQQ.YX)
  return(Transformed)
}


