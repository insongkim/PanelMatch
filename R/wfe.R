wfe <- function (formula, data, treat = "treat.name",
                 unit.index, time.index = NULL, method = "unit",
                 qoi = "ate", estimator = NULL, C.it = NULL,
                 hetero.se = TRUE, auto.se = TRUE,
                 White = TRUE, White.alpha = 0.05,
                 verbose = TRUE, unbiased.se = FALSE, unweighted = FALSE,
                 store.wdm = FALSE, maxdev.did= NULL,
                 tol = sqrt(.Machine$double.eps)){


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
        if (verbose)
            cat(" \nMissing values are removed\n")
    }

    data$y <- y
    
    ## Creating dummies variables for White test in the end
    X <- as.data.frame(x[,-1])
    p <- ncol(X)

    ## e <- environment()
    ## save(file = "temp.RData", list = ls(), env = e)
    
### --------------------------------------------------------
### Warnings
### --------------------------------------------------------

    ## Warning for missing unit & time index
    if (missing(unit.index))
        stop("'unit.index' or index for strata should be provided")

    if (is.null(time.index) & method == "time")
        stop("'time.index' should be provided")

    ## Warning for methods
    if(method=="time" && !is.null(estimator) && estimator == "fd")
        stop("First Difference is not compatible with 'time' method: set method == 'unit'")

    if(method=="time" && !is.null(estimator) && estimator == "did")
        stop("Difference-in-Differences is not compatible with 'time' method: set method == 'unit'")

    if(method=="time" && !is.null(estimator) && estimator == "Mdid")
        stop("Match-Difference-in-Differences is not compatible with 'time' method: set method == 'unit'")

    if(is.null(time.index) && !is.null(estimator) && estimator == "fd")
        stop("First Difference cannot calculate when 'time.index' is missing")

    ## Warning for C.it
    if (!is.null(C.it)){
        Cit <- data[,C.it]
        if (!is.numeric(Cit) && length(Cit)!= tn.row)
            stop("'C.it' must be a numeric vector with length equal to number of observations")
        if ( sum(Cit < 0) > 0 )
            stop("'C.it' must be a non-negative numeric vector")
    }

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
    ## treat should be 0,1 where 1 indicates treatment

    ## Warning for maxdev.did
    if(!is.null(maxdev.did) && maxdev.did < 0){
        stop("Warning: maxdev.did should be a positive numeric value")
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

    if (verbose)
        cat(" \nNumber of unique", method, "is", unit.number, "\n")

    ## e <- environment()
    ## save(file = "temp.RData", list = ls(), env = e)
### -------------------------------------------------------- 

    ## order data by unit index

    tmp <- cbind(X, data$u.index, data$t.index) 
    tmp <- tmp[order(tmp[,(p+1)], tmp[,(p+2)]),]
    X <- as.data.frame(tmp[,-((p+1):(p+2))])
    colnames(X) <- colnames(x)[-1]

    data <- data[order(data$u.index, data$t.index),]
    y <- data$y[order(data$u.index, data$t.index)]

    ## e <- environment()
    ## save(file = "temp.RData", list = ls(), env = e)

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

    if (is.null(estimator)) {
        est <- "NULL"
    }
    if (!is.null(estimator) && estimator == "fd"){
        est <- "FD (First-Difference)"
    }
    if (!is.null(estimator) && estimator == "did"){
        est <- "DID (Difference-in-Differences)"
    }
    if (!is.null(estimator) && estimator == "Mdid"){
        est <- "DID (Difference-in-Differences) with Matching on Pre-treatment Outcome"
    }

    if (unweighted == TRUE){
        causal <- "Unweighted (Standard) Fixed Effect"
    }

    
### Weights calculation

    ## One-way Weighted FE
    if (is.null(estimator) || estimator=="fd") {

        if (verbose) {
            cat("\nWeight calculation started ")
            flush.console()
        }
        
        ## Standard Fixed effect
        if (unweighted == TRUE) {
            data$W.it <- rep(1, nrow(data))
            W <- matrix(1,nrow=length(uniq.t), ncol=length(uniq.u))
        } else {
            
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
        }
        
        
        if (verbose) {
            cat("Weight calculation done \n")
            flush.console()
        }
        
        if(( (method=="unit") && (qoi == "ate") && (!is.null(estimator) && estimator == "fd")) | ((method == "unit") && (qoi =="att") && (!is.null(estimator) && estimator == "fd"))) {
            nz.obs <- sum(as.numeric(data$W.it !=0))
            if (verbose)
                cat("\nTotal number observations with non-zero weight:", nz.obs,"\n")
            flush.console()
        }
        
        


### Demean based on the weights
        
        wdm.Data <- mf.sorted
        
        
        ## Add Weight and unit index
        wdm.Data$W <- data$W.it
        wdm.Data$unit <- data$u.index
        wdm.Data$time <- data$t.index
        wdm.Data <- as.data.frame(wdm.Data)

        ## data for traditional fixed effect
        dm.Data <- wdm.Data

        ## excluding zero weights data for weighted fixed effect
        
        nc <- ncol(wdm.Data)

        ## unique number of units and time with positive weights
        if (method == "unit"){
            unit.number.pw <- length(unique(wdm.Data$unit))
        } else if (method == "time"){
            unit.number.pw <- length(unique(wdm.Data$time))
        } else {
            stop("method should be either unit or time")
        }


        ## ## Creating the matrix for final analysis
        ## Data.dm <- as.data.frame(matrix(NA, ncol = nc-3, nrow = nrow(wdm.Data)))
        ## Data.wdm <- as.data.frame(matrix(NA, ncol = nc-3, nrow = nrow(wdm.Data)))

        ## Creating the matrix for final analysis
        Data.dm <- dm.Data
        Data.wdm <- wdm.Data


        ## colume names
        colnames(Data.dm) <- colnames(wdm.Data)[1:nc]
        colnames(Data.wdm) <- colnames(wdm.Data)[1:nc]

        
        for (k in 1:(nc-3)) {
            ## in C
            if (method == "unit"){
                w.wdemean <- WWDemean(wdm.Data[,k], wdm.Data$W, wdm.Data$unit, unit.number.pw, nrow(wdm.Data))
                demean <- Demean(dm.Data[,k], dm.Data$unit, unit.number, nrow(dm.Data))
                Data.wdm[,k] <- as.vector(w.wdemean)
                Data.dm[,k] <- as.vector(demean)
            }
            if (method == "time"){
                w.wdemean <- WWDemean(wdm.Data[,k], wdm.Data$W, wdm.Data$time, unit.number.pw, nrow(wdm.Data))
                demean <- Demean(dm.Data[,k], dm.Data$time, unit.number, nrow(dm.Data))
                Data.wdm[,k] <- as.vector(w.wdemean)
                Data.dm[,k] <- as.vector(demean)
            }

        }

        ## save weighted demeaned dataframe
        if (store.wdm == TRUE){
            Y.wdm <- Data.wdm[,1]
            X.wdm <- Data.wdm[,(2:(nc-3))]
        } else {
            Y.wdm <- NULL
            X.wdm <- NULL
        }
        

        
        ## change formula without intercept
        a <- unlist(strsplit(as.character(formula), "~"))
        formula.ni <- as.formula(paste(a[2], "~ -1 + ",  a[3]))
        ## print(formula.ni)
        
        ## final regression on weighted demeaned data
        fit.final <- lm(formula.ni, data = Data.wdm)
        fit.ols <- lm(formula.ni, data = Data.dm)

        ## brute force matrix calculation
        ## V <- as.matrix(Data.wdm[,2:(nc-3)])
        ## coef2 <- ginv(crossprod(V, V))%*% crossprod(V, Data.wdm[,1])
        
        
        ## residuals
        
        ## u.tilde <- sqrt(wdm.Data$W)*resid(fit.final) # NT x 1
        u.tilde <- resid(fit.final) # NT x 1
        
        ## alternatively
        ## u.tilde <- sqrt(wdm.Data$W) * (Data.wdm[,1] - as.matrix(Data.wdm[,2:(nc-3)])%*%c(fit.final$coef))
        u.hat <- as.matrix(resid(fit.ols))  # NT x 1

        ## saving results
        coef.wls <- fit.final$coef
        coef.ols <- fit.ols$coef
        d.f <- fit.final$df - J.u
        ## sigma2 <- (sum((sqrt(wdm.Data$W)*resid(fit.final))^2))/d.f
        sigma2 <- (sum((resid(fit.final))^2))/d.f

        ## alternatively
        ## sigma2.a <- sum(u.tilde^2)/d.f
        ## cat("compare:", sigma2, sigma2.a, "\n") # same...
        
        
### Robust Standard Errors

        ## contructing de(weighted)-meaned X matrix

        ## in case only one covariate
        if (p == 1){ # one covariate case
            X.tilde <- as.matrix(Data.wdm[,2])
            X.hat <- as.matrix(Data.dm[,2])
        } else {
            X.tilde <- as.matrix(Data.wdm[,2:(1+p)])
            X.hat <- as.matrix(Data.dm[,2:(1+p)])
        }


        ## unit vector (the length should be same as nrow(X.tilde) and length(u.tilde)
        wdm.unit <- as.vector(Data.wdm$unit)
        
### (Robust) standard errors

        ## cat("3: Robust error calculation started\n")
        

        ginv.XX.tilde <- ginv(crossprod(X.tilde, X.tilde))
        ginv.XX.hat <- ginv(crossprod(X.hat, X.hat))

        diag.ee.tilde <- c(u.tilde^2)
        diag.ee.hat <- c(u.hat^2)

        
        if ((hetero.se == TRUE) & (auto.se == TRUE)) {# Default is Arellano
            std.error <- "Heteroscedastic / Autocorrelation Robust Standard Error"

            ## 1. arbitrary autocorrelation as well as heteroskedasticity (Eq 12)

            ## degrees of freedom adjustment
            df.adjust <- 1/(nrow(X.tilde)) * ((nrow(X.tilde)-1)/(nrow(X.tilde)-J.u-p)) * (J.u/(J.u-1))
            
            Omega.hat.HAC <- OmegaHatHAC(nrow(X.tilde), p, wdm.unit, J.u, X.tilde, u.tilde)
            Omega.hat.HAC <- matrix(Omega.hat.HAC, nrow = p, ncol = p)
            ## Omega.hat.HAC <- (1/(nrow(X.tilde)))* Omega.hat.HAC # without degree of freedom adjustment 
            Omega.hat.HAC <- df.adjust * Omega.hat.HAC
            
            Omega.hat.fe.HAC <- OmegaHatHAC(nrow(X.hat), p, wdm.unit, J.u, X.hat, u.hat)
            Omega.hat.fe.HAC <- matrix(Omega.hat.fe.HAC, nrow = p, ncol = p)
            ## Omega.hat.fe.HAC <- (1/(nrow(X.hat))) * Omega.hat.fe.HAC # without degree of freedom adjustment
            Omega.hat.fe.HAC <- df.adjust * Omega.hat.fe.HAC

            ## Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Omega.hat.HAC %*% (nrow(X.tilde) * ginv.XX.tilde)
            ## Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (nrow(X.hat) * ginv.XX.hat)

            Psi.hat.wfe <- (J.u*ginv.XX.tilde) %*% Omega.hat.HAC %*% (J.u*ginv.XX.tilde)
            Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (J.u*ginv.XX.hat)
            
        } else if ( (hetero.se == TRUE) & (auto.se == FALSE) ) {# independence across observations but heteroskedasticity
            std.error <- "Heteroscedastic Robust Standard Error"

            ## 2. independence across observations but heteroskedasticity (Eq 11)
            
            ## Omega.hat.HC <- OmegaHatHC(nrow(X.tilde), p, wdm.unit, J.u, X.tilde, u.tilde)
            ## Omega.hat.HC <- matrix(Omega.hat.HC, nrow = p, ncol = p)
            ## ## same as the following matrix multiplication but slower
            ## ## Omega.hat.he <- t(X.tilde) %*% diag(c(u.tilde)^2, nrow=nrow(X.tilde)) %*% X.tilde
            ## Omega.hat.HC <- (1/(nrow(X.tilde) - J.u - p)) * Omega.hat.HC


            ## Omega.hat.HC <- (1/(nrow(X.tilde) - J.u - p))*(crossprod((X.tilde*diag.ee.tilde), X.tilde)) 
            ## Omega.hat.fe.HC <- (1/(nrow(X.hat) - J.u - p))*(crossprod((X.hat*diag.ee.hat), X.hat))  

            Omega.hat.HC <- (1/J.u)*(crossprod((X.tilde*diag.ee.tilde), X.tilde)) 
            Omega.hat.fe.HC <- (1/J.u)*(crossprod((X.hat*diag.ee.hat), X.hat))  
            
            
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
                    ## Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Sigma_HRFE %*% (nrow(X.tilde) * ginv.XX.tilde)
                    Psi.hat.wfe <- (J.u*ginv.XX.tilde) %*% Sigma_HRFE %*% (J.u*ginv.XX.tilde)

                    Sigma2_HRFE <- ((J.t-1)/(J.t-2))*(Omega.hat.fe.HC - (1/(J.t-1))*B2.hat)
                    ## Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Sigma2_HRFE %*% (nrow(X.hat) * ginv.XX.hat)
                    Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Sigma2_HRFE %*% (J.u*ginv.XX.hat)
                    
                    
                } else { 
                    stop ("unbiased.se == TRUE is allowed only when panel is balanced")
                }
            }
            
            ## Psi.hat.wfe <- (nrow(X.tilde) * ginv.XX.tilde) %*% Omega.hat.HC %*% (nrow(X.tilde) * ginv.XX.tilde)
            ## Psi.hat.fe <- (nrow(X.hat) * ginv.XX.hat) %*% Omega.hat.fe.HC %*% (nrow(X.hat) * ginv.XX.hat)

            Psi.hat.wfe <- (J.u*ginv.XX.tilde) %*% Omega.hat.HC %*% (J.u*ginv.XX.tilde)
            Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Omega.hat.fe.HC %*% (J.u*ginv.XX.hat)
            

        } else if ( (hetero.se == FALSE) & (auto.se == FALSE) ) {# indepdence and homoskedasticity

            stop("standard errors with independence and homoskedasticity is not supported")

            ## std.error <- "Homoskedastic Standard Error"
            
            ## ## Psi.hat.wfe <- nrow(X.tilde) * (sigma2 * ginv.XX.tilde)
            ## Psi.hat.wfe <- J.u * (sigma2 * ginv.XX.tilde)            
            
            ## ## Psi.hat.fe <- nrow(X.hat) * vcov(fit.ols)
            ## Psi.hat.fe <- J.u * vcov(fit.ols)            

        } else if ( (hetero.se == FALSE) & (auto.se == TRUE) ) {# Kiefer

            stop ("robust standard errors with autocorrelation and homoskedasiticy is not supported")

        }

        ## var.cov <- Psi.hat.wfe * (1/nrow(X.tilde))
        ## var.cov.fe <- Psi.hat.fe *(1/nrow(X.hat))

        var.cov <- Psi.hat.wfe * (1/J.u)
        var.cov.fe <- Psi.hat.fe *(1/J.u)
        
### traditional one way fixed effect results

        
        ## if (verbose) {
        ##   cat("Traditional one-way fixed effect\n")
        ##   print(summary(fit.ols))
        ##   cat("Robust Standard errors for Standard FE \n")
        ##   print(sqrt(diag(var.cov.fe)))
        ##   flush.console()
        
        ## }
        

### White (1980) Test: Theorem 4

        diag.ee <- c(u.hat) * c(u.tilde)
        
        Lambda.hat1 <-  1/((nrow(X.hat) - J.u - p))* (crossprod((X.hat*diag.ee), X.tilde))  
        Lambda.hat2 <-  1/((nrow(X.tilde) - J.u - p))* (crossprod((X.tilde*diag.ee), X.hat))  


        Phi.hat <- Psi.hat.wfe + Psi.hat.fe - (nrow(X.hat)*ginv.XX.hat) %*% Lambda.hat1 %*% (nrow(X.tilde)*ginv.XX.tilde) - (nrow(X.tilde)*ginv.XX.tilde) %*% Lambda.hat2 %*% (nrow(X.hat)*ginv.XX.hat)

        rm(Lambda.hat1, Lambda.hat2)
        gc()

        ## White test: null hypothesis is ``no misspecification''
        white.stat <- nrow(X.hat) * t(coef.ols - coef.wls) %*% ginv(Phi.hat) %*% (coef.ols - coef.wls)
        test.null <- pchisq(as.numeric(white.stat), df=p, lower.tail=F) < White.alpha
        white.p <- pchisq(as.numeric(white.stat), df=p, lower.tail=F)

        
        flush.console()
        
        
        ## ## compare with sandwich standard error (essentially same)

        ## ## cat("sandwich package regression se HC:", print(sqrt(diag((vcovHC(fit.traditional, type="HC")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC0:", print(sqrt(diag((vcovHC(fit.traditional, type="HC0")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC1:", print(sqrt(diag((vcovHC(fit.traditional, type="HC1")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC2:", print(sqrt(diag((vcovHC(fit.traditional, type="HC2")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC3:", print(sqrt(diag((vcovHC(fit.traditional, type="HC3")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC4:", print(sqrt(diag((vcovHC(fit.traditional, type="HC4")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC4m:", print(sqrt(diag((vcovHC(fit.traditional, type="HC4m")[1:p,1:p])))), "\n")
        ## ## cat("sandwich package regression se HC5:", print(sqrt(diag((vcovHC(fit.traditional, type="HC5")[1:p,1:p])))), "\n")
        

        ## Creating a weight verctor
        ## original index
        idx <- paste(orig.unit.idx, orig.time.idx, sep="_")
        a <- units
        b <- times
        Wv <- as.vector(W) # as vector
        a1 <- rep(as.character(a$unit), each=nrow(W))
        b1 <- rep(as.character(b$time), ncol(W))
        idxall <- paste(a1, b1, sep="_")
        idxall.sub <- idxall[which(idxall %in% idx)]
        W.it <- Wv[which(idxall %in% idx)]
        u.sub <- unlist(lapply(idxall.sub,
                               function(x) strsplit(x, "_")[[1]][1]))
        t.sub <- unlist(lapply(idxall.sub,
                               function(x) strsplit(x, "_")[[1]][2]))
        cmd <- paste("W1 <- data.frame(", unit.index, "= u.sub)", sep="")
        eval(parse(text=cmd))
        if(is.null(time.index)){
            W1$obs.idx <- t.sub
        } else {
            cmd2 <- paste("W1$", time.index, " <- t.sub", sep="")
            eval(parse(text=cmd2))
        }
        W1$W.it <- W.it

        ## ensuring the order reflects the original idx
        mf$W.it <- W.it[match(idxall.sub, idx)]
        u.orig <- unlist(lapply(idx,
                                function(x) strsplit(x, "_")[[1]][1]))
        t.orig <- unlist(lapply(idx,
                                function(x) strsplit(x, "_")[[1]][2]))
        cmd <- paste("D <- data.frame(", unit.index, "= u.orig)", sep="")
        eval(parse(text=cmd))
        if(is.null(time.index)){
            D$obs.idx <- t.orig
        } else {
            cmd2 <- paste("D$", time.index, " <- t.orig", sep="")
            eval(parse(text=cmd2))
        }
        
        mf <- cbind(D, mf)
        
        Num.nonzero <- length(which(!W1$weights==0))
        
###Saving results


        z <- list(coefficients = coef.wls,
                  x = x[,-1,drop=FALSE],
                  y = y,
                  mf = mf,
                  call = wfe.call,
                  vcov = var.cov,
                  se = sqrt(diag(var.cov)),
                  sigma = sqrt(sigma2),
                  df = d.f,
                  residuals = y - (x[,-1,drop=FALSE] %*% coef.wls),
                  W = W1,
                  Num.nonzero = Num.nonzero,
                  uniq.n.units = J.u,
                  units = units,
                  times = times,
                  method = method,
                  causal = causal,
                  est = est,
                  std.error = std.error,
                  White.pvalue = white.p,
                  White.alpha = White.alpha,
                  White.stat = white.stat,
                  White.test = test.null,
                  Y.wdm = Y.wdm,
                  X.wdm = X.wdm)
        class(z) <- "wfe"
        z

### ********************************************************
### Two-way Weighted Fixed Effects
### ********************************************************
    } else {
        ## if (verbose) {
        ##     did <- CalDID(data$u.index, data$t.index, data$TR, data$C.it,
        ##                   y, tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, verbose)
        ##     cat("\nMulti-period DID estimate with no covariate adjustments is", did ,"\n")
        ##     flush.console()
        ## }
        
        ## Differences-in-difference
        if(( (method=="unit") & (qoi == "ate") & (!is.null(estimator) & estimator == "did")) |
           ( (method == "unit") & (qoi =="att") & (!is.null(estimator) & estimator == "did")) |
           ( (method=="unit") & (qoi == "ate") & (!is.null(estimator) & estimator == "Mdid")) |
           ( (method == "unit") & (qoi =="att") & (!is.null(estimator) & estimator == "Mdid"))       
           ) {

            method <- "Weighted Two-way"
            ## Standard Fixed effect
            if (unweighted == TRUE) {
                data$W.it <- rep(1, nrow(data))
                W <- matrix(1,nrow=length(uniq.t), ncol=length(uniq.u))
            } else {
                if (verbose) {
                    cat("\nWeight calculation started ")
                    flush.console()
                }
                if(estimator == "Mdid"){
                    if(is.null(maxdev.did)){
                        maxdev.did <- -1
                        if (verbose) {
                            cat(": Nearest Neighbor Matching\n")
                            flush.console()
                        }
                        
                    } else {
                        if (verbose) {
                            cat(": Matching on Pre-Treatment Outcome Within Maximum Deviation", maxdev.did,"\n")
                            flush.console()
                        }
                        
                        maxdev.did <- as.numeric(maxdev.did)
                    }
                    WDiD <- GenWeightsMDID(data$u.index, data$t.index, data$TR, data$C.it, y, maxdev.did,
                                           tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, verbose)

                } else {
                    WDiD <- GenWeightsDID(data$u.index, data$t.index, data$TR, data$C.it,
                                          tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, verbose)
                }
                
                W <- matrix(WDiD, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)            
                data$W.it <- VectorizeC(as.matrix(W), data$t.index, data$u.index, tn.row)
                
                if (verbose) { 
                    cat("\nWeight calculation done \n")
                    flush.console()
                }
                
            }
            ## e <- environment()
            ## save(file = "temp.RData", list = ls(), env = e)
            
            ## creating index for sparse dummy matrix

            u <- as.matrix(table(data$u.index))
            
            Udummy.i <- seq(1:sum(u))
            Udummy.j <- c()

            for (j in 1:length(uniq.u)) {
                Udummy.j <- c(Udummy.j, rep(j, u[j,1]))
            }
            Udummy  <- sparseMatrix(x=1, i=Udummy.i, j=Udummy.j)

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
                Tdummy <- rBind(Tdummy,  Diagonal(x = t[j,], n=length(uniq.t)))
            }

            ## when panel is unbalanced there will be rows of zeros
            zero <- which(apply(Tdummy, 1, mean) == 0)
            if(length(zero) > 0) {
                Tdummy <- Tdummy[-zero,]
            }

            ## if (verbose)
            ##   cat("\n Dummy creation done \n")
            ## flush.console()


#########################################################################
### Projection for standard twoway FE

            ## if (verbose)
            ##   cat("\n Standard FE Projection Started \n")
            ## flush.console()

            ## this step takes time
            P1 <- Udummy %*% tcrossprod(Diagonal(x=1/as.vector(table(data$u.index))), Udummy)
            
            ## e <- environment()
            ## save(file = "test.RData", list = ls(), env = e)
            
            Q1 <- Diagonal(x=1, n=nrow(Udummy)) - P1
            

            ## Q: not too sparse
            Q <- Q1 %*% Tdummy
            Q.QQginv <- Q %*% ginv(as.matrix(crossprod(Q)))

            X <- as.matrix(X)
            Y <- as.matrix(data$y)
            
            YX <- cbind(Y,X)
            
            Data.2wdm <- as.data.frame(as.matrix(YX - P1%*%YX - Q.QQginv %*% crossprod(Q,YX)))
            colnames(Data.2wdm) <- colnames(mf.sorted)

            rm(YX, P1, Q.QQginv, Q)
            gc()

            ## if (verbose)
            ##   cat("\n Standard FE Projection done \n")
            ## flush.console()
            
            ## e <- environment()
            ## save(file = "temp2.RData", list = ls(), env = e)
            
            
            
            a <- unlist(strsplit(as.character(formula), "~"))
            formula.ni <- as.formula(paste(a[2], "~ -1 + ",  a[3]))


            ## final regression on 2way demeaned data
            fit.ols <- lm(formula.ni, data = Data.2wdm)
            
            coef.ols <- fit.ols$coef
            resid.ols <- resid(fit.ols)

            u.hat <- as.matrix(resid.ols)
            X.hat <- as.matrix(Data.2wdm[,-1])
            rm(Data.2wdm)
            gc()
            


            
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
            if (verbose)
                cat("\nTotal number observations with non-zero weight:", nz.obs,"\n")
            flush.console()
            

            
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
            
            ## if (verbose)
            ##   cat("\n Udummy done\n")
            ## flush.console()


            ## 2. Time dummies

            if (length(which(as(apply(Tdummy, 2, sum), "sparseVector")==0)) > 0) {
                zero <- which(as(apply(Tdummy, 2, sum), "sparseVector") == 0)
                n.zero <- length(zero)
                Tdummy <- Tdummy[,-zero]
                ## delete last column
                ## last <- ncol(Tdummy) 
                ## Tdummy <- Tdummy[,-last]
                ## adding a column of 1000's for numerical stability of ginv
                Tdummy <- cBind(Tdummy, rep(1000, nrow(Tdummy)))
                gc()
            } else {
                last <- ncol(Tdummy)
                Tdummy <- Tdummy[,-last] # for identification exclude the last year dummy
                ## adding a column of 1's for numerical stability of ginv
                Tdummy <- cBind(Tdummy, rep(1000, nrow(Tdummy)))
                gc()
            }

            ## if (verbose)
            ##   cat("\n Tdummy done\n")
            ## flush.console()

            ## number of columns for Unit/Time dummy matrix
            n.Udummy <- ncol(Udummy)
            n.Tdummy <- ncol(Tdummy)

            
            ## combining unit dummy matrix and X matrix
            D <- Matrix(cBind(Udummy, Tdummy))
            ## Note: Tdummy part does not have zero columns, but Udummy part has it
            ## will be addressed thie issue below by n.zero
            
            ## final number of dummies
            fn.dummies <- ncol(D) # final number of dummies



### Projection onto Complex-plane

            
            ## if (verbose)
            ##   cat("\n Calculation for Projection Matrix Started \n")
            ## flush.console()

            ## sqrt of weights (imaginary numbers)
            Im.sqrt <- function(weights) {
                if (weights >= 0) # real number
                    im.w <- complex(real=sqrt(weights), imaginary =0)
                if (weights < 0) # imaginary number
                    im.w <- complex(real=0, imaginary = sqrt(-weights))
                invisible(im.w)
            }

            ## vector of sqrt(W.it)
            w.sqrt <- sapply(data$W.it[nz.index], Im.sqrt)
            

#########################################################################

            
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

            
### moving the first column of D1 to D2
            
            
            ## ## index for "zero dummy columns" after deleting zero-weights observations
            ## u.zero <- try(which(as(apply(Udummy, 2, sum), "sparseVector") == 0), silent=TRUE)
            ## if ((class(u.zero) == "try-error") & (White == TRUE)) {
            ##   stop ("Insufficient memory. White = FALSE option is needed")
            ## }
            
            ## rm(Udummy)
            ## gc()
            
            ## ## unit index for which weights sum to zero
            ## sum.0.u.index <- which(round(apply(W, 2, sum), digits=10) == 0)

            ## ## unit index for which weights do not sum to zero
            ## sum.no0.u.index <- which(round(apply(W, 2, sum), digits=10) !=0)

            ## ## D1 index: unit index with 1)0-columns after deleting zero weights observation should be excluded, and 2) weights sum to zero
            ## temp <- 1:J.u
            ## ## temp2 <- 1:(J.u+J.t-1)

            ## ## D1.index <- which((!temp %in% u.zero) & (!temp %in% sum.0.u.index))
            ## D1.index <- which((!temp %in% u.zero))
            
            
            ## ## attaching part of D1 to D2 (units weights sum to zero)
            ## ## unit index that should be moved from D1 to D2
            
            ## ## if (length(which(!temp%in%u.zero & temp%in%sum.0.u.index)) > 0) {
            ## ##   toD2.index <- which(!temp%in%u.zero & temp%in%sum.0.u.index)
            ## ##   D2.index <- c(toD2.index, (n.Udummy+1):(n.Udummy+n.Tdummy))
            ## ## } else {
            ## ##   D2.index <- c((n.Udummy+1):(n.Udummy+n.Tdummy))
            ## ## }
            
            ## D2.index <- c((n.Udummy+1):(n.Udummy+n.Tdummy))
            
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
            
            ## D1.starL[[1]] <- drop0(D.starL[[1]][,D1.index])
            ## D1.starL[[2]] <- drop0(D.starL[[2]][,D1.index])
            ## D2.starL[[1]] <- drop0(D.starL[[1]][,D2.index])
            ## D2.starL[[2]] <- drop0(D.starL[[2]][,D2.index])
            ## cat("Number of columns for D1", length(D1.index), "\n")
            ## cat("Number of columns for D2", length(D2.index), "\n")

            D1.starL[[1]] <- drop0(D.starL[[1]][,1:n.Udummy])
            D1.starL[[2]] <- drop0(D.starL[[2]][,1:n.Udummy])
            D2.starL[[1]] <- drop0(D.starL[[1]][,((n.Udummy+1):(n.Udummy+n.Tdummy))])
            D2.starL[[2]] <- drop0(D.starL[[2]][,((n.Udummy+1):(n.Udummy+n.Tdummy))])
            ## cat("Number of columns for D1", n.Udummy, "\n")
            ## cat("Number of columns for D2", n.Tdummy, "\n")


            

            ## D1.index.full <- 1:ncol(D.starL[[1]]) %in% D1.index
            ## D2.index.full <- 1:ncol(D.starL[[1]]) %in% D2.index
            ## D1.starL[[1]] <- drop0(D.starL[[1]][,D1.index.full])
            ## D1.starL[[2]] <- drop0(D.starL[[2]][,D1.index.full])
            ## D2.starL[[1]] <- drop0(D.starL[[1]][,D2.index.full])
            ## D2.starL[[2]] <- drop0(D.starL[[2]][,D2.index.full])

            
            rm(D.starL)
            gc()
#########################################################################

            
            
            ## sum of sqrt(W.it) across years for each dyad
            general.inv <- function(weight, tol = tol){
                if(abs(weight) < tol) {
                    out <- 0
                } else {
                    out <- 1/weight
                }
                out
            }
            
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

            ## if (verbose) {
            ##   cat("\n P1 created\n")
            ##   flush.console()      
            ## }
            
            Q1L <- list()
            Q1L[[1]] <- drop0(Diagonal(x=1, n=nrow(P1L[[1]])) - P1L[[1]])
            Q1L[[2]] <- drop0(Diagonal(x=0, n=nrow(P1L[[1]])) - P1L[[2]])

            ## if (verbose) {
            ##   cat("\n Q1 created\n")
            ##   flush.console()
            ## }
            
            
            Q <- try(Sparse_compMatrixMultiply(Q1L[[1]], Q1L[[2]], D2.starL[[1]], D2.starL[[2]]), silent = TRUE)
            if ((class(Q) == "try-error") & (White == TRUE)) {
                stop ("Insufficient memory. White = FALSE option is needed")
            }
            rm(Q1L, D2.starL)
            gc()

            
            Q.matrix <- matrix(complex(real=as.matrix(Q[[1]]), imaginary=as.matrix(Q[[2]])), nrow=nrow(Q[[1]]))


            QQ.inv <- list()
            QQ.inv[[1]] <- drop0(Matrix(Re(ginv(crossprod(Q.matrix)))))
            QQ.inv[[2]] <- drop0(Matrix(Im(ginv(crossprod(Q.matrix)))))
            rm(Q.matrix)
            gc()

            PL <- list()

            Q.QQinv <- Sparse_compMatrixMultiply(Q[[1]], Q[[2]], QQ.inv[[1]], QQ.inv[[2]]) 
            rm(QQ.inv)
            gc()
            
            ## if (verbose) {
            ##   cat("\n Q.QQinv created\n")
            ##   flush.console()
            ## }

#########################################################################
### Fast Projection in R


            YX.starL <- list()
            YX.starL[[1]] <- cBind(y.starL[[1]], x.starL[[1]])
            YX.starL[[2]] <- cBind(y.starL[[2]], x.starL[[2]])
            rm(y.starL, x.starL)
            gc()
            
            
            P1.YX <- Sparse_compMatrixMultiply(P1L[[1]], P1L[[2]], YX.starL[[1]], YX.starL[[2]])
            rm(P1L)
            gc()


            Q.YX <- Sparse_compMatrix_crossprod(Q[[1]], Q[[2]], YX.starL[[1]], YX.starL[[2]])
            rm(Q)
            gc()
            
            QQQQ.YX <- Sparse_compMatrixMultiply(Q.QQinv[[1]], Q.QQinv[[2]], Q.YX[[1]], Q.YX[[2]])
            
            ## cat("dimension of P1.YX:", dim(P1.YX[[1]]), "\n")
            ## cat("dimension of QQQQ.YX:", dim(QQQQ.YX[[1]]), "\n")

            ## TransformedL <- list()
            ## TransformedL[[1]] <- YX.starL[[1]] - P1.YX[[1]] - P2.YX[[1]]
            ## TransformedL[[2]] <- YX.starL[[2]] - P1.YX[[2]] - P2.YX[[2]]

            Transformed <- matrix(complex(real=as.vector(YX.starL[[1]] - P1.YX[[1]] - QQQQ.YX[[1]]), imaginary=as.vector(YX.starL[[2]] - P1.YX[[2]] - QQQQ.YX[[2]])), nrow=tn.row)
            rm(YX.starL, P1.YX, QQQQ.YX)
            
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

            
            ginv.XX.tilde <- ginv(crossprod(X.tilde))
            betaT <- ginv.XX.tilde%*% crossprod(X.tilde, y.tilde)
            if (length(betaT) == 1) {
                colnames(betaT) <- a[3]
            }
            ## print(betaT)
            coef.wls <- matrix(as.double(Re(betaT)))
            rownames(coef.wls) <- colnames(X.tilde)
            

            ## e <- environment()
            ## save(file = "temp.RData", list = ls(), env = e)
            

            ## heteroskedasticity robust standard errors

            ## if (verbose) {
            ##   cat("\nStd.error calculation start\n")
            ##   Sys.time()
            ##   flush.console()
            ## }
            

            ## #######################################################################            
            ## (Robust) standard errors (GMM asymptotic variance for wfe)
            ## calculating GMM standard errors
            ## #######################################################################

            ## weighted residuals
            e.tilde <- (y.tilde - X.tilde %*% betaT)

            ## true residuals
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
            

            ## check residuals      
            ## print(cbind(as.matrix(true.resid), e.tilde, data$W.it))
            ## print(cbind(sum(true.resid^2), sum(e.tilde^2)))


            ## diag.ee.tilde <- diag(tcrossprod(e.tilde,e.tilde))
            ## diag.resid <- as.vector(resid * resid)
            diag.ee.tilde <- as.vector(e.tilde * e.tilde)
            
#########################################################################

            ## cat("dimension of X.tilde:", dim(X.tilde), "\n")
            
            ## XX.hat <- crossprod(X.hat, X.hat)
            ginv.XX.hat <- ginv(crossprod(X.hat, X.hat))
            ## d.f <- length(y.tilde) - n.Udummy - n.Tdummy - dim(X.tilde)[2]
            d.f <- length(y.tilde)

            ## cat("Sum of squared residuals:", sum(resid^2), "\n")
            sigma2 <- as.double(Re(sum(resid^2)/d.f))
            
            ## cat("sigma2", sigma2, "\n")


            ## e <- environment()
            ## save(file = "temp.RData", list = ls(), env = e)

            ## two-way WFE robust standard errors calculation

            if ((hetero.se == TRUE) & (auto.se == TRUE)){

                ## 1. arbitrary autocorrelation as well as heteroskedasticity (Eq 12)
                std.error <- "Heteroscedastic / Autocorrelation Robust Standard Error"
                ## stop ("Robust standard errors with autocorrelation is currently not supported")

                ## Remove observations with zero weights
                zero.ind <- which(data$W.it==0)
                if(length(zero.ind)>0){
                    data <- data[-zero.ind, ]
                }
                Nnonzero <- nrow(data)
                
                ## Demean data
                ## -----------------------------------------------------
                ## 2way demean variables 
                ## -----------------------------------------------------

                x.vars <- colnames(x)
                x.vars <- x.vars[-grep("Intercept", x.vars)]
                variables <- c("y", x.vars)

                DemeanedMatrix <- matrix(NA, nrow=nrow(data), ncol=length(variables))
                colnames(DemeanedMatrix) <- variables
                year.counts <- as.numeric(table(data$unit))
                unit.counts <- as.numeric(table(data$time))
                obs.counts <- nrow(data)
                
                for(k in 1:length(variables)){
                    v <- variables[k]
                    cmd1 <- paste("demean.unit <- tapply(data$", v, ", as.factor(data$u.index), mean, na.rm=T)", sep="")
                    cmd2 <- paste("demean.time <- tapply(data$", v, ", as.factor(data$t.index), mean, na.rm=T)", sep="")
                    cmd3 <- paste("demean.all <- mean(data$", v, ", na.rm=T)", sep="")
                    cmd4 <- paste("demean.units <- demean.unit[data$u.index]", sep="")
                    cmd5 <- paste("demean.times <- demean.time[data$t.index]", sep="")
                    cmd6 <- paste("demean.alls <- rep(demean.all, times=obs.counts)", sep="")
                    cmd7 <- paste("DemeanedMatrix[,k] <- data$", v, "- demean.units - demean.times + demean.alls", sep="")

                    eval(parse(text=cmd1))
                    eval(parse(text=cmd2))
                    eval(parse(text=cmd3))
                    eval(parse(text=cmd4))
                    eval(parse(text=cmd5))
                    eval(parse(text=cmd6))
                    eval(parse(text=cmd7))
                }
                ## verify: should return the same results (checked!)
                ## lm(y~ as.factor(u.index) + as.factor(t.index) + tr+x1+x2, data=data)
                ## lm(DemeanedMatrix[,1]~ -1 + DemeanedMatrix[,-1])

                ## standard error calculation
                n.units <- length(unique(data$u.index))
                U <- matrix(0, nrow=length(x.vars), ncol=length(x.vars))
                V <- matrix(0, nrow=length(x.vars), ncol=length(x.vars))
                Beta <- as.matrix(coef.wls)
                
                for(g in 1:length(n.units)){
                    Y.dm <- DemeanedMatrix[which(data$u.index==g),1]
                    X.dm <- DemeanedMatrix[which(data$u.index==g),-1]
                    W.diag <- diag(data$W.it[data$u.index==g])

                    U.i <- t(X.dm) %*% W.diag %*% X.dm
                    U <- U + U.i
                    V.i.tmp <- t(X.dm) %*% W.diag %*% (Y.dm-X.dm %*% Beta) %*% t(Y.dm-X.dm %*% Beta) %*% W.diag %*% X.dm
                    V.i <- V.i.tmp %*% t(V.i.tmp)
                    V <- V + V.i
                }

                ## asymptotic variance using Methods of Moments
                inv.U <- solve(1/Nnonzero * U)
                V <- 1/Nnonzero * V
                Psi.hat.wfe <- inv.U %*% V %*% inv.U

                ## -----------------------------------------------------
                ## vcov matrix for FE for White statistics calculation
                ## -----------------------------------------------------


                Omega.hat.fe.HAC <- OmegaHatHAC(nrow(X.hat), ncol(X.hat), data$u.index, J.u, X.hat, u.hat)
                Omega.hat.fe.HAC <- matrix(Omega.hat.fe.HAC, nrow = ncol(X.hat), ncol = ncol(X.hat))
                ## Omega.hat.fe.HAC <- (1/(nrow(X.hat)-J.u-J.t-p)) * Omega.hat.fe.HAC
                Omega.hat.fe.HAC <- (1/J.u) * Omega.hat.fe.HAC

                ## Psi.hat.fe <- (nrow(X.hat)*ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (nrow(X.hat)*ginv.XX.hat)
                Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (J.u*ginv.XX.hat)
                ## garbage collection
                rm(Omega.hat.fe.HAC)

                
                ## -----------------------------------------------------
                ## old code 
                ## -----------------------------------------------------
                
                ## ## degrees of freedom adjustment
                ## ## cat("degrees of freedom:", J.u, J.t, p, "\n")
                ## df.adjust <- 1/(nrow(X.tilde)) * ((nrow(X.tilde)-1)/(nrow(X.tilde)-J.u-J.t-p+1)) * (J.u/(J.u-1))
                
                ## Omega.hat.HAC <- as.double(comp_OmegaHAC(c(X.tilde), e.tilde, c(X.tilde), e.tilde, dim(X.tilde)[1], dim(X.tilde)[2], data$u.index, J.u))
                ## Omega.hat.HAC <- matrix(Omega.hat.HAC, nrow=ncol(X.tilde), ncol=ncol(X.tilde), byrow=T)
                ## ## Omega.hat.HAC <- (1/(nrow(X.tilde)-J.u-J.t-p+1))* Omega.hat.HAC
                ## ## Omega.hat.HAC <- (1/(nrow(X.tilde)))* Omega.hat.HAC 
                ## ## Omega.hat.HAC <- df.adjust * Omega.hat.HAC
                ## Omega.hat.HAC <- (1/J.u) * Omega.hat.HAC


                ## ## check positive definiteness of Omega.hat.HAC
                ## if ( sum(as.numeric(eigen(Omega.hat.HAC)$values < 0)) > 0 ) {
                ##     ## cat ("*** Omega.hat is not positive definite ***\n")                    
                ##     stop ("*** Omega.hat is not positive definite ***")

                ## }
                
                ## Omega.hat.fe.HAC <- OmegaHatHAC(nrow(X.hat), ncol(X.hat), data$u.index, J.u, X.hat, u.hat)
                ## Omega.hat.fe.HAC <- matrix(Omega.hat.fe.HAC, nrow = ncol(X.hat), ncol = ncol(X.hat))
                ## ## Omega.hat.fe.HAC <- (1/(nrow(X.hat)-J.u-J.t-p)) * Omega.hat.fe.HAC
                ## Omega.hat.fe.HAC <- (1/J.u) * Omega.hat.fe.HAC
                
                ## ## Psi.hat.wfe <- ((nrow(X.tilde))*ginv.XX.tilde) %*% Omega.hat.HAC %*% ((nrow(X.tilde))*ginv.XX.tilde)
                ## ## Psi.hat.fe <- (nrow(X.hat)*ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (nrow(X.hat)*ginv.XX.hat)
                ## ## Psi.hat.wfe <- (J.u*ginv.XX.tilde) %*% Omega.hat.HAC %*% (J.u*ginv.XX.tilde)
                ## Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Omega.hat.fe.HAC %*% (J.u*ginv.XX.hat)
                
                ## ## garbage collection
                ## rm(Omega.hat.HAC, Omega.hat.fe.HAC)
                ## gc()

            } else if ( (hetero.se == TRUE) & (auto.se == FALSE)) {
                stop("Please set hetero.se == TRUE & auto.se == TRUE when you run two-way FE")
                
                ## 2. independence across observations but heteroskedasticity (Eq 11)
                
                std.error <- "Heteroscedastic Robust Standard Error"

                Omega.hat.HC <- as.double(comp_OmegaHC(c(X.tilde), e.tilde, c(X.tilde), e.tilde, dim(X.tilde)[1], dim(X.tilde)[2], data$u.index, J.u))
                Omega.hat.HC <- matrix(Omega.hat.HC, nrow=ncol(X.tilde), ncol=ncol(X.tilde), byrow=T)
                ## Omega.hat.HC <- (1/(nrow(X.tilde)-J.u-J.t-p+1))* Omega.hat.HC
                Omega.hat.HC <- (1/J.u)* Omega.hat.HC                

                ## Psi.hat.wfe <- ((nrow(X.tilde))*ginv.XX.tilde) %*% Omega.hat.HC %*% ((nrow(X.tilde))*ginv.XX.tilde)
                Psi.hat.wfe <- (J.u*ginv.XX.tilde) %*% Omega.hat.HC %*% (J.u*ginv.XX.tilde)
                
                ## ## alternatively calculation (don't need to invert)
                ## Psi.hat.wfe2 <- (length(y.tilde)*ginv.XX.tilde) %*% ( (1/length(y.tilde)) * (crossprod((X.tilde*diag.ee.tilde), X.tilde)) ) %*% ((length(y.tilde))*ginv.XX.tilde)

                
                ## Omega.hat for FE

                Omega.hat.fe.he <- OmegaHatHC(nrow(X.hat), ncol(X.hat), data$u.index, J.u, X.hat, u.hat)
                Omega.hat.fe.he <- matrix(Omega.hat.fe.he, nrow = ncol(X.hat), ncol = ncol(X.hat))
                ## Omega.hat.fe.he <- (1/(nrow(X.hat)-J.u-J.t-p+1)) * Omega.hat.fe.he
                Omega.hat.fe.he <- (1/J.u) * Omega.hat.fe.he                

                ## Psi.hat.fe <- (nrow(X.hat)*ginv.XX.hat) %*% Omega.hat.fe.he %*% (nrow(X.hat)*ginv.XX.hat)
                Psi.hat.fe <- (J.u*ginv.XX.hat) %*% Omega.hat.fe.he %*% (J.u*ginv.XX.hat)                

                ## Same as the following matrix multiplication
                ## Psi.hat.fe <- (solve(XX.hat) %*% (1/d.f *(t(X.hat) %*% diag(diag(u.hat %*% t(u.hat))) %*% X.hat)) %*% solve(XX.hat))

                ## garbage collection
                rm(Omega.hat.HC, Omega.hat.fe.he)
                gc()
                

            } else if ( (hetero.se == FALSE) & (auto.se == FALSE) ) {# indepdence and homoskedasticity

                stop("Please set hetero.se == TRUE & auto.se == TRUE when you run two-way FE")

                ## std.error <- "Homoskedastic Standard Error"
                
                ## ## Psi.hat.wfe <- sigma2 * ( (length(y.tilde)*ginv.XX.tilde) %*% ((1/(length(y.tilde)-J.u-J.t-p+1))*(crossprod((X.tilde*data$W.it[nz.index]), X.tilde))) %*% (length(y.tilde)*ginv.XX.tilde) ) 
                ## ## Psi.hat.fe <- (nrow(X.hat)) * vcov(fit.ols)

                ## Psi.hat.wfe <- sigma2 * ( (J.u*ginv.XX.tilde) %*% ((1/J.u)*(crossprod((X.tilde*data$W.it[nz.index]), X.tilde))) %*% (J.u*ginv.XX.tilde) ) 
                ## Psi.hat.fe <- (J.u) * vcov(fit.ols)
                
            } else if ( (hetero.se == FALSE) & (auto.se == TRUE) ) {# Kiefer
                stop ("Robust standard errors with autocorrelation and homoskedasiticy is not supported")
            }


            ## vcov of wfe model
            ## vcov.wfe <- Psi.hat.wfe * (1/nrow(X.tilde))
            vcov.wfe <- Psi.hat.wfe * (1/Nnonzero)            
            ## cat("dimension of vcov:", dim(vcov.wfe), "\n")
            se.did <- as.double(Re(sqrt(diag(vcov.wfe))))
            
            ## ## check vcov of wfe model
            ## vcov.wfe2 <- Psi.hat.wfe2 * (1/nrow(X.tilde))
            ## se.did2 <- as.double(Re(sqrt(diag(vcov.wfe2))))

            ## cat("\nStd.errors for wfe:", se.did, se.did2, "\n")
            
            
            ## vcov of standard fe model (note:already divided by J.u)
            ## var.cov.fe <- Psi.hat.fe * (1/nrow(X.hat))
            var.cov.fe <- Psi.hat.fe * (1/J.u)            
            se.ols <- sqrt(diag(var.cov.fe))


            
            ## if (verbose) {
            ##   cat("\nStd.error calculation done")
            ##   flush.console()
            ## }

            
            ## e <- environment()
            ## save(file = "temp.RData", list = ls(), env = e)
            

### traditional one way fixed effect results

            
            ## if (verbose) {
            ##   cat("Traditional two-way fixed effect\n")
            ##   print(summary(fit.ols))
            ##   cat("Robust Standard errors for Standard FE \n")
            ##   print(se.ols)
            ##   flush.console()
            ## }
            
            

### White (1980) Test: Theorem 4

            if (White == TRUE){
                
                diag.ee <- c(u.hat) * c(e.tilde)
                
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
            


            ## Creating a weight verctor
            ## original index
            idx <- paste(orig.unit.idx, orig.time.idx, sep="_")
            a <- units
            b <- times
            Wv <- as.vector(W) # as vector
            a1 <- rep(as.character(a$unit), each=nrow(W))
            b1 <- rep(as.character(b$time), ncol(W))
            idxall <- paste(a1, b1, sep="_")
            idxall.sub <- idxall[which(idxall %in% idx)]
            W.it <- Wv[which(idxall %in% idx)]
            u.sub <- unlist(lapply(idxall.sub,
                                   function(x) strsplit(x, "_")[[1]][1]))
            t.sub <- unlist(lapply(idxall.sub,
                                   function(x) strsplit(x, "_")[[1]][2]))

            cmd <- paste("W1 <- data.frame(", unit.index, "= u.sub)", sep="")
            eval(parse(text=cmd))
            if(is.null(time.index)){
                W1$obs.idx <- t.sub
            } else {
                cmd2 <- paste("W1$", time.index, " <- t.sub", sep="")
                eval(parse(text=cmd2))
            }
            W1$W.it <- W.it

            ## ensuring the order reflects the original idx
            mf$W.it <- W.it[match(idxall.sub, idx)]
            u.orig <- unlist(lapply(idx,
                                    function(x) strsplit(x, "_")[[1]][1]))
            t.orig <- unlist(lapply(idx,
                                    function(x) strsplit(x, "_")[[1]][2]))

            cmd <- paste("D <- data.frame(", unit.index, "= u.orig)", sep="")
            eval(parse(text=cmd))
            if(is.null(time.index)){
                D$obs.idx <- t.orig
            } else {
                cmd2 <- paste("D$", time.index, " <- t.orig", sep="")
                eval(parse(text=cmd2))
            }

            mf <- cbind(D, mf)
            Num.nonzero <- length(which(!W1$weights==0))

            
            
### Saving results

            z <- list(coefficients = coef.wls,
                      x = x[,-1,drop=FALSE],
                      y = y,
                      mf = mf,
                      call = wfe.call,
                      vcov = as.numeric(vcov.wfe),
                      se = se.did,
                      sigma = try(sqrt(sigma2)),
                      df = d.f,
                      residuals = y - (x[,-1,drop=FALSE] %*% coef.wls),
                      W = W1,
                      Num.nonzero = Num.nonzero,
                      units = units,
                      times = times,
                      method = method,
                      causal = causal,
                      est = est,
                      std.error = std.error,
                      White.pvalue = white.p,
                      White.alpha = White.alpha,
                      White.stat = white.stat,
                      White.test = test.null,
                      Y.wdm = Y.wdm,
                      X.wdm = X.wdm)

            class(z) <- "wfedid"
            z
        }

    }
}




### print wfe class

print.wfe <- function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nStd.Err:\n")
    print(x$se)
}

summary.wfe <- function(object, signif.stars = getOption("show.signif.stars"),...){
    se <- object$se
    sigma <- object$sigma
    df <- object$df
    tval <- try(coef(object) / se)
    TAB <- cbind(Estimate = coef(object),
                 Std.Error = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df = object$df))
    res <- list(call = object$call,
                coefficients = TAB,
                sigma = object$sigma,
                df = object$df,
                W = object$W,
                Num.nonzero = object$Num.nonzero,
                units = object$units,
                times = object$times,
                residuals = object$residuals,
                method = object$method,
                causal = object$causal,
                estimator = object$est,
                std.error = object$std.error,
                White.pvalue = object$White.pvalue,
                White.alpha = object$White.alpha,
                White.stat = object$White.stat,
                White.test = object$White.test,
                Y = object$y,
                X = object$x,
                Y.wdm = object$Y.wdm,
                X.wdm = object$X.wdm              
                )
    class(res) <- "summary.wfe"
    res
}

print.summary.wfe <- function(x, ...){
    cat("\nMethod:", x$method, "Fixed Effects\n")
    cat("\nQuantity of Interest:", x$causal)
    cat("\nEstimator:", x$estimator)
    cat("\nStandard Error:", x$std.error)
    cat("\n")
    cat("\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Coefficients:\n")
    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
    cat("\nResidual standard error:", format(signif(x$sigma,
                                                    4)), "on", x$df, "degrees of freedom")
    cat("\nWhite statistics for functional misspecification:", x$White.stat, "with Pvalue=", x$White.pvalue)
    cat("\nReject the null of NO misspecification:", x$White.test)
    cat("\n")
}



### print wfedid class

print.wfedid <- function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nStd.Err:\n")
    print(x$se)
}

summary.wfedid <- function(object, signif.stars = getOption("show.signif.stars"),...){
    coef <- object$coefficients
    se <- object$se
    sigma <- object$sigma
    df <- object$df
    tval <- coef(object) / se
    TAB <- cbind(Estimate = coef,
                 Std.Error = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df = object$df))
    res <- list(call = object$call,
                coefficients = TAB,
                sigma = object$sigma,
                df = object$df,
                W = object$W,
                Num.nonzero = object$Num.nonzero,
                residuals = object$residuals,
                method = object$method,
                causal = object$causal,
                estimator = object$est,
                std.error = object$std.error,
                White.pvalue = object$White.pvalue,
                White.alpha = object$White.alpha,
                White.stat = object$White.stat,
                White.test = object$White.test,
                Y = object$y,
                X = object$x,
                Y.wdm = object$Y.wdm,
                X.wdm = object$X.wdm
                )
    class(res) <- "summary.wfedid"
    res
}

print.summary.wfedid <- function(x, ...){
    cat("\nMethod:", x$method, "Fixed Effects\n")
    cat("\nQuantity of Interest:", x$causal)
    cat("\nEstimator:", x$estimator)
    cat("\nStandard Error:", x$std.error)
    cat("\n")
    cat("\n")
    cat("Call:\n")
    cat("Coefficients:\n")
    print(x$call)
    cat("\n")
    colnames(x$coefficients) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
    cat("\nResidual standard error:", format(signif(x$sigma,
                                                    4)), "on", x$df, "degrees of freedom")
    if (!is.null(x$White.stat)){
        cat("\nWhite statistics for functional misspecification:", x$White.stat, "with Pvalue=", x$White.pvalue)
        cat("\nReject the null of NO misspecification:", x$White.test)
    }
    cat("\n")
}



