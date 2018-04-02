
Index <- function(index_name, uniq_index_name, len_u_index, len_data) {

  return(.C("Index", as.integer(index_name), as.integer(uniq_index_name), as.integer(len_u_index),
            as.integer(len_data),
            result = integer(len_data))$result)

}



VectorizeC <- function(W, time.index, dyad.index, n.row) {

  return(.C("VectorizeC", as.double(W), as.integer(nrow(W)), as.integer(ncol(W)),
            as.integer(time.index), as.integer(dyad.index), as.integer(n.row),
            results = double(n.row))$results)

}



Transform <- function(y, treat, pscore) {

  return(.C("Transform", as.double(y), as.integer(length(y)), as.integer(treat),
            as.double(pscore), ytrans = double(length(y)))$ytrans)

}



GenTime <- function(unit_index, len_data, len_u_index) {
  return(.C("GenTime", as.integer(unit_index), as.integer(len_data),
            as.integer(len_u_index), time_index = double(len_data))$time_index)
}



GenWeightsUnit <- function(unit_index, time_index, tr, C_it, len_data, len_u_index, len_t_index, ate, att, size, verbose) {
  return(.C("GenWeightsUnit", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att), as.integer(verbose),
            weight = double(len_u_index*len_t_index))$weight)
}


GenWeightsTime <- function(time_index, unit_index, tr, C_it, len_data, len_t_index, len_u_index, ate, att, size, verbose) {
  return(.C("GenWeightsTime", as.integer(time_index), as.integer(unit_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_t_index), as.integer(len_u_index),
            as.integer(ate), as.integer(att), as.integer(verbose),
            weight = double(len_u_index*len_t_index))$weight)
}





GenWeightsFD <- function(unit_index, time_index, tr, C_it, len_data, len_u_index, len_t_index, ate, att, verbose) {
  return(.C("GenWeightsFD", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att), as.integer(verbose),
            weightfd = double(len_u_index*len_t_index))$weightfd)
}


CalDID <- function(unit_index, time_index, tr, C_it, y, len_data, len_u_index, len_t_index, ate, att, verbose) {
  return(.C("CalDID", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.double(y), as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att), as.integer(verbose),
            did = double(1))$did)
}



DemeanDID <- function(var_name, weight, unit_index, time_index, len_u_index, len_t_index, len_data) {
  return(.C("DemeanDID", as.double(var_name), as.double(weight), 
            as.integer(unit_index), as.integer(time_index),
            as.integer(len_u_index), as.integer(len_t_index), as.integer(len_data),
            DemeanDID = double(len_data))$DemeanDID)
}


GenWeightsMDID <- function(unit_index, time_index, tr, C_it, y, maxdev.did, len_data,
                           len_u_index, len_t_index, ate, att, verbose) {
    return(.C("GenWeightsMDID", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
              as.double(y), as.double(maxdev.did), as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
              as.integer(ate), as.integer(att), as.integer(verbose),
              weightmdid = double(len_u_index*len_t_index))$weightmdid)
}



GenWeightsDID <- function(unit_index, time_index, tr, C_it, len_data, len_u_index, len_t_index, ate, att, verbose) {
  return(.C("GenWeightsDID", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att), as.integer(verbose),
            weightdid = double(len_u_index*len_t_index))$weightdid)
}

Demean <- function(var_name, index, len_index, len_data) {
  return(.C("Demean", as.double(var_name), as.integer(index), as.integer(len_index),
            as.integer(len_data), demean = double(len_data))$demean)
}




WDemean <- function(var_name, weight, index, len_index, len_data) {
  return(.C("WDemean", as.double(var_name), as.double(weight), as.integer(index), as.integer(len_index),
            as.integer(len_data), Wdemean = double(len_data))$Wdemean)
}




WWDemean <- function(var_name, weight, index, len_index, len_data) {
  return(.C("WWDemean", as.double(var_name), as.double(weight), as.integer(index), as.integer(len_index),
            as.integer(len_data), WWdemean = double(len_data))$WWdemean)
}



TwayDemean <- function(var_name, unit_index, time_index, len_u_index, len_t_index, len_data) {
  return(.C("TwayDemean", as.double(var_name), as.integer(unit_index), as.integer(time_index),
            as.integer(len_u_index), as.integer(len_t_index), as.integer(len_data),
            TwayDemean = double(len_data))$TwayDemean)
}



MDummy <-  function(index, len_index, len_data) {
  return(.C("MDummy", as.integer(index), as.integer(len_index), as.integer(len_data),
            dummy = integer(len_data*len_index))$dummy)
}




XXiSum <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde) {
  return(.C("XXiSum", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde),
            result = double(n_cov*n_cov))$result)
}



XWXiSum <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, weights) {
  return(.C("XWXiSum", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(weights),
            result = double(n_cov*n_cov))$result)
}


OmegaHatHAC <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, u.tilde) {
  return(.C("OmegaHatHAC", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(u.tilde),
            Omega_hat_HAC = double(n_cov*n_cov))$Omega_hat_HAC)
}


OmegaHatHC <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, u.tilde) {
  return(.C("OmegaHatHC", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(u.tilde),
            Omega_hat_HC = double(n_cov*n_cov))$Omega_hat_HC)
}




OmegaDiDHAC <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, u.tilde, W) {
  return(.C("OmegaDiDHAC", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(u.tilde), as.double(W),
            Omega_DiD_HAC = double(n_cov*n_cov))$Omega_DiD_HAC)
}


OmegaDiDHAC2 <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, u.tilde, W) {
  return(.C("OmegaDiDHAC2", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(u.tilde), as.double(W),
            Omega_DiD_HAC = double(n_cov*n_cov))$Omega_DiD_HAC)
}



OmegaDiDHAC <-  function(len_data, n_cov, unit_index, len_uniq_unit_index, X.tilde, u.tilde, W) {
  return(.C("OmegaDiDHAC", as.integer(len_data), as.integer(n_cov),
            as.integer(unit_index), as.integer(len_uniq_unit_index),
            as.double(X.tilde), as.double(u.tilde), as.double(W),
            Omega_DiD_HAC = double(n_cov*n_cov))$Omega_DiD_HAC)
}




ProjectionM <-  function(Q_QQinv, Q, P1.first, P1.second,
                         P1.third, P1.fourth, P1.fifth,
                         Y, X, len_data, n_col, n_var,
                         n_p1, n_p2, n_p3, n_p4){
  return(.C("ProjectionM", as.complex(Q_QQinv), as.complex(Q),
            as.complex(P1.first), as.complex(P1.second),
            as.complex(P1.third), as.complex(P1.fourth), as.complex(P1.fifth),
            as.complex(Y), as.complex(X),
            as.integer(len_data), as.integer(n_col),
            as.integer(n_var), as.integer(n_p1),
            as.integer(n_p2), as.integer(n_p3), as.integer(n_p4),
            Projection = complex(len_data*n_var))$Projection)
}



comp_OmegaHAC <-  function(X_1, u_1, X_2, u_2, len_data,
                           n_cov, unit_index, len_unit){
  return(.C("comp_OmegaHAC", as.complex(X_1), as.complex(u_1),
            as.complex(X_2), as.complex(u_2), as.integer(len_data),
            as.integer(n_cov), as.integer(unit_index), as.integer(len_unit),
            OmegaHAC = complex(n_cov*n_cov))$OmegaHAC)
}



comp_OmegaHC <-  function(X_1, u_1, X_2, u_2, len_data,
                          n_cov, unit_index, len_unit){
  return(.C("comp_OmegaHC", as.complex(X_1), as.complex(u_1),
            as.complex(X_2), as.complex(u_2), as.integer(len_data),
            as.integer(n_cov), as.integer(unit_index), as.integer(len_unit),
            OmegaHC = complex(n_cov*n_cov))$OmegaHC)
}



LamdaDID1 <-  function(len_Xtrow, len_Xhrow, Tunit_index, len_uniq_Tunit_index,
                       Hunit_index, len_uniq_Hunit_index,
                       X.tilde, len_Xtcol, u.tilde,
                       X.hat, len_Xhcol, u.hat, W) {
  return(.C("LamdaDID1", as.integer(len_Xtrow), as.integer(len_Xhrow),
            as.integer(Tunit_index), as.integer(len_uniq_Tunit_index),
            as.integer(Hunit_index), as.integer(len_uniq_Hunit_index),
            as.double(X.tilde), as.integer(len_Xtcol), as.double(u.tilde),
            as.double(X.hat), as.integer(len_Xhcol), as.double(u.hat), as.double(W),
            LamdaDID1 = double(len_Xhcol*len_Xtcol))$LamdaDID1)
}



LamdaDID2 <-  function(len_Xtrow, len_Xhrow, Tunit_index, len_uniq_Tunit_index,
                       Hunit_index, len_uniq_Hunit_index,
                       X.tilde, len_Xtcol, u.tilde,
                       X.hat, len_Xhcol, u.hat, W) {
  return(.C("LamdaDID2", as.integer(len_Xtrow), as.integer(len_Xhrow),
            as.integer(Tunit_index), as.integer(len_uniq_Tunit_index),
            as.integer(Hunit_index), as.integer(len_uniq_Hunit_index),
            as.double(X.tilde), as.integer(len_Xtcol), as.double(u.tilde),
            as.double(X.hat), as.integer(len_Xhcol), as.double(u.hat), as.double(W),
            LamdaDID2 = double(len_Xtcol*len_Xhcol))$LamdaDID2)
}





Wdemean <- function (x, w, unit.index, time.index, data) {
  data$u.index <- Index(data[,unit.index])
  data$t.index <- Index(data[,time.index])
  uniq.unit <- unique(data$u.index)
  uniq.time <- unique(data$t.index)
  wdemean.x <- c()
  wdemean.x <- new.tr.x <- matrix(NA, nrow = length(uniq.time), ncol = length(uniq.unit))
  for (i in 1:length(uniq.unit)){
    sub.i <- data[data$u.index == uniq.unit[i], ]
    nr <- nrow(sub.i)
    denom <- as.numeric(sum(sub.i[,w]))
    x.star <- as.numeric((sub.i[,x]%*%sub.i[,w])/(denom))
    tr.x <- as.vector(sub.i[,x] - rep(x.star, nr))
    wdemean.x[,i] <- tr.x
      # sqrt(w)* (weighted demean) for SE
      for (j in 1:length(uniq.time)){
       new.tr.x[j,i] <- sqrt(sub.i[,w][j])*(sub.i[,x][j] - x.star)
      }
  }
  w.demean <- as.vector(wdemean.x)
  w.tr.x <- as.vector(new.tr.x)
  list(w.demeaned = wdemean.x, new.tr.x = w.tr.x)
}




## a function that checks memory usage

.ls.objects <- function (pos = 1, pattern, order.by = "Size", decreasing=TRUE, head = TRUE, n = 10) {
    # based on postings by Petr Pikal and David Hinds to the r-help list in 2004
    # modified by: Dirk Eddelbuettel (http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session)
    # I then gave it a few tweaks (show size as megabytes and use defaults that I like)
    # a data frame of the objects and their associated storage needs.
    napply <- function(names, fn) sapply(names, function(x)
                                                   fn(get(x, pos = pos)))
      names <- ls(pos = pos, pattern = pattern)
      obj.class <- napply(names, function(x) as.character(class(x))[1])
      obj.mode <- napply(names, mode)
      obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
      obj.size <- napply(names, object.size) / 10^6 # megabytes
      obj.dim <- t(napply(names, function(x)
                                      as.numeric(dim(x))[1:2]))
      vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
      obj.dim[vec, 1] <- napply(names, length)[vec]
      out <- data.frame(obj.type, obj.size, obj.dim)
      names(out) <- c("Type", "Size", "Rows", "Columns")
      out <- out[order(out[[order.by]], decreasing=decreasing), ]
      if (head)
            out <- head(out, n)
      out
  }


## ### wfe objects for xtable

## xtable.wfe <- function(x,caption=NULL,label=NULL,align=NULL,
##                       digits=NULL,display=NULL,...) {
##   return(xtable.summary.wfe(summary(x),caption=caption,label=label,
##                            align=align, digits=digits,display=display))
## }

## xtable.summary.wfe <- function(x,caption=NULL,label=NULL,align=NULL,
##                               digits=NULL,display=NULL,...) {
##   x <- data.frame(x$coef,check.names=FALSE)

##   class(x) <- c("xtable","data.frame")
##   caption(x) <- caption
##   label(x) <- label
##   align(x) <- switch(1+is.null(align),align,c("r","r","r","r","r"))
##   digits(x) <- switch(1+is.null(digits),digits,c(0,4,4,2,4))
##   display(x) <- switch(1+is.null(display),display,c("s","f","f","f","f"))
##   return(x)
## }

