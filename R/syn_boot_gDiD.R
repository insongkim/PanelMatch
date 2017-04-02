syn_boot_gDiD <- function(d, unit.id = "ccode", time.id = "year", ITER = 10, qoi = "ATE", L = 4, F = 0, 
                          covariate, weights_att = "weights_att", weights_atc = "weights_atc",
                          dit_att = "dit_att", dit_atc = "dit_atc",
                          treatment, dependent) {
  
  coefs <- rep(NA, ITER) 
  if (qoi == "ATT") {
    o.coef <- syn_DID(L = L, F, unit.id = unit.id, treatment = treatment, time.id = time.id, dependent = dependent, 
                      covariate = covariate, d = d, qoi = "ATT")
    d <- syn_DID_weights(L = L, F, time.id = time.id, qoi = "ATT",
                         unit.id = unit.id,
                         treatment = treatment, covariate = covariate, dependent, d)
    d <- na.omit(d[c(unit.id, time.id, treatment, dependent,
                     weights_att, dit_att)])
    for (k in 1:ITER) {  
      # make new data
      clusters <- unique(d[, unit.id])
      units <- sample(clusters, size = length(clusters), replace=T)
      # create bootstap sample with sapply
      df.bs <- lapply(units, function(x) which(d[,unit.id]==x))
      d.sub1 <- d[unlist(df.bs),]
      colnames(d.sub1)[3:4] <- c("treatment", "dv")
      att.new <- (sum(1 * d.sub1[d.sub1$treatment== 1, ]$dv * d.sub1[d.sub1$treatment== 1, ]$weights_att) -
                    sum((d.sub1[d.sub1$treatment== 0, ]$treatment+1)*d.sub1[d.sub1$treatment== 0, ]$dv*d.sub1[d.sub1$treatment== 0, ]$weights_att))/sum(d.sub1$dit_att)
      coefs[k] <- att.new
    }
    return(list("o.coef" = o.coef, "boots" = coefs))
  } else {
    o.coef <- syn_DID(L = L, F, unit.id = unit.id, treatment = treatment, time.id = time.id, dependent = dependent, 
                      covariate = covariate, d = d, qoi = "ATE")
    d <- syn_DID_weights(L = L, F, time.id = time.id, qoi = "ATE",
                         unit.id = unit.id,
                         treatment = treatment, covariate = covariate, dependent, d)
    d <- na.omit(d[c(unit.id, time.id, treatment, dependent, weights_atc,
                     weights_att, dit_att, dit_atc)])
    for (k in 1:ITER) {
      # make new data
      clusters <- unique(d[, unit.id])
      units <- sample(clusters, size = length(clusters), replace=T)
      # create bootstap sample with sapply
      df.bs <- lapply(units, function(x) which(d[,unit.id]==x))
      d.sub1 <- d[unlist(df.bs),]
      colnames(d.sub1)[3:4] <- c("treatment", "dv")
      att.new <- (sum(1 * d.sub1[d.sub1$treatment== 1, ]$dv * d.sub1[d.sub1$treatment== 1, ]$weights_att) -
                    sum((d.sub1[d.sub1$treatment== 0, ]$treatment+1)*d.sub1[d.sub1$treatment== 0, ]$dv*d.sub1[d.sub1$treatment== 0, ]$weights_att))/sum(d.sub1$dit_att)
      atc.new <- -((sum(1 * d.sub1[d.sub1$treatment== 0, ]$dv * d.sub1[d.sub1$treatment== 0, ]$weights_atc) -
                      sum(1 * d.sub1[d.sub1$treatment== 1, ]$dv * d.sub1[d.sub1$treatment== 1, ]$weights_atc))/sum(d.sub1$dit_atc))
      coefs[k] <- (att.new*sum(d.sub1$dit_att) + atc.new*sum(d.sub1$dit_atc))/(sum(d.sub1$dit_att) + sum(d.sub1$dit_atc))
      
    }
    return(list("o.coef" = o.coef, "boots" = coefs))
  }
  
}
