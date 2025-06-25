#
# Sourcing this R file contains h.fct function used in mph.fit().
#
# Create an environment to store all f-functions
  constraint <- new.env(parent = emptyenv())
#
#
########################  BEGIN define constraint  #############################

###----------------------------------------------------###
###        Constraint function for the ME model        ###
###----------------------------------------------------###

constraint$ME <- function(p, f.fam = NULL, score = NULL, n.dim = NULL, n.cat = NULL) {
  n.alpha <- n.dim - 1      # df for alpha parameters

  # Build the design matrix Delta_1
  X.mat <- build.design.mat_alt(model = "ME", n.dim = n.dim, n.cat = n.cat, score = score)
  Delta_1 <- X.mat[, 1:n.alpha]
  constraint.vec <- t(Delta_1) %*% p
  return(as.matrix(constraint.vec))
}




###----------------------------------------------------###
###        Constraint function for the VE model        ###
###----------------------------------------------------###

constraint$VE <- function(p, f.fam = NULL, score = NULL, n.dim = NULL, n.cat = NULL) {
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  
  ## compute marginal variances
  var.vec <- numeric(n.dim)
  for (h in seq_len(n.dim)) {
    var.vec[h] <- t(score[idx.mat[ , h]]^2) %*% p - (t(score[idx.mat[ , h]]) %*% p)^2
  }
  constraint.vec <- diff(var.vec)
  return(as.matrix(constraint.vec))
}




###----------------------------------------------------###
###        Constraint function for the CE model        ###
###----------------------------------------------------###

constraint$CE <- function(p, f.fam = NULL, score = NULL, n.dim = NULL, n.cat = NULL) {
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  
  ## compute marginal means and variances
  mean.vec <- var.vec <- numeric(n.dim)
  for (h in seq_len(n.dim)) {
    mean.vec[h] <- t(score[idx.mat[ , h]]) %*% p
    var.vec[h] <- t(score[idx.mat[ , h]]^2) %*% p - (t(score[idx.mat[ , h]]) %*% p)^2
  }
  if (any(var.vec < .Machine$double.eps))
    stop("Zero (or near-zero) marginal variance — correlations undefined.")
  
  ## compute correlations
  pair.mat <- t(combn(n.dim, 2))
  n.pairs <- nrow(pair.mat)
  corr.vec <- numeric(n.pairs)
  for (k in seq_len(n.pairs)) {
    s <- pair.mat[k, 1];    t <- pair.mat[k, 2]
    cov <- t(score[idx.mat[ , s]] * score[idx.mat[ , t]]) %*% p - (mean.vec[s] * mean.vec[t])
    corr.vec[k] <- cov / sqrt(var.vec[s] * var.vec[t])
  }
  constraint.vec <- diff(corr.vec)
  return(as.matrix(constraint.vec))
}




###------------------------------------------------------###
###        Constraint function for the ME_2 model        ###
###------------------------------------------------------###

constraint$ME_2 <- function(p, f.fam = NULL, score = NULL, n.dim = NULL, n.cat = NULL) {
  # Build the design matrix Delta_V2
  X.mat <- build.design.mat_alt(model = "ME_2", n.dim = n.dim, n.cat = n.cat, score = score)
  constraint.vec <- t(X.mat) %*% p
  return(as.matrix(constraint.vec))
}

##########################  END define constraint  #############################