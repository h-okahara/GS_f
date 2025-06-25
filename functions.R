#
# Sourcing this R file (> source("functions.R") ) results in the creation 
# of the following 13 functions:
#
# Calc.F.vec, build.design.mat, build.design.mat_alt
# Calc.M.moments, Calc.variances, Calc.covariances, Calc.correlations,
# get.constraints, get.f_function, Calc.theta.vec, build.covariance,
# generate.3dimNormal, simulate.3dimNormal
#
###############  BEGIN Subroutines for the GS[f] and S model  ##################

###--------------------------------------###
###    Calculate F.vec : F(pi / pi^S)    ###
###--------------------------------------###

## INPUT:
# p      : Numeric vector, joint cell probabilities.
# F.func : Vectorised function that returns F(x) for numeric vector x (F = f').
# n.dim  : Integer, dimensionality of the contingency table.
# n.cat  : Integer, number of categories of the contingency table.

## OUTPUT:
# An n.cat^n.dim array whose element in i = (i_1,...,i_T) equals pi_i / pi_i^S,
# where pi^S is the complete symmetry structure defined in Sec. 2.1.

Calc.F.vec <- function(p, F.func = NULL, n.dim = NULL, n.cat = NULL) {
  ## Input validation
  if (is.null(n.dim)) stop("Argument 'n.dim' is missing – specify the table dimensionality (T).")
  if (is.null(n.cat)) stop("Argument 'n.cat' is missing – specify the number of categories (r).")
  if (length(p) != n.cat^n.dim) stop("Length of 'p' must be n.cat^n.dim.")
  
  ## Generate all cell configurations and compute pi^S
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))   # n.cat^n.dim × n.dim matrix
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  sym.positions <- apply(idx.mat, 1, function(v) paste(sort(v), collapse = "-"))
  p.sym  <- ave(p, sym.positions, FUN = mean)  # pi^S
  
  ## Compute F(pi / pi^S)
  p.ratio <- p / p.sym
  F.vec <- F.func(p.ratio)
  return(as.matrix(F.vec))
}




###-----------------------------------------------###
###    Build Design Matrix for the GS[f] model    ###
###     - build.design.mat     for the Eq.(2)     ###
###     - build.design.mat_alt for the Eq.(11)    ###
###-----------------------------------------------###

# INPUT:
# model : Character string that chooses which sub-model to build.
# n.dim : Integer, dimensionality of the contingency table.
# n.cat : Integer, number of categories of the contingency table.
# score : A vector of scores attached to each category of a table.

# OUTPUT:
# Design Matrix X (r^T x n.df), where df is number of df the model.

build.design.mat <- function(model = NULL, n.dim = NULL, n.cat = NULL, score = NULL) 
  {
  ## Input validation
  if (is.null(model)) stop("Argument 'model' is missing - choose from 'S', 'ME', 'ME_2', 'GS_f', 'ELS_f', 'LS_f'.")
  if (is.null(n.dim)) stop("Argument 'n.dim' is missing – specify the table dimensionality (T).")
  if (is.null(n.cat)) stop("Argument 'n.cat' is missing – specify the number of categories (r).")
  if (is.null(score)) score <- seq_len(n.cat)   # default 1:r
  n.cells <- n.cat^n.dim
  
  ## Lexicographic cell index matrix
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  if (nrow(idx.mat) != n.cells) stop("nrow of 'idx.mat' must be n.cat^n.dim.")
  
  # Construct X_S for gamma' parameters
  sym.positions <- apply(idx.mat, 1, function(v) paste(sort(v), collapse = "-"))
  unique.sym.positions <- unique(sym.positions)
  L_rT <- length(unique.sym.positions)
  X_S <- matrix(0, nrow = n.cells, ncol = L_rT)
  for (k in seq_len(L_rT)) {
    X_S[, k] <- as.integer(sym.positions == unique.sym.positions[k])
  }
  if (model == "S") return(X_S)
  
  ## Construct columns (1,...,d_2) for alpha' and B parameters
  ## Identifiability constraints: alpha_T = 0, beta_{TT} = 0, beta_{(T-1)T} = 0.
  n.alpha <- n.dim - 1      # df for alpha parameters
  n.beta_diag <- n.dim - 1  # df for diagonal of matrix B
  pair.mat <- t(combn(n.dim, 2))
  n.beta_off <- nrow(pair.mat) - 1  # df for off-diagonal of matrix B
  alpha.mat <- matrix(0, nrow = n.cells, ncol = n.alpha)
  beta_diag.mat <- matrix(0, nrow = n.cells, ncol = n.beta_diag)
  beta_off.mat <- matrix(0, nrow = n.cells, ncol = n.beta_off)
  
  # Columns for alpha_s (s = 1,...,T-1)
  for (h in seq_len(n.alpha)) {
    alpha.mat[ , h] <- score[idx.mat[ , h]]
  }
  if (model == "LS_f") return(cbind(alpha.mat, X_S))
  if (model == "ME") return(alpha.mat)
  
  # Columns for beta_{ss} (s = 1,...,T-1)
  for (h in seq_len(n.beta_diag)) {
    beta_diag.mat[ , h] <- score[idx.mat[ , h]]^2
  }
  if (model == "ELS_f") return(cbind(alpha.mat, beta_diag.mat, X_S))
  
  # Columns for beta_{st} (s < t, (s,t) != (T-1, T))
  for (k in seq_len(n.beta_off)) {
    s <- pair.mat[k, 1];    t <- pair.mat[k, 2]
    beta_off.mat[ , k] <- score[idx.mat[ , s]] * score[idx.mat[ , t]]
  }
  if (model == "GS_f") return(cbind(alpha.mat, beta_diag.mat, beta_off.mat, X_S))
  if (model == "ME_2") return(cbind(alpha.mat, beta_diag.mat, beta_off.mat))
}


build.design.mat_alt <- function(model = NULL, n.dim = NULL, n.cat = NULL, score = NULL) 
  {
  ## Input validation
  if (is.null(model)) stop("Argument 'model' is missing - choose from 'S', 'ME', 'ME_2', 'GS_f', 'ELS_f', 'LS_f'.")
  if (is.null(n.dim)) stop("Argument 'n.dim' is missing – specify the table dimensionality (T).")
  if (is.null(n.cat)) stop("Argument 'n.cat' is missing – specify the number of categories (r).")
  if (is.null(score)) score <- seq_len(n.cat)   # default 1:r
  n.cells <- n.cat^n.dim
  
  ## Lexicographic cell index matrix
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  if (nrow(idx.mat) != n.cells) stop("nrow of 'idx.mat' must be n.cat^n.dim.")
  
  # Construct X_S for gamma' parameters
  sym.positions <- apply(idx.mat, 1, function(v) paste(sort(v), collapse = "-"))
  unique.sym.positions <- unique(sym.positions)
  L_rT <- length(unique.sym.positions)
  X_S <- matrix(0, nrow = n.cells, ncol = L_rT)
  for (k in seq_len(L_rT)) {
    X_S[, k] <- as.integer(sym.positions == unique.sym.positions[k])
  }
  if (model == "S") return(X_S)
    
  ## Construct columns x_1,...,x_{d_2} for alpha' and B parameters
  ## Identifiability constraints: alpha_T = 0, beta_{TT} = 0, beta_{(T-1)T} = 0.
  n.alpha <- n.dim - 1      # df for alpha parameters
  n.beta_diag <- n.dim - 1  # df for diagonal of matrix B
  pair.mat <- t(combn(n.dim, 2))
  n.beta_off <- nrow(pair.mat) - 1  # df for off-diagonal of matrix B
  alpha.mat <- matrix(0, nrow = n.cells, ncol = n.alpha)
  beta_diag.mat <- matrix(0, nrow = n.cells, ncol = n.beta_diag)
  beta_off.mat <- matrix(0, nrow = n.cells, ncol = n.beta_off)
  
  # Columns for alpha_s (s = 1,...,T-1)
  for (h in seq_len(n.alpha)) {
    alpha.mat[ , h] <- score[idx.mat[ , h]] - score[idx.mat[ , h + 1]]
  }
  if (model == "LS_f") return(cbind(alpha.mat, X_S))
  if (model == "ME") return(alpha.mat)
  
  # Columns for beta_{ss} (s = 1,...,T-1)
  for (h in seq_len(n.beta_diag)) {
    beta_diag.mat[ , h] <- score[idx.mat[ , h]]^2 - score[idx.mat[ , h + 1]]^2
  }
  if (model == "ELS_f") return(cbind(alpha.mat, beta_diag.mat, X_S))
  
  # Columns for beta_{st} (s < t, (s,t) != (T-1, T))
  for (k in seq_len(n.beta_off)) {
    s <- pair.mat[k, 1];    t <- pair.mat[k, 2]
    u <- pair.mat[k+1, 1];  v <- pair.mat[k+1, 2]
    beta_off.mat[ , k] <- score[idx.mat[ , s]] * score[idx.mat[ , t]] - score[idx.mat[ , u]] * score[idx.mat[ , v]]
  }
  if (model == "GS_f") return(cbind(alpha.mat, beta_diag.mat, beta_off.mat, X_S))
  if (model == "ME_2") return(cbind(alpha.mat, beta_diag.mat, beta_off.mat))
}


################  END Subroutines for the GS[f] and S model  ###################


##############  BEGIN Subroutines for the ME,VE,CE and ME_2 models  ############

###----------------------###
###    Calc.M.moments    ###
###----------------------###

## INPUT:
# df      : A data.frame whose rows correspond to the cells of the table.
# p       : A vector of length r^T holding each cell probability. 
# score   : A vector of scores attached to each category of a table.
# n.dim   : Integer, dimensionality of the contingency table.
# M.order : Positive integer m indicating that the routine returns the m‑th
#           marginal moments E[X_t^m] for t = 1,…,n.dim.

## OUTPUT:
# A vactor that contains each marginal mean E[X_i^m] for i=1,...,n.dim.

Calc.M.moments <- function(df, p, n.dim = NULL, score = NULL, M.order = 1) {
  ## Input validation
  if (is.null(n.dim)) {
    n.dim <- length(grep("^Var", names(df))) 
  }
  if (is.null(score)) {
    n.cat <- nlevels(df[["Var1"]])
    score <- seq_len(n.cat)         # 1,2,…,n_cat for every dimension
  }
  
  ## Build marginalizing matrix of frequencies and probabilities for this dimension
  M.moment.vec <- numeric(n.dim)
  for (t in seq_len(n.dim)) {
    col.name <- paste0("Var", t)
    M.freq <- M.fct(df[[col.name]])
    M.prob <- as.vector(M.freq %*% p)
    M.moment.vec[t] <- sum(M.prob * score^M.order)
  }
  return(M.moment.vec)
}




###----------------------###
###    Calc.variances    ###
###----------------------###

## INPUT:
# df    : A data.frame whose rows correspond to the cells of the table.
# p     : A vector of length r^T holding each cell probability. 
# n.dim : Integer, dimensionality of the contingency table.
# score : A vector of scores attached to each category of a table.

## OUTPUT:
# A vactor that contains each marginal variances Var[X_i] for i=1,...,n.dim.

Calc.variances <- function(df, p, n.dim = NULL, score = NULL) {
  ## Input validation
  if (is.null(n.dim)) {
    n.dim <- length(grep("^Var", names(df))) 
  }
  if (is.null(score)) {
    n.cat <- nlevels(df[["Var1"]])
    score <- seq_len(n.cat)         # 1,2,…,n_cat for every dimension
  }
  
  M.variance.vec <- numeric(n.dim)
  M.moments1 <- Calc.M.moments(df, p, n.dim, score, M.order = 1)
  M.moments2 <- Calc.M.moments(df, p, n.dim, score, M.order = 2)
  M.variance.vec <- M.moments2 - M.moments1^2
  return(M.variance.vec)
}




###------------------------###
###    Calc.covariances    ###
###------------------------###

## INPUT:
# df    : A data.frame whose rows correspond to the cells of the table.
# p     : A vector of length r^T holding each cell probability. 
# n.dim : Integer, dimensionality of the contingency table.
# score : A vector of scores attached to each category of a table.

## OUTPUT:
# A vactor that contains each covariance Cov[X_i, X_j] for i<j.

Calc.covariances <- function(df, p, n.dim = NULL, score = NULL) {
  ## Input validation
  if (is.null(n.dim)) {
    n.dim <- length(grep("^Var", names(df)))
  }
  if (is.null(score)) {
    n_cat <- nlevels(df[["Var1"]])
    score <- seq_len(n_cat)         # 1,2,…,n_cat for every dimension
  }
  
  ## Build marginalizing matrix of frequencies and probabilities for this pair.
  M.mean.vec <- Calc.M.moments(df, p, n.dim = n.dim, score = score, M.order = 1)
  covariance.vec <- numeric(choose(n.dim, 2))
  idx <- 1
  for (i in 1:(n.dim - 1)) {
    for (j in (i + 1):n.dim) {
      col.name_i <- paste0("Var", i)
      col.name_j <- paste0("Var", j)
      pair_ij <- interaction(df[[col.name_i]], df[[col.name_j]], drop = TRUE)
      M.freq <- M.fct(pair_ij)
      M.prob <- as.vector(M.freq %*% p)
      
      # Compute covariance vector
      M.prod.score <- as.vector(outer(score, score))
      covariance.vec[idx] <- sum(M.prob * M.prod.score) - M.mean.vec[i] * M.mean.vec[j]
      idx <- idx + 1
    }
  }
  return(covariance.vec)
}




###-------------------------###
###    Calc.correlations    ###
###-------------------------###

## INPUT:
# df    : A data.frame whose rows correspond to the cells of the table.
# p     : A vector of length r^T holding each cell probability. 
# n.dim : Integer, dimensionality of the contingency table.
# score : A vector of scores attached to each category of a table.

## OUTPUT:
# A vactor that contains each correlations Cor[X_i X_j] for i < j.

Calc.correlations <- function(df, p, n.dim = NULL, score = NULL) {
  ## Input validation
  if (is.null(n.dim)) n.dim <- length(grep("^Var", names(df)))
  if (is.null(score)) {
    n.cat <- nlevels(df[["Var1"]])
    score <- seq_len(n.cat)   # 1,2,…,n_cat for every dimension
  }
  
  M.variance.vec <- Calc.variances(df, p, n.dim = n.dim, score = score)
  covariance.vec <- Calc.covariances(df, p, n.dim = n.dim, score = score)
  correlation.vec <- numeric(length(covariance.vec))
  idx     <- 1
  for (i in 1:(n.dim - 1)) {
    for (j in (i + 1):n.dim) {
      normalization <- sqrt(M.variance.vec[i] * M.variance.vec[j])
      correlation.vec[idx] <- if (normalization > 0) covariance.vec[idx] / normalization else NA_real_
      idx <- idx + 1
    }
  }
  return(correlation.vec)
}

###############  END Subroutines for the ME,VE,CE and ME_2 models  #############


###############  BEGIN Subroutines or Complementary functions  #################

###--------------------------------------###
###    Get constraints for each model    ###
###--------------------------------------###

## INPUT:
# model : Character string — one of "GS_f", "ME_2", "ME", "VE", or "CE".
# f.fam : An f-function stored in the 'f_functions.R'.
#          Only required for "GS_f"; if NULL KL-divergence is used by default.
# score : A vector of scores attached to each category of a table.
# n.dim : Integer, dimensionality of the contingency table.
# n.cat : Integer, number of categories of the contingency table.

## OUTPUT:
# A function with signature 'p' that calls the corresponding constraint in `constraints.R`.

get.constraints <- function(model = NULL, f.fam = NULL, 
                            score = NULL, n.dim = NULL, n.cat = NULL) 
  {
  ## Input validation
  if (is.null(model)) stop("Argument 'model' is missing - choose from 'ME_2', 'ME', 'VE', 'CE'.")
  if ((model == "GS_f" || model == "S") && is.null(f.fam)) {
    stop("Argument 'f.fam' is missing - choose from 'kl', 'revkl', 'pearson', 'hellinger', 'power', 'alpha'.") 
  }
  if (is.null(n.dim)) stop("Argument 'n.dim' is missing – specify the table dimensionality (T).")
  if (is.null(n.cat)) stop("Argument 'n.cat' is missing – specify the number of categories (r).")
  if (is.null(score)) score <- seq_len(n.cat)   # default 1:r
  if (!exists(model, envir = constraint))
    stop("Unknown model '", model, "'. Check the spell of 'model' and constraints.R environment.")
  h.fct <- get(model, envir = constraint)
  
  ## Create closures for mph.fit
  force(f.fam); force(score); force(n.dim); force(n.cat)
  
  function(p) { # 'mph.fit' sees only this signature 'p'
    h.fct(p, f.fam = f.fam, score = score, n.dim = n.dim, n.cat = n.cat)
  } 
}




###----------------------###
###    get.f_function    ###
###----------------------###

## INPUT:
# name   : Character string indicating the f-divergence family.
#          - one of "kl", "revkl", "pearson", "hellinger", "power", or "alpha".
# lambda : Numeric scalar, parameter for the power divergence family.
# alpha  : Numeric scalar, parameter for the alpha divergence family.

## OUTPUT:
# A list with two vectorised components
#   - $f.func – the f-function itself (f)
#   - $F.func – its first derivative (F = f′)

get.f_function <- function(name = c("kl", "revkl", "pearson", "hellinger", "power", "alpha"), 
                           lambda = NULL, alpha = NULL) {
  name <- match.arg(name)
  switch(name,
         power = {
           if (is.null(lambda)) stop("Specify 'lambda' for the power divergence.")
           f_divergence$power(lambda)
         },
         alpha = {
           if (is.null(alpha)) stop("Specify 'alpha' for the alpha divergence.")
           f_divergence$alpha(alpha)
         },
         {
           get(name, envir = f_divergence)
         }
  )
}




###-----------------------------------------------------###
###    Calculate kernels within potential parameters    ###
###-----------------------------------------------------###

## INPUT:
# model     : Character string that chooses which sub-model to calculate.
# theta.vec : mph.fit()$beta (coefficient vector).
# lambda    : A hyperparameter of power divergence.
# n.dim     : Integer, dimensionality of the contingency table.
# n.cat     : Integer, number of categories of the contingency table.
# score     : A vector of scores attached to each category of a table.

## OUTPUT:
# Data frame of cell position and the plug-in estimator θ_i within each 'model'.

Calc.theta.vec <- function(model = NULL, theta.vec, lambda = 0, 
                           n.dim = NULL, n.cat = NULL, score = NULL) 
  {
  ## Input validation
  if (is.null(model)) stop("Argument 'model' is missing - choose from 'GS_f', 'ELS_f', 'LS_f'.")
  if (is.null(n.dim)) stop("Argument 'n.dim' is missing – specify the table dimensionality (T).")
  if (is.null(n.cat)) stop("Argument 'n.cat' is missing – specify the number of categories (r).")
  if (is.null(score)) stop("Argument 'score' is missing – specify the same score as used in the goodness-of-fit test above.")
  
  ## Construct a vector α and a symmetric matrix B
  n.alpha <- n.dim - 1
  n.beta_diag <- n.dim - 1
  n.beta_off <- choose(n.dim, 2) - 1
  alpha.hat <- theta.vec[seq_len(n.alpha)]
  beta_diag.hat <- theta.vec[n.alpha + seq_len(n.beta_diag)]
  beta_off.hat  <- theta.vec[n.alpha + n.beta_diag + seq_len(n.beta_off)]
  
  alpha.hat <- c(alpha.hat, 0)
  B.hat <- matrix(0, n.dim, n.dim)
  diag(B.hat)[seq_len(n.dim - 1)] <- beta_diag.hat
  pair.mat <- t(combn(n.dim, 2))
  for (k in seq_len(n.beta_off)) {
    s <- pair.mat[k, 1];  t <- pair.mat[k, 2]
    B.hat[s, t] <- B.hat[t, s] <- beta_off.hat[k]/2
  }
  
  ## Lexicographic cell index matrix
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  sym.positions <- apply(idx.mat, 1, function(v) paste(sort(v), collapse = "-"))
  n.sym <- ave(rep(1, n.cells), sym.positions, FUN = length)  # |D(i)|
  U.mat <- matrix(score[idx.mat], nrow = n.cells, ncol = n.dim) # score matrix (n.cells × n.dim)
  
  ## Compute potential parameters for specific divergence
  lin.term  <- U.mat %*% alpha.hat                  # u_i^T α
  quad.term <- rowSums((U.mat %*% B.hat) * U.mat)   # u_i^T B u_i
  kernel.vec <- lin.term + quad.term
  
  theta.hat <- if (lambda == 0) {
    exp(kernel.vec)                       # KL-divergence
  } else {
    (lambda / n.sym^lambda) * kernel.vec  # power divergence (λ≠0)
  }
  print(data.frame(idx.mat, theta.hat = theta.hat, row.names = NULL, check.names = FALSE))
  return(theta.hat)
}




###-------------------------------###
###    Build covariance matrix    ###
###-------------------------------###

## INPUT:
# var  : A vector of variances.
# corr : A vector of pairwise correlations in lexicographic order.

## OUTPUT:
# A covariance matrix.

build.covariance <- function(var, corr) {
  n.dim <- length(var)
  cov.mat <- diag(1, n.dim)
  pairs <- t(combn(n.dim, 2))
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k,1];  j <- pairs[k,2]
    cov.mat[i,j] <- cov.mat[j,i] <- corr[k]
  }
  D <- diag(sqrt(var), n.dim)
  cov.mat <- D %*% cov.mat %*% D
  return(cov.mat)
}


###----------------------------------------------###
###    Generate true cell probabilities          ###
###    from 3-dimensional normal distribution    ###
###----------------------------------------------###

## INPUT:
# mean.vec   : Numeric vector of length 3 (marginal means μ).
# covariance : A 3 × 3 matrix (covariance matrix Σ).
# n.cat      : Integer, number of categories (must be 4 or 6).
# cuts.list  : Optional. A list of length 3. Each element is a numeric.
#              vector of cut-points for that dimension (length n.cat + 1),
#              including -Inf and Inf.
#              If NULL, cuts are generated based on mean.vec and covariance.
# n.cluster  : #CPU cores to use for parallel integration.

## OUTPUT:
# Numeric vector containing the true probabilities.

generate.3dimNormal <- function(mean.vec = NULL, covariance = NULL,
                                n.cat = NULL, cuts.list = NULL,
                                n.cluster = parallel::detectCores() - 2)
  {
  n.dim <- 3 # Hardcoded for this specific version
  
  ## Input validation
  if (is.null(n.cat) || !is.element(n.cat, c(4, 6))) {
    stop("'n.cat' must be 4 or 6 for this function.")
  }
  if (is.null(mean.vec) || length(mean.vec) != n.dim) {
    stop(paste0("'mean.vec' must be a numeric vector of length ", n.dim, "."))
  }
  if (is.null(covariance) || !all(dim(covariance) == c(n.dim, n.dim))) {
    stop(paste0("'covariance' must be a ", n.dim, "x", n.dim, " matrix."))
  }
  if (any(diag(covariance) < 0)) stop("Variances cannot be negative.")
  
  
  ## Build cut-points and intervals for each dimension
  if (is.null(cuts.list)) {
    cuts.list <- vector("list", n.dim)
    sigma.vec <- sqrt(diag(covariance))
    
    # Make intervals for specific categories n.cat = 4 or 6
    make.intervals <- function(mu, sigma, current.n.cat){
      if (current.n.cat == 4) {
        width_1  <- 0.6 * sigma
        cuts <- sort(unique(c(-Inf, mu - width_1, mu, mu + width_1, Inf)))
      } else {
        width_1  <- 0.6 * sigma
        width_2  <- 1.2 * sigma
        cuts <- sort(unique(c(-Inf, mu - width_2, mu - width_1, mu, mu + width_1, mu + width_2, Inf)))
      }
      return(cuts)
    }
    for (d in 1:n.dim) {
      cuts.list[[d]] <- make.intervals(mean.vec[1], sigma.vec[1], n.cat)
    }
  } else { # Validate provided 'cuts.list'
    if (!is.list(cuts.list) || length(cuts.list) != n.dim) {
      stop(paste0("'cuts.list' must be a list of length ", n.dim, "."))
    }
    for (d in 1:n.dim) {
      if (!is.numeric(cuts.list[[d]]) || length(cuts.list[[d]]) != (n.cat + 1)) {
        stop(paste0("'cuts.list[[", d, "]]' must be a numeric vector of length 'n.cat+1': ", n.cat+1, "."))
      }
    }
  }
  
  ## Lexicographic cell index matrix
  n.cells <- n.cat^n.dim
  idx.mat <- arrayInd(seq_len(n.cells), .dim = rep(n.cat, n.dim))
  idx.mat <- idx.mat[do.call(order, as.data.frame(idx.mat)), ]
  if (nrow(idx.mat) != n.cells) stop("nrow of 'idx.mat' must be n.cat^n.dim.")
  
  ## Setting for parallel processing
  if (n.cluster > 1) {
    cl <- parallel::makeCluster(n.cluster)
    doSNOW::registerDoSNOW(cl)
    
    # Export functions and variables to each node
    clusterEvalQ(cl, {
      library(mvtnorm)
      library(cubature)
    })
    clusterExport(cl, c("mean.vec", "covariance", "n.dim", "n.cat", "cuts.list", "idx.mat"),
                  envir = environment())
    
    # Display a progress bar
    pb <- progress::progress_bar$new(
      format = "  Integrating cells [:bar] :percent eta: :eta (:current/:total)",
      total = n.cells, clear = FALSE, width = 80)
    progress <- function(n_prog) pb$tick()
    opts <- list(progress = progress)
    
    #######################  BEGIN Parallel Processing  ########################
    cell.probabilities <- foreach::foreach(i = 1:n.cells, .combine = "c", .options.snow = opts) %dopar% {
      cell.indices <- as.numeric(idx.mat[i, ])
      idx1 <- cell.indices[1]
      idx2 <- cell.indices[2]
      idx3 <- cell.indices[3]
      lower.limits <- c(cuts.list[[1]][idx1], cuts.list[[2]][idx2], cuts.list[[3]][idx3])
      upper.limits <- c(cuts.list[[1]][idx1 + 1], cuts.list[[2]][idx2 + 1], cuts.list[[3]][idx3 + 1])
      
      # Numerical integral
      integral.result <- tryCatch({
        cubature::adaptIntegrate(
          f          = function(x) mvtnorm::dmvnorm(x, mean = mean.vec, sigma = covariance, log = FALSE),
          lowerLimit = lower.limits,
          upperLimit = upper.limits
        )$integral
      }, error = function(e) {
        warning(paste("Integration error for cell (indices:", paste(cell.indices, collapse=","),
                      "),  Error:", e$message), call.=FALSE, immediate.=TRUE)
        return(NA_real_)
      })
      
      if (!is.na(integral.result) && integral.result < 0) {
        stop("integral.result must be positive.")
      }
      return(integral.result)
    }
    parallel::stopCluster(cl)
    if (exists("pb") && inherits(pb, "progress_bar")) try(pb$terminate(), silent=TRUE)
    ########################  END Parallel Processing  #########################
  } else {
    #######################  BEGIN Sequential Processing  ######################
    cell.probabilities <- numeric(n.cells)
    pb <- progress::progress_bar$new(
      format = "  Integrating cells (sequential) [:bar] :percent eta: :eta (:current/:total)",
      total = n.cells, clear = FALSE, width = 80)
  
    for (i in 1:n.cells) {
      cell.indices <- as.numeric(idx.mat[i, ])
      idx1 <- cell.indices[1]
      idx2 <- cell.indices[2]
      idx3 <- cell.indices[3]
      lower.limits <- c(cuts.list[[1]][idx1], cuts.list[[2]][idx2], cuts.list[[3]][idx3])
      upper.limits <- c(cuts.list[[1]][idx1 + 1], cuts.list[[2]][idx2 + 1], cuts.list[[3]][idx3 + 1])

      integral.result <- tryCatch({
        cubature::adaptIntegrate(
          f          = function(x) mvtnorm::dmvnorm(x, mean = mean.vec, sigma = covariance, log = FALSE),
          lowerLimit = lower.limits,
          upperLimit = upper.limits
        )$integral
      }, error = function(e) {
        warning(paste("Integration error for cell (indices:", paste(cell.indices, collapse=","),
                      "),  Error:", e$message), call.=FALSE, immediate.=TRUE)
        return(NA_real_)
      })
      
      if (!is.na(integral.result) && integral.result < 0) {
        stop("integral.result must be positive.")
      }
      cell.probabilities[i] <- integral.result
    }
    if (exists("pb") && inherits(pb, "progress_bar")) try(pb$terminate(), silent=TRUE)
    ########################  END Sequential Processing  #######################
  }
  
  if (any(is.na(cell.probabilities))) {
    stop(paste(sum(is.na(cell.probabilities)), "cell probabilities are NA due to integration errors."), immediate.=TRUE)
  }
  total.prob <- sum(cell.probabilities, na.rm = TRUE) # total probability
  
  if (abs(1-total.prob) > 1e-3 && total.prob > 0) {
    warning(paste("Total probability mass is", formatC(total.prob), "- significantly different from 1. Normalizing probabilities."), immediate. = TRUE)
  } 
  
  cell.prob <- as.matrix(cell.probabilities / total.prob)
  return(cell.prob)
}





###-----------------------------------------------------###
###    Simulation of fitting for normal distribution    ###
###-----------------------------------------------------###

## INPUT:
# p.true    : Numeric vector containing the true probabilities.
# model     : Character string that chooses which sub-model to calculate.
# f.fam     : An f-function stored in the 'f_functions.R'.
# n.cat     : Integer, number of categories of the contingency table.
# score     : A vector of scores attached to each category of a table.
# level     : Significance level (Type I error rate).
# n.freq    : Number of frequency of each contingency table.
# n.rep     : Number of iterations.
# n.cluster : #CPU cores to use for parallel integration.

## OUTPUT:
# Empirical power (i.e. proportion of rejections) for each specified 'model'.
# at the given significance level.

simulate.3dimNormal <- function(p.true, model = c("S", "GS_f", "ELS_f", "LS_f"), 
                                f.fam = NULL, n.cat = NULL, score = NULL, 
                                level = 0.05, n.freq = 1e3, n.rep = 1e3,
                                n.cluster = parallel::detectCores() - 2) 
  {
  n.dim <- 3 # Hardcoded for this specific version
  
  ## Input validation
  if (is.null(n.cat) || !is.element(n.cat, c(4, 6))) {
    stop("'n.cat' must be 4 or 6 for this function.")
  }

  ## Setting for parallel processing
  if (n.cluster > 1) {
    ## Setting for parallel processing
    cl <- makeCluster(n.cluster)
    registerDoSNOW(cl)
    clusterSetRNGStream(cl, 20)
    
    # Export functions and variables to each node
    clusterEvalQ(cl, {
      #library(foreach)
      source('libraries.R')
      source('mph.Rcode.R')
      source('functions.R')
    })
    clusterExport(cl, c("p.true", "n.freq", "n.dim", "n.cat", "score",
                        "f.fam", "model"),
                  envir = environment())
    
    # Display a progress bar
    pb <- progress::progress_bar$new(
      format = "  Fitting process [:bar] :percent eta: :eta (:current/:total)",
      total = n.rep, clear = FALSE, width = 80)
    progress <- function(n_prog) pb$tick()
    opts <- list(progress = progress)
    
    #######################  BEGIN Parallel Processing  ########################
    rej.mat <- foreach(rep = seq_len(n.rep), .combine = cbind, .options.snow = opts) %dopar% {
      tryCatch({
        y <- as.vector(rmultinom(1, n.freq, p.true))
        
        ## Get G^2 for each 'model'
        n.rejection <- numeric(length(model))
        names(n.rejection) <- model
        for (m in model) {
          if (m == "S") {
            L.fct <- identity
            X     <- build.design.mat(model = "S", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "GS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "GS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "ELS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "ELS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "LS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "LS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else {
            stop("model: ", m, " is not implemented")
          }
          
          m.result <- mph.fit(y = y, L.fct = L.fct, X = X)
          p.value <- 1 - pchisq(m.result$Gsq, m.result$df)
          n.rejection[m] <- (p.value < level) * 1
        }
        n.rejection
      }, error = function(e) {
        warning("In rep = ", rep, ", Error: ", conditionMessage(e))
        n.rejection[m] <- NA_real_
        n.rejection
      })
    }
    parallel::stopCluster(cl)
    if (exists("pb") && inherits(pb, "progress_bar")) try(pb$terminate(), silent=TRUE)
    ########################  END Parallel Processing  #########################
  } else {
    #######################  BEGIN Sequential Processing  ######################
    rej.mat <- matrix(NA, nrow = length(model), ncol = n.rep)
    pb <- progress::progress_bar$new(
      format = "  Integrating cells (sequential) [:bar] :percent eta: :eta (:current/:total)",
      total = n.rep, clear = FALSE, width = 80)
    
    for (i in 1:n.rep) {
      rej.mat[ ,i] <- tryCatch({
        y <- as.vector(rmultinom(1, n.freq, p.true))
        
        ## Get G^2 for each 'model'
        n.rejection <- numeric(length(model))
        names(n.rejection) <- model
        for (m in model) {
          if (m == "S") {
            L.fct <- identity
            X     <- build.design.mat(model = "S", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "GS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "GS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "ELS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "ELS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else if (m == "LS_f") {
            L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
            X     <- build.design.mat(model = "LS_f", n.dim = n.dim, n.cat = n.cat, score = score)
          } else {
            stop("model: ", m, " is not implemented")
          }
          
          m.result <- mph.fit(y = y, L.fct = L.fct, X = X)
          p.value <- 1 - pchisq(m.result$Gsq, m.result$df)
          n.rejection[m] <- (p.value < level) * 1
        }
        n.rejection
      }, error = function(e) {
        warning("In rep = ", rep, ", Error: ", conditionMessage(e))
        n.rejection[m] <- NA_real_
        n.rejection
      })
      print(rej.mat[ ,i])
    }
    if (exists("pb") && inherits(pb, "progress_bar")) try(pb$terminate(), silent=TRUE)
    ########################  END Sequential Processing  #######################
  }
  
  ## Output empirical power
  rownames(rej.mat) <- model
  power <- rowMeans(rej.mat)
  return(power)
}

################  END Subroutines or Complementary functions  ##################