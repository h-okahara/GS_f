#
# Sourcing this R file contains f-functions used in main.R.
# 
# Create an environment to store all f-functions
  f_divergence <- new.env(parent = emptyenv())
#
#
########################  BEGIN define f divergence  ###########################

###------------------------------------------------###
###        Kullback–Leibler (KL) divergence        ###
###------------------------------------------------###

f_divergence$kl <- list(
  f.func = function(x) x * log(x) - x + 1,
  F.func = function(x) log(x)
  )


###-------------------------------------###
###        Reverse KL divergence        ###
###-------------------------------------###

f_divergence$revkl <- list(
  f.func = function(x) -log(x) + x - 1,
  F.func = function(x) -1 / x + 1
  )


###---------------------------------------------###
###        Pearsonian's chi^2 divergence        ###
###---------------------------------------------###

f_divergence$pearson <- list(
  f.func = function(x) (x - 1)^2 / 2,
  F.func = function(x) x - 1
  )


###------------------------------------###
###        Hellinger divergence        ###
###------------------------------------###

f_divergence$hellinger <- list(
  f.func = function(x) 2 * (sqrt(x) - 1)^2,
  F.func = function(x) 2 - 2 / sqrt(x)
  )


###---------------------------------------###
###        Power divergence family        ###
###---------------------------------------###

f_divergence$power <- function(lambda = NULL) {
  force(lambda)
  if (lambda == 0) return(f_divergence$kl)
  if (lambda == -1) return(f_divergence$revkl)
  
  list(
    f.func = function(x)
      (x^(lambda + 1) - x) / (lambda * (lambda + 1)) - (x - 1) / (lambda + 1),
    F.func = function(x) (x^lambda - 1) / lambda
  )
}


###-------------------------------------------###
###        alpha divergence divergence        ###
###-------------------------------------------###

f_divergence$alpha <- function(alpha = NULL) {
  force(alpha)
  if (alpha ==  1) return(f_divergence$kl)
  if (alpha == -1) return(f_divergence$revkl)
  
  const <- 4 / (1 - alpha^2)
  a  <- (1 - alpha) / 2
  b  <- (1 + alpha) / 2
  
  force(const); force(b); force(b) 
  list(
    f.func = function(x) const * (a + b * x - x^b),
    F.func = function(x) 2 * (x^((alpha - 1) / 2) - 1) / (1 - alpha)
  )
}

#########################  END define f_divergence  ############################