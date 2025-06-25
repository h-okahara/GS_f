#
# Sourcing this R file (> source("libraries.R") ) contains libraries used in main.R.
#
###########################  BEGIN import libraries  ###########################

library(Rsolnp)   # Non-linear optimisation for constrained problems
library(MASS)     # Statistical functions
library(foreach)  # Looping construct that works with parallel backends
library(doSNOW)   # ‘snow’ backend for foreach
library(parallel) # Parallel infrastructure
library(progress) # Progress bars in parallel processing

############################  END import libraries  ############################