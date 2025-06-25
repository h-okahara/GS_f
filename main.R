#
# Sourcing this R file contains 5 parts.
#
# import & setting, mph.fit, parameter estimation, simulation
#
######################  BEGIN import & setting  ################################

source("libraries.R")
source("mph.Rcode.R")
source("functions.R")
source("database.R")
source("f_divergence.R")
source("constraints.R")

# Preparation
table <- database$party
n.dim <- length(dim(table)) # dimensionality
n.cat <- dim(table)[1]      # the number of category
score <- c(1:n.cat)         # equal interval score

df <- as.data.frame.table(table,
                          base = list(as.character(1:max(dim(table)))),
                          responseName = "value")
var.col <- paste0("Var", seq_len(n.dim))
df.lex   <- df[do.call(order, df[var.col]), ] # lexicographic re-order
rownames(df.lex) <- NULL 

# Choose an model to test from "GS_f", "ELS_f" or "LS_f"
# Choose an f divegence
lambda <- 0     # KL-divergence
# lambda <- 1     # χ-squared divergence
# lambda <- -1/2  # Hellinger distance
f.fam <- get.f_function(name = "power", lambda = lambda)

######################  END import & setting  ##################################


###########################  BEGIN mph.fit  ####################################

# the S model
X_S <- build.design.mat(model = "S", n.dim = n.dim, n.cat = n.cat, score = score)
S.result <- mph.fit(y = df.lex$value, L.fct = identity, X = X_S)
mph.summary(S.result, cell.stats = n.dim)


# the GS[f] model
L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
X.mat <- build.design.mat(model = "GS_f", n.dim = n.dim, n.cat = n.cat, score = score)
GS_f.result <- mph.fit(y = df.lex$value, L.fct = L.fct, X = X.mat)
mph.summary(GS_f.result, cell.stats = n.dim)

# ME_2 model
h.fct <- get.constraints(model = "ME_2", score = score, n.dim = n.dim, n.cat = n.cat)
ME_2.result <- mph.fit(y = df.lex$value, h.fct = h.fct)
mph.summary(ME_2.result, cell.stats = n.dim)

# ME model
h.fct <- get.constraints(model = "ME", score = score, n.dim = n.dim, n.cat = n.cat)
ME.result <- mph.fit(y = df.lex$value, h.fct = h.fct)
mph.summary(ME.result, cell.stats = n.dim)

# VE model
h.fct <- get.constraints(model = "VE", score = score, n.dim = n.dim, n.cat = n.cat)
VE.result <- mph.fit(y = df.lex$value, h.fct = h.fct)
mph.summary(VE.result, cell.stats = n.dim)

# CE model
h.fct <- get.constraints(model = "CE", score = score, n.dim = n.dim, n.cat = n.cat)
CE.result <- mph.fit(y = df.lex$value, h.fct = h.fct)
mph.summary(CE.result, cell.stats = n.dim)

###########################  END mph.fit  ######################################



########################  BEGIN parameter estimation  ##########################

# Compute plug-in estimater of potential parameters of the 'model'
model <- "GS_f"
L.fct <- function(p) Calc.F.vec(p, F.func = f.fam$F.func, n.dim = n.dim, n.cat = n.cat)
X.mat <- build.design.mat(model = model, n.dim = n.dim, n.cat = n.cat, score = score)
GS_f.result <- mph.fit(y = df.lex$value, L.fct = L.fct, X = X.mat)
theta.hat <- Calc.theta.vec(model = model, theta.vec = GS_f.result$beta, lambda = lambda, 
                            n.dim = n.dim, n.cat = n.cat, score = score)

theta.hat[2] - theta.hat[4] # (1,1,2) - (1,2,1)
theta.hat[11] - theta.hat[13] # (2,1,2) - (2,2,1)
theta.hat[6] - theta.hat[8] # (1,2,3) - (1,3,2)
theta.hat[25] - theta.hat[21] # (3,3,1) - (3,1,3)
theta.hat[3] / theta.hat[7] # (1,1,3) - (1,3,1)
theta.hat[15] - theta.hat[17] # (2,2,3) - (2,3,2)
theta.hat[26] - theta.hat[24] # (3,3,2) - (3,2,3)

#########################  END parameter estimation  ###########################


###########################  BEGIN simulation  #################################

scenarios <- list(
  list(name  = "Scenario 1", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0), var = c(1, 1, 1),     corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 2", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0), var = c(1, 1.2, 1.4), corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 3", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0), var = c(1, 1, 1),     corr = c(0.2, 0.3, 0.4)),
  list(name  = "Scenario 4", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0), var = c(1, 1.2, 1.4), corr = c(0.2, 0.3, 0.4)),
  list(name  = "Scenario 5", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0.1), var = c(1, 1, 1),     corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 6", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0.1), var = c(1, 1.2, 1.4), corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 7", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0.1), var = c(1, 1, 1),     corr = c(0.2, 0.3, 0.4)),
  list(name  = "Scenario 8", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, 0, 0.1), var = c(1, 1.2, 1.4), corr = c(0.2, 0.3, 0.4)),
  list(name  = "Scenario 9",  n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, -0.1, 0.1),  var = c(1, 1, 1),     corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 10", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, -0.1, 0.1),  var = c(1, 1.2, 1.4), corr = c(0.2, 0.2, 0.2)),
  list(name  = "Scenario 11", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, -0.1, 0.1),  var = c(1, 1, 1),     corr = c(0.2, 0.3, 0.4)),
  list(name  = "Scenario 12", n.dim = 3, n.cat = 4, score = 1:4, mean = c(0, -0.1, 0.1),  var = c(1, 1.2, 1.4), corr = c(0.2, 0.3, 0.4))
)

# Choose an f divegence
lambda <- 0     # KL-divergence
# lambda <- 1     # χ-squared divergence
# lambda <- -1/2  # Hellinger distance
f.fam <- get.f_function(name = "power", lambda = lambda)

# Scenario 1
for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  print(scenario$name)
  covariance <- build.covariance(scenario$var, scenario$corr)
  p.true <- generate.3dimNormal(mean.vec = scenario$mean, covariance = covariance, n.cat = scenario$n.cat)
  power <- simulate.3dimNormal(p.true, model = c("S", "GS_f", "ELS_f", "LS_f"), 
                               f.fam = f.fam, n.cat = scenario$n.cat, score = scenario$score, 
                               level = 0.05, n.freq = 10000, n.rep = 10000)
  print(power)
}

############################  END simulation  ##################################