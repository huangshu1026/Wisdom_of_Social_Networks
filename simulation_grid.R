if (!require('invgamma')) install.packages('invgamma'); library('invgamma') # rinvgamma
if (!require('MASS')) install.packages('MASS'); library('MASS') # mvrnorm
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc') # is.positive.definite
if (!require('MBESS')) install.packages('MBESS'); library('MBESS') # cor2cov
# set up the path and functions 
setwd("Document/WOC_SN_otherdata/ABS/")
source("BayesianNetModel.R")
source("CommunityDetection.R")
source("ExternalEstimation.R")
source("GlobalFunctions.R")

##### step 0: prepare the parameters #####
K <- seq(2, 10, by = 2) # number of types
N <- seq(20, 100, by = 20) # number of judges
## parameter of target value
target.mu <- 0
target.var <- 1
## hyper parameters for judges characteristics
mean.mu <- c(0, 0.5, 1)
sigma.mu <- c(0.5, 1, 2)
alpha.sigma2 <- c(0.5, 1, 2)
beta.sigma2 <- 1
alpha.rhoxy <- c(0.5, 1, 2)
beta.rhoxy <- 1
## connection probabilities
beta.net <- c(0.5, 0.7, 0.9)
gamma.net <- c(0.1, 0.2, 0.3)
## tolerance of connection probabilities
beta.net.hat <- c(0.4, 0.6, 0.8, 1)
gamma.net.hat <- c(0, 0.1, 0.2, 0.3)
## set up the correlations within type and across type
rho.delta0 <- c(-0.1, 0, 0.1, 0.2)
rho.delta1 <- c(0.3, 0.5, 0.7, 0.9)
## estimated correlation of judgments: need to estimate additionally
par <- expand.grid(K, N, target.mu, target.var, mean.mu, 
                   sigma.mu, alpha.sigma2, beta.sigma2, alpha.rhoxy, beta.rhoxy, 
                   beta.net, gamma.net, beta.net.hat, gamma.net.hat, rho.delta0, 
                   rho.delta1)
##### step 1: run the simulation #####
cor.within.hat.vec <- c()
cor.across.hat.vec <- c()
mse.allcase <- c()
for(i in 20:nrow(par)){
  parameters <- unlist(par[i, ])
  ##### step 0: set up parameter values #####
  K <- parameters[1]
  N <- parameters[2]
  target.mu <- parameters[3]
  target.var <- parameters[4]
  mean.mu <- parameters[5]
  sigma.mu <- parameters[6]
  alpha.sigma2 <- parameters[7]
  beta.sigma2 <- parameters[8]
  alpha.rhoxy <- parameters[9]
  beta.rhoxy <- parameters[10]
  beta.net <- parameters[11]
  gamma.net <- parameters[12]
  beta.net.hat <- parameters[13]
  gamma.net.hat <- parameters[14]
  rho.delta0 <- parameters[15]
  rho.delta1 <- parameters[16]
  num.A <- 100
  n <- 1000
  nDP <- 10000
  ##### step 1: set up the type information #####
  # determine the type vector 
  set.seed(123)
  type.vec <- sample(1:K, N, replace = TRUE)
  type.mat <- TypeVecToMat(type.vec)
  # type characteristics 
  type.mu <- rnorm(K, mean = mean.mu, sd = sigma.mu)
  type.var <- rinvgamma(K, shape = alpha.sigma2, scale = beta.sigma2)
  type.rhoxy <- rbeta(K, alpha.rhoxy, beta.rhoxy)
  ##### step 2: generate social networks #####
  A <- GenNet(beta.net, gamma.net, type.mat)
  ##### step 3: generate judgments #####
  # parameter for judges 
  judge.mu <- type.mu[type.vec]
  judge.var <- type.var[type.vec]
  judge.rhoxy <- type.rhoxy[type.vec]
  judge.cor <- type.mat * rho.delta1 + (1 - type.mat) * rho.delta0
  # generate judgments 
  ## assuming normal distribution 
  judge.env.Sigma <- GetCov(judge.cor, judge.rhoxy, judge.var, target.var)
  judge.env.Mu <- c(judge.mu, target.mu)
  obs1 <- mvrnorm(1000, judge.env.Mu, judge.env.Sigma)
  X1 <- obs1[, 1:N]
  y1 <- obs1[, N+1]
  ##### step 4: calcualte weights #####
  # estimating correlations within types and across types 
  rho.hat <- EstCorrelation(type.vec, beta.net, gamma.net, judge.mu, judge.env.Sigma[1:N, 1:N], num.A)
  cor.within.hat <- rho.hat$rho.within.hat
  cor.across.hat <- rho.hat$rho.across.hat
  cor.within.hat.vec[i] <- cor.within.hat
  cor.across.hat.vec[i] <- cor.across.hat
  # Bayesian Model 
  Delta.Postprob <- EstPost_MCMC(A, n, nDP, alpha = K, beta = beta.net.hat, gamma = gamma.net.hat)
  weights.bayes <- Weight_Bayes(matrix(unlist(Delta.Postprob$type), nrow = length(Delta.Postprob$logpost), ncol = N),
                                unlist(Delta.Postprob$logpost), 
                                cor.within.hat, cor.across.hat)
  # modularity maximization 
  est.type.modmax <- ModularityMaximum(A)
  weights.modmax <- Weight_Modularity(est.type.modmax, cor.within.hat, cor.across.hat)
  ##### step 5: comparing the simple average and the weighted average #####
  SA <- rowMeans(X1)
  WA.bayes <- X1 %*% weights.bayes
  WA.modmax <- X1 %*% weights.modmax
  agg <- cbind(SA, WA.bayes, WA.modmax)
  # compute the MSE 
  mse <- colMeans((agg - y1)^2)
  test1 <- t.test(agg[ , 1], agg[ , 2], paired = T)
  test2 <- t.test(agg[ , 1], agg[ , 3], paired = T)
  results <- c(parameters, cor.within.hat, cor.across.hat, mse, test1$p.value, test2$p.value) 
  write.table(t(results), "simulation_20210217.csv", row.names = F, col.names = F, append = T, sep = ",")
  mse.allcase <- rbind(mse.allcase, mse)
  #if(i %% 100 == 0){print(i)}
  print(i)
}
