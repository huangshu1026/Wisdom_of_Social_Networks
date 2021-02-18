##### Estimating empirical correlations within types and across types #####
EstCorrelation <- function(type.vec, beta.net, gamma.net, judge.mu, judge.Sigma, num.A){
  # number of judges 
  M <- length(type.vec)
  # network generation 
  type.mat <- TypeVecToMat(type.vec)
  A.array <- array(unlist(sapply(1:num.A, function(x) GenNet(beta.net, gamma.net, type.mat))), dim = c(M, M, num.A))
  # generate judgments 
  X <- mvrnorm(num.A, judge.mu, judge.Sigma)
  # estimating correlation within types and across types 
  est.rho <- matrix(unlist(sapply(1:num.A, function(x) GetRho(A.array[ , , x], X[x, ]))), ncol = 2, nrow = num.A, byrow = T)
  rho.hat <- colMeans(est.rho)
  return(list(rho.within.hat = rho.hat[1], rho.across.hat = rho.hat[2]))
}
# for each network and each round of judgments 
GetRho <- function(A, X1){
  M <- length(X1)
  pair.connect <- c()
  pair.unconnect <- c()
  for(i in 1:M){
    neighbor <- which(A[i, ] == 1)
    stranger <- which(A[i, ] == 0)
    pair.connect <- rbind(pair.connect, cbind(rep(X1[i], length(neighbor)), X1[neighbor]))
    pair.unconnect <- rbind(pair.unconnect, cbind(rep(X1[i], length(stranger)), X1[stranger]))
  }
  rho.within <- cor(pair.connect[ , 1], pair.connect[ , 2])
  rho.across <- cor(pair.unconnect[ , 1], pair.unconnect[ , 2])
  return(c(rho.within, rho.across))
}



