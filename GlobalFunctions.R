# additional functions 
TypeVecToMat <- function(type.vec){
  K <- unique(type.vec)
  N <- length(type.vec)
  type.mat <- matrix(0, N, N)
  for(i in 1:length(K)){
    index.withintype <- which(type.vec == K[i])
    type.mat[index.withintype, index.withintype] <- 1
  }
  return(type.mat)
}
GenNet <- function(beta.net, gamma.net, type.mat){
  mat0.beta <- type.mat * beta.net
  mat0.gamma <- (1-type.mat) * gamma.net
  mat0 <- mat0.beta + mat0.gamma
  N <- ncol(type.mat)
  mat1 <- matrix(runif(N*N, 0, 1), N, N)
  mat1[lower.tri(mat1)] <- t(mat1)[lower.tri(t(mat1))]
  A <- matrix(as.numeric(mat1 < mat0), N, N, byrow = FALSE)
  diag(A) <- 0
  return(A)
}
ApproxPositiveDefinite <- function(matrix){
  if(is.positive.definite(matrix) == FALSE){
    matrix <- nearPD(matrix)$mat
  } 
  return(as.matrix(matrix))
}
GetCov <- function(judge.cor, judge.rhoxy, judge.var, target.var){
  judge.cormat <- rbind(cbind(judge.cor, judge.rhoxy), c(judge.rhoxy, 1))
  judge.sdvec <- sqrt(c(judge.var, target.var))
  judge.covmat <- diag(judge.sdvec) %*% judge.cormat %*% diag(judge.sdvec)
  judge.covmat[lower.tri(judge.covmat)] <- t(judge.covmat)[lower.tri(t(judge.covmat))]
  judge.Sigma <- ApproxPositiveDefinite(judge.covmat)
  return(judge.Sigma)
}
