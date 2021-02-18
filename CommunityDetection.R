if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc') # is.positive.definite
if (!require('Matrix')) install.packages('Matrix'); library('Matrix') # nearPD
if (!require('quadprog')) install.packages('quadprog'); library('quadprog') # solve.QP
####################################################################
######## Functions of modularity maximization algorithm ############
################### citation: Newman, PNAS, 2006 ###################
####################################################################
# community detection 
## Given the adjancecy matrix: A
FirstTwoPartition <- function(A){
  M <- nrow(A)
  degree <- rowSums(A)
  m <- sum(A)/2
  # get the modularity matrix 
  B <- matrix(0, M, M)
  for(i in 1:M){
    for(j in 1:M){
      B[i, j] <- A[i, j] - degree[i]*degree[j]/(2*m)
    }
  }
  # divide into two groups at the first level
  lambda <- eigen(B)$values
  u <- eigen(B)$vectors[ , which.max(lambda)]
  type.assign <- rep(-1, M)
  type.assign[which(u >= 0)] <- 1
  mod <- (t(type.assign) %*% B %*% type.assign)/(4*m)
  return(list(mod = mod, type.assign = type.assign, B = B, m = m))
} # return mod, type.assign, B, and m
FurtherPartition <- function(Bij, subindex, m){
  ng <- length(subindex)
  Bij.g <- Bij[subindex, subindex]
  delta <- diag(rowSums(Bij.g))
  Bij.g <- Bij.g - delta
  # furture partition
  lambda <- eigen(Bij.g)$values
  u <- eigen(Bij.g)$vectors[ , which.max(lambda)]
  type.assign <- rep(-1, ng)
  type.assign[which(u >= 0)] <- 1
  delta.mod <- (t(type.assign) %*% Bij.g %*% type.assign)/(4*m)
  return(list(delta.mod = delta.mod, type.assign = type.assign))
} # return delta.mod and type.assign
ModularityMaximum <- function(A){
  # step 1
  step1 <- FirstTwoPartition(A)
  mod.first <- step1$mod
  type.first <- step1$type.assign
  Bij <- step1$B
  m <- step1$m
  
  # future steps 
  type <- rep(1, nrow(A))
  type[which(type.first == -1)] <- 2
  undivided <- c()
  for(r in 1:nrow(A)){
    labels <- setdiff(unique(type), undivided)
    if(length(labels) == 0){
      break
    } else {
      for(i in 1:length(labels)){
        subindex <- which(type == labels[i])
        if(length(subindex) > 1){
          partition <- FurtherPartition(Bij, subindex, m)
          deltamod <- partition$delta.mod
          if(deltamod <= 0){
            undivided <- c(undivided, labels[i])
          } else {
            typeassign <- partition$type.assign
            type[subindex[which(typeassign == -1)]] <- max(labels) + 1
          }
        }
      }
    }
  }
  
  return(type)
}
## return type of all agents 
# weighted average 
Weight_Modularity <- function(type, cor.within, cor.across){
  M <- length(type)
  K <- unique(type)
  cor.mat <- matrix(NA, M, M)
  for(k in 1:length(K)){
    index.within <- which(type == K[k])
    index.across <- c(1:M)[-index.within]
    cor.mat[index.within, index.within] <- cor.within
    cor.mat[index.within, index.across] <- cor.across
    cor.mat[index.across, index.within] <- cor.across
  }
  weights <- optw(cor.mat)
  return(weights)
}

optw <- function(cor.mat){
  M <- ncol(cor.mat)
  mat <- cor.mat
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights 
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M), diag(M))
  b <- c(1, rep(0, M))
  d <- rep(1, M)
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1POS!")
  })
  return(weights)
}