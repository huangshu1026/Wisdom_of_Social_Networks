#########################################################################################
######### Bayesian Network Model ##########
#########################################################################################
## 3 ideas: 
#### (1): all combinations from MOD as starting points 
#### (2): save the path of random walk 
#### (3): search larger posterior probability: setting new samples with large posterior as new starting point
# given number of types 
# alpha  <- 2
# connection probability within types: beta 
# connection probability across types: gamma 

## Bayesian Model 
DPlogDensity <- function(type, alpha){
  sorttype <- sort(type)
  type.label <- min(sorttype)
  logprob <- 0
  for(i in 2:length(sorttype)){
    if(sorttype[i] > type.label){
      logprob <- logprob + log(alpha/(alpha + i - 1))
      type.label <- sorttype[i]
    } else {
      nx <- length(which(sorttype[1:(i-1)] == sorttype[i]))
      logprob <- logprob + log(nx/(alpha + i - 1))
    }
  }
  return(logprob)
}
transTmat <- function(type){
  tmat <- matrix(0, length(type), length(type))
  type.num <- max(type)
  for(i in 1:type.num){
    index <- which(type == i)
    tmat[index, index] <- 1
  }
  diag(tmat) <- 1
  return(tmat)
}
llik <- function(A, tmat, beta, gamma){
  M <- nrow(tmat)
  prob <- 0
  for(i in 1:(M-1)){
    for(j in (i+1):M){
      prob <- prob + log(((tmat[i,j]*beta+(1-tmat[i,j])*gamma)^(A[i,j])) * ((1-tmat[i,j]*beta-(1-tmat[i,j])*gamma)^(1-A[i, j])))
    }
  }
  return(prob)
}
## Dirichlet process to generate random partitions 
MultiUrnGenType <- function(M, alpha){
  type <- c(1)
  next.label <- 2
  for (i in 1:(M-1)) {
    if (runif(1,0,1) < alpha/(alpha + i)) {
      # Add a new ball color
      type <- c(type, next.label)
      next.label <- next.label + 1
    } else {
      # Pick out a ball from the urn, and add back a
      # ball of the same color.
      select.label <- type[sample(1:length(type), 1)]
      type <- c(type, select.label)
    }
  }
  return(type)
}
typefind.small <- function(num, M, alpha){
  n <- ceiling(num/length(alpha))
  typeconf <- c()
  for(i in 1:length(alpha)){
    typeconf0 <- t(sapply(1:n, function(x) MultiUrnGenType(M, alpha[i])))
    typeconf <- rbind(typeconf, unique(typeconf0))
  }
  typeconf <- unique(typeconf)
  return(typeconf)
}
## multiple starting points 
MultiStartTypes <- function(StartTypes){
  typelist <- StartTypes
  type.num <- unique(StartTypes)
  if(length(type.num) == 1){
    return(typelist)
  } else {
    for(i in 2:length(type.num)){
      combine <- combn(type.num, i)
      for(j in 1:ncol(combine)){
        index <- which(StartTypes %in% combine[ , j])
        new <- StartTypes
        new[index] <- combine[1, j]
        typelist <- rbind(typelist, new)
      }
    }
    return(typelist)
  }
}
RandomWalk <- function(type.start, n){
  if(max(type.start) == 1){
    num.sample <- n
    type.rw <- c()
    ###### initialization ######
    type.old <- type.start
    num.type <- min(max(type.old), length(type.old)-1)
    ###### random walk ######
    for(t in 1:num.sample){
      # select a person
      personprime <- sample(1:M, size = 1, replace = F)
      # select a type
      groupprime <- sample(c(1:(num.type+1))[-type.old[personprime]], size = 1, replace = F)
      # generate
      typeprime <- type.old
      typeprime[personprime] <- groupprime
      type.rw <- rbind(type.rw, typeprime)
      # calculate ratio 
      logp0.prime <- DPlogDensity(typeprime, alpha)
      loglik.prime <- llik(A, transTmat(typeprime), beta, gamma)
      logp0.old <- DPlogDensity(type.old, alpha)
      loglik.old <- llik(A, transTmat(type.old), beta, gamma)
      ratio <- exp(logp0.prime + loglik.prime - logp0.old - loglik.old)
      ## original MH: accept or reject
      #u <- runif(1, 0, 1)
      #if(u <= ratio){
      #  type.new <- typeprime
      #} else {
      #  type.new <- type.old
      #}
      ## select a larger posterior probabilities 
      if(ratio >= 1){
        type.new <- typeprime
      } else {
        type.new <- type.old
      }
      type.old <- type.new
      num.type <- min(max(type.new), length(type.new)-1)
    }
  } else {
    num.sample <- n/nrow(type.start)
    type.rw <- c()
    for(i in 1:nrow(type.start)){
      ###### initialization ######
      type.old <- type.start[i, ]
      num.type <- min(max(type.old), length(type.old)-1)
      ###### random walk ######
      for(t in 1:num.sample){
        # select a person
        personprime <- sample(1:M, size = 1, replace = F)
        # select a type
        groupprime <- sample(c(1:(num.type+1))[-type.old[personprime]], size = 1, replace = F)
        # generate
        typeprime <- type.old
        typeprime[personprime] <- groupprime
        type.rw <- rbind(type.rw, typeprime)
        # calculate ratio 
        logp0.prime <- DPlogDensity(typeprime, alpha)
        loglik.prime <- llik(A, transTmat(typeprime), beta, gamma)
        logp0.old <- DPlogDensity(type.old, alpha)
        loglik.old <- llik(A, transTmat(type.old), beta, gamma)
        ratio <- exp(logp0.prime + loglik.prime - logp0.old - loglik.old)
        ## original MH: accept or reject
        u <- runif(1, 0, 1)
        if(u <= ratio){
          type.new <- typeprime
        } else {
          type.new <- type.old
        }
        ## select a larger posterior probabilities 
        #if(ratio >= 1){
        #  type.new <- typeprime
        #} else {
        #  type.new <- type.old
        #}
        type.old <- type.new
        num.type <- min(max(type.new), length(type.new)-1)
      }
    }
  }
  return(type.rw)
}
## remove the duplicated types 
GetCluster <- function(typevector){
  typeindex <- unique(typevector)
  cluster <- sapply(1:length(typeindex), function(x) paste(as.character(sort(which(typevector == typeindex[x]))), collapse = ""))
  return(cluster)
}
MHsearch <- function(type.start.piecerow, n, M, alpha, beta, gamma){
  type.rw <- c()
  ###### initialization ######
  type.old <- type.start.piecerow
  num.type <- min(length(unique(type.old)), length(type.old)-1)
  type.name <- unique(type.old)
  ###### random walk ######
  for(t in 1:n){
    # select a person
    personprime <- sample(1:M, size = 1, replace = F)
    # select a type
    groupprime <- sample(c(type.name, max(type.name)+1)[-type.old[personprime]], size = 1, replace = F)
    # generate
    typeprime <- type.old
    typeprime[personprime] <- groupprime
    type.rw <- rbind(type.rw, typeprime)
    # calculate ratio 
    logp0.prime <- DPlogDensity(typeprime, alpha)
    loglik.prime <- max(llik(A, transTmat(typeprime), beta, gamma), -100)
    logp0.old <- DPlogDensity(type.old, alpha)
    loglik.old <- max(llik(A, transTmat(type.old), beta, gamma), -100)
    ratio <- exp(logp0.prime + loglik.prime - logp0.old - loglik.old)
    ## original MH: accept or reject
    u <- runif(1, 0, 1)
    if(u <= ratio){
      type.new <- typeprime
      #type.path <- rbind(type.path, typeprime)
    } else {
      type.new <- type.old
    }
    type.old <- type.new
    num.type <- min(length(unique(type.new)), length(type.new)-1)
    type.name <- unique(type.new)
  }
  type.rw <- unique(type.rw)
  return(type.rw)
}
##### using MCMC method to generate samples #####
EstPost_MCMC <- function(A, n, nDP, alpha, beta, gamma){
  # n is the sample size starting from each starttype vector
  M <- ncol(A)
  ## Initialization: multiple starting points 
  type.start <- MultiStartTypes(ModularityMaximum(A))
  ## Iterations: based on Metropolis-Hastings algorithm 
  num.sample <- n/nrow(type.start) # control the total number of samples is only 400
  #threshold <- 0.5
  type.rw <- apply(type.start, 1, function(x) MHsearch(x, num.sample, M, alpha, beta, gamma))
  type.rw.mat <- c()
  for(i in 1:length(type.rw)){
    type.rw.mat <- rbind(type.rw.mat, type.rw[[i]])
  }
  type.rw <- unique(type.rw.mat)
  ## Dirichlet process the extend the random type configurations 
  type.dp <- typefind.small(nDP, M, alpha = c(seq(0.05, 1, by = 0.05), seq(2, 100, by = 1)))
  ## results: need to remove redandunt types
  typeMCMC <- rbind(type.start, type.rw, type.dp)
  ## remove the redundant type configurations
  typeMCMC <- t(sapply(1:nrow(typeMCMC), function(x) rank(typeMCMC[x, ], ties.method = "min")))
  typeMCMC <- typeMCMC[!duplicated(typeMCMC), ]
  #clusters <- sapply(1:nrow(typeMCMC), function(x) GetCluster(typeMCMC[x, ]))
  #cluster.results <- do.call(rbind,clusters)
  #rownames(cluster.results) <- as.character(1:nrow(typeMCMC))
  #cluster.results <- unique(cluster.results)
  #typeMCMC <- typeMCMC[as.numeric(rownames(cluster.results)), ]
  # calculate probabilities 
  logprior.BAY <- sapply(1:nrow(typeMCMC), function(x) DPlogDensity(typeMCMC[x, ], alpha))
  loglik.BAY <- sapply(1:nrow(typeMCMC), function(x) llik(A, transTmat(typeMCMC[x, ]), beta, gamma))
  logpost.BAY <- logprior.BAY + loglik.BAY
  return(list(type = typeMCMC, logpost = logpost.BAY))
}
transTypetoCormat <- function(type.vec, cor.within, cor.across){
  M <- length(type.vec)
  K <- unique(type.vec)
  cor.mat <- matrix(NA, M, M)
  for(k in 1:length(K)){
    index.within <- which(type.vec == K[k])
    index.across <- c(1:M)[-index.within]
    cor.mat[index.within, index.within] <- cor.within
    cor.mat[index.within, index.across] <- cor.across
    cor.mat[index.across, index.within] <- cor.across
  }
  return(cor.mat)
}
Weight_Bayes <- function(types, logposts, cor.within, cor.across){
  M <- ncol(types)
  cor.array <- array(unlist(apply(types, 1, function(x) transTypetoCormat(x, cor.within, cor.across))), dim = c(M, M, nrow(types)))
  exp.cor.mat <- matrix(0, M, M)
  # update the logposts
  logposts.new <- logposts
  logposts.new[which(logposts.new == -Inf)] <- NA
  logposts.new <- logposts.new - mean(logposts.new[which(is.na(logposts.new) == FALSE)])
  posts <- exp(logposts.new)/sum(na.omit(exp(logposts.new)))
  for(i in which(is.na(posts) == FALSE)){
    exp.cor.mat <- exp.cor.mat + cor.array[ , , i] * posts[i]
  }
  diag(exp.cor.mat) <- 1
  weights <- optw(exp.cor.mat)
  return(weights)
}

#Delta.Postprob$type, Delta.Postprob$logpost, beta.net.hat, gamma.net.hat, cor.within.hat, cor.across.hat