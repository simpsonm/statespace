## functions for replicating previous work done on this data

library(dlm)
library(MCMCpack)
library(coda)

## function for creating the dlm model object to feed into MLE algorithm
## definitely works for model 3b... not sure about others
mymodMLE <- function(parm, other){
  V <- exp(parm[1])
  order <- other[1]
  s <- other[2]
  q <- other[3]
  ##ev <- other[4]
  W1 <- exp(parm[2:(order+1)])
  W2 <- exp(parm[(order+2):length(parm)])
  mod <- dlmModPoly(order, dV=V, dW=W1) + dlmModTrig(s=s, q=q, dV=0, dW=W2)
  return(mod)
}

## function for creating the dlm model object for non MLE purposes
## definitely works for model 3b... not sure about others
mymod <- function(parm, other){
  V <- parm[1]
  order <- other[1]
  s <- other[2]
  q <- other[3]
  ##ev <- other[4]
  W1 <- parm[2:(order+1)]
  W2 <- parm[(order+2):length(parm)]
  mod <- dlmModPoly(order, dV=V, dW=W1) + dlmModTrig(s=s, q=q, dV=0, dW=W2)
  return(mod)
}

## sample from the posterior of model 3b
postsam <- function(n, data, prior, start){
  ay <- prior[1]
  by <- prior[2]
  ath1 <- prior[3]
  bth1 <- prior[4]
  ath2 <- prior[5]
  bth2 <- prior[6]
  ath3 <- prior[7]
  bth3 <- prior[8]
  parm <- start
  T <- length(data)
  Va <- ay+T/2
  W1a <- ath1 + T/2
  W2a <- ath2 + T/2
  W3a <- ath3 + 3*T
  out <- mcmc(matrix(0,ncol=4, nrow=n))
  colnames(out) <- c("V", "Wmu", "Wbeta", "Wseas")
  for(i in 1:n){
    mod <- mymod(parm, c(2,24,3))
    FF <- mod$FF
    GG <- mod$GG
    filt <- dlmFilter(data, mod)
    theta <- dlmBSample(filt)
    SSY <- sum((data - FF%*%t(theta[-1,]))^2)
    SSTH <- apply((t(theta[-1,]) - GG%*%t(theta[-(T+1),]))^2,1,sum)
    Vp <- rgamma(1, Va, rate=by+0.5*SSY)  # note everything on precision scale here
    W1p <- rgamma(1, W1a, rate=bth1 + 0.5*SSTH[1])
    W2p <- rgamma(1, W2a, rate=bth2 + 0.5*SSTH[2])
    W3p <- rgamma(1, W3a, rate=bth3 + 0.5*sum(SSTH[-c(1,2)]))
    parm <- c(1/Vp, 1/W1p, 1/W2p, 1/W3p)
    out[i,] <- parm
  }
  return(out)
}

