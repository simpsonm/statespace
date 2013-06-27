library(dlm)
library(MCMCpack)
library(coda)


## function for creating the dlm model object to feed into MLE algorithm
## note: won't check inputs - make sure q is specified correctly
## note2: won't allow you to order the evolution variances if there are
## greater than 1 but less than q of them.
mymodMLE <- function(parm, other){
  V <- exp(parm[1])
  order <- other[1]
  s <- other[2]
  q <- other[3]
  ##ev <- other[4]
  W1 <- exp(parm[2:(order+1)])
  W2 <- exp(parm[(order+2):length(parm)])
  nW2 <- length(W2)
  repW2 <- q/nW2
  W2star <- rep(W2, each=2*repW2)
  mod <- dlmModPoly(order, dV=V, dW=W1) + dlmModTrig(s=s, q=q, dV=0, dW=W2star)
  return(mod)
}

## function for creating the dlm model object for non MLE purposes
## note: won't check inputs - make sure q is specified correctly
## note2: won't allow you to order the evolution variances if there are
## greater than 1 but less than q of them.
mymod <- function(parm, other){
  V <- parm[1]
  order <- other[1]
  s <- other[2]
  q <- other[3]
  ##ev <- other[4]
  W1 <- parm[2:(order+1)]
  W2 <- parm[(order+2):length(parm)]
  nW2 <- length(W2)
  repW2 <- q/nW2
  W2star <- rep(W2, each=2*repW2)
  mod <- dlmModPoly(order, dV=V, dW=W1) + dlmModTrig(s=s, q=q, dV=0, dW=W2star)
  return(mod)
}

## sample from the posterior of model 3b, IG prior
postsamIG <- function(n, data, prior, start){
  av <- prior[1]
  bv <- prior[2]
  am <- prior[3]
  bm <- prior[4]
  ab <- prior[5]
  bb <- prior[6]
  as <- prior[7]
  bs <- prior[8]
  parm <- start
  T <- length(data)
  Va <- av+T/2
  W1a <- am + T/2
  W2a <- ab + T/2
  W3a <- as + 3*T
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
    Vp <- rinvgamma(1, Va, bv+0.5*SSY)
    W1p <- rinvgamma(1, W1a, bm + 0.5*SSTH[1])
    W2p <- rinvgamma(1, W2a, bb + 0.5*SSTH[2])
    W3p <- rinvgamma(1, W3a, bs + 0.5*sum(SSTH[-c(1,2)]))
    parm <- c(Vp, W1p, W2p, W3p)
    out[i,] <- parm
  }
  return(out)
}


## sample from the posterior of model 3b, uniform prior, rejection steps
postsamUR <- function(n, data, prior, start){
  KV <- prior[1]
  Km <- prior[2]
  Kb <- prior[3]
  Ks <- prior[4]
  parm <- start
  T <- length(data)
  out <- mcmc(matrix(0,ncol=4, nrow=n))
  colnames(out) <- c("V", "Wmu", "Wbeta", "Wseas")
  for(i in 1:n){
    print(i)
    mod <- mymod(parm, c(2,24,3))
    FF <- mod$FF
    GG <- mod$GG
    filt <- dlmFilter(data, mod)
    theta <- dlmBSample(filt)
    SSY <- sum((data - FF%*%t(theta[-1,]))^2)
    SSTH <- apply((t(theta[-1,]) - GG%*%t(theta[-(T+1),]))^2,1,sum)

    rej <- TRUE
    while(rej){
      Vp <- rinvgamma(1, Va, 0.5*SSY)
      if(Vp<KV)
          rej <- FALSE
    }

    rej <- TRUE
    while(rej){
      W1p <- rinvgamma(1, W1a, 0.5*SSTH[1])
      if(W1p<Km)
          rej <- FALSE
    }

    rej <- TRUE
    while(rej){
      W2p <- rinvgamma(1, W2a, 0.5*SSTH[2])
      if(W2p<Kb)
          rej <- FALSE
    }

    rej <- TRUE
    while(rej){
      W3p <- rinvgamma(1, W3a, 0.5*sum(SSTH[-c(1,2)]))
      if(W3p<Ks)
          rej <- FALSE
    }
    parm <- c(Vp, W1p, W2p, W3p)
    out[i,] <- parm
  }
  return(out)
}

## sample from the posterior of model 3b, uniform prior, metropolis steps
postsamUM <- function(n, data, prior, start){
  KV <- prior[1]
  Km <- prior[2]
  Kb <- prior[3]
  Ks <- prior[4]
  parm <- start
  T <- length(data)
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

    Vprop <- rinvgamma(1, Va, 0.5*SSY)
    if(Vprop<KV)
        Vp <- Vprop

    W1prop <- rinvgamma(1, W1a, 0.5*SSTH[1])
    if(W1prop<Km)
        W1p <- W1prop

    W2prop <- rinvgamma(1, W2a, 0.5*SSTH[2])
    if(W2prop<Kb)
        W2p <- W2prop

    W3prop <- rinvgamma(1, W3a, 0.5*sum(SSTH[-c(1,2)]))
    if(W3prop<Ks)
        W3p <- W3prop

    parm <- c(Vp, W1p, W2p, W3p)
    out[i,] <- parm
  }
  return(out)
}


## sample from the posterior of model 3b, half-T prior
postsamHT <- function(n, data, prior, start){
  vv <- prior[1]
  vm <- prior[2]
  vb <- prior[3]
  vs <- prior[4]
  Av <- prior[5]
  Am <- prior[6]
  Ab <- prior[7]
  As <- prior[8]

  parm <- start
  V <- parm[1]
  Wm <- parm[2]
  Wb <- parm[3]
  Ws <- parm[4]

  aV <- (T+vv)/2
  aWm <- (T+vm)/2
  aWb <- (T+vb)/2
  aWs <- 3*T + vs/2

  T <- length(data)
  out <- mcmc(matrix(0,ncol=4, nrow=n))
  colnames(out) <- c("V", "Wmu", "Wbeta", "Wseas")
  for(i in 1:n){
    print(c(i,round(parm,10)))
    mod <- mymod(parm, c(2,24,3))
    FF <- mod$FF
    GG <- mod$GG
    filt <- dlmFilter(data, mod)
    theta <- dlmBSample(filt)
    SSY <- sum((data - FF%*%t(theta[-1,]))^2)
    SSTH <- apply((t(theta[-1,]) - GG%*%t(theta[-(T+1),]))^2,1,sum)

    av <- rinvgamma(1, 1/2, 1/Av^2 + vv/V )
    V <- rinvgamma(1, aV, vv/av + 0.5*SSY )

    am <- rinvgamma(1, 1/2, 1/Am^2 + vm/Wm )
    Wm <- rinvgamma(1, aWm, vm/am + 0.5*SSTH[1] )

    ab <- rinvgamma(1, 1/2, 1/Ab^2 + vb/Wb )
    Wb <- rinvgamma(1, aWb, vb/ab + 0.5*SSTH[2] )

    as <- rinvgamma(1, 1/2, 1/As^2 + vs/Ws )
    Ws <- rinvgamma(1, aWs, vs/as + 0.5*sum(SSTH[3:8]) )

    parm <- c(V, Wm, Wb, Ws)
    out[i,] <- parm
  }
  return(out)
}


## Simulates from a local level model
llsim <- function(T, V, W, m0, C0){
  out <- rep(0,T)
  theta <- rnorm(1, m0, sqrt(C0))
  for(t in 1:T){
    theta <- theta + rnorm(1, 0, sd=sqrt(W))
    out[t] <- theta + rnorm(1, 0, sd=sqrt(V))
  }
  return(out)
}


statesam <- function(n, start, dat, prior){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=2))
  colnames(out) <- c("V","W")
  vv <- prior[1]
  vw <- prior[2]
  Av <- prior[3]
  Aw <- prior[4]
  for(i in 1:n){
    print(c(i,V,W))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    av <- rinvgamma(1, 1/2, 1/(Av^2) + vv/V)
    V <- rinvgamma(1,(vv + T)/2, vv/av + 0.5*sum((dat - theta[-1])^2))
    aw <- rinvgamma(1, 1/2, 1/(Aw^2) + vw/W)
    W <- rinvgamma(1, (vw + T)/2, vw/aw + 0.5*sum((theta[-1]-theta[-(T+1)])^2))
    out[i,] <- c(V,W)
  }
    return(out)
}
