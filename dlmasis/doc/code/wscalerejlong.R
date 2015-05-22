## A set of functions for sampling from the posterior density of the
## local level model using the "wrongly-scaled" parameterizations

## log density of the t location scale family
## x: point to evaluate density at
## mu: location parameter
## sig2: scale parameter sigma^2
## df: degrees of freedom
dtproptilde <- function(x, mu, sig2, df){
  sig <- sqrt(sig2)
  z <- (x-mu)/sig
  out <- dt(z, df, log=TRUE) - log(sig)
  return(out)
}

## simulate from the t location scale family
## n: sample size
## mu: location parameter
## sig2: scale parameter sigma^2
## df: degrees of freedom
rtproptilde <- function(n, mu, sig2, df){
  sig <- sqrt(sig2)
  temp <- rt(n, df=df)
  out <- mu + temp*sig
  return(out)
}

## log of the conditional posterior density of W (V) given V (W), gamma.tilde (psi.tilde)
## VW: point to evaluate the density at
## a, b, cc, avw: parameters of the density
logpiVWtilde <- function(VW, a, b, cc, avw){
  out <-  - avw*VW - a*exp(-VW) + b*exp(-VW/2) - cc*exp(VW)
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
## VW: point to evaluate the density at
## a, b, cc, avw: parameters of the density
logpiVWprimetilde <- function(VW, a, b, cc, avw){
  out <-  - avw + a*exp(-VW) - (b/2)*exp(-VW/2) - cc*exp(VW)
  return(out)
}

## checks whether the conditional posterior of V (W) is log concave
## a, b, cc: parameters of the density
## eps: set > 0 to allow for numerical error
logcontilde <- function(a, b, cc, eps=.01){
  out <- (b <= 0)
  if(!out)
      out <- (a > (b^2/(cc*16))*(1/4 - (b/cc)^2/16^3) + eps)
  return(out)
}

## difference between log conditional posterior of W (V) and the proposal density
## VW: point to evaluate the difference of the densities at
## a, b, cc, avw: parameters of the true density
## mu, sig2, df: parameters of the proposal density
logpirejtilde <- function(VW, a, b, cc, avw, mu, sig2, df){
  out <- logpiVWtilde(VW, a, b, cc, avw) - dtproptilde(VW, mu, sig2, df)
  return(out)
}

## draw W from its conditional posterior given V and gamma.tilde
## dat: time series
## gam: gamma.tilde - wrongly scaled disturbances
## aw, bw: hyperparameters of W's IG prior
Wgamitertilde <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  Wa <- aw + T/2
  Wb <- bw + V*sum(gam[-1]^2)/2
  W <- rinvgamma(1, Wa, Wb)
  return(W)
}

## draw V from its conditional posterior given W and psi.tilde
## dat: time series
## psi: psi.tilde - wrongly scaled errors
## av, bv: hyperparameters of V's IG prior
Vpsiitertilde <- function(dat, psi, W, av, bv){
  T <- length(dat)
  Va <- av + T/2
  Vb <- bv + W*sum(psi[-1]^2)/2
  V <- rinvgamma(1, Va, Vb)
  return(V)
}

## draw V from its conditional posterior given W and gamma.tilde
## dat: time series
## gam: gamma.tilde - wrongly scaled disturbances
## av, bv: hyperparameters of V's IG prior
Vgamitertilde <- function(dat, gam, W, av, bv){
  T <- length(dat)
  gam0 <- gam[1]
  cgams <- cumsum(gam[-1])
  a <- bv + sum((dat - gam0)^2)/2
  b <- sum((dat-gam0)*cgams)
  cc <- cgams[T]^2/(2*W)
  V <- VWrejtilde(a,b,cc,av)
  return(V)
}

## draw W from its conditional posterior given V and psi.tilde
## dat: time series
## psi: psi.tilde - wrongly scaled errors
## aw, bw: hyperparameters of W's IG prior
Wpsiitertilde <- function(dat, psi, V, aw, bw){
  T <- length(dat)
  ys <- c(psi[1], dat)
  psis <- c(0,psi[-1])
  Ly <- ys[-1]-ys[-(T+1)]
  Lpsi <- psis[-1]-psis[-(T+1)]
  a <- bw + sum(Ly^2)/2
  b <- sum(Ly*Lpsi)
  cc <- sum(Lpsi^2)/(2*V)
  W <- VWrejtilde(a,b,cc,aw)
  return(W)
}

## rejection sampler for drawing V (W) given W (V) and gamma.tilde (psi.tilde)
## a, b, cc, avw: parameters of V's (W's) conditional posterior
## df: degrees of freedom for the proposal; 10 seems to work well
VWrejtilde <- function(a, b, cc, avw, df=10){
  ## check for log concavity of target density
  lcon <- logcontilde(a, b, cc)
  lcon <- FALSE
  adrej <- lcon
  ## find mode of target
  mu <- uniroot(logpiVWprimetilde, c(-10^2, 10^2), a=a, b=b, cc=cc, avw=avw)$root
  if(lcon){
    ## if target is log concave, use adaptive rejection sampling
    ## using -10*mu, mu, and 10*mu as starting poitns for the hull
    try(VW <- ars(n=1, logpiVWtilde, logpiVWprimetilde, x=c(-10*mu, mu, 10*mu),
                  a=a, b=b, cc=cc, avw=avw))
    if(VW==0){
      ## error handling - if working properly, this should never happen
      ## but if there is an error, drops out to do the t approximation
      adrej <- FALSE
    }
  }
  if(!adrej){
    ## if target is not log concave, use a t approxmation rejection sampler on log scale
    sig2 <- - 1 /( -a*exp(-mu) +(b/4)*exp(-mu/2) - cc*exp(mu) )
    M <- optimize(logpirejtilde, c(-10^2,10^2), maximum=TRUE, a, b, cc, avw, mu, sig2, df)$objective
    rej <- TRUE
    rejit <- 1
    while(rej){
      VW <- rtproptilde(1, mu, sig2, df)
      R <- logpirejtilde(VW, a, b, cc, avw, mu, sig2, df) - M
      u <- runif(1,0,1)
      if(log(u)<R){
        rej <- FALSE
      }
      rejit <- rejit + 1
      if(rejit%%1000000==0){
        ## for troubleshooting stuck rejection samplers
        print(paste("Rejection sampler appears to be stuck for VW, ",
                    rejit, " iterations", sep=""))
        print(paste("Last value of W in the chain: ", VW, sep=""))
        print(paste("Values of a b and c:", a, b, cc, sep="  "))
        print(paste("T = ", T, sep=""))
        print(paste("Prior value of alphaVW:", avw, sep="  "))
        print(paste("Attempted df: ", df, sep=""))
        print(paste("Attempted location and scale:", mu, sig2, sep="  "))
      }
    }
  }
  return(c(exp(VW), lcon, adrej))
}

## wrongly-scaled error sampler: samples V and W conditional on the wrongly scaled errors
## n: posterior sample size
## start: initial values for V and W
## av, bv, aw, bw: hyperparameters for the priors of V and W
## m0, C0: hyperparameters for the prior on theta_0
werrorsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n)  
  for(i in 1:n){
    ## draw psi using MCFA
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    psi <- c(psi[1], psi[-1]*sqrt(V/W))
    ## draw V from easy full conditional
    V <- Vpsiitertilde(dat, psi, W, av, bv)
    ## draw W from hard full conditional
    Wout <- Wpsiitertilde(dat, psi, V,  aw, bw)
    W <- Wout[1]
    ## transform back to theta
    thetat <- dat - sqrt(W)*psi[-1]
    theta <- c(psi[1], thetat)
    out[i,] <- c(NA,NA,Wout[2:3], V, W)
  }
  return(out)
}


## wrongly-scaled disturbance sampler: samples V and W conditional on the 
## wrongly scaled disturbances
## n: posterior sample size
## start: initial values for V and W
## av, bv, aw, bw: hyperparameters for the priors of V and W
## m0, C0: hyperparameters for the prior on theta_0
wdistsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n)  
  for(i in 1:n){
    ## draw theta using MCFA and transform to gamma
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gamt <- (theta[-1] - theta[-(T+1)])/sqrt(V)
    gam <- c(theta[1], gamt)
    ## draw V from hard full conditional
    Vout <- Vgamitertilde(dat, gam, W, av, bv)
    V <- Vout[1]
    ## draw W from easy full conditional
    W <- Wgamitertilde(dat, gam, V, aw, bw)
    ## transform back to theta
    thetat <- gam[1] + cumsum(gam[-1])*sqrt(V)
    theta <- c(gam[1], thetat)
    out[i,] <- c(Vout[2:3],NA,NA,V,W)
  }
  return(out)
}

