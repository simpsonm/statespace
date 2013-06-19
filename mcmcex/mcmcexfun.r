## A set of functions for simulated from and fitting local level models

library(dlm)
library(coda)
library(MCMCpack)
library(ars)

## scaled disturbance sampler: samples V and W conditional on system disturbances
## (plus the initial state, theta_0)
distsam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  for(i in 1:n){
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - gam0 - sqrt(W)*cgam )^2 )/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    sb <- 1*(b>0) - 1*(b<0)
    mn <- optimize(distlogW, c(0,10^10), maximum=TRUE, a, b, a2)$maximum
    if(b^2 < (a2 + 1)^3*(1 - 2/3*sb)*32/9/b2 + .1){
      ## If true, density isn't log concave. + .1 to make sure ARS
      ## algorithm doesn't fail on near non-log-concave cases
      propvar <- - 1 /( (a2 + 1)*mn^(-2) - b*mn^(-3/2)/4 - 2*b2*mn^(-3) )
      d <- optimize(propM, c(1,10^10), maximum=FALSE,  a, b, a2, mn, propvar)
      M <- d$objective
      df <- d$minimum
      rej <- TRUE
      while(rej){
        prop <- rtprop(1, mn, propvar, df)
        if(prop>0){
          R <- distrej(prop, a, b, a2, mn, propvar, df) - M
          u <- runif(1,0,1)
          if(log(u)<R){
            W <- prop
            rej <- FALSE
          }
        }
      }
    }
    else{
      W <- ars(n=1, distlogW, distlogWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a2=a2)
    }
    A <- B/sqrt(W)
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    out[i,] <- c(theta,V,W)
  }
  return(out)
}

## compute M in logtarget - logprop <= M
propM <- function(df, a, b, a2, mn, propvar){
  M <- optimize(distrej, c(0,10^10), maximum=TRUE, a, b, a2, mn, propvar, df)
  return(M$objective)
}

## log density of the t location scale family
dtprop <- function(df, x, mn, var){
  sd <- sqrt(var)
  z <- (x-mn)/sd
  out <- dt(z, df, log=TRUE) - log(sd)
  return(out)
}

## simulate from the t location scale family
rtprop <- function(n, mn, var, df){
  sd <- sqrt(var)
  temp <- rt(n, df=df)
  out <- mn + temp*sd
  return(out)
}

## log of the conditional posterior density of W given V, gamma, data
distlogW <- function(W, a, b, a2){
  out <- -a*W + b*sqrt(W) - (a2 + 1)*log(W) - b2/W
  return(out)
}

## first derivative of log of the conditional posterior of W
distlogWprime <- function(W, a, b, a2){
  out <- -a + b/2/sqrt(W) - (a2 + 1)/W + b2/(W^2)
  return(out)
}

## difference between log conditional posterior of W and the proposal density
distrej <- function(W, a, b, a2, mn, propvar, df){
  out <- distlogW(W, a, b, a2) - dtprop(W, mn, propvar, df)
  return(out)
}


## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  for(i in 1:n){
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-T1])^2)/2)
    out[i,] <- c(theta,V,W)
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


