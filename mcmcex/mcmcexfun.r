## A set of functions for simulated from and fitting local level models

library(dlm)
library(coda)
library(MCMCpack)
library(ars)

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

## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, time=FALSE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")

  if(time){
    times <- matrix(0, nrow=n, ncol=3)
    colnames(times) <- c("states", "V", "W")
    itime <- rep(0,4)
  }

  for(i in 1:n){
    if(time){
      itime[1] <- proc.time()[3]
    }

    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)

    if(time){
      itime[2] <- proc.time()[3]
    }

    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)

    if(time){
      itime[3] <- proc.time()[3]
    }

    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-(T+1)])^2)/2)

    if(time){
      itime[4] <- proc.time()[3]
      times[i,] <- itime[-1] - itime[-4]
    }

    out[i,] <- c(theta,V,W)
  }

  if(time){
    return(list(out, times))
  }
  else{
    return(out)
  }
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, time=FALSE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")

  if(time){
    times <- matrix(0, nrow=n, ncol=4)
    itime <- rep(0,5)
    rejsam <- rep(0,n)
  }

  for(i in 1:n){
    if(time){
      itime[1] <- proc.time()[1]
    }

    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    psi <- ytild + A%*%theta

    if(time){
      itime[2] <- proc.time()[1]
    }

    Wa <- a2 + T/2
    Wb <- b2 + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)

    if(time){
      itime[3] <- proc.time()[1]
    }

    psi0 <- psi[1]
    psiLT <- c(0,psi[-1])
    Lpsi <- psiLT[-1] - psiLT[-(T+1)]
    ys <- c(psi0, dat)
    Ly <- ys[-1] - ys[-(T+1)]
    a <- sum(Lpsi^2)/2/W
    b <- sum(Lpsi*Ly)/W
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a1, b12=b1)$maximum
    if(logcon(b, a1, b1)){
      rejsam[i] <- "ARS"
      V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a1, b12=b1)
    }
    else{
      rejsam[i] <- "TPROP"
      propvar <- - 1 /( (a1 + 1)*mn^(-2) - b*mn^(-3/2)/4 - 2*b1*mn^(-3) )
      d <- optimize(propM, c(1,10^10), maximum=FALSE,  a=a, b=b, a12=a1, b12=b1, mn=mn, propvar=propvar)
      M <- d$objective
      df <- d$minimum
      rej <- TRUE
      while(rej){
        prop <- rtprop(1, mn, propvar, df)
        if(prop>0){
          R <- logpirej(prop, a, b, a1, b1, mn, propvar, df) - M
          u <- runif(1,0,1)
          if(log(u)<R){
            V <- prop
            rej <- FALSE
          }
        }
      }
    }

    if(time){
      itime[4] <- proc.time()[1]
    }

    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    theta <- solve(A)%*%(psi - ytild)

    if(time){
      itime[5] <- proc.time()[1]
      times[i,] <- c(itime[-1] - itime[-4])
    }

    out[i,] <- c(theta,V,W)
  }

  if(time){
    timeout <- data.frame(states=times[,1] + times[,4], V=times[,3], W=times[,2], rejsam=rejsam)
    return(list(out, timeout))
  }
  else{
    return(out)
  }
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, time=FALSE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B

  if(time){
    times <- matrix(0, nrow=n, ncol=4)
    itime <- rep(0,5)
    rejsam <- rep(0,n)
  }

  for(i in 1:n){
    if(time){
      itime[1] <- proc.time()[1]
    }

    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]

    if(time){
      itime[2] <- proc.time()[1]
    }

    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)

    if(time){
      itime[3] <- proc.time()[1]
    }

    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a2, b12=b2)$maximum
    if(logcon(b, a2, b2)){
      rejsam[i] <- "ARS"
      W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a2, b12=b2)
    }
    else{
      rejsam[i] <- "TPROP"
      propvar <- - 1 /( (a2 + 1)*mn^(-2) - b*mn^(-3/2)/4 - 2*b2*mn^(-3) )
      d <- optimize(propM, c(1,10^10), maximum=FALSE,  a=a, b=b, a12=a2, b12=b2, mn=mn, propvar=propvar)
      M <- d$objective
      df <- d$minimum
      rej <- TRUE
      while(rej){
        prop <- rtprop(1, mn, propvar, df)
        if(prop>0){
          R <- logpirej(prop, a, b, a2, b2, mn, propvar, df) - M
          u <- runif(1,0,1)
          if(log(u)<R){
            W <- prop
            rej <- FALSE
          }
        }
      }
    }

    if(time){
      itime[4] <- proc.time()[1]
    }

    A <- B/sqrt(W)
    A[1,1] <- 1
    theta <- solve(A)%*%gam

    if(time){
      itime[5] <- proc.time()[1]
      times[i,] <- c(itime[-1] - itime[-4])
    }

    out[i,] <- c(theta,V,W)
  }

  if(time){
    timeout <- data.frame(states=times[,1] + times[,4], V=times[,2], W=times[,3], rejsam=rejsam)
    return(list(out, timeout))
  }
  else{
    return(out)
  }
}

## returns TRUE if target density of log-concave, FALSE otherwise
logcon <- function(b, a12, b12, eps=0.2){
  sb <- 1*(b>0) - 1*(b<0)  ##sign(b)
  RHS <- (a12 + 1)^3*(1 - 2/3*sb)*32/9/b12 + eps
  ## + eps to make sure ARS algorithm doesn't fail on
  ## near non-log-concave cases
  LHS <- b^2
  out <- (LHS > RHS)
  return(out)
}
## compute M in logtarget - logprop <= M
propM <- function(df, a, b, a12, b12, mn, propvar){
  M <- optimize(logpirej, c(0,10^10), maximum=TRUE, a, b, a12, b12, mn, propvar, df)
  return(M$objective)
}

## log density of the t location scale family
dtprop <- function(x, mn, var, df){
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

## log of the conditional posterior density of W (V) given V (W), gamma (psi), data
logpiVW <- function(VW, a, b, a12, b12){
  out <- -a*VW + b*sqrt(VW) - (a12 + 1)*log(VW) - b12/VW
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprime <- function(VW, a, b, a12, b12){
  out <- -a + b/2/sqrt(VW) - (a12 + 1)/VW + b12/(VW^2)
  return(out)
}

## difference between log conditional posterior of W (V) and the proposal density
logpirej <- function(VW, a, b, a12, b12, mn, propvar, df){
  out <- logpiVW(VW, a, b, a12, b12) - dtprop(VW, mn, propvar, df)
  return(out)
}




