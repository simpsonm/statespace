## log density of the t location scale family
dtproptilde <- function(x, mn, var, df){
  sd <- sqrt(var)
  z <- (x-mn)/sd
  out <- dt(z, df, log=TRUE) - log(sd)
  return(out)
}

## simulate from the t location scale family
rtproptilde <- function(n, mn, var, df){
  sd <- sqrt(var)
  temp <- rt(n, df=df)
  out <- mn + temp*sd
  return(out)
}

## log of the conditional posterior density of W (V) given V (W), gamma (psi), tildeta
logpiVWtilde <- function(VW, a, b, cc, avw){
  out <-  - avw*VW - a*exp(-VW) + b*exp(-VW/2) - cc*exp(VW)
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprimetilde <- function(VW, a, b, cc, avw){
  out <-  - avw + a*exp(-VW) - (b/2)*exp(-VW/2) - cc*exp(VW)
  return(out)
}

logcontilde <- function(a, b, cc, eps=.01){
  out <- (b <= 0)
  if(!out)
      out <- (a > (b^2/(cc*16))*(1/4 - (b/cc)^2/16^3) + eps)
  return(out)
}

## difference between log conditional posterior of W (V) and the proposal density
logpirejtilde <- function(VW, a, b, cc, avw, mn, propvar, df){
  out <- logpiVWtilde(VW, a, b, cc, avw) - dtproptilde(VW, mn, propvar, df)
  return(out)
}

Wgamitertilde <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  Wa <- aw + T/2
  Wb <- bw + V*sum(gam[-1]^2)/2
  W <- rinvgamma(1, Wa, Wb)
  return(W)
}

Vpsiitertilde <- function(dat, psi, W, av, bv){
  T <- length(dat)
  Va <- av + T/2
  Vb <- bv + W*sum(psi[-1]^2)/2
  V <- rinvgamma(1, Va, Vb)
  return(V)
}

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


VWrejtilde <- function(a, b, cc, avw){
  mn <- uniroot(logpiVWprimetilde, c(-10^2, 10^2), a=a, b=b, cc=cc, avw=avw)$root
  lcon <- logcontilde(a, b, cc)
  lcon <- FALSE
  adrej <- lcon
  if(lcon){
    try(VW <- ars(n=1, logpiVWtilde, logpiVWprimetilde, x=c(mn-5*mn*2, mn, 5*mn*2),
                  a=a, b=b, cc=cc, avw=avw))
    if(VW==0){
      adrej <- FALSE
    }
  }
  if(!adrej){ 
    propvar <- - 1 /( -a*exp(-mn) +(b/4)*exp(-mn/2) - cc*exp(mn) )
    df <- 10
    M <- optimize(logpirejtilde, c(-10^2,10^2), maximum=TRUE, a, b, cc, avw, mn, propvar, df)$objective
    rej <- TRUE
    rejit <- 1
    while(rej){
      VW <- rtproptilde(1, mn, propvar, df)
      R <- logpirejtilde(VW, a, b, cc, avw, mn, propvar, df) - M
      u <- runif(1,0,1)
      if(log(u)<R){
        rej <- FALSE
      }
      rejit <- rejit + 1
      if(rejit%%1000000==0){
        print(paste("Rejection sampler appears to be stuck for VW, ",
                    rejit, " iterations", sep=""))
        print(paste("Last value of W in the chain: ", VW, sep=""))
        print(paste("Values of a b and c:", a, b, cc, sep="  "))
        print(paste("T = ", T, sep=""))
        print(paste("Prior value of alphaVW:", avw, sep="  "))
        print(paste("Attempted df: ", df, sep=""))
        print(paste("Attempted mean and variance:", mn, propvar, sep="  "))
      }
    }
  }
  return(c(exp(VW), lcon, adrej))
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
werrorsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- data.frame(matrix(0, nrow=n, ncol=T+9))
  colnames(out) <- c("logconV", "adrejV", "logconW",
                     "adrejW", "kernel", "stime",
                     "V", "W", paste("theta",0:T,sep=""))
  for(i in 1:n){
    ptma <- proc.time()
    psi <- awolpssmooth(dat, V, W, m0, C0)
    psi <- c(psi[1], psi[-1]*sqrt(V/W))
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vpsiitertilde(dat, psi, W, av, bv)
    Wout <- Wpsiitertilde(dat, psi, V,  aw, bw)
    W <- Wout[1]
    thetat <- dat - sqrt(W)*psi[-1]
    theta <- c(psi[1], thetat)
    out[i,] <- c(NA,NA,Wout[2:3], NA, smoothtime, V, W, theta)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
wdistsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- data.frame(matrix(0, nrow=n, ncol=T+9))
  colnames(out) <- c("logconV", "adrejV", "logconW",
                     "adrejW", "kernel", "stime",
                     "V", "W", paste("theta",0:T,sep=""))
  for(i in 1:n){
    ptma <- proc.time()
    theta <- awolthsmooth(dat, V, W, m0, C0)
    gamt <- (theta[-1] - theta[-(T+1)])/sqrt(V)
    gam <- c(theta[1], gamt)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    Vout <- Vgamitertilde(dat, gam, W, av, bv)
    V <- Vout[1]
    W <- Wgamitertilde(dat, gam, V, aw, bw)
    thetat <- gam[1] + cumsum(gam[-1])*sqrt(V)
    theta <- c(gam[1], thetat)
    out[i,] <- c(Vout[2:3],NA,NA,NA,smoothtime,V,W,theta)
  }
  return(out)
}

