## compute M in logtarget - logprop <= M
propMda <- function(df, a, b, c, avw, mn, propvar){
  M <- optimize(logpirejda, c(0,10^10), maximum=TRUE, a, b, c, avw, mn, propvar, df)
  return(M$objective)
}

## log density of the t location scale family
dtpropda <- function(x, mn, var, df){
  sd <- sqrt(var)
  z <- (x-mn)/sd
  out <- dt(z, df, log=TRUE) - log(sd)
  return(out)
}

## simulate from the t location scale family
rtpropda <- function(n, mn, var, df){
  sd <- sqrt(var)
  temp <- rt(n, df=df)
  out <- mn + temp*sd
  return(out)
}

## log of the conditional posterior density of W (V) given V (W), gamma (psi), data
logpiVWda <- function(VW, a, b, c, avw){
  out <-  - (avw + 1)*log(VW) - a/VW + b/sqrt(VW) - c*VW
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprimeda <- function(VW, a, b, c, avw){
  out <-  - (avw + 1)/VW +a/(VW^2) -b/(2*VW^(3/2)) - c
  return(out)
}

## difference between log conditional posterior of W (V) and the proposal density
logpirejda <- function(VW, a, b, c, avw, mn, propvar, df){
  out <- logpiVWda(VW, a, b, c, avw) - dtpropda(VW, mn, propvar, df)
  return(out)
}

Wgamiterda <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  Wa <- aw + T/2
  Wb <- bw + V*sum(gam[-1]^2)/2
  W <- rinvgamma(1, Wa, Wb)
  return(W)
}

Vpsiiterda <- function(dat, psi, W, av, bv){
  T <- length(dat)
  Va <- av + T/2
  Vb <- bv + W*sum(psi[-1]^2)/2
  V <- rinvgamma(1, Va, Vb)
  return(V)
}

Vgamiterda <- function(dat, gam, W, av, bv){
  T <- length(dat)
  gam0 <- gam[1]
  cgams <- cumsum(gam[-1])
  a <- bv + sum((dat - gam0)^2)/2
  b <- sum((dat-gam0)*cgams)
  cc <- cgams[T]^2/(2*W)
  V <- VWrejda(a,b,cc,av)
  return(V)
}

Wpsiiterda <- function(dat, psi, V, aw, bw){
  T <- length(dat)
  ys <- c(dat,psi[1])
  psis <- c(psi[-1],0)
  Ly <- ys[-1]-ys[-(T+1)]
  Lpsi <- psis[-1]-psis[-(T+1)]
  a <- bw + sum(Ly^2)/2
  b <- sum(Ly*Lpsi)
  cc <- sum(Lpsi^2)/(2*V)
  W <- VWrejda(a,b,cc,aw)
  return(W)
}



VWrejda <- function(a, b, cc, avw){
  mn <- optimize(logpiVWda, c(0,10^10), maximum=TRUE, a=a, b=b, c=cc, avw=avw)$maximum
  propvar <- - 1 /( (avw + 1)*mn^(-2) - 2*a*mn^(-3) + (3/4)*b*mn^(-5/2) )
  d <- optimize(propMda, c(1,10^10), maximum=FALSE,  a=a, b=b, c=cc, avw=avw, 
                mn=mn, propvar=propvar)
  M <- d$objective
  df <- d$minimum
  rej <- TRUE
  rejit <- 1
  while(rej){
    VW <- rtpropda(1, mn, propvar, df)
    if(VW>0){
      R <- logpirejda(VW, a, b, cc, avw, mn, propvar, df) - M
      u <- runif(1,0,1)
      if(log(u)<R){
        rej <- FALSE
      }
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

  return(VW)
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsamda <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W","logconV",
                     "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    psit <- psitrans(dat, theta, V)
    psi <- psit*sqrt(V)/sqrt(W)
    psi[1] <- psit[1]
    V <- Vpsiiterda(dat, psi, W, av, bv)
    W <- Wpsiiterda(dat, psi, V,  aw, bw)
    psit <- psi*sqrt(W)/sqrt(V)
    psit[1] <- psi[1]
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(theta,V,W,c(NA,NA),c(NA,NA), NA)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsamda <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W",
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    gamt <- gamtrans(theta, W)
    gam <- gamt*sqrt(W)/sqrt(V)
    gam[1] <- gamt[1]
    V <- Vgamiterda(dat, gam, W, av, bv)
    W <- Wgamiterda(dat, gam, V, aw, bw)
    gamt <- gam*sqrt(V)/sqrt(W)
    gamt[1] <- gam[1]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(theta,V,W,c(NA,NA),c(NA,NA), NA)
  }
  return(out)
}

