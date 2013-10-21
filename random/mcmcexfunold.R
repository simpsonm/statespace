library(dlm)
library(coda)
library(MCMCpack)
library(ars)

sdsmooth <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, jac=TRUE, simp=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  for(i in 1:n){
    ##print(c(i, V, W))
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
    k <- (!jac)*T/2 + a2 + 1
    mn <- optimize(hSD, c(0,10000), maximum=TRUE, a, b, k, V, dat, gam0, cgam, simp)$maximum
    propvar <- - 1 /( (k)*mn^(-2) -b*mn^(-3/2)/2-2*b2*mn^(-3) )
    M <- optimize(RSD, c(0,100000), maximum=TRUE, a, b, k, V, dat, gam0, cgam, simp, mn, propvar)
    M <- M$objective
    rej <- TRUE
    while(rej){
      prop <- rcauchy(1, mn, sqrt(propvar))
      if(prop>0){
        R <- RSD(prop, a, b, k, V, dat, gam0, cgam, simp, mn, propvar) - M
        u <- runif(1,0,1)
        if(log(u)<R){
          W <- prop
          rej <- FALSE
        }
      }
    }
    A <- B/sqrt(W)
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    out[i,] <- c(theta,V,W)
  }
  return(out)
}

hSD <- function(W, a, b, k, V, dat, gam0, cgam, simp){
  if(simp){
    out <- -a*W + b*sqrt(W) - (k)*log(W) - b2/W
  }
  else{
    out <- -k*log(W) - b2/W - sum((dat - gam0 - sqrt(W)*cgam)^2)/2/V
  }
  return(out)
}

RSD <- function(W, a, b, k, V, dat, gam0, cgam, simp, mn, propvar){
  out <- hSD(W, a, b, k, V, dat, gam0, cgam, simp) - dcauchy(W, mn, sqrt(propvar), log=TRUE)
  return(out)
}


notranssmooth <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  T1 <- T+1
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  for(i in 1:n){
    print(c(i, V, W))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/W
    A[1,1] <- 1
    gam <- A%*%theta
    Va <- a1 + T/2
    Vb <- b1 + sum((dat - gam[1] - W*cumsum(gam[-1]))^2)/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cumsum(gam[-1])^2 )/2/V
    b <- sum( (dat - gam[-1]) * cumsum( gam[-1]) )/V
    mn <- optimize(hNT, c(0,10000), maximum=TRUE, V, a, b, T, a2, b2)$maximum
    propvar <- - 1 /( (T/2 + a2 + 1)/mn^2 - 2*a - 2*b2/mn^3)
    M <- optimize(RNT, c(0,100000), maximum=TRUE, V, a, b, T, a2, b2, mn, propvar)
    M <- M$objective
    rej <- TRUE
    while(rej){
      prop <- rcauchy(1, mn, sqrt(propvar))
      if(prop>0){
        R <- RNT(prop, V, a, b, T, a2, b2, mn, propvar) - M
        u <- runif(1,0,1)
        if(log(u)<R){
          W <- prop
          rej <- FALSE
        }
      }
    }
    print(W)
    A <- B/W
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    out[i,] <- c(theta,V,W)
  }
  return(out)
}

hNT <- function(W, V, a, b, T, a2, b2){
  out <- -a*W^2 + b*W - (T/2 + a2 + 1)*log(W) - b2/W
  return(out)
}

RNT <- function(W, V, a, b, T, a2, b2, mn, propvar){
  out <- hNT(W, V, a, b, T, a2, b2) - dcauchy(W, mn, sqrt(propvar), log=TRUE)
  return(out)
}

exhNT <- function(W, V, a, b, T, a2){
  out <- hNT(W, V, a, b, T, a2)
  return(exp(out))
}

hpNT <- function(W, V, a, b, T, a2){
  out <- -2*a*W + b - (T/2 + a2 + 1)/W
  return(out)
}

distFFBS <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, trans=TRUE){
  T <- length(dat)
  T1 <- T+1
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  A <- solve(B)
  for(i in 1:n){
    mod <- dlmModPoly(1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    gam <- distLLSample(filt)
    theta <- A%*%gam
    if(trans==TRUE){
      V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
      W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-T1])^2)/2)
    }
    else{
      Va <- a1 + T/2
      Vb <- b1 + 0.5 * sum( ( dat - gam[1] - W*cumsum(gam[-1]) )^2 )
      V <- rinvgamma(1, Va, Vb)
      WR <- wrej(gam, dat, V, a2, b2)
      W <- WR[1]
    }
    out[i,] <- c(theta,V,W)
  }
  return(out)

}

##provides one sample from the (unscaled) disturbances
distLLSample <- function(filt){
  dat <- filt$y
  T <- length(dat)
  v <- dat - filt$f
  F <- rep(0,T)
  P <- F
  K <- F
  mod <- filt$mod
  V <- mod$V
  W <- mod$W
  ##extract info from Kalman filter
  for(t in 1:T){
    U.R <- filt$U.R[[t]]
    D.R <- filt$D.R[t,]
    P[t] <-  U.R^2*D.R^2
    F[t] <- P[t] + V
    K[t] <- P[t]/F[t]
  }
  ##construct means and variances of disturbances
  r <- rep(0,T+1)
  e <- rep(0,T)
  D <- e
  N <- r
  VHu <- e
  for(t in T:1){
    e[t] <- v[t]/F[t] - K[t]*r[t+1]
    r[t] <- e[t] + r[t+1]
    D[t] <-  1/F[t] + N[t+1]*K[t]^2
    N[t] <- D[t]-2*N[t+1]*K[t] + N[t+1]
  }
  mns <- W*r
  VHu <- W*(1-W*N[1:T])
  ##Used for covariances of disturbances
  L <- 1 - K
  Ls <- matrix(0,T,T)
  for(t in 1:T){
    for(s in 1:t){
      if(t==s){
        Ls[t,s] <- 1
      }
      else{
        Ls[t,s] <- prod(L[(t-1):(s)])
      }
    }
  }
  ##construct covariance matrix of disturbances
  VSig <- diag(VHu)
  for(t in 2:T){
    for(s in 1:(t-1)){
      VSig[t,s] <- -W^2*N[t]*(1-K[t-1])*Ls[t-1,s]
      VSig[s,t] <- VSig[t,s]
    }
  }
  CK <- chol(VSig)
  ##sample from the T disturbances
  gam <- mns[1:T] + t(CK)%*%rnorm(T)
  mn0 <- T/V*mean(dat-cumsum(gam)) + mod$m0/mod$C0
  vr0 <- 1/(1/mod$C0 + T/V)
  mn0 <- vr0*mn0
  out <- c(rnorm(1,mn0,sqrt(vr0)),gam)
  return(out)
}

## FFBS: samples V and W conditional on theta in Gibbs
simsmooth <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  T1 <- T+1
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


## uses the rejection sampler with cauchy proposal
dissmoothR <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, Cau=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta", 0:T, sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  rejs <- 0
  for(i in 1:n){
    print(c(i,V,W, rejs))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/W
    A[1,1] <- 1
    gam <- A%*%theta
    Va <- a1 + T/2
    Vb <- b1 + 0.5 * sum( ( dat - gam[1] - W*cumsum(gam[-1]) )^2 )
    V <- rinvgamma(1, Va, Vb)
    if((Cau) & (rejs < 100)){
      WR <- wrejC(gam, dat, V, a2, b2)
    }
    else{
      WR <- wrej(gam, dat, V, a2, b2)
    }
    W <- WR[1]
    rejs <- WR[2]
    out[i,] <- c(theta, V, W)
  }
  return(out)
}



## uses the adaptive rejection sampler - but breaks down
dissmoothAR <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+3))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  x <- start[2]*c(.25, 1.75)
  for(i in 1:n){
    print(c(i,V,W,x))
    mod <- dlmModPoly(1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/W
    A[1,1] <- 1
    gam <- A%*%theta
    Va <- a1 + T/2
    Vb <- b1 + 0.5 * sum( ( dat - gam[1] - W*cumsum(gam[-1]) )^2 )
    V <- rinvgamma(1, Va, Vb)

    ##generate W
    a <- sum( cumsum(gam[-1])^2 )/2/V
    b <- sum( (dat - gam[1]) * cumsum(gam[-1]) )/V
    hp <- hprimefun(x, a, b, a2, b2, T)
    while(hp[1]<0){
      x[1] <- x[1]/2
      hp <- hprimefun(x, a, b, a2, b2, T)
    }
    while(hp[2]>0){
      x[2] <- x[2]*2
      hp <- hprimefun(x, a, b, a2, b2, T)
    }
    h <- hfun(x, a, b, a2, b2, T)
    z <- NULL
    k <- 2
    for(j in 1:(k-1)){
      zj <- ( h[j+1]-h[j] - x[j+1]*hp[j+1] + x[j]*hp[j] )/( hp[j] - hp[j+1] )
      z <- c(z,zj)
    }
    rej <- TRUE
    while(rej){
      WXstar <- ssam(z, x, h, hp)
      W <- WXstar[1]
      u <- runif(1,0,1)
      ukstar <- ukfun(W, z, x, h, hp)
      lkstar <- lkfun(W, x, h, hp)
      if( u <= exp( lkstar - ukstar) ){
        rej <- FALSE
      }
      else{
        hstar <- hfun(W, a, b, a2, b2, T)
        if(u <= exp(hstar - ukstar)){
          rej <- FALSE
        }
        else{
          hpstar <- hprimefun(W, a, b, a2, b2, T)
          x <- c(x,W)
          h <- c(h,hstar)
          hp <- c(hp,hpstar)
          idx <- order(x)
          x <- x[idx]
          h <- h[idx]
          hp <- hp[idx]
          k <- k+1
          z <- NULL
          for(j in 1:(k-1)){
            zj <- ( h[j+1]-h[j] - x[j+1]*hp[j+1] + x[j]*hp[j] )/( hp[j] - hp[j+1] )
            z <- c(z,zj)
          }

        }
      }
    }
    x <- WXstar[c(2,3)]
    A <- B/W
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    out[i,] <- c(theta,V, W)
  }
  return(out)
}

## uses metropolis with cauchy proposal for W
dissmoothM <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=2))
  colnames(out) <- c("V","W")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  for(i in 1:n){
    ##print(c(i,V,W))
    mod <- dlmModPoly(1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/W
    A[1,1] <- 1
    gam <- A%*%theta
    Va <- a1 + T/2
    Vb <- b1 + 0.5 * sum((dat - gam[1] - W*cumsum(gam[-1]))^2)
    V <- rinvgamma(1, Va, Vb)
    W <- wmet(gam, dat, V, W, a2, b2)
    out[i,] <- c(V, W)
  }
  return(out)
}

hfun <- function(W, a, b, a2, b2, T){
  I <- length(W)
  out <- rep(0,I)
  for(i in 1:I){
    out[i] <- -a*W[i]^2 + b*W[i] + (T/2 - a2 - 2)*log(W[i]) - b2/W[i]
  }
  return(out)
}

hprimefun <- function(W, a, b, a2, b2, T){
  I <- length(W)
  out <- rep(0,I)
  for(i in 1:I){
    out[i] <- -2*a*W[i] + b + (T/2 - a2 - 2)/W[i] + b2/W[i]^2
  }
  return(out)
}

ukfun <- function(W, z, x, h, hp){
  I <- length(W)
  out <- rep(0,I)
  for(i in 1:I){
    j <- sum(W[i] > z) + 1
    out[i] <- h[j] + (W[i]-x[j])*hp[j]
  }
  return(out)
}

lkfun <- function(W, x, h, hp){
  I <- length(W)
  out <- rep(0,I)
  for(i in 1:I){
    j <- sum(W[i]>x)
    if(j == 0 | j == length(x)){
      out[i] <- -Inf
    }
    else{
      out[i] <- ( (x[j+1]-W[i])*h[j] + (W[i]-x[j])*h[j+1] )/( x[j+1] - x[j] )
    }
  }
  return(out)
}

ssam <- function(z, x, h, hp){
  zz <- c(0,z,Inf)
  con <- rep(0,length(x))
  lcon <- con
  for(j in 1:length(x)){
    if(hp[j]>0){
      lcon[j] <- h[j] - x[j]*hp[j] - log(hp[j]) + log( exp(zz[j+1]*hp[j]) - exp(zz[j]*hp[j]) )
    }
    else{
      lcon[j] <- h[j] - x[j]*hp[j] - log(-hp[j]) + log( exp(zz[j]*hp[j]) - exp(zz[j+1]*hp[j]) )
    }
    con[j] <- exp(lcon[j])
  }
  p <- runif(1,0,1)
  sam <- sinv(p, lcon, zz, x, h, hp)
  xnew <- c(sinv(0.25, lcon, zz, x, h, hp), sinv(0.75, lcon, zz, x, h, hp))
  return(c(sam, xnew))
}

sinv <- function(p, lcon, zz, x, h, hp){
  lc1 <- lcon[1]
  lcumcon <- c(-Inf, lc1 + log(cumsum(exp(lcon-lc1))) )
  ltcon <- lcumcon[length(lcumcon)]
  i <- sum(ltcon + log(p) >= lcumcon)
  out <- log( (exp(ltcon + log(p)) - exp(lcumcon[i]))*hp[i]*exp(x[i]*hp[i]-h[i]) + exp(zz[i]*hp[i]))/hp[i]
  return(out)
}


##log(( exp(ltcon)*p - exp(lcumcon[i]) )*hp[i]*exp(x[i]*hp[i] - h[i])+exp(zz[i]*hp[i]))/hp[i]

##sinv(.01, lcon, zz, x, h, hp)


llsim <- function(T, V, W, m0, C0){
  out <- rep(0,T)
  theta <- rnorm(1, m0, sqrt(C0))
  for(t in 1:T){
    theta <- theta + rnorm(1, 0, sd=sqrt(W))
    out[t] <- theta + rnorm(1, 0, sd=sqrt(V))
  }
  return(out)
}

rejfun <- function(W, a, b, a2, b2, T, rm, rv, df){
  out <- hfun(W, a, b, a2, b2, T) - dt((W-rm)/sqrt(rv), df, log=TRUE)
  return(out)
}

rejfunC <- function(W, a, b, a2, b2, T, rm, rsd){
  out <- hfun(W, a, b, a2, b2, T) - dcauchy(W, rm, rsd, log=TRUE)
  return(out)
}


wrej <- function(gam, dat, V, a2, b2){
  T <- length(dat)
  a <- sum(cumsum(gam[-1])^2)/2/V
  b <- sum( (dat - gam[1]) * cumsum(gam[-1]))/V
  rm <- optimize(hfun, c(0,100000), a, b, a2, b2, T, maximum=TRUE)$maximum
  rv <- 1/(2*a + (T/2-a2-1)/rm^2 + 2*b2/rm^3)
  df <- c(1:49, seq(50, 95, by=5))
  ##df <- 1
  M <- rep(0,length(df))
  for(i in 1:length(df)){
    Mtemp <- optimize(rejfun, c(0,100000), a, b, a2, b2, T, rm, rv, df[i], maximum=TRUE)
    M[i] <- Mtemp$objective
  }
  idx <- which.min(M)
  dfst <- df[idx]
  Mst <- M[idx]
  rej <- TRUE
  cnt <- 0
  while(rej){
    cnt <- cnt + 1
    prop <- rt(1, dfst)*sqrt(rv) + rm
    if(prop>0){
      u <- runif(1,0,1)
      R <- hfun(prop, a, b, a2, b2, T) - dt((prop-rm)/sqrt(rv), dfst, log=TRUE) - Mst
      if(u <= exp(R)){
        rej <- FALSE
      }
    }
  }
  return(c(prop, cnt))
}

wrejC <- function(gam, dat, V, a2, b2){
  T <- length(dat)
  a <- sum( cumsum(gam[-1])^2 )/2/V
  b <- sum( (dat - gam[1]) * cumsum(gam[-1]) )/V
  rm <- optimize(hfun, c(0,10000000), a, b, a2, b2, T, maximum=TRUE)$maximum
  rv <- 1/(2*a + (T/2-a2-1)/rm^2 + 2*b2/rm^3)
  rsd <- sqrt(rv)
  M <- optimize(rejfunC, c(0,10000000), a, b, a2, b2, T, rm, rsd, maximum=TRUE)
  M <- M$objective
  rej <- TRUE
  cnt <- 0
  while(rej){
    cnt <- cnt + 1
    prop <- rcauchy(1,rm,rsd)
    if(prop>0){
      u <- runif(1,0,1)
      R <- hfun(prop, a, b, a2, b2, T) - dcauchy(prop, rm, rsd, log=TRUE) - M
      hfun(prop, a, b, a2, b2, T) - dcauchy(prop, rm, rsd, log=TRUE)
      if(u <= exp(R)){
        rej <- FALSE
      }
    }
  }
  return(c(prop, cnt))
}


wmet <- function(gam, dat, V, W, a2, b2){
  T <- length(dat)
  a <- sum(cumsum(gam[-1])^2)/2/V
  b <- sum( (dat - gam[1]) * cumsum(gam[-1]))/V
  rm <- optimize(hfun, c(0,100000), a, b, a2, b2, T, maximum=TRUE)$maximum
  rv <- 1/(2*a + (T/2-a2-1)/rm^2 + 2*b2/rm^3)
  rsd <- sqrt(rv)
  ##df <- c(1:49, seq(50, 95, by=5), seq(100, 150, by=10), Inf)
  df <- 1
  M <- rep(0,length(df))
  for(i in 1:length(df)){
    Mtemp <- optimize(rejfun, c(0,100000), a, b, a2, b2, T, rm, rv, df[i], maximum=TRUE)
    M[i] <- Mtemp$objective
  }
  idx <- which.min(M)
  dfst <- df[idx]
  Mst <- M[idx]
  prop <- 0
  while(prop <= 0){
    prop <- rt(1, dfst)*rsd + rm
  }
  R1 <- hfun(prop, a, b, a2, b2, T) - dt((prop-rm)/rsd, dfst, log=TRUE)
  R2 <- hfun(W, a, b, a2, b2, T) - dt((W-rm)/rsd, dfst, log=TRUE)
  R <- R1-R2
  u <- runif(1,0,1)
  if(u<exp(R)){
    out <- prop
  }
  else{
    out <- W
  }
  return(out)
}





