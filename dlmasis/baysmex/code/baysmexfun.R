library(ars)
source("mcfa.R")


naivegibbs <- function(n,start, dat, au, bu, av, bv, aw, bw, V0, W0, m0, U0, beta0, B){
  J <- ncol(dat)
  N <- nrow(dat)
  U <- start[1]
  V <- start[2:(J+1)]
  W <- start[(J+2):(2*J+1)]
  beta <- start[(2*J+2):(3*J+1)]
  mu <- start[(3*J+2):(3*J+2+N)]
  theta <- matrix(start[(3*J+2+N+1):(3*J + (J+1)*(N+1)+1)], byrow=FALSE, ncol=J)
  out <- matrix(0, nrow=n, ncol=3*J + (J+1)*(N+1) + 1)
  thetanames <- NULL
  for(j in 1:J){
    thetanames <- c(thetanames, paste("theta", j, 0:N, sep="."))
  }
  colnames(out) <- c("U", paste("V", 1:J, sep="."), paste("W", 1:J, sep="."),
                     paste("beta", 1:J, sep="."), paste("mu", 0:N, sep="."),
                     thetanames)
  aut <- au + N/2
  avt <- av + N/2
  awt <- aw + N/2
  but <- bu
  bvt <- bv
  bwt <- bw
  for(i in 1:n){
    mu <- mcfamu(W, W0, U, U0, m0, theta, beta)
    but <- bu + sum(diff(mu)^2)/2
    U <- 1/rgamma(1, shape=aut, rate=but)
    beta <- rnorm(J)
    for(j in 1:J){
      thetaj <- theta[,j]
      Obeta <- sum(thetaj[-(N+1)]^2)/W[j] + 1/B
      obeta <- sum((thetaj[-1]-mu[-1])*thetaj[-(N+1)])/W[j] + beta0/B
      Obetainv <- 1/Obeta
      beta[j] <- Obetainv*obeta + sqrt(Obetainv)*beta[j]
      theta[,j] <- mcfatheta(dat[,j], V[j], W[j], mu, W0[j], beta[j])
      thetaj <- theta[,j]
      bvt[j] <- bv[j] + sum((dat[,j] - thetaj[-1])^2)/2
      bwt[j] <- bw[j] + sum((thetaj[-1] - mu[-1] - thetaj[-(N+1)]*beta[j])^2)/2
      V[j] <- 1/rgamma(1,shape=avt[j], rate=bvt[j])
      W[j] <- 1/rgamma(1,shape=awt[j], rate=bwt[j])
    }
    out[i,] <- c(U,V,W,beta,mu,theta)
  }
  return(out)
}


disterrorinter <- function(n,start, dat, au, bu, av, bv, aw, bw, V0, W0, m0, U0, beta0, B){
  J <- ncol(dat)
  N <- nrow(dat)
  U <- start[1]
  V <- start[2:(J+1)]
  W <- start[(J+2):(2*J+1)]
  beta <- start[(2*J+2):(3*J+1)]
  mu <- start[(3*J+2):(3*J+2+N)]
  theta <- matrix(start[(3*J+2+N+1):(3*J + (J+1)*(N+1)+1)], byrow=FALSE, ncol=J)
  out <- matrix(0, nrow=n, ncol=3*J + (J+1)*(N+1) + 1)
  thetanames <- NULL
  for(j in 1:J){
    thetanames <- c(thetanames, paste("theta", j, 0:N, sep="."))
  }
  colnames(out) <- c("U", paste("V", 1:J, sep="."), paste("W", 1:J, sep="."),
                     paste("beta", 1:J, sep="."), paste("mu", 0:N, sep="."),
                     thetanames)
  aut <- au + N/2
  avt <- av + N/2
  awt <- aw + N/2
  for(i in 1:n){
    mu <- mcfamu(W, W0, U, U0, m0, theta, beta)
    but <- bu + sum(diff(mu)^2)/2
    U <- 1/rgamma(1, shape=aut, rate=but)
    beta <- rnorm(J)
    for(j in 1:J){
      Obeta <- sum(theta[-(N+1),j]^2)/W[j] + 1/B
      obeta <- sum((theta[-1,j]-mu[-1])*theta[-(N+1),j])/W[j] + beta0/B
      Obetainv <- 1/Obeta
      beta[j] <- Obetainv*obeta + sqrt(Obetainv)*beta[j]
      ##theta[,j] <- mcfatheta(dat[,j], V[j], W[j], mu, W0[j], beta[j])
      ##bvt[j] <- bv[j] + sum((dat[,j] - theta[-1,j])^2)/2
      ##bwt[j] <- bw[j] + sum((theta[-1,j] - mu[-1] - theta[-(N+1),j]*beta[j])^2)/2
      ##V[j] <- 1/rgamma(1,shape=avt[j], rate=bvt[j])
      ##W[j] <- 1/rgamma(1,shape=awt[j], rate=bwt[j])
    }
    out[i,] <- c(U,V,W,beta,mu,theta)
  }
  return(out)

}

## returns TRUE if target density of log-concave, FALSE otherwise
logcon <- function(b, avw, bvw, eps=.1){
  RHS <- sqrt((avw + 1)^3/bvw)*4/3*sqrt(2/3) + eps
  ## + eps to make sure ARS algorithm doesn't fail on
  ## near non-log-concave cases
  out <- (b > RHS)
  return(out)
}

logconls <- function(a, b, avw, bvw, eps=.1){
  RHS <- 16*(bvw/27*a^3)^(1/4) - eps
  out <- (b < RHS)
  return(out)
}

## compute M in logtarget - logprop <= M
propM <- function(df, a, b, avw, bvw, mn, propvar){
  M <- optimize(logpirej, c(-10^2,10^2), maximum=TRUE, a, b, avw, bvw, mn, propvar, df)
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

## log of the conditional posterior density of logW (logV) given V (W), gamma (psi), data
logpilVW <- function(x, a, b, avw, bvw){
  out <- -a*exp(x) + b*exp(x/2) - avw*x - bvw*exp(-x)
  if(exp(x) == Inf)
      out <- -Inf
  return(out)
}

## first derivative of log of the conditional posterior of logW (logV)
logpilVWprime <- function(x, a, b, avw, bvw){
  out <- -a*exp(x) + b*exp(x/2)/2 - avw + bvw*exp(-x)
  return(out)
}

## log of the conditional posterior density of W (V) given V (W), gamma (psi), data
logpiVW <- function(VW, a, b, avw, bvw){
  out <- -a*VW + b*sqrt(VW) - (avw + 1)*log(VW) - bvw/VW
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprime <- function(VW, a, b, avw, bvw){
  out <- -a + b/2/sqrt(VW) - (avw + 1)/VW + bvw/(VW^2)
  return(out)
}

## difference between log conditional posterior of logW (logV) and the proposal density
logpirej <- function(VW, a, b, avw, bvw, mn, propvar, df){
  out <- logpilVW(VW, a, b, avw, bvw) - dtprop(VW, mn, propvar, df)
  return(out)
}

## transforms theta to gamma
gamtrans <- function(theta, W){
  T <- length(theta) - 1
  gam1 <- (theta[-1]-theta[-(T+1)])/sqrt(W)
  gam <- c(theta[1], gam1)
  return(gam)
}

## transforms theta to psi
psitrans <- function(dat, theta, V){
  psi1 <- (dat - theta[-1])/sqrt(V)
  psi <- c(theta[1], psi1)
  return(psi)
}

## transforms gamma to psi
psigamtrans <- function(dat, gam, V, W){
  diff <- dat-gam[1]
  cgam <- cumsum(gam[-1])
  psi1 <- (diff - sqrt(W)*cgam)/sqrt(V)
  psi <- c(gam[1], psi1)
  return(psi)
}

## transforms psi to gamma
gampsitrans <- function(dat, psi, V, W){
  T <- length(dat)
  gam2 <- (dat[-1]-dat[-T] - sqrt(V)*(psi[-c(1:2)] - psi[-c(1,T+1)]))/sqrt(W)
  gam <- c(psi[1], (dat[1] - sqrt(V)*psi[2] - psi[1])/sqrt(W), gam2)
  return(gam)
}


## transforms gamma to theta
thetagamtrans <- function(gam, W){
  gam0 <- gam[1]
  cgam <- cumsum(gam[-1])
  theta <- c(gam0, gam0 + sqrt(W)*cgam)
  return(theta)
}

## transforms psi to theta
thetapsitrans <- function(dat, psi, V){
  theta1 <- dat - sqrt(V)*psi[-1]
  theta <- c(psi[1], theta1)
  return(theta)
}

## samples V,W conditional on theta
VWthetaiter <- function(dat, theta, av, aw, bv, bw){
  T <- length(dat)
  V <- rinvgamma(1, av + T/2, bv + sum((dat-theta[-1])^2)/2)
  W <- rinvgamma(1, aw + T/2, bw + sum((theta[-1]-theta[-(T+1)])^2)/2)
  return(c(V,W))
}

VWrejiter <- function(a, b, avw, bvw){
  mn <- optimize(logpilVW, c(-10^2,10^2), maximum=TRUE, a=a, b=b, avw=avw, bvw=bvw)$maximum
  propvar <- - 1 /( -a*exp(mn) + b*exp(mn/2)/4 - bvw*exp(-mn) )
  df <- 1
  rej <- TRUE
  rejit <- 1
  M <- optimize(logpirej, c(-10^2, 10^2), a=a, b=b, avw=avw, bvw=bvw, mn=mn, propvar=propvar, df=df, maximum=TRUE)$objective
  while(rej){
    prop <- rtprop(1, mn, propvar, df)
    R <- logpirej(prop, a, b, avw, bvw, mn, propvar, df) - M
    u <- runif(1,0,1)
    if(log(u)<R){
      W <- exp(prop)
      rej <- FALSE
    }
    rejit <- rejit + 1
    if(rejit == 100000)
        rej <- FALSE
  }
  return(W)
}

## samples W conditional on V,gamma
Wgamiter <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2/V
  b <- sum( (dat - gam0) * cgam )/V
  lcon <- logcon(b, aw, bw)
  adrej <- lcon
  if(lcon){
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=aw, bvw=bw)$maximum
    try(W <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(mn/2, mn, mn*2),
                      lb=TRUE, xlb=0, a=a, b=b, avw=aw, bvw=bw))
    if(W==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    W <- VWrejiter(a, b, aw, bw)
  }
  return(c(W, lcon, adrej))
  }

## samples V conditional on W,gamma
Vgamiter <- function(dat, gam, W, av, bv){
  T <- length(dat)
  theta <- thetagamtrans(gam, W)
  Va <- av + T/2
  Vb <- bv + sum( (dat - theta[-1])^2 )/2
  V <- rinvgamma(1, Va, Vb)
  return(V)
}


## samples V conditional on W,psi 
Vpsiiter <- function(dat, psi, W, av, bv){
  T <- length(dat)
  psi0 <- psi[1]
  psiLT <- c(0,psi[-1])
  Lpsi <- psiLT[-1] - psiLT[-(T+1)]
  ys <- c(psi0, dat)
  Ly <- ys[-1] - ys[-(T+1)]
  a <- sum(Lpsi^2)/2/W
  b <- sum(Lpsi*Ly)/W
  lcon <- logcon(b, av, bv)
  adrej <- lcon
  if(lcon){
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=av, bvw=bv)$maximum
    try(V <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(mn/2, mn, mn*2),
                 lb=TRUE, xlb=0, a=a, b=b, avw=av, bvw=bv))
    if(V==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    V <- VWrejiter(a, b, av, bv)
  }
  return(c(V, lcon, adrej))
}


## samples W conditional on V,psi
Wpsiiter <- function(dat, psi, V,  aw, bw){
  T <- length(dat)
  theta <- thetapsitrans(dat, psi, V)
  theta <- c(theta)
  Wa <- aw + T/2
  Wb <- bw + sum( (theta[-1] - theta[-(T+1)])^2 )/2
  W <- rinvgamma(1, Wa, Wb)
  return(W)
}

randkerniter <- function(dat, V, W, theta, av, aw, bv, bw, m0=0, C0=10^7, probs=c(1/3, 1/3, 1/3)){
  kernels <- c("state", "dist", "error")
  kernel <- sample(kernels, 1, prob=probs)
  if(kernel=="state"){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    rejsW <- c(NA,NA)
    rejsV <- c(NA,NA)
  }
  if(kernel=="error"){
    ptma <- proc.time()
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    rejsW <- c(NA,NA)
    theta <- thetapsitrans(dat, psi, V)
  }
  if(kernel=="dist"){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    rejsV <- c(NA,NA)
    theta <- thetagamtrans(gam, W)
  }
  out <- list(theta=theta, V=V, W=W, rejsV=rejsV, rejsW=rejsW, kernel=kernel, smoothtime=smoothtime)
  return(out)
}

randkernsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7, probs=c(1/3, 1/3, 1/3)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)  
  for(i in 1:n){
    iter <- randkerniter(dat, V, W, theta, av, aw, bv, bw, m0, C0, probs)
    theta <- iter$theta
    V <- iter$V
    W <- iter$W
    kernel <- iter$kernel
    rejsV <- iter$rejsV
    rejsW <- iter$rejsW
    smoothtime <- iter$smoothtime
    out[i,-5] <- c(rejsV, rejsW, smoothtime, V, W, theta)
    out[i,5] <- kernel
  }
  return(out)
}

partialcissam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)  
  for(i in 1:n){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    gam <- gamtrans(theta, W)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(rep(NA,2), Wout[2:3], NA, smoothtime, V, W, theta)
  }
  return(out)
}


fullcissam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  Va <- av + T/2
  Wa <- aw + T/2
  out <- samoutsetup(n, T)  
  for(i in 1:n){
    ptma <- proc.time()
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    theta <- thetapsitrans(dat, psi, V)
    Vb <- bv + sum((dat-theta[-1])^2)/2
    V <- rinvgamma(1, Va, Vb)
    Wb <- bw + sum((theta[-1]-theta[-(T+1)])^2)/2
    W <- rinvgamma(1, Wa, Wb) 
    gam <- gamtrans(theta, W)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(Vout[2:3], Wout[2:3], NA, smoothtime, V, W, theta)
  }
  return(out)
}
