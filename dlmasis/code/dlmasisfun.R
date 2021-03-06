## A set of functions for simulating from and fitting local level models
library(ars)
library(plyr)
source("../code/mcfa.R")
source("../code/wscalerej.R")

fullsim <- function(samplers, simdata, n, burn, parallel){
  out <- ddply(samplers, .(sampler), samsim, simdata=simdata, n=n,
               burn=burn, parallel=parallel)
  return(out)
}


## simulates from a given sampler for each dataset and for multiple
## chains, and returns summary info on the first chain.
samsim <- function(samplers, simdata, n, burn, parallel){
  sampler <- samplers$sampler[1]
  print(sampler)
  sam <- ddply(simdata, .(V.T, W.T, T.T), samwrap, .parallel=parallel,
               n=n, samp=sampler)
  samnam <- paste(sampler, "SAM.RData", sep="")
  colnam <- grep("(V.T|W.T|T.T|V|W|time|stime|logconW|adrejW|logconWls|adrejWls|logconV|adrejV|logconVls|adrejVls|kernel)$", colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samshort)
  out <- ddply(sam, .(V.T, W.T, T.T), samsummary,
               .parallel=parallel, dat=simdata, burn=burn,
               sampler=sampler)
  if(sampler=="trialt"){
    postcors <- ddply(sam, .(V.T, W.T, T.T), postcor, .parallel=parallel,
                           dat=simdata, burn=burn)
    save(postcors, file="postcors.RData")
    rm(postcors)
  }
  rm(sam)
  save(out, file=paste(sampler, "OUT.RData", sep=""))
  print(paste(sampler, "finished", sep=" "))
  return(out)
}

postcor <- function(sam, dat, burn){
  V.T <- sam$V.T[1]
  W.T <- sam$W.T[1]
  T.T <- sam$T.T[1]
  V <- sam$V[-c(1:burn)]
  W <- sam$W[-c(1:burn)]
  theta0s <- sam[-c(1:burn),grep("theta", colnames(sam))]
  theta0s <- theta0s[,1:(T.T+1)]
  bv <- apply(t((yt - t(theta0s[,-1]))^2), 1, sum)
  bw <- apply( t((theta0s[,-1] - theta0s[,-(T.T+1)])^2 ), 1, sum)
  thetas <- theta0s[,-1]
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  cgam <- t(apply(gammas, 1, cumsum))
  agam <- apply(cgam^2, 1, sum)/(2*V)
  bgam <- apply( (yt-theta0)*cgam, 1, sum)/V
  psis <- (data - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  ys <- cbind(theta0, matrix(yt, ncol=length(yt), nrow=length(W), byrow=TRUE))
  Ly <- ys[,-1] - ys[,-(T.T+1)]
  psiLT <- cbind(0,psis)
  Lpsi <- psiLT[,-1] - psiLT[,-(T.T+1)]
  apsi <- apply(Lpsi^2, 1, sum)/2/W
  bpsi <- apply(Lpsi*Ly, 1, sum)/W
  VWcor <- cor(V,W)
  Vbv <- cor(V, bv)
  Wbv <- cor(W, bv)
  Vbw <- cor(V, bw)
  Wbw <- cor(W, bw)
  Wagam <- cor(W, agam)
  Wbgam <- cor(W, bgam)
  Vapsi <- cor(V, apsi)
  Vbpsi <- cor(V, bpsi)
  out <- data.frame(V=V.T, W=W.T, T=T.T)
  out$VW <- VWcor
  out$Vbv <- Vbv
  out$Wbv <- Wbv
  out$Vbw <- Vbw
  out$Wbw <- Wbw
  out$Wagam <- Wagam
  out$Wbgam <- Wbgam
  out$Vapsi <- Vapsi
  out$Vbpsi <- Vbpsi
  rownames(out)[1] <- ""
  return(out)
}

## Finds autocorrelation and effective sample size info from a
## sample from a given sampler
samsummary <- function(sam, dat, burn, sampler){
  V.T <- sam$V.T[1]
  W.T <- sam$W.T[1]
  T.T <- sam$T.T[1]
  V <- sam$V[-c(1:burn)]
  W <- sam$W[-c(1:burn)]
  theta0s <- sam[-c(1:burn),grep("theta", colnames(sam))]
  theta0s <- theta0s[,1:(T.T + 1)]
  thetas <- theta0s[,-1]
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
###time per 1000 iterations
  time <- sam$time[1]/length(sam$time)*1000 
  stime <- mean(sam$stime)*1000
  logconV <- mean(sam$logconV, na.rm=TRUE)
  adrejV <- mean(sam$adrejV, na.rm=TRUE)
  logconW <- mean(sam$logconW, na.rm=TRUE)
  adrejW <- mean(sam$adrejW, na.rm=TRUE)
  statkern <- mean(sam$kernel=="state", na.rm=TRUE)
  distkern <- mean(sam$kernel=="dist", na.rm=TRUE)
  errorkern <- 1 - (statkern + distkern)
  init <- data.frame(time=time, stime=stime,
                     logconV=logconV,     adrejV=adrejV,
                     logconW=logconW,     adrejW=adrejW,
                     statkern=statkern, distkern=distkern,
                     errorkern=errorkern)
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  thetaAC <- apply(thetas, 2, corfun)
  gammaAC <- apply(gammas, 2, corfun)
  psiAC <- apply(psis, 2, corfun)
  thetaES <- apply(thetas, 2, effectiveSize)
  gammaES <- apply(gammas, 2, effectiveSize)
  psiES <- apply(psis, 2, effectiveSize)
  theta0.AC <- corfun(theta0)
  theta1.AC <- thetaAC[1]
  thetaT4.AC <- thetaAC[ceiling(T.T/4)]
  thetaT2.AC <- thetaAC[T.T/2]
  theta3T4.AC <- thetaAC[floor(3*T.T/4)]
  thetaT.AC <- thetaAC[T.T]
  gamma1.AC <- gammaAC[1]
  gammaT4.AC <- gammaAC[ceiling(T.T/4)]
  gammaT2.AC <- gammaAC[T.T/2]
  gamma3T4.AC <- gammaAC[floor(3*T.T/4)]
  gammaT.AC <- gammaAC[T.T]
  psi1.AC <- psiAC[1]
  psiT4.AC <- psiAC[ceiling(T.T/4)]
  psiT2.AC <- psiAC[T.T/2]
  psi3T4.AC <- psiAC[floor(3*T.T/4)]
  psiT.AC <- psiAC[T.T]
  theta.ACmax <- thetaAC[which.max(abs(thetaAC))]
  gamma.ACmax <- gammaAC[which.max(abs(gammaAC))]
  psi.ACmax <- psiAC[which.max(abs(psiAC))]
  theta.ACavg <- mean(abs(thetaAC))
  gamma.ACavg <- mean(abs(gammaAC))
  psi.ACavg <- mean(abs(psiAC))
  theta0.ES <- effectiveSize(theta0)
  theta1.ES <- thetaES[1]
  thetaT4.ES <- thetaES[ceiling(T.T/4)]
  thetaT2.ES <- thetaES[T.T/2]
  theta3T4.ES <- thetaES[floor(3*T.T/4)]
  thetaT.ES <- thetaES[T.T]
  gamma1.ES <- gammaES[1]
  gammaT4.ES <- gammaES[ceiling(T.T/4)]
  gammaT2.ES <- gammaES[T.T/2]
  gamma3T4.ES <- gammaES[floor(3*T.T/4)]
  gammaT.ES <- gammaES[T.T]
  psi1.ES <- psiES[1]
  psiT4.ES <- psiES[ceiling(T.T/4)]
  psiT2.ES <- psiES[T.T/2]
  psi3T4.ES <- psiES[floor(3*T.T/4)]
  psiT.ES <- psiES[T.T]
  theta.ESmin <- thetaES[which.min(thetaES)]
  gamma.ESmin <- gammaES[which.min(gammaES)]
  psi.ESmin <- psiES[which.min(psiES)]
  theta.ESavg <- mean(thetaES)
  gamma.ESavg <- mean(gammaES)
  psi.ESavg <- mean(psiES)
  V.AC <- corfun(V)
  W.AC <- corfun(W)
  V.ES <- effectiveSize(V)
  W.ES <- effectiveSize(W)
  
  out <- cbind(init, V.AC, W.AC, theta0.AC,
               theta1.AC, thetaT4.AC, thetaT2.AC, theta3T4.AC,
               thetaT.AC, theta.ACmax,theta.ACavg,
               gamma1.AC, gammaT4.AC, gammaT2.AC, gamma3T4.AC,
               gammaT.AC, gamma.ACmax, gamma.ACavg,
               psi1.AC, psiT4.AC, psiT2.AC, psi3T4.AC, psiT.AC,
               psi.ACmax, psi.ACavg,
               V.ES, W.ES, theta0.ES,
               theta1.ES, thetaT4.ES, thetaT2.ES, theta3T4.ES,
               thetaT.ES, theta.ESmin, theta.ESavg,
               gamma1.ES, gammaT4.ES, gammaT2.ES, gamma3T4.ES,
               gammaT.ES, gamma.ESmin, gamma.ESavg,
               psi1.ES, psiT4.ES, psiT2.ES, psi3T4.ES,
               psiT.ES, psi.ESmin, psi.ESavg)
  rownames(out) <- ""
  return(out)
}

## function for finding the first order autocorrelation of a TS
corfun <- function(x){
  acf(x, lag.max=1, plot=FALSE)[[1]][2]
}

## used to quickly create the simulation grid, used because we don't
## want a full grid of starting values, just a grid of dispersed points
## plus one point at the true values
dfun <- function(M, simgrid){
  data.frame(M, simgrid)
}

## Wrapper for quickly simulating from all samplers
samwrap <- function(par, n, samp){
  dat <- par$y[order(par$t)]
  T <- length(dat)
  start <- c(par$V.T[1], par$W.T[1])
  av <- par$av[1]
  aw <- par$aw[1]
  bv <- par$bv[1]
  bw <- par$bw[1]
  m0 <- par$m0[1]
  C0 <- par$C0[1]
  print(paste(c(samp, " T=", T, " V=", start[1], " W=", start[2]), collapse=""))
  if(samp=="state")
    time <- system.time(out <- statesam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="dist")
    time <- system.time(out <- distsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="error")
    time <- system.time(out <- errorsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="sdint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, TRUE, 1))
  if(samp=="seint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, TRUE, 2))
  if(samp=="deint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, TRUE, 3))
  if(samp=="triint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, c(TRUE, TRUE), 4))
  if(samp=="sdalt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, FALSE, 1))
  if(samp=="sealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, FALSE, 2))
  if(samp=="dealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, FALSE, 3))
  if(samp=="trialt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, m0, C0, c(FALSE, FALSE), 4))
  if(samp=="sdkern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, m0, C0, c(1/2, 1/2, 0)))
  if(samp=="sekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, m0, C0, c(1/2, 0, 1/2)))
  if(samp=="dekern")
   time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, m0, C0, c(0, 1/2, 1/2)))
  if(samp=="trikern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="partialcis")
    time <- system.time(out <- partialcissam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="stnstate")
    time <- system.time(out <- stnstatesam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="distda")
    time <- system.time(out <- distsamda(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="errorda")
    time <- system.time(out <- errorsamda(n, start, dat, av, aw, bv, bw, m0, C0))

  outdat <- data.frame(out)
  outdat$time <- time[1]
  cols <- ncol(outdat)
  outdat <- outdat[,c( cols, 1:(cols-1) )]
  print(paste(c(samp, " T=", T, " V=", start[1], " W=", start[2], " FINISHED"), collapse=""))
  return(outdat)
  
}


## A wrapper for quickly simulating from each of the interweaving/
## alternating samplers.
samwrapper <- function(n, start, dat, av, aw, bv, bw, m0, C0, inter, samp){
  if(samp==1){
    out <- statedistinter(n, start, dat, av, aw, bv, bw, m0, C0, inter)
  }
  if(samp==2){
    out <- stateerrorinter(n, start, dat, av, aw, bv, bw, m0, C0, inter)
  }
  if(samp==3){
    out <- disterrorinter(n, start, dat, av, aw, bv, bw, m0, C0, inter)
  }
  if(samp==4){
    out <- tripleinter(n, start, dat, av, aw, bv, bw, m0, C0, c(inter, inter))
  }
  return(out)
}

## Wrapper for llsim for use with expand.grid and ddply
lldsim <- function(df, m0, C0){
  T <- df$T.T[1]
  V <- df$V.T[1]
  W <- df$W.T[1]
  out <- llsim(T, V, W, m0, C0)
  return(data.frame(t=1:T, y=out))
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

samoutsetup <- function(n, T){
  out <- data.frame(matrix(0, nrow=n, ncol=T+9))
  colnames(out) <- c("logconV", "adrejV", 
                     "logconW", "adrejW", 
                     "kernel", "stime", "V", "W",
                     paste("theta",0:T,sep=""))
  return(out)
}

## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
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
    out[i,] <- c(NA,NA,NA,NA,NA,smoothtime,V,W,theta)
  }
  return(out)
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,NA,NA,NA, smoothtime, V, W, theta)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(NA,NA,rejsW, NA, smoothtime, V, W, theta)
  }
  return(out)
}

## state + dist interveaving or alternating sampler
statedistinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampler
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled disturbance sampler
    if(!inter){
      ptma <- proc.time()
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, W)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      gam <- gamtrans(theta, W)
    }
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(NA,NA,rejsW, NA, smoothtime, V, W, theta)
  }
  return(out)
}

## state + error interveaving/alternating sampler
stateerrorinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampnler
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled error sampler
    if(!inter){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      psi <- psitrans(dat, theta, V)
    }
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,NA,NA,NA, smoothtime, V, W, theta)
  }
  return(out)
}

## dist + error interveaving/alternating sampler
disterrorinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## scaled disturbance sampler
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## scaled error sampler
    if(!inter){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      psi <- psigamtrans(dat, gam, V, W)
    }
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,rejsW, NA, smoothtime, V, W, theta)
  }
    return(out)
}

## state + dist + error interweaving/alternating sampler
tripleinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, m0=0, C0=10^7, inter=c(TRUE, TRUE)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampler
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled disturbance sampler
    if(!inter[1]){
      ptma <- proc.time()
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, W)
      ptmb <- proc.time()
      smoothtime <- smoothtime + ptmb[3]-ptma[3]
    }
    else{
      gam <- gamtrans(theta, W)
    }
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## scaled error sampler
    if(!inter[2]){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- smoothtime + ptmb[3]-ptma[3]
    }
    else{
      psi <- psigamtrans(dat, gam, V, W)
    }
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,rejsW, NA, smoothtime, V, W, theta)
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
