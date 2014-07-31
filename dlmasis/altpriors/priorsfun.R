## A set of functions for simulating from and fitting local level models
library(dlm)
library(coda)
library(MCMCpack)
library(ars)
library(plyr)
library(GeneralizedHyperbolic)
source("../code/mcfa.R")

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
  colnam <- grep("(V.T|W.T|T.T|V|W|sigv|sigw|time|kernel)$",
                 colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samshort)
  out <- ddply(sam, .(V.T, W.T, T.T), samsummary,
               .parallel=parallel, dat=simdata, burn=burn,
               sampler=sampler)
#  if(sampler=="trialt"){
#    postcors <- ddply(sam, .(V.T, W.T, T.T), postcor, .parallel=parallel,
#                           dat=simdata, burn=burn)
#    save(postcors, file="postcors.RData")
#    rm(postcors)
#  }
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
  sigv <- sam$sigv[-c(1:burn)]
  sigw <- sam$sigw[-c(1:burn)]
  theta0s <- sam[-c(1:burn),grep("theta", colnames(sam))]
  theta0s <- theta0s[,1:(T.T+1)]
  thetas <- theta0s[,-1]
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  VWcor <- cor(V,W)
  svswcor <- cor(sigv,sigw)
  Vth0cor <- cor(V,theta0)
  Wth0cor <- cor(W,theta0)
  Vthcors <- cor(V, thetas)
  Wthcors <- cor(W, thetas)
  Vgacors <- cor(V, gammas)
  Wgacors <- cor(W, gammas)
  Vpscors <- cor(V, psis)
  Wpscors <- cor(W, psis)
  Vthmaxcor <- Vthcors[1,which.max(abs(Vthcors))]
  Wthmaxcor <- Wthcors[1,which.max(abs(Wthcors))]
  Vgamaxcor <- Vgacors[1,which.max(abs(Vgacors))]
  Wgamaxcor <- Wgacors[1,which.max(abs(Wgacors))]
  Vpsmaxcor <- Vpscors[1,which.max(abs(Vpscors))]
  Wpsmaxcor <- Wpscors[1,which.max(abs(Wpscors))]
  Vthavgcor <- mean(abs(Vthcors[1,]))
  Wthavgcor <- mean(abs(Wthcors[1,]))
  Vgaavgcor <- mean(abs(Vgacors[1,]))
  Wgaavgcor <- mean(abs(Wgacors[1,]))
  Vpsavgcor <- mean(abs(Vpscors[1,]))
  Wpsavgcor <- mean(abs(Wpscors[1,]))

  out <- data.frame(VW=VWcor, svsw=svswcor, Vtheta0=Vth0cor, Wtheta0=Wth0cor,
                    Vtheta1=Vthcors[1,1], Wtheta1=Wthcors[1,1],
                    VthetaT=Vthcors[1,T.T], WthetaT=Wthcors[1,T.T],
                    Vgamma1=Vgacors[1,1], Wgamma1=Wgacors[1,1],
                    VgammaT=Vgacors[1,T.T], WgammaT=Wgacors[1,T.T],
                    Vpsi1=Vpscors[1,1], Wpsi1=Wpscors[1,1],
                    VpsiT=Vpscors[1,T.T], WpsiT=Wpscors[1,T.T],
                    VthetaMax=Vthmaxcor, WthetaMax=Wthmaxcor,
                    VgammaMax=Vgamaxcor, WgammaMax=Wgamaxcor,
                    VpsiMax=Vpsmaxcor, WpsiMax=Wpsmaxcor,
                    VthetaAvg=Vthavgcor, WthetaAvg=Wthavgcor,
                    VgammaAvg=Vgaavgcor, WgammaAvg=Wgaavgcor,
                    VpsiAvg=Vpsavgcor, WpsiAvg=Wpsavgcor,
                    VthT4 = Vthcors[1, ceiling(T.T/4)],
                    VthT2 = Vthcors[1, T.T/2],
                    Vth3T4 = Vthcors[1, floor(3*T.T/4)],
                    WthT4 = Wthcors[1, ceiling(T.T/4)],
                    WthT2 = Wthcors[1, T.T/2],
                    Wth3T4 = Wthcors[1, floor(3*T.T/4)],
                    VgaT4 = Vgacors[1, ceiling(T.T/4)],
                    VgaT2 = Vgacors[1, T.T/2],
                    Vga3T4 = Vgacors[1, floor(3*T.T/4)],
                    WgaT4 = Wgacors[1, ceiling(T.T/4)],
                    WgaT2 = Wgacors[1, T.T/2],
                    Wga3T4 = Wgacors[1, floor(3*T.T/4)],
                    VpsT4 = Vpscors[1, ceiling(T.T/4)],
                    VpsT2 = Vpscors[1, T.T/2],
                    Vps3T4 = Vpscors[1, floor(3*T.T/4)],
                    WpsT4 = Wpscors[1, ceiling(T.T/4)],
                    WpsT2 = Wpscors[1, T.T/2],
                    Wps3T4 = Wpscors[1, floor(3*T.T/4)])
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
  sigv <- sam$sigv[-c(1:burn)]
  sigw <- sam$sigw[-c(1:burn)]
  theta0s <- sam[-c(1:burn),grep("theta", colnames(sam))]
  theta0s <- theta0s[,1:(T.T + 1)]
  thetas <- theta0s[,-1]
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  stime <- mean(sam$stime)*1000
  time <- sam$time[1]/length(sam$time)*1000
  logconV <- mean(sam$logconV, na.rm=TRUE)
  adrejV <- mean(sam$adrejV, na.rm=TRUE)
  logconW <- mean(sam$logconW, na.rm=TRUE)
  adrejW <- mean(sam$adrejW, na.rm=TRUE)
  statkern <- mean(sam$kernel=="state", na.rm=TRUE)
  distkern <- mean(sam$kernel=="dist", na.rm=TRUE)
  errorkern <- mean(sam$kernel=="error", na.rm=TRUE)
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
  sigv.AC <- corfun(sigv)
  sigw.AC <- corfun(sigw)
  sigv.ES <- effectiveSize(sigv)
  sigw.ES <- effectiveSize(sigw)
  
  out <- cbind(init, V.AC, W.AC, sigv.AC, sigw.AC, theta0.AC,
               theta1.AC, thetaT4.AC, thetaT2.AC, theta3T4.AC,
               thetaT.AC, theta.ACmax,theta.ACavg,
               gamma1.AC, gammaT4.AC, gammaT2.AC, gamma3T4.AC,
               gammaT.AC, gamma.ACmax, gamma.ACavg,
               psi1.AC, psiT4.AC, psiT2.AC, psi3T4.AC, psiT.AC,
               psi.ACmax, psi.ACavg,
               V.ES, W.ES, sigv.ES, sigw.ES, theta0.ES,
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
  Qv <- par$V.T[1]
  Qw <- par$W.T[1]
  m0 <- par$m0[1]
  C0 <- par$C0[1]
  start <- c(Qv, Qw)
  print(paste(c("sampler: ", samp, "  T=", T, "  V=", Qv, "  W=", Qw), collapse=""))
  if(samp=="state")
    time <- system.time(out <- statesam(n, start, dat, Qv, Qw, m0, C0))
  if(samp=="dist")
    time <- system.time(out <- distsam(n, start, dat, Qv, Qw, m0, C0))
  if(samp=="error")
    time <- system.time(out <- errorsam(n, start, dat, Qv, Qw, m0, C0))
  if(samp=="sdint")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, TRUE, 1))
  if(samp=="seint")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, TRUE, 2))
  if(samp=="deint")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, TRUE, 3))
  if(samp=="triint")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, c(TRUE, TRUE), 4))
  if(samp=="sdalt")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, FALSE, 1))
  if(samp=="sealt")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, FALSE, 2))
  if(samp=="dealt")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, FALSE, 3))
  if(samp=="trialt")
    time <- system.time(out <- samwrapper(n, start, dat, Qv, Qw, m0, C0, c(FALSE, FALSE), 4))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, Qv, Qw, m0, C0))
    outdat <- data.frame(out)
  outdat$time <- time[1]
  cols <- ncol(outdat)
  outdat <- outdat[,c( cols, 1:(cols-1) )]
  print(paste(c(samp, " T=", T, " V=", start[1], " W=", start[2], " FINISHED"), collapse=""))
  return(outdat)
}


## A wrapper for quickly simulating from each of the interweaving/
## alternating samplers.
samwrapper <- function(n, start, dat, Qv, Qw, m0, C0, inter, samp){
  if(samp==1){
    out <- statedistinter(n, start, dat, Qv, Qw, m0, C0, inter)
  }
  if(samp==2){
    out <- stateerrorinter(n, start, dat, Qv, Qw, m0, C0, inter)
  }
  if(samp==3){
    out <- disterrorinter(n, start, dat, Qv, Qw, m0, C0, inter)
  }
  if(samp==4){
    out <- tripleinter(n, start, dat, Qv, Qw, m0, C0, c(inter, inter))
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
  out <- data.frame(matrix(0, nrow=n, ncol=T+11))
  colnames(out) <- c("logconV", "adrejV", 
                     "logconW", "adrejW", 
                     "kernel", "stime", "V", "W",
                     "sigv", "sigw",
                     paste("theta",0:T,sep=""))
  return(out)
}

## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VWiter <- VWthetaiter(dat, theta, Qv, Qw)
    V <- VWiter[1]
    W <- VWiter[2]
    sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    out[i,] <- c(NA,NA,NA,NA,NA,smoothtime,V,W,sigv,sigw,theta)
  }
  return(out)
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsam <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    Vout <- Vpsiiter(dat, psi, W, Qv)
    V <- Vout[1]
    sigv <- Vout[2]
    W <- Wpsiiter(dat, psi, sigv, Qw)
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    theta <- thetapsitrans(dat, psi, sigv)
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, sigv, sigw, theta)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsam <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
  for(i in 1:n){
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, sigw)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vgamiter(dat, gam, sigw, Qv)
    Wout <- Wgamiter(dat, gam, V, Qw) 
    W <- Wout[1]
    sigw <- Wout[2]
    theta <- thetagamtrans(gam, sigw)
    sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, sigv, sigw, theta)

  }
  return(out)
}

## state + dist interveaving or alternating sampler
statedistinter <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7, inter=TRUE){
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
    VWiter <- VWthetaiter(dat, theta, Qv, Qw)
    V <- VWiter[1]
    W <- VWiter[2]
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    ## scaled disturbance sampler
    if(!inter){
      ptma <- proc.time()
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, sigw)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      gam <- gamtrans(theta, sigw)
    }
    V <- Vgamiter(dat, gam, sigw, Qv)
    Wout <- Wgamiter(dat, gam, V, Qw)
    W <- Wout[1]
    sigw <- Wout[2]
    theta <- thetagamtrans(gam, sigw)
    sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, sigv, sigw, theta)
  }
  return(out)
}

## state + error interveaving/alternating sampler
stateerrorinter <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7, inter=TRUE){
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
    VWiter <- VWthetaiter(dat, theta, Qv, Qw)
    V <- VWiter[1]
    W <- VWiter[2]
    sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
    ## scaled error sampler
    if(!inter){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      psi <- psitrans(dat, theta, sigv)
    }
    Vout <- Vpsiiter(dat, psi, W, Qv)
    V <- Vout[1]
    sigv <- Vout[2]
    W <- Wpsiiter(dat, psi, sigv,  Qw)
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    theta <- thetapsitrans(dat, psi, sigv)
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, sigv, sigw, theta)
  }
  return(out)
}

## dist + error interveaving/alternating sampler
disterrorinter <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## scaled disturbance sampler
    ptma <- proc.time()
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, sigw)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    V <- Vgamiter(dat, gam, sigw, Qv)
    Wout <- Wgamiter(dat, gam, V, Qw)
    W <- Wout[1]
    sigw <- Wout[2]
    ## scaled error sampler
    if(!inter){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- ptmb[3]-ptma[3] + smoothtime
    }
    else{
      sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
      theta <- thetagamtrans(gam, sigw)
      psi <- psitrans(dat, theta, sigv)
    }
    Vout <- Vpsiiter(dat, psi, W, Qv)
    V <- Vout[1]
    sigv <- Vout[2]
    W <- Wpsiiter(dat, psi, sigv, Qw)
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    theta <- thetapsitrans(dat, psi, sigv)
    out[i,] <- c(rep(NA,5), smoothtime, V, W, sigv, sigw, theta)
  }
    return(out)
}

## state + dist + error interweaving/alternating sampler
tripleinter <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7, inter=c(TRUE, TRUE)){
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
    VWiter <- VWthetaiter(dat, theta, Qv, Qw)
    V <- VWiter[1]
    W <- VWiter[2]
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    ## scaled disturbance sampler
    if(!inter[1]){
      ptma <- proc.time()
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, sigw)
      ptmb <- proc.time()
      smoothtime <- smoothtime + ptmb[3]-ptma[3]
    }
    else{
      gam <- gamtrans(theta, sigw)
    }
    V <- Vgamiter(dat, gam, sigw, Qv)
    Wout <- Wgamiter(dat, gam, V, Qw)
    W <- Wout[1]
    sigw <- Wout[2]
    ## scaled error sampler
    if(!inter[2]){
      ptma <- proc.time()
      psi <- mcfapssmooth(dat, V, W, m0, C0)
      ptmb <- proc.time()
      smoothtime <- smoothtime + ptmb[3]-ptma[3]
    }
    else{
      sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
      theta <- thetagamtrans(gam, sigw)
      psi <- psitrans(dat, theta, sigv)
    }
    Vout <- Vpsiiter(dat, psi, W, Qv)
    V <- Vout[1]
    sigv <- Vout[2]
    W <- Wpsiiter(dat, psi, sigv,  Qw)
    theta <- thetapsitrans(dat, psi, sigv)
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    out[i,] <- c(rep(NA,5), smoothtime, V, W, sigv, sigw, theta)
  }
    return(out)
}


## transforms theta to gamma
gamtrans <- function(theta, sigw){
  T <- length(theta) - 1
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  A <- B/sigw
  A[1,1] <- 1
  gam <- A%*%matrix(theta, ncol=1)
  return(gam)
}

## transforms theta to psi
psitrans <- function(dat, theta, sigv){
  T <- length(dat)
  A <- diag(-1/sigv, T+1 )
  A[1,1] <- 1
  ytild <- c(0, dat/sigv)
  psi <- ytild + A%*%matrix(theta, ncol=1)
  return(psi)
}

## transforms gamma to theta
thetagamtrans <- function(gam, sigw){
  T <- length(gam) - 1
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  A <- B/sigw
  A[1,1] <- 1
  theta <- solve(A)%*%matrix(gam, ncol=1)
  return(theta)
}

## transforms psi to theta
thetapsitrans <- function(dat, psi, sigv){
  T <- length(dat)
  A <- diag(-1/sigv, T+1 )
  A[1,1] <- 1
  ytild <- c(0, dat/sigv)
  theta <- solve(A)%*%(matrix(psi, ncol=1) - matrix(ytild, ncol=1))
  return(theta)
}

rgig2 <- function(n, a, b, p, eps=.0001){
  ## transform into a GIG RV that won't have numerical probs if
  ## a or b are near zero
  ## note: a is the coef of 1/x, b is the coef of x in dgig(x)
  if(a < eps | b < eps){
    g <- sqrt(b/a)
    x <- rgig(n, a*g, b/g, p)/g
  }
  else{
    x <- rgig(n, a, b, p)
  }
  return(x)
}

## samples V,W conditional on theta
VWthetaiter <- function(dat, theta, Qv, Qw){
  T <- length(dat)
  av <- sum((dat-theta[-1])^2)
  bv <- 1/Qv
  p <- -(T-1)/2
  aw <- sum((theta[-1]-theta[-(T+1)])^2)
  bw <- 1/Qw
  V <- rgig2(1, av, bv, p)
  W <- rgig2(1, aw, bw, p)
  return(c(V,W))
}
## samples W conditional on V,gamma
Wgamiter <- function(dat, gam, V, Qw){
  T <- length(dat)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  Sighatinv <- 1/Qw + sum(cgam^2)/V
  Sighat <- 1/Sighatinv
  muhat <- Sighat*sum( (dat - gam0)*cgam )/V
  sigw <- rnorm(1, muhat, sqrt(Sighat))
  W <- sigw^2
  return(c(W,sigw))
}

## samples V conditional on sigw,gamma
Vgamiter <- function(dat, gam, sigw, Qv){
  T <- length(dat)
  theta <- thetagamtrans(gam, sigw)
  av <- sum((dat-theta[-1])^2)
  bv <- 1/Qv
  p <- -(T-1)/2
  V <- rgig2(1, av, bv, p)
  return(V)
}

## samples V conditional on theta
Vthetaiter <- function(dat, theta,  Qv){
  T <- length(dat)
  av <- sum((dat-theta[-1])^2)
  bv <- 1/Qv
  p <- -(T-1)/2
  V <- rgig2(1, av, bv, p)
  return(V)
}


## samples V conditional on W,psi 
Vpsiiter <- function(dat, psi, W, Qv){
  T <- length(dat)
  ys <- c(psi[1],dat)
  Dys <- ys[-1] - ys[-(T+1)]
  psistar <- c(0, psi[-1])
  Dpsistar <- psistar[-1] - psistar[-(T+1)]
  Sighatinv <- 1/Qv + 1/W*sum(Dpsistar^2)
  Sighat <- 1/Sighatinv
  muhat <- Sighat/W*sum(Dys*Dpsistar)
  sigv <- rnorm(1,muhat,sqrt(Sighat))
  V <- sigv^2
  return(c(V,sigv))
}


## samples W conditional on sigv,psi
Wpsiiter <- function(dat, psi, sigv, Qw){
  T <- length(dat)
  theta <- thetapsitrans(dat, psi, sigv)
  theta <- c(theta)
  aw <- sum((theta[-1]-theta[-(T+1)])^2)
  bw <- 1/Qw
  p <- -(T-1)/2
  W <- rgig2(1, aw, bw, p) 
  return(W)
}

## samples W conditional on theta
Wthetaiter <- function(dat, theta, Qw){
  T <- length(dat)
  aw <- sum((theta[-1]-theta[-(T+1)])^2)
  bw <- 1/Qw
  p <- -(T-1)/2
  W <- rgig2(1, aw, bw, p) 
  return(W)
}

fullcissam <- function(n, start, dat, Qv, Qw, m0=0, C0=10^7){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  p <- -(T-1)/2
  out <- samoutsetup(n, T)  
  for(i in 1:n){
    ptma <- proc.time()
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    ## V step
    Vout <- Vpsiiter(dat, psi, W, Qv)
    V <- Vout[1]
    sigv <- Vout[2]
    theta <- thetapsitrans(dat, psi, sigv)
    av <- sum((dat-theta[-1])^2)
    bv <- 1/Qv
    V <- rgig2(1, av, bv, p)
    sigv <- (2*rbinom(1,1,1/2)-1)*sqrt(V)
    ## W step
    aw <- sum((theta[-1]-theta[-(T+1)])^2)
    bw <- 1/Qw
    W <- rgig2(1, aw, bw, p) 
    sigw <- (2*rbinom(1,1,1/2)-1)*sqrt(W)
    gam <- gamtrans(theta, sigw)
    Wout <- Wgamiter(dat, gam, V, Qw)
    W <- Wout[1]
    sigw <- Wout[2]
    theta <- thetagamtrans(gam, sigw)
    out[i,] <- c(rep(NA,5), smoothtime, V, W, sigv, sigw, theta)
  }
  return(out)
}
