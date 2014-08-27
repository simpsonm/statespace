## A set of functions for simulating from and fitting local level models

library(dlm)
library(coda)
library(MCMCpack)
library(ars)
library(plyr)
source("mcfa.R") ## code for performing MCFA
source("wscalerej.R") ## code for wrongly scaled samplers

## simulates from the posterior given a data set and a list of given samplers.
## samplers: vector of sampler names
## simdata: data frame containing data and hyperparameters
## n: posterior sample size
## burn: burn in size
## parallel: if TRUE, will use DOMC to parallelize sampling
## returns a dataframe of summary infor for all samplers
fullsim <- function(samplers, simdata, n, burn, parallel){
  out <- ddply(samplers, .(sampler), samsim, simdata=simdata, n=n,
               burn=burn, parallel=parallel)
  return(out)
}


## simulates from a given sampler for each dataset and for multiple
## chains, and returns summary info on the first chain.
## called by fullsim
samsim <- function(samplers, simdata, n, burn, parallel){
  sampler <- samplers$sampler[1]
  print(sampler)
  ## Simulate from the posterior for each dataset using a given sampler
  ## V.T, W.T and T.T denotes the true values of V, W and T used to simulate the series
  sam <- ddply(simdata, .(V.T, W.T, T.T), samwrap, .parallel=parallel,
               n=n, samp=sampler)
  samnam <- paste(sampler, "SAM.RData", sep="")
  colnam <- grep("(V.T|W.T|T.T|V|W|time|stime|logconW|adrejW|logconWls|adrejWls|logconV|adrejV|logconVls|adrejVls|kernel)$", colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samshort)
  ## summarize the posterior information for that sampler
  out <- ddply(sam, .(V.T, W.T, T.T), samsummary,
               .parallel=parallel, dat=simdata, burn=burn,
               sampler=sampler)
  ## Obtain posterior correlation information using only one of the samplers
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

## Finds the posterior correlation between various quantities
## called by samsim
## sam: posterior sample
## dat: data frame containing each of the time series
## burn: burn in size
postcor <- function(sam, dat, burn){
  ## True values of V, W, and T used to simulate the series
  V.T <- sam$V.T[1]
  W.T <- sam$W.T[1]
  T.T <- sam$T.T[1]
  V <- sam$V[-c(1:burn)]
  W <- sam$W[-c(1:burn)]
  theta0s <- sam[-c(1:burn),grep("theta", colnames(sam))] 
  theta0s <- theta0s[,1:(T.T+1)] ## the thetas including theta_0
  thetas <- theta0s[,-1] ## thetas not including theta_0
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  ## compute correlations between V, W, and various quantities
  VWcor <- cor(V,W)
  Vth0cor <- cor(V,theta0)
  Wth0cor <- cor(W,theta0)
  Vthcors <- cor(V, thetas)
  Wthcors <- cor(W, thetas)
  Vgacors <- cor(V, gammas)
  Wgacors <- cor(W, gammas)
  Vpscors <- cor(V, psis)
  Wpscors <- cor(W, psis)
  ## Find the maximum correlation between, e.g., V and any of the theta's
  Vthmaxcor <- Vthcors[1,which.max(abs(Vthcors))]
  Wthmaxcor <- Wthcors[1,which.max(abs(Wthcors))]
  Vgamaxcor <- Vgacors[1,which.max(abs(Vgacors))]
  Wgamaxcor <- Wgacors[1,which.max(abs(Wgacors))]
  Vpsmaxcor <- Vpscors[1,which.max(abs(Vpscors))]
  Wpsmaxcor <- Wpscors[1,which.max(abs(Wpscors))]
  ## Find the average correlation between, e.g., V and any of the theta's
  Vthavgcor <- mean(abs(Vthcors[1,]))
  Wthavgcor <- mean(abs(Wthcors[1,]))
  Vgaavgcor <- mean(abs(Vgacors[1,]))
  Wgaavgcor <- mean(abs(Wgacors[1,]))
  Vpsavgcor <- mean(abs(Vpscors[1,]))
  Wpsavgcor <- mean(abs(Wpscors[1,]))

  out <- data.frame(VW=VWcor, Vtheta0=Vth0cor, Wtheta0=Wth0cor,
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
## called by samsim
## sam: posterior sample
## dat: dataframe containing each simulated time series
## burn: burn in size
## sampler: name of the sampler used to generate the posterior draws
samsummary <- function(sam, dat, burn, sampler){
  ## True values of V, W, and T
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
  ## time per 1000 iterations
  time <- sam$time[1]/length(sam$time)*1000
  ## tracks when relevant densities are log concave and when
  ## adaptive rejection sampling is used
  logconV <- mean(sam$logconV, na.rm=TRUE)
  adrejV <- mean(sam$adrejV, na.rm=TRUE)
  logconW <- mean(sam$logconW, na.rm=TRUE)
  adrejW <- mean(sam$adrejW, na.rm=TRUE)
  init <- data.frame(time=time,
                     logconV=logconV,     adrejV=adrejV,
                     logconW=logconW,     adrejW=adrejW)
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  ## Find autocorrelation and effective sample size for various quantities
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

## convenience function for finding the first order
## autocorrelation of a time series
corfun <- function(x){
  acf(x, lag.max=1, plot=FALSE)[[1]][2]
}

## Wrapper for quickly simulating from all samplers
## Called by samsim
## par: a list containing relevant parameters, including the data
## n: desired sample size (including burn in)
## samp: name of the desired sampler
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
  if(samp=="trikern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="partialcis")
    time <- system.time(out <- partialcissam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="stnstate")
    time <- system.time(out <- stnstatesam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="wdist")
    time <- system.time(out <- wdistsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="werror")
    time <- system.time(out <- werrorsam(n, start, dat, av, aw, bv, bw, m0, C0))

  outdat <- data.frame(out)
  outdat$time <- time[1]
  cols <- ncol(outdat)
  outdat <- outdat[,c( cols, 1:(cols-1) )]
  print(paste(c(samp, " T=", T, " V=", start[1], " W=", start[2], " FINISHED"), collapse=""))
  return(outdat)
}


## A wrapper for quickly simulating from each of the interweaving/
## alternating samplers.
## Called by samwrap
## n: desired sample size
## start: starting values for the chain
## dat: time series of data
## av, aw, bv, bw: hyperparameters for the independent IG priors for V, W
## m0, C0: hyperparmaeters for Gaussian prior on theta_0
## inter: if TRUE interweave, else alternate
## samp: number indicating desired sample: SD=1, SE=2, DE=3 and Triple=4.
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
## df: initialized data frame of V, W, and T.
## m0, C0: prior for simulating theta_0
lldsim <- function(df, m0, C0){
  T <- df$T.T[1]
  V <- df$V.T[1]
  W <- df$W.T[1]
  out <- llsim(T, V, W, m0, C0)
  return(data.frame(t=1:T, y=out))
}


## Simulates from a local level model
## T: length of the time series
## V: observation variance
## W: system variance
## m0: mean of Gaussian prior on theta_0
## C0: variance of Gaussian prior on theta_0
llsim <- function(T, V, W, m0, C0){
  out <- rep(0,T)
  theta <- rnorm(1, m0, sqrt(C0))
  for(t in 1:T){
    theta <- theta + rnorm(1, 0, sd=sqrt(W))
    out[t] <- theta + rnorm(1, 0, sd=sqrt(V))
  }
  return(out)
}

## Function for quickly setting up the output data frame for a given sampler
## Called by each individual sampler
samoutsetup <- function(n, T){
  out <- data.frame(matrix(0, nrow=n, ncol=T+7))
  colnames(out) <- c("logconV", "adrejV", "logconW", "adrejW", 
                     "V", "W", paste("theta",0:T,sep=""))
  return(out)
}

## state sampler: samples V and W conditional on states
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
statesam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## draw theta using MCFA
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ## draw (V,W) from independent IG distributions
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    out[i,] <- c(NA,NA,NA,NA,V,W,theta)
  }
  return(out)
}

## scaled error sampler: samples V and W conditional on the scaled observation errors
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
errorsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## draw psi using MCFA
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ## draw V from complicated full conditional
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    ## draw W from IG full conditional
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    ## transform back to theta for storage purposes
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,NA,NA, V, W, theta)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
distsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## draw theta using MCFA then transform to gamma
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ## draw V from IG full conditional
    V <- Vgamiter(dat, gam, W, av, bv)
    ## draw W from complicated full conditional
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## transform back to theta for storage purposes
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(NA,NA,rejsW, V, W, theta)
  }
  return(out)
}

## state + dist interveaving or alternating sampler
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
## inter: TRUE if interweaving, FALSE if alternating
statedistinter <- function(n, start, dat, av, aw, bv, bw, m0, C0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampler
    ## draw theta using MCFA
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ## draw (V,W) from independent IG distributions
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    if(!inter){
      ## alternate to gamma
      ## draw theta using MCFA and transform to gamma
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, W)
    }
    else{
      ## interweave to gamma
      ## transform from theta to gamma
      gam <- gamtrans(theta, W)
    }
    ## scaled disturbance sampler
    ## draw V from IG full conditional
    V <- Vgamiter(dat, gam, W, av, bv)
    ## draw W from complicated full conditional
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## transform back to theta for storage purposes
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(NA,NA,rejsW, V, W, theta)
  }
  return(out)
}

## state + error interveaving/alternating sampler
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
## inter: TRUE if interweaving, FALSE if alternating
stateerrorinter <- function(n, start, dat, av, aw, bv, bw, m0, C0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampnler
    ## draw theta using MCFA
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ## draw (V,W) from independent IG full conditionals
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    if(!inter){
      ## alternate to psi
      ## draw psi using MCFA
      psi <- mcfapssmooth(dat, V, W, m0, C0)
    }
    else{
      ## interweave to psi
      ## transform to psi
      psi <- psitrans(dat, theta, V)
    }
    ## scaled error sampler
    ## draw psi from complicated full conditional
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    ## draw W from IG full conditional
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    ## transform back to theta for storage
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,NA,NA, V, W, theta)
  }
  return(out)
}

## dist + error interveaving/alternating sampler
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
## inter: TRUE if interweaving, FALSE if alternating
disterrorinter <- function(n, start, dat, av, aw, bv, bw, m0, C0=, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## scaled disturbance sampler
    ## draw theta using MCFA and transform to gamma
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ## draw V from IG full conditional
    V <- Vgamiter(dat, gam, W, av, bv)
    ## draw W from complicated full conditional
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    if(!inter){
      ## alternate to psi
      ## draw psi using MCFA
      psi <- mcfapssmooth(dat, V, W, m0, C0)
    }
    else{
      ## interweave to psi
      ## transform from gamma to psi
      psi <- psigamtrans(dat, gam, V, W)
    }
    ## scaled error sampler
    ## draw V from complicated full conditional
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    ## draw W from IG full conditional
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    ## transform back to theta for storage
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,rejsW, V, W, theta)
  }
  return(out)
}

## state + dist + error interweaving/alternating sampler
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
## inter: vector of two elements,
##   (TRUE, TRUE) if interweaving, (FALSE, FALSE) if alternating
tripleinter <- function(n, start, dat, av, aw, bv, bw, m0, C0, inter=c(TRUE, TRUE)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ## state sampler
    ## draw theta using MCFA
    theta <- mcfathsmooth(dat, V, W, m0, C0)
    ## draw (V,W) from independent IG full conditionals
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    if(!inter[1]){
      ## alternate to gamma
      ## draw theta using MCFA and transform to gamma
      theta <- mcfathsmooth(dat, V, W, m0, C0)
      gam <- gamtrans(theta, W)
    }
    else{
      ## interweave to gamma
      ## transform from theta to gamma
      gam <- gamtrans(theta, W)
    }
    ## scaled disturbance sampler
    ## draw V from IG full conditional
    V <- Vgamiter(dat, gam, W, av, bv)
    ## draw W from complicated full conditional
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    if(!inter[2]){
      ## alternate to psi
      ## draw psi using MCFA
      psi <- mcfapssmooth(dat, V, W, m0, C0)
    }
    else{
      ## interweave to psi
      ## transform from gamma to psi
      psi <- psigamtrans(dat, gam, V, W)
    }
    ## scaled error sampler
    ## draw V from complicated full conditional
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    ## draw W from IG full conditional
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    ## transform back to theta for storage
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(rejsV,rejsW, V, W, theta)
  }
    return(out)
}

## CIS interweaving sampler
## n: posterior sample size, including burn in
## start: starting values for (V,W)
## av, aw, bv, bw: hyperparameters for independent IG priors on V and W
## m0, C0: hperparameters for Gaussian prior on theta_0=gamma_0=psi_0
## inter: vector of two elements,
##   (TRUE, TRUE) if interweaving, (FALSE, FALSE) if alternating
fullcissam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  Va <- av + T/2
  Wa <- aw + T/2
  out <- samoutsetup(n, T)  
  for(i in 1:n){
    ## V step
    ## draw psi using MCFA
    psi <- mcfapssmooth(dat, V, W, m0, C0)
    ## draw V from complicated full conditional
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    ## transform from psi to theta
    theta <- thetapsitrans(dat, psi, V)
    ## draw V from IG full conditional
    Vb <- bv + sum((dat-theta[-1])^2)/2
    V <- rinvgamma(1, Va, Vb)
    ## W step
    ## draw W from IG full conditional
    Wb <- bw + sum((theta[-1]-theta[-(T+1)])^2)/2
    W <- rinvgamma(1, Wa, Wb)
    ## transform from theta to gamma
    gam <- gamtrans(theta, W)
    ## draw W from complicated full conditional
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    ## transform back to theta for storage
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(Vout[2:3], Wout[2:3], V, W, theta)
  }
  return(out)
}

## returns TRUE if target density is log-concave, FALSE otherwise
## used in sampling from complicated full conditional using SDs or SEs
## if TRUE is returned, adaptive rejection sampling will be used
## b, avw, and bvw are hyperparameters
## eps allows for numerical error in computing condition
logcon <- function(b, avw, bvw, eps=.1){
  RHS <- sqrt((avw + 1)^3/bvw)*4/3*sqrt(2/3) + eps
  out <- (b > RHS)
  return(out)
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
## used in t-approximation for rejection sampling
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

## rejection sampler for drawing V (W) conditional on W (V), psi (gamma) and the data
## uses t approximation for drawing log(V) ( log(W) ) then transforms
## a, b, avw, bvw: hyperparameters for the target density
VWrejiter <- function(a, b, avw, bvw){
  ## set the location and scale of t approximation
  mn <- optimize(logpilVW, c(-10^2,10^2), maximum=TRUE, a=a, b=b, avw=avw, bvw=bvw)$maximum
  propvar <- - 1 /( -a*exp(mn) + b*exp(mn/2)/4 - bvw*exp(-mn) )
  ## degrees of freedom always set to 1
  df <- 1
  rej <- TRUE
  rejit <- 1 ## rejection iteration counter
  ## Bound the approximation error, logtarget - logprop <= M
  M <- optimize(logpirej, c(-10^2, 10^2), a=a, b=b, avw=avw, bvw=bvw, mn=mn, propvar=propvar, df=df, maximum=TRUE)$objective
  ## Sample proposal points until one is accepted
  ## (This step can sometimes be slow)
  while(rej){
    prop <- rtprop(1, mn, propvar, df)
    R <- logpirej(prop, a, b, avw, bvw, mn, propvar, df) - M
    u <- runif(1,0,1)
    if(log(u)<R){
      W <- exp(prop)
      rej <- FALSE
    }
    rejit <- rejit + 1
  }
  return(W)
}

## samples W conditional on V, gamma from its complicated full conditional
## dat: time series
## gam: scaled disturbances, including gamma_0
## aw, bw: hyperparmeters of W's IG prior
Wgamiter <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  ## compute parameters of target density
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2/V
  b <- sum( (dat - gam0) * cgam )/V
  ## check for log concavity of target density
  lcon <- logcon(b, aw, bw)
  adrej <- lcon
  if(lcon){
    ## if target is log concave, use adaptive rejection sampling
    ## using mode/2, mode, and 2*mode as starting points for hull
    md <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=aw, bvw=bw)$maximum
    try(W <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(md/2, md, md*2),
                      lb=TRUE, xlb=0, a=a, b=b, avw=aw, bvw=bw))
    if(W==0){
      ## error handling - if working properly, this should never happen
      ## but if there is an error, drops out to do the t approximation
      adrej <- FALSE
    }
  }
  if(!adrej){
    ## if target is not log concave, use a t approxmation rejection sampler on log scale
    W <- VWrejiter(a, b, aw, bw)
  }
  return(c(W, lcon, adrej))
  }

## samples V conditional on W, gamma from its IG full conditional
## dat: time series
## gam: scaled disturbances, including gamma_0
## av, bv: hperparameters of V's IG prior
Vgamiter <- function(dat, gam, W, av, bv){
  T <- length(dat)
  theta <- thetagamtrans(gam, W)
  Va <- av + T/2
  Vb <- bv + sum( (dat - theta[-1])^2 )/2
  V <- rinvgamma(1, Va, Vb)
  return(V)
}


## samples V conditional on W, psi from its complicated full conditional
## dat: time series
## psi: scaled errors, including psi_0
## av, bv: hyperparmeters of V's IG prior
Vpsiiter <- function(dat, psi, W, av, bv){
  T <- length(dat)
  ## compute parameters of target density
  psi0 <- psi[1]
  psiLT <- c(0,psi[-1])
  Lpsi <- psiLT[-1] - psiLT[-(T+1)]
  ys <- c(psi0, dat)
  Ly <- ys[-1] - ys[-(T+1)]
  a <- sum(Lpsi^2)/2/W
  b <- sum(Lpsi*Ly)/W
  ## check for log concavity of target density
  lcon <- logcon(b, av, bv)
  adrej <- lcon
  if(lcon){
    ## if target is log concave, use adaptive rejection sampling
    ## using mode/2, mode, and 2*mode as starting points for hull
    md <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=av, bvw=bv)$maximum
    try(V <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(md/2, md, md*2),
                 lb=TRUE, xlb=0, a=a, b=b, avw=av, bvw=bv))
    if(V==0){
      ## error handling - if working properly, this should never happen
      ## but if there is an error, drops out to do the t approximation
      adrej <- FALSE
    }
  }
  if(!adrej){
    ## if target is not log concave, use a t approxmation rejection sampler on log scale
    V <- VWrejiter(a, b, av, bv)
  }
  return(c(V, lcon, adrej))
}


## samples W conditional on V, psi from its IG full conditional
## dat: time series
## psi: scaled errors, including psi_0
## aw, bw: hperparameters of W's IG prior
Wpsiiter <- function(dat, psi, V,  aw, bw){
  T <- length(dat)
  theta <- thetapsitrans(dat, psi, V)
  theta <- c(theta)
  Wa <- aw + T/2
  Wb <- bw + sum( (theta[-1] - theta[-(T+1)])^2 )/2
  W <- rinvgamma(1, Wa, Wb)
  return(W)
}


