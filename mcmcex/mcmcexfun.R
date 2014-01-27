## A set of functions for simulated from and fitting local level models
library(dlm)
library(coda)
library(MCMCpack)
library(ars)
library(plyr)



chains.diag <- function(ns, samplers, simdata, av, aw, parallel, mix=FALSE){
  n <- ns$n[1]
  G.D.full <- ddply(samplers, .(sampler, iter), fullsam.diag, simdata=simdata,
                    n=n, av=av, aw=aw, parallel=parallel, .parallel=parallel, mix=mix)
  return(G.D.full)
}


fullsam.diag <- function(samplers, simdata, n, av, aw, parallel=parallel, mix=FALSE){
  sampler <- samplers$sampler[1]
  sam <- ddply(simdata, .(V.T, W.T, T.T, V.S, W.S, ch), samwrapstart,
               n=n, av=av, aw=aw, samp=sampler)
  samdiag <- ddply(sam, .(V.T, W.T, T.T), sam.diag, .parallel=parallel,
                   parallel=parallel, mix=mix)
  return(samdiag)
}

sam.diag <- function(sam, parallel, mix=FALSE){
  T <- sam$T.T[1]
  namid <- grep("(V.T|W.T|T.T|time)", colnames(sam))
  sam.par <- sam[,-namid]
  sam.par <- sam.par[,1:(T+1+2+2+1)]
  if(!mix){
    sam.list <- dlply(sam.par, .(ch, V.S, W.S), function(x){
      namid <- grep("(V.S|W.S|ch)", colnames(x))
      x.par <- x[,-namid]
      return(mcmc(x.par))
    }, .parallel=parallel)
    GD <- gelman.diag(mcmc.list(sam.list))
    G.D <- data.frame(G.D.M=GD[[2]], G.D.V=GD[[1]][1,1], G.D.W=GD[[1]][2,1])
  }
  else{
    ns <- seq(from=100, to=dim(sam.par)[1]/5, by= 100)
    G.D <- NULL
    for(n in ns){
      sam.list <- dlply(sam.par, .(ch, V.S, W.S), function(x,n){
        namid <- grep("(V.S|W.S|ch)", colnames(x))
        x.par <- x[1:n,-namid]
        return(mcmc(x.par))
      }, .parallel=parallel, n=n)
      GD <- gelman.diag(mcmc.list(sam.list))
      G.D <- rbind(G.D, data.frame(n2=n, G.D.M=GD[[2]], G.D.V=GD[[1]][1,1], G.D.W=GD[[1]][2,1]))
    }
  }
  return(G.D)
}


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
  colnam <- grep("(V.T|W.T|T.T|V|W|time|logconW|adrejW|logconV|adrejV|kernel)$",
                 colnames(sam))
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
  thetas <- theta0s[,-1]
  theta0 <- theta0s[,1]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  VWcor <- cor(V,W)
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
  time <- sam$time[1]/length(sam$time)*1000 #time per 1000 iterations
  logconV <- mean(sam$logconV, na.rm=TRUE)
  adrejV <- mean(sam$adrejV, na.rm=TRUE)
  logconW <- mean(sam$logconW, na.rm=TRUE)
  adrejW <- mean(sam$adrejW, na.rm=TRUE)
  statkern <- mean(sam$kernel=="state", na.rm=TRUE)
  distkern <- mean(sam$kernel=="dist", na.rm=TRUE)
  errorkern <- 1 - (statkern + distkern)
  init <- data.frame(time=time, logconV=logconV, adrejV=adrejV,
                     logconW=logconW, adrejW=adrejW, statkern=statkern,
                     distkern=distkern, errorkern=errorkern)
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

## Wrapper for quickly simulating from all samplers w/ multiple
## chains at diff starting values
samwrapstart <- function(par, n, av, aw, samp){
  dat <- par$y[order(par$t)]
  T <- length(dat)
  start <- c(par$V.S[1]*par$V.T[1], par$W.S[1]*par$W.T[1])
  bv <- (av-1)*par$V.T[1]
  bw <- (aw-1)*par$W.T[1]
  if(samp=="state")
    time <- system.time(out <- statesam(n, start, dat, av, aw, bv, bw))
  if(samp=="dist")
    time <- system.time(out <- distsam(n, start, dat, av, aw, bv, bw))
  if(samp=="error")
    time <- system.time(out <- errorsam(n, start, dat, av, aw, bv, bw))
  if(samp=="sdint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 1))
  if(samp=="seint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 2))
  if(samp=="deint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 3))
  if(samp=="triint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(TRUE, TRUE), 4))
  if(samp=="sdalt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 1))
  if(samp=="sealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 2))
  if(samp=="dealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 3))
  if(samp=="trialt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(FALSE, FALSE), 4))
  if(samp=="sdkern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 1/2, 0)))
  if(samp=="sekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 0, 1/2)))
  if(samp=="dekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(0, 1/2, 1/2)))
  if(samp=="trikern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw))
  if(samp=="partialcis")
    time <- system.time(out <- partialcissam(n, start, dat, av, aw, bv, bw))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, av, aw, bv, bw))
  outdat <- data.frame(out)
  outdat$time <- time[1]
  return(outdat[,c(T+9, (T+2):(T+8), 1:(T+1))])
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
  if(samp=="state")
    time <- system.time(out <- statesam(n, start, dat, av, aw, bv, bw))
  if(samp=="dist")
    time <- system.time(out <- distsam(n, start, dat, av, aw, bv, bw))
  if(samp=="error")
    time <- system.time(out <- errorsam(n, start, dat, av, aw, bv, bw))
  if(samp=="sdint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 1))
  if(samp=="seint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 2))
  if(samp=="deint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 3))
  if(samp=="triint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(TRUE, TRUE), 4))
  if(samp=="sdalt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 1))
  if(samp=="sealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 2))
  if(samp=="dealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 3))
  if(samp=="trialt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(FALSE, FALSE), 4))
  if(samp=="sdkern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 1/2, 0)))
  if(samp=="sekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 0, 1/2)))
  if(samp=="dekern")
   time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(0, 1/2, 1/2)))
  if(samp=="trikern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw))
  if(samp=="partialcis")
    time <- system.time(out <- partialcissam(n, start, dat, av, aw, bv, bw))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, av, aw, bv, bw))
  if(samp=="stnstate")
    time <- system.time(out <- stnstatesam(n, start, dat, av, aw, bv, bw))
  outdat <- data.frame(out)
  outdat$time <- time[1]
  return(outdat[,c(T+9, (T+2):(T+8), 1:(T+1))])
}


## A wrapper for quickly simulating from each of the interweaving/
## alternating samplers.
samwrapper <- function(n, start, dat, av, aw, bv, bw, inter, samp){
  if(samp==1){
    out <- statedistinter(n, start, dat, av, aw, bv, bw, inter)
  }
  if(samp==2){
    out <- stateerrorinter(n, start, dat, av, aw, bv, bw, inter)
  }
  if(samp==3){
    out <- disterrorinter(n, start, dat, av, aw, bv, bw, inter)
  }
  if(samp==4){
    out <- tripleinter(n, start, dat, av, aw, bv, bw, c(inter, inter))
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

## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""), "V", "W", "logconV",
                     "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    out[i,] <- c(theta,V,W,c(NA,NA,NA,NA), NA)
  }
  return(out)
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W","logconV",
                     "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(theta,V,W,rejsV,c(NA,NA), NA)
  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W",
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    gam <- gamtrans(theta, W)
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(theta,V,W,c(NA,NA),rejsW, NA)
  }
  return(out)
}

## state + dist interveaving or alternating sampler
statedistinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    ## state sampler
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled disturbance sampler
    if(!inter){
      theta <- thetaiter(dat, V, W)
    }
    gam <- gamtrans(theta, W)
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(theta,V,W,c(NA,NA),rejsW, NA)
  }
  return(out)
}

## state + error interveaving/alternating sampler
stateerrorinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    ## state sampler
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled error sampler
    if(!inter){
      theta <- thetaiter(dat, V, W)
    }
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(theta,V,W,rejsV,c(NA,NA), NA)
  }
  return(out)
}

## dist + error interveaving/alternating sampler
disterrorinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    ## scaled disturbance sampler
    theta <- thetaiter(dat, V, W)
    gam <- gamtrans(theta, W)
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## scaled error sampler
    if(!inter){
      theta <- thetaiter(dat, V, W)
    }
    else{
      theta <- thetagamtrans(gam, W)
    }
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(theta,V,W, rejsV,rejsW, NA)
  }
    return(out)
}

## state + dist + error interweaving/alternating sampler
tripleinter <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, inter=c(TRUE, TRUE)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W",
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    ## state sampler
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    ## scaled disturbance sampler
    if(!inter[1]){
      theta <- thetaiter(dat, V, W)
    }
    gam <- gamtrans(theta, W)
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    ## scaled error sampler
    if(!inter[2]){
      theta <- thetaiter(dat, V, W)
    }
    else{
      theta <- thetagamtrans(gam, W)
    }
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    theta <- thetapsitrans(dat, psi, V)
    out[i,] <- c(theta,V,W,rejsV,rejsW, NA)
  }
    return(out)
}

## returns TRUE if target density of log-concave, FALSE otherwise
logcon <- function(b, avw, bvw, eps=.1){
  sb <- 1*(b>0) - 1*(b<0)  ##sign(b)
  RHS <- (avw + 1)^3*(1 - 2/3*sb)*32/9/bvw + eps
  ## + eps to make sure ARS algorithm doesn't fail on
  ## near non-log-concave cases
  LHS <- b^2
  out <- (LHS > RHS)
  return(out)
}
## compute M in logtarget - logprop <= M
propM <- function(df, a, b, avw, bvw, mn, propvar){
  M <- optimize(logpirej, c(0,10^10), maximum=TRUE, a, b, avw, bvw, mn, propvar, df)
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
logpiVW <- function(VW, a, b, avw, bvw){
  out <- -a*VW + b*sqrt(VW) - (avw + 1)*log(VW) - bvw/VW
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprime <- function(VW, a, b, avw, bvw){
  out <- -a + b/2/sqrt(VW) - (avw + 1)/VW + bvw/(VW^2)
  return(out)
}

## difference between log conditional posterior of W (V) and the proposal density
logpirej <- function(VW, a, b, avw, bvw, mn, propvar, df){
  out <- logpiVW(VW, a, b, avw, bvw) - dtprop(VW, mn, propvar, df)
  return(out)
}



## samples theta conditional on V and W
thetaiter <- function(dat, V, W){
  mod <- dlmModPoly(order=1, dV=V, dW=W)
  filt <- dlmFilter(dat, mod)
  theta <- dlmBSample(filt)
  return(theta)
}

## transforms theta to gamma
gamtrans <- function(theta, W){
  T <- length(theta) - 1
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  A <- B/sqrt(W)
  A[1,1] <- 1
  gam <- A%*%matrix(theta, ncol=1)
  return(gam)
}

## transforms theta to psi
psitrans <- function(dat, theta, V){
  T <- length(dat)
  A <- diag(-1/sqrt(V), T+1 )
  A[1,1] <- 1
  ytild <- c(0, dat/sqrt(V))
  psi <- ytild + A%*%matrix(theta, ncol=1)
  return(psi)
}

## transforms gamma to theta
thetagamtrans <- function(gam, W){
  T <- length(gam) - 1
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  A <- B/sqrt(W)
  A[1,1] <- 1
  theta <- solve(A)%*%matrix(gam, ncol=1)
  return(theta)
}

## transforms psi to theta
thetapsitrans <- function(dat, psi, V){
  T <- length(dat)
  A <- diag(-1/sqrt(V), T+1 )
  A[1,1] <- 1
  ytild <- c(0, dat/sqrt(V))
  theta <- solve(A)%*%(matrix(psi, ncol=1) - matrix(ytild, ncol=1))
  return(theta)
}

## samples V,W conditional on theta
VWthetaiter <- function(dat, theta, av, aw, bv, bw){
  T <- length(dat)
  V <- rinvgamma(1, av + T/2, bv + sum((dat-theta[-1])^2)/2)
  W <- rinvgamma(1, aw + T/2, bw + sum((theta[-1]-theta[-(T+1)])^2)/2)
  return(c(V,W))
}
## samples W conditional on V,gamma
Wgamiter <- function(dat, gam, V, aw, bw){
  T <- length(dat)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2/V
  b <- sum( (dat - gam0) * cgam )/V
  mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=aw, bvw=bw)$maximum
  adrej <- logcon(b, aw, bw)
  adrej1 <- adrej
  if(adrej){
    try(W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2),
                      lb=TRUE, xlb=0, a=a, b=b, avw=aw, bvw=bw))
    if(W==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    propvar <- - 1 /( (aw + 1)*mn^(-2) - b*mn^(-3/2)/4 - 2*bw*mn^(-3) )
    d <- optimize(propM, c(1,10^10), maximum=FALSE,  a=a, b=b, avw=aw, bvw=bw,
                  mn=mn, propvar=propvar)
    M <- d$objective
    df <- d$minimum
    rej <- TRUE
    rejit <- 1
    while(rej){
      prop <- rtprop(1, mn, propvar, df)
      if(prop>0){
        R <- logpirej(prop, a, b, aw, bw, mn, propvar, df) - M
        u <- runif(1,0,1)
        if(log(u)<R){
          W <- prop
          rej <- FALSE
        }
      }
      rejit <- rejit + 1
      if(rejit%%1000000==0){
        print(paste("Rejection sampler appears to be stuck for W|V, ",
                    rejit, " iterations", sep=""))
        print(paste("Last value of V in the chain: ", V, sep=""))
        print(paste("Values of a and b:", a, b, sep="  "))
        print(paste("T = ", T, sep=""))
        print(paste("Prior values of alphaV and betaV:", aw, bw, sep="  "))
        print(paste("Attempted df: ", df, sep=""))
        print(paste("Attempted mean and variance:", mn, propvar, sep="  "))
      }
    }
  }
  return(c(W, adrej1, adrej))
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
  mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, avw=av, bvw=bv)$maximum
  adrej <- logcon(b, av, bv, eps=0)
  adrej1 <- adrej
  if(adrej){
    try(V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2),
                 lb=TRUE, xlb=0, a=a, b=b, avw=av, bvw=bv))
    if(V==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    propvar <- - 1 /( (av + 1)*mn^(-2) - b*mn^(-3/2)/4 - 2*bv*mn^(-3) )
    d <- optimize(propM, c(1,10^10), maximum=FALSE,  a=a, b=b, avw=av, bvw=bv,
                  mn=mn, propvar=propvar)
    M <- d$objective
    df <- d$minimum
    rej <- TRUE
    rejit <- 1
    while(rej){
      prop <- rtprop(1, mn, propvar, df)
      if(prop>0){
        R <- logpirej(prop, a, b, av, bv, mn, propvar, df) - M
        u <- runif(1,0,1)
        if(log(u)<R){
          V <- prop
          rej <- FALSE
        }
      }
      rejit <- rejit + 1
      if(rejit%%1000000==0){
        print(paste("Rejection sampler appears to be stuck for V|W, ",
              rejit, " iterations.", sep=""))
        print(paste("T = ", T, sep=""))
        print(paste("Last value of W in the chain: ", W, sep=""))
        print(paste("Values of a and b:", a, b, sep="  "))
        print(paste("Prior values of alphaV and betaV:", av, bv, sep="  "))
        print(paste("Attempted df: ", df, sep=""))
        print(paste("Attempted mean and variance:", mn, propvar, sep="  "))
      }
    }
  }
  return(c(V, adrej1, adrej))
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

randkerniter <- function(dat, V, W, theta, av, aw, bv, bw, probs=c(1/3, 1/3, 1/3)){
  kernels <- c("state", "dist", "error")
  kernel <- sample(kernels, 1, prob=probs)
  if(kernel=="state"){
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    rejsW <- c(NA,NA)
    rejsV <- c(NA,NA)
  }
  if(kernel=="error"){
    theta <- thetaiter(dat, V, W)
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    rejsV <- Vout[2:3]
    W <- Wpsiiter(dat, psi, V,  aw, bw)
    rejsW <- c(NA,NA)
    theta <- thetapsitrans(dat, psi, V)
  }
  if(kernel=="dist"){
    theta <- thetaiter(dat, V, W)
    gam <- gamtrans(theta, W)
    V <- Vgamiter(dat, gam, W, av, bv)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    rejsW <- Wout[2:3]
    rejsV <- c(NA,NA)
    theta <- thetagamtrans(gam, W)
  }
  out <- list(theta=theta, V=V, W=W, rejsV=rejsV, rejsW=rejsW, kernel=kernel)
  return(out)
}

randkernsam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0, probs=c(1/3, 1/3, 1/3)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- data.frame(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                                      "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    iter <- randkerniter(dat, V, W, theta, av, aw, bv, bw, probs)
    theta <- iter$theta
    V <- iter$V
    W <- iter$W
    kernel <- iter$kernel
    rejsV <- iter$rejsV
    rejsW <- iter$rejsW
    out[i,1:(T+7)] <- c(theta,V,W, rejsV, rejsW)
    out[i,T+8] <- kernel
  }
  return(out)
}

partialcissam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- data.frame(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                                      "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    VWiter <- VWthetaiter(dat, theta, av, aw, bv, bw)
    V <- VWiter[1]
    W <- VWiter[2]
    gam <- gamtrans(theta, W)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(theta,V,W, rep(NA,5))
  }
  return(out)
}


fullcissam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- data.frame(matrix(0, nrow=n, ncol=T+8))
  colnames(out) <- colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", 
                                      "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    theta <- thetaiter(dat, V, W)
    V <- rinvgamma(1, av + T/2, bv + sum((dat-theta[-1])^2)/2)
    psi <- psitrans(dat, theta, V)
    Vout <- Vpsiiter(dat, psi, W, av, bv)
    V <- Vout[1]
    theta <- thetapsitrans(dat, psi, V)
    W <- rinvgamma(1, aw + T/2, bw + sum((theta[-1]-theta[-(T+1)])^2)/2)
    gam <- gamtrans(theta, W)
    Wout <- Wgamiter(dat, gam, V, aw, bw)
    W <- Wout[1]
    theta <- thetagamtrans(gam, W)
    out[i,] <- c(theta,V,W, rep(NA,5))
  }
  return(out)
}
