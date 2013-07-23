## A set of functions for simulated from and fitting local level models
library(dlm)
library(coda)
library(MCMCpack)
library(ars)
library(plyr)

fullsim <- function(samplers, simdata, n, burn, a1, a2){
  parallel <- require(doMC, quietly=TRUE)
  if(parallel){
    registerDoMC()
  }
  out <- ddply(samplers, .(sams), samsim, .parallel=parallel,
               simdata=simdata, n=n, burn=burn, a1=a1, a2=a2)
  return(out)
}


## simulates from a given sampler for each dataset and for multiple
## chains, and returns summary info on the first chain.
samsim <- function(sampler, simdata, n, burn, a1, a2){
  sampler <- sampler$sams[1]
  print(sampler)
  parallel <- require(doMC, quietly=TRUE)
  sam <- ddply(simdata, .(V.T, W.T, T.T, ch), samwrap, .parallel=parallel,
               n=n, a1=a1, a2=a2, samp=sampler)
  samnam <- paste(sampler, "SAM.RData", sep="")
  colnam <- grep("(V.T|W.T|T.T|ch|V|W|time|theta(0|1|10|100|1000)$)",
                 colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samsamshort)
  out <- ddply(sam[sam$ch==1,], .(V.T, W.T, T.T), samsummary,
               .parallel=parallel, dat=simdata[simdata$ch==1,], burn=burn,
               sampler=sampler)
  if(sampler=="trialt"){
    posteriorcors <- ddply(sam[sam$ch==1,], .(V.T, W.T, T.T), postcor, .parallel=parallel,
                           dat=simdata[simdata$ch==1,], burn=burn)
    save(posteriorcors, file="postcors.RData")
  }

  rm(sam)
  return(out)
}

postcor <- function(sam, dat, burn){
  V.T <- sam$V.T[1]
  W.T <- sam$W.T[1]
  T.T <- sam$T.T[1]
  V <- sam$V[-c(1:burn)]
  W <- sam$W[-c(1:burn)]
  theta0s <- sam[,grep("theta", colnames(sam))]
  theta0s <- theta0s[-c(1:burn),1:(T.T+1)]
  thetas <- theta0s[-c(1:burn),-1]
  theta0 <- sam$theta0[-c(1:burn)]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  gammas <- (theta0s[,-1] - theta0s[,-(T.T+1)])/sqrt(W)
  colnames(gammas) <- paste("gamma", 1:T.T, sep="")
  psis <- (matrix(data, ncol=1) - thetas)/sqrt(V)
  colnames(psis) <- paste("psi", 1:T.T, sep="")
  VWcor <- cor(V,W)
  Vth0cor <- cor(V,theta)
  Wth0cor <- cor(W,theta)
  Vthcors <- cor(V, thetas)
  Wthcors <- cor(W, thetas)
  Vgacors <- cor(V, gammas)
  Wgacors <- cor(W, gammas)
  Vgacors <- cor(V, psis)
  Wgacors <- cor(W, psis)
  Vthmaxcor <- Vthcors[1,which.max(abs(Vthcors))]
  Wthmaxcor <- Wthcors[1,which.max(abs(Wthcors))]
  Vgamaxcor <- Vgacors[1,which.max(abs(Vgacors))]
  Wgamaxcor <- Wgacors[1,which.max(abs(Wgacors))]
  Vpsmaxcor <- Vpscors[1,which.max(abs(Vpscors))]
  Wpsmaxcor <- Wpscors[1,which.max(abs(Wpscors))]
  out <- data.frame(VW=VWcor, Vtheta0=Vth0cor, Wtheta0=Wth0cor, Vtheta1=Vthcors[1,1],
                    Wtheta1=Wthcors[1,1], VthetaT=Vthcors[1,T.T], WthetaT=Wthcors[1,T.T],
                    Vgamma1=Vgacors[1,1], Wgamma1=Wgacors[1,1], VgammaT=Vgacors[1,T.T],
                    WgammaT=Wgacors[1,T.T], Vpsi1=Vpscors[1,1], Wpsi1=Wpscors[1,1],
                    VpsiT=Vpscors[1,T.T], WpsiT=Wpscors[1,T.T], Vtheta=Vthmaxcor,
                    Wtheta=Wthmaxcor, Vgamma=Vgamaxcor, Wgamma=Wgamaxcor, Vpsi=Vpsmaxcor,
                    Wpsi=Wpsmaxcor)
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
  theta0s <- sam[,grep("theta", colnames(sam))]
  theta0s <- theta0s[-c(1:burn),1:(T.T+1)]
  thetas <- theta0s[-c(1:burn),-1]
  theta0 <- sam$theta0[-c(1:burn)]
  data <- dat$y[dat$V.T==V.T & dat$W.T==W.T & dat$T.T==T.T]
  time <- mean(sam$time[-c(1:burn)])
  init <- data.frame(sampler=sampler, time=time)
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
  thetaT.AC <- thetaAC[T.T]
  gamma1.AC <- gammaAC[1]
  gammaT.AC <- gammaAC[T.T]
  psi1.AC <- psiAC[1]
  psiT.AC <- psiAC[T.T]
  theta.AC <- thetaAC[which.max(abs(thetaAC))]
  gamma.AC <- gammaAC[which.max(abs(gammaAC))]
  psi.AC <- psiAC[which.max(abs(psiAC))]
  theta0.ES <- effectiveSize(theta0)
  theta1.ES <- thetaES[1]
  thetaT.ES <- thetaES[T.T]
  gamma1.ES <- gammaES[1]
  gammaT.ES <- gammaES[T.T]
  psi1.ES <- psiES[1]
  psiT.ES <- psiES[T.T]
  theta.ES <- thetaES[which.min(thetaES)]
  gamma.ES <- gammaES[which.min(gammaES)]
  psi.ES <- psiES[which.min(psiES)]
  V.AC <- corfun(V)
  W.AC <- corfun(W)
  V.ES <- effectiveSize(V)
  W.ES <- effectiveSize(W)
  out <- cbind(init, V.AC, W.AC, theta0.AC, theta1.AC, thetaT.AC,
               theta.AC, gamma1.AC, gammaT.AC, gamma.AC, psi1.AC, psiT.AC,
               psi.AC, V.ES, W.ES, theta0.ES, theta1.ES, thetaT.ES,
               theta.ES, gamma1.ES, gammaT.ES, gamma.ES, psi1.ES,
               psiT.ES, psi.ES)
  rownames(out) <- ""
  return(out)
}

## function for finding the first order autocorrelation of a TS
corfun <- function(x){
  acf(x, lag.max=1, plot=FALSE)[[1]][2]
}

## Wrapper for quickly simulating from all samplers
samwrap <- function(par, n, a1, a2, samp){
  dat <- par$y[order(par$t)]
  M <- c(1, 1/100, 100)
  ch <- par$ch[1]
  start <- c(par$V[1], par$W[1])*M[ch]
  b1 <- (a1-1)*par$V[1]
  b2 <- (a2-1)*par$W[1]
  if(samp=="state")
      out <- statesam(n, start, dat, a1, a2, b1, b2)
  if(samp=="dist")
      out <- distsam(n, start, dat, a1, a2, b1, b2)
  if(samp=="error")
      out <- errorsam(n, start, dat, a1, a2, b1, b2)
  if(samp=="sdint")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, TRUE, 1)
  if(samp=="seint")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, TRUE, 2)
  if(samp=="deint")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, TRUE, 3)
  if(samp=="triint")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, c(TRUE, TRUE), 4)
  if(samp=="sdalt")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, FALSE, 1)
  if(samp=="sealt")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, FALSE, 2)
  if(samp=="dealt")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, FALSE, 3)
  if(samp=="trialt")
      out <- samwrapper(n, start, dat, a1, a2, b1, b2, c(FALSE, FALSE), 4)
  return(data.frame(out))
}


## A wrapper for quickly simulating from each of the interweaving/
## alternating samplers.
samwrapper <- function(n, start, dat, a1, a2, b1, b2, inter, samp){
  if(samp==1){
    out <- statedistinter(n, start, dat, a1, a2, b1, b2, inter)
  }
  if(samp==2){
    out <- stateerrorinter(n, start, dat, a1, a2, b1, b2, inter)
  }
  if(samp==3){
    out <- disterrorinter(n, start, dat, a1, a2, b1, b2, inter)
  }
  if(samp==4){
    out <- tripleinter(n, start, dat, a1, a2, b1, b2, c(inter, inter))
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

## state sampler: samples V and W conditional on states
statesam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")


  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-(T+1)])^2)/2)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W, time)
  }
  return(out)
}

## scaled error sampler: samples V and W conditional on the scaled observation
## errors (plus the initial state, theta_0)
errorsam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")

  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    psi <- ytild + A%*%theta
    Wa <- a2 + T/2
    Wb <- b2 + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
    psi0 <- psi[1]
    psiLT <- c(0,psi[-1])
    Lpsi <- psiLT[-1] - psiLT[-(T+1)]
    ys <- c(psi0, dat)
    Ly <- ys[-1] - ys[-(T+1)]
    a <- sum(Lpsi^2)/2/W
    b <- sum(Lpsi*Ly)/W
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a1, b12=b1)$maximum
    if(logcon(b, a1, b1)){
      V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a1, b12=b1)
    }
    else{
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

    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    theta <- solve(A)%*%(psi - ytild)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W,time)

  }
  return(out)
}


## scaled disturbance sampler: samples V and W conditional on the scaled
## system disturbances (plus the initial state, theta_0)
distsam <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B

  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a2, b12=b2)$maximum
    if(logcon(b, a2, b2)){
      W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a2, b12=b2)
    }
    else{
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
    A <- B/sqrt(W)
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W, time)
  }
  return(out)
}

## state + dist interveaving or alternating sampler
statedistinter <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B

  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    ## sample states
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    ## sampler V,W conditional on states
    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-(T+1)])^2)/2)
    ## If not interweaving, sample states again
    if(!inter){
      mod <- dlmModPoly(order=1, dV=V, dW=W)
      filt <- dlmFilter(dat, mod)
      theta <- dlmBSample(filt)
    }
    ## convert states to disturbances and sample (V,W) conditional on disturbances
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a2, b12=b2)$maximum
    if(logcon(b, a2, b2)){
      W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a2, b12=b2)
    }
    else{
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
    ## convert disturbances back to states and save
    A <- B/sqrt(W)
    A[1,1] <- 1
    theta <- solve(A)%*%gam
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W,time)
  }
  return(out)
}

## state + error interveaving/alternating sampler
stateerrorinter <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W","time")
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    ## sample states
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    ## sample (V,W) conditional on states
    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-(T+1)])^2)/2)
    ## if not interweaving, sample states again
    if(!inter){
      mod <- dlmModPoly(order=1, dV=V, dW=W)
      filt <- dlmFilter(dat, mod)
      theta <- dlmBSample(filt)
    }
    ## convert states to errors and sample (V,W) conditional on errors
    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    psi <- ytild + A%*%theta
    Wa <- a2 + T/2
    Wb <- b2 + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
    psi0 <- psi[1]
    psiLT <- c(0,psi[-1])
    Lpsi <- psiLT[-1] - psiLT[-(T+1)]
    ys <- c(psi0, dat)
    Ly <- ys[-1] - ys[-(T+1)]
    a <- sum(Lpsi^2)/2/W
    b <- sum(Lpsi*Ly)/W
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a1, b12=b1)$maximum
    if(logcon(b, a1, b1)){
      V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a1, b12=b1)
    }
    else{
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
    ## convert errors to states and save
    A <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    theta <- solve(A)%*%(psi - ytild)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W,time)
  }
  return(out)
}

## dist + error interveaving/alternating sampler
disterrorinter <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, inter=TRUE){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    ## sample states and convert to disturbances
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    ## sample (V,W) conditional on disturbances
    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a2, b12=b2)$maximum
    if(logcon(b, a2, b2)){
      W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a2, b12=b2)
    }
    else{
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
    ## if not interweaving, sample the states. Otherwise convert disturbances
    ## to states
    if(!inter){
      mod <- dlmModPoly(order=1, dV=V, dW=W)
      filt <- dlmFilter(dat, mod)
      theta <- dlmBSample(filt)
    }
    else{
      A <- B/sqrt(W)
      A[1,1] <- 1
      theta <- solve(A)%*%gam
    }
    ## convert states to errors and sample (V,W) conditional on errors
    Ap <- diag(-1/sqrt(V), T+1 )
    Ap[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    psi <- ytild + Ap%*%theta
    Wa <- a2 + T/2
    Wb <- b2 + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
    psi0 <- psi[1]
    psiLT <- c(0,psi[-1])
    Lpsi <- psiLT[-1] - psiLT[-(T+1)]
    ys <- c(psi0, dat)
    Ly <- ys[-1] - ys[-(T+1)]
    a <- sum(Lpsi^2)/2/W
    b <- sum(Lpsi*Ly)/W
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a1, b12=b1)$maximum
    if(logcon(b, a1, b1)){
      V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a1, b12=b1)
    }
    else{
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
    ## convert errors to states and save
    Ap <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    theta <- solve(Ap)%*%(psi - ytild)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W,time)
  }
    return(out)
}

## state + dist + error interweaving/alternating sampler
tripleinter <- function(n, start, dat, a1=0, a2=0, b1=0, b2=0, inter=c(TRUE, TRUE)){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W","time")
  B <- diag(1,T+1,T+1)
  C <- diag(-1,T,T)
  B <- rbind(rep(0,T+1),cbind(C,rep(0,T))) + B

  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    ## sample states
    mod <- dlmModPoly(order=1, dV=V, dW=W)
    filt <- dlmFilter(dat, mod)
    theta <- dlmBSample(filt)
    ## sample (V,W) conditional on states
    V <- rinvgamma(1, a1 + T/2, b1 + sum((dat-theta[-1])^2)/2)
    W <- rinvgamma(1, a2 + T/2, b2 + sum((theta[-1]-theta[-(T+1)])^2)/2)
    ## If not interweaving, sample the states again
    if(!inter[1]){
      mod <- dlmModPoly(order=1, dV=V, dW=W)
      filt <- dlmFilter(dat, mod)
      theta <- dlmBSample(filt)
    }
    ## Convert states to disturbances
    A <- B/sqrt(W)
    A[1,1] <- 1
    gam <- A%*%theta
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    ## Sample (V,W) conditional on disturbances
    Va <- a1 + T/2
    Vb <- b1 + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
    a <- sum( cgam^2 )/2/V
    b <- sum( (dat - gam0) * cgam ) /V
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a2, b12=b2)$maximum
    if(logcon(b, a2, b2)){
      W <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a2, b12=b2)
    }
    else{
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
    ## If not interweaving, sample the states.
    ## Otherwise, convert disturbances to states
    if(!inter[2]){
      mod <- dlmModPoly(order=1, dV=V, dW=W)
      filt <- dlmFilter(dat, mod)
      theta <- dlmBSample(filt)
    }
    else{
      A <- B/sqrt(W)
      A[1,1] <- 1
      theta <- solve(A)%*%gam
    }
    ## Convert states to errors and sample (V,W) conditional on errors
    Ap <- diag(-1/sqrt(V), T+1 )
    Ap[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    psi <- ytild + Ap%*%theta
    Wa <- a2 + T/2
    Wb <- b2 + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
    psi0 <- psi[1]
    psiLT <- c(0,psi[-1])
    Lpsi <- psiLT[-1] - psiLT[-(T+1)]
    ys <- c(psi0, dat)
    Ly <- ys[-1] - ys[-(T+1)]
    a <- sum(Lpsi^2)/2/W
    b <- sum(Lpsi*Ly)/W
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, a12=a1, b12=b1)$maximum
    if(logcon(b, a1, b1)){
      V <- ars(n=1, logpiVW, logpiVWprime, x=c(mn/2, mn, mn*2), lb=TRUE, xlb=0, a=a, b=b, a12=a1, b12=b1)
    }
    else{
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
    ## Convert errors to states and save
    Ap <- diag(-1/sqrt(V), T+1 )
    A[1,1] <- 1
    ytild <- c(0, dat/sqrt(V))
    theta <- solve(Ap)%*%(psi - ytild)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W,time)
  }
    return(out)
}

## returns TRUE if target density of log-concave, FALSE otherwise
logcon <- function(b, a12, b12, eps=10){
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




