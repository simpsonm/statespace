lpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  out <- (ayou + T/2)*( lvw +  log(a + 2*b*exp(-lvw/2) + c*exp(-lvw))) + ame*lvw + bet*exp(-lvw)
  return(-out)
}

gradlpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]
  g <- a*exp(lvw) + 2*b*exp(lvw/2) + c
  gp <- a*exp(lvw) + b*exp(lvw/2)
  out <- - (ayou + T/2) * gp / g - ame + bet*exp(-lvw)
  return(out)
}

lvwmetrop <- function(old, par, std2){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]
  gradroot <- uniroot(gradlpr, c(-100, 100), par=par)
  mn <- gradroot$root
  g <- a*exp(mn) + 2*b*exp(mn/2) + c
  gp <- a*exp(mn) + b*exp(mn/2)
  gpp <- a*exp(mn) + b*exp(mn/2)/2
  prec <- abs((ayou + T/2)*(g*gpp - gp^2)/gpp^2 + bet*exp(-mn))
  std <- sqrt(1/prec)*std2
  prop <- rnorm(1, mn, std)
  num1 <- lpr(prop, par)
  denom1 <- lpr(old, par)
  num2 <- dnorm(old, mn, std, log=TRUE)
  denom2 <- dnorm(prop, mn, std, log=TRUE)
  u <- runif(1,0,1)
  out <- old
  if(log(u)< num1 - denom1 + num2 - denom2)
    out <- prop
  return(out)
}


lvwrwmetrop <- function(old, par, std){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  prop <- rnorm(1, old, std)
  num1 <- lpr(prop, par)
  denom1 <- lpr(old, par)
  u <- runif(1,0,1)
  out <- old
  if(log(u)< num1 - denom1)
    out <- prop
  return(out)
}


VWgamiter <- function(dat, gam, av, aw, bv, bw, Vold, Wold, std){
  T <- length(dat)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2
  b <- -sum( (dat - gam0) * cgam )/2
  c <- bv + sum( (dat - gam0)^2 )/2
  par <- c(a, b, c, aw, av, bw, T)
  lold <- log(Wold)
  lnew <- lvwmetrop(lold, par, std)
  W <- exp(lnew)
  if(lold==lnew){
    V <- Vold
  }
  else{
    theta <- thetagamtrans(gam, W)
    Va <- av + T/2
    Vb <- bv + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
  }
  return(c(V,W))
}

VWpsiiter <- function(dat, psi, av, aw, bv, bw, Vold, Wold, std){
  T <- length(dat)
  psi0 <- psi[1]
  psiLT <- c(0,psi[-1])
  Lpsi <- psiLT[-1] - psiLT[-(T+1)]
  ys <- c(psi0, dat)
  Ly <- ys[-1] - ys[-(T+1)]
  a <- sum(Lpsi^2)/2
  b <- -sum(Lpsi*Ly)/2
  c <- bw + sum( Ly^2 )/2
  par <- c(a, b, c, av, aw, bv, T)
  lold <- log(Vold)
  lnew <- lvwmetrop(lold, par, std)
  V <- exp(lnew)
  if(lold==lnew){
    W <- Wold
  }
  else{
    theta <- thetapsitrans(dat, psi, V)
    theta <- c(theta)
    Wa <- aw + T/2
    Wb <- bw + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
  }
  return(c(V,W))
}

distjointsam <- function(n, start, dat, av, aw, bv, bw){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")
  V.T <- bv/(av-1)
  W.T <- bw/(aw-1)
  R <- W.T/V.T
  R.T <- (W.T/V.T)^(min(R^(1/8),1))
  std <- 1/R.T
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, V, W)
    gam <- gamtrans(theta, W)
    VWout <- VWgamiter(dat, gam, av, aw, bv, bw, V, W, std)
    V <- VWout[1]
    W <- VWout[2]
    theta <- thetagamtrans(gam, W)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta, V, W, time)
  }
  return(out)
}

errorjointsam <- function(n, start, dat, av, aw, bv, bw){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+4))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time")
  V.T <- bv/(av-1)
  W.T <- bw/(aw-1)
  R <- W.T/V.T
  R.T <- (W.T/V.T)^(min(R^(1/15),0.8))
  std <- R.T
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, V, W)
    psi <- psitrans(dat, theta, V)
    VWout <- VWpsiiter(dat, psi, av, aw, bv, bw, V, W, std)
    V <- VWout[1]
    W <- VWout[2]
    theta <- thetapsitrans(dat, psi, V)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta, V, W, time)
  }
  return(out)
}

metwrap <- function(par, n, samp){
  dat <- par$y[order(par$t)]
  T <- length(dat)
  V.T <- par$V.T[1]
  W.T <- par$W.T[1]
  av <- par$av[1]
  aw <- par$aw[1]
  bv <- par$bv[1]
  bw <- par$bw[1]
  start <- c(V.T,W.T)

  if(samp=="error")
    out <- errorjointsam(n, start, dat, av, aw, bv, bw)
  if(samp=="dist")
    out <- distjointsam(n, start, dat, av, aw, bv, bw)
  
  return(data.frame(out[,c(T+4, T+2, T+3, 1:(T+1))]))
}

## simulates from a given sampler for each dataset and for multiple
## chains, and returns summary info on the first chain.
metsim <- function(samplers, simdata, n, burn, parallel){
  sampler <- samplers$sampler[1]
  print(sampler)
  sam <- ddply(simdata, .(V.T, W.T, T.T), metwrap, .parallel=parallel,
               n=n, samp=sampler)
  samnam <- paste(sampler, "SAM.RData", sep="")
  colnam <- grep("(V.T|W.T|T.T|V|W|time)$",
                 colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samshort)
  out <- ddply(sam, .(V.T, W.T, T.T), metsummary,
               .parallel=parallel, dat=simdata, burn=burn,
               sampler=sampler)
  rm(sam)
  save(out, file=paste(sampler, "OUT.RData", sep=""))
  print(paste(sampler, "finished", sep=" "))
  return(out)
}

metfullsim <- function(samplers, simdata, n, burn, parallel){
  out <- ddply(samplers, .(sampler), metsim, simdata=simdata, n=n,
               burn=burn, parallel=parallel)
  return(out)
}

metsummary <- function(sam, dat, burn, sampler){
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
  time <- mean(sam$time)
  init <- data.frame(time=time)
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

