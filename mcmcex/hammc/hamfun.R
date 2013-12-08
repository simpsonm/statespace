# SIMPLE IMPLEMENTATION OF HAMILTONIAN MONTE CARLO.
#
# Radford M. Neal, 2010.
#
# This program appears in Figure 2 of "MCMC using Hamiltonian dynamics",
# to appear in the Handbook of Markov Chain Monte Carlo.
#
# The arguments to the HMC function are as follows:
#
#   U          A function to evaluate minus the log of the density of the
#              distribution to be sampled, plus any constant - ie, the
#              "potential energy".
#
#   grad_U     A function to evaluate the gradient of U.
#
#   epsilon    The stepsize to use for the leapfrog steps.
#
#   L          The number of leapfrog steps to do to propose a new state.
#
#   current_q  The current state (position variables only).
#
# Momentum variables are sampled from independent standard normal
# distributions within this function.  The value return is the vector
# of new position variables (equal to current_q if the endpoint of the
# trajectory was rejected).
#
# This function was written for illustrative purposes.  More elaborate
# implementations of this basic HMC method and various variants of HMC
# are available from my web page, http://www.cs.utoronto.ca/~radford/

#### Note: I modified the function

HMC = function (U, grad_U, epsilon, L, current_q, par)
{
  q = current_q
  p = rnorm(length(q),0,1)  # independent standard normal variates
  current_p = p

  # Make a half step for momentum at the beginning

  p = p - epsilon * grad_U(q, par) / 2

  # Alternate full steps for position and momentum

  for (i in 1:L)
  {
    # Make a full step for the position

    q = q + epsilon * p

    # Make a full step for the momentum, except at end of trajectory

    if (i!=L) p = p - epsilon * grad_U(q, par)
  }

  # Make a half step for momentum at the end.

  p = p - epsilon * grad_U(q, par) / 2

  # Negate momentum at end of trajectory to make the proposal symmetric

  p = -p

  # Evaluate potential and kinetic energies at start and end of trajectory

  current_U = U(current_q, par)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q, par)
  proposed_K = sum(p^2) / 2

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position

  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K))
  {
    return (c(q,1))  # accept
  }
  else
  {
    return (c(current_q,0))  # reject
  }
}

lpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  out <- (ayou + T/2)*log(a*exp(lvw) + 2*b*exp(lvw/2) + c) + ame*lvw + bet*exp(-lvw)
  
  return(out)
}

lpr2 <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  out <- (ayou + T/2)*( lvw +  log(a + 2*b*exp(-lvw/2) + c*exp(-lvw))) + ame*lvw + bet*exp(-lvw)
  
  return(out)
}




gradlpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  out <- (ayou + T/2) * (a*exp(lvw) + b*exp(lvw/2)) / (a*exp(lvw) + 2*b*exp(lvw/2) + c) +
    ame - bet*exp(-lvw)
  
  return(out)
}

VWgamiter <- function(dat, gam, av, aw, bv, bw, Vold, Wold, eps, L){
  T <- length(dat)

  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2
  b <- -sum( (dat - gam0) * cgam )/2
  c <- bv + sum( (dat - gam0)^2 )/2
  par <- c(a, b, c, aw, av, bw, T)
  lold <- log(Wold)
  epsprime <- runif(1, 0.8*eps, 1.2*eps)
  Wout <- HMC(lpr, gradlpr, epsprime, L, lold, par)
  W <- exp(Wout[1])
  acc <-  Wout[2]
  if(acc==0){
    V <- Vold
  }
  else{
    theta <- thetagamtrans(gam, W)
    Va <- av + T/2
    Vb <- bv + sum( (dat - theta[-1])^2 )/2
    V <- rinvgamma(1, Va, Vb)
  }
  return(c(V,W,acc))
}

VWpsiiter <- function(dat, psi, av, aw, bv, bw, Vold, Wold, eps, L){
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
  epsprime <- runif(1, 0.8*eps, 1.2*eps)
  Vout <- HMC(lpr, gradlpr, epsprime, L, lold, par)
  V <- exp(Vout[1])
  acc <- Vout[2]
  if(acc==0){
    W <- Wold
  }
  else{
    theta <- thetapsitrans(dat, psi, V)
    theta <- c(theta)
    Wa <- aw + T/2
    Wb <- bw + sum( (theta[-1] - theta[-(T+1)])^2 )/2
    W <- rinvgamma(1, Wa, Wb)
  }
  return(c(V,W,acc))
}

distjointsam <- function(n, start, dat, av, aw, bv, bw, eps, L){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+5))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time", "acc")
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, V, W)
    gam <- gamtrans(theta, W)
    VWout <- VWgamiter(dat, gam, av, aw, bv, bw, V, W, eps, L)
    V <- VWout[1]
    W <- VWout[2]
    acc <- VWout[3]
    theta <- thetagamtrans(gam, W)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta, V, W, time, acc)
  }
  return(out)
}

errorjointsam <- function(n, start, dat, av, aw, bv, bw, eps, L){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- mcmc(matrix(0, nrow=n, ncol=T+5))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time", "acc")
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, V, W)
    psi <- psitrans(dat, theta, V)
    VWout <- VWpsiiter(dat, psi, av, aw, bv, bw, V, W, eps, L)
    V <- VWout[1]
    W <- VWout[2]
    acc <- VWout[3]
    theta <- thetapsitrans(dat, psi, V)
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta, V, W, time, acc)
  }
  return(out)
}

hamwrap <- function(par, n, samp){
  dat <- par$y[order(par$t)]
  T <- length(dat)
  V.T <- par$V.T[1]
  W.T <- par$W.T[1]
  av <- par$av[1]
  aw <- par$aw[1]
  bv <- par$bv[1]
  bw <- par$bw[1]
  eps <- 0.01349 * min(V.T/W.T, W.T/V.T)
  L <- 100 * max(V.T/W.T, W.T/V.T)
  start <- c(V.T,W.T)

  if(samp=="error")
    out <- errorjointsam(n, start, dat, av, aw, bv, bw, eps, L)
  if(samp=="dist")
    out <- distjointsam(n, start, dat, av, aw, bv, bw, eps, L)
  
  return(data.frame(out[,c(T+4, T+2, T+3, T+5, 1:(T+1))]))
}

## simulates from a given sampler for each dataset and for multiple
## chains, and returns summary info on the first chain.
hamsim <- function(samplers, simdata, n, burn, parallel){
  sampler <- samplers$sampler[1]
  print(sampler)
  sam <- ddply(simdata, .(V.T, W.T, T.T), hamwrap, .parallel=parallel,
               n=n, samp=sampler)
  samnam <- paste(sampler, "SAM.RData", sep="")
  colnam <- grep("(V.T|W.T|T.T|V|W|time|acc)$",
                 colnames(sam))
  samshort <- sam[,colnam]
  save(samshort, file=samnam)
  rm(samshort)
  out <- ddply(sam, .(V.T, W.T, T.T), hamsummary,
               .parallel=parallel, dat=simdata, burn=burn,
               sampler=sampler)
  rm(sam)
  save(out, file=paste(sampler, "OUT.RData", sep=""))
  print(paste(sampler, "finished", sep=" "))
  return(out)
}

hamfullsim <- function(samplers, simdata, n, burn, parallel){
  out <- ddply(samplers, .(sampler), hamsim, simdata=simdata, n=n,
               burn=burn, parallel=parallel)
  return(out)
}

hamsummary <- function(sam, dat, burn, sampler){
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
  acc <- mean(sam$acc)
  init <- data.frame(time=time, accept=acc)
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

