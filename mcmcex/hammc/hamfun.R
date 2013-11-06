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
  Wout <- HMC(lpr, gradlpr, eps, L, lold, par)
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
  b <- sum(Lpsi*Ly)/2
  c <- bv + sum( Ly^2 )/2
  par <- c(a, b, c, av, aw, bv, T)
  lold <- log(Vold)
  Vout <- HMC(lpr, gradlpr, eps, L, lold, par)
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
