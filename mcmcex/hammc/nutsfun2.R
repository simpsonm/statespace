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

  out <- (ayou + T/2) * (a*exp(lvw) + b*exp(lvw/2)) / (a*exp(lvw) + 2*b*exp(lvw/2) + c) +
    ame - bet*exp(-lvw)
  
  return(-out)
}

f <- function(lvw, par){
  logprob <- lpr(lvw, par)
  grad <- gradlpr(lvw, par)
  out <- list(logprob=logprob, grad=grad)
  return(out)
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

distjointsam <- function(M, Madapt, start, dat, av, aw, bv, bw, delta=0.6){
  T <- length(dat)
  V <- start[1]
  W <- start[2]

  theta <- thetaiter(dat, V, W)
  gam <- gamtrans(theta, W)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  a <- sum( cgam^2 )/2
  b <- -sum( (dat - gam0) * cgam )/2
  c <- bv + sum( (dat - gam0)^2 )/2
  par <- c(a, b, c, aw, av, bw, T)

  samples <- matrix(0, nrow=M, ncol=T+5)
  fout <- f(W, par)
  logp <- fout$logp
  grad <- fout$grad
  colnames(samples) <- c(paste("theta",0:T,sep=""),"V","W", "time", "epsilon")

  ## Choose a reasonable first epsilon by a simple heuristic.
  epsilon <- find_reasonable_epsilon(W, grad, logp, f, par)

  ## Initialize dual averaging algorithm.
  epsilonbar <- 1
  Hbar <- 0

  ## Parameters to the dual averaging algorithm.
  gamma <- 0.05
  t0 <- 10
  kappa <- 0.75
  mu <- log(10*epsilon)
  samples[1,] <- c(theta, V, W, 0, epsilon)

  for(m in 2:(M+Madapt)){
    Vold <- V
    Wold <- W
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, Vold, Wold)
    gam <- gamtrans(theta, Wold)
    cgam <- cumsum(gam[-1])
    gam0 <- gam[1]
    a <- sum( cgam^2 )/2
    b <- -sum( (dat - gam0) * cgam )/2
    c <- bv + sum( (dat - gam0)^2 )/2
    par <- c(a, b, c, aw, av, bw, T)
    VW0 <- log(Wold)
        
    ## Resample momenta.
    r0 <- rnorm(1, 0, 1);
    ## Joint log-probability of theta and momenta r.
    joint <- logp - 0.5 * sum(r0^2)    
    ## Resample u ~ uniform([0, exp(joint)]).
    ## Equivalent to (log(u) - joint) ~ exponential(1).
    logu <- joint - rexp(1)
    ## Initialize tree.
    VWminus <- VW0
    VWplus <- VW0
    rminus <- r0
    rplus <- r0
    gradminus <- grad
    gradplus <- grad
    ## Initial height j = 0.
    j <- 0
    ## If all else fails, the next sample is the previous sample.
    VW <- VW0
    ##% Initially the only valid point is the initial point.
    n <- 1
    
    ##% Main loop---keep going until the criterion s == 0.
    s <- 1
    while(s == 1){
      ## % Choose a direction. -1=backwards, 1=forwards.
      v <- 2*(runif(1,0,1) < 0.5)-1;
      ##% Double the size of the tree.
      if(v == -1){
        BTout <- build_tree(VWminus, rminus, gradminus, logu, v, j, epsilon, f, par, joint)
        VWminus <- BTout[[1]]
        rminus <- BTout[[2]]
        gradminus <- BTout[[3]]
        VWprime <- BTout[[7]]
        gradprime <- BTout[[8]]
        logpprime <- BTout[[9]]
        nprime <- BTout[[10]]
        sprime <- BTout[[11]]
        alpha <- BTout[[12]]
        nalpha <- BTout[[13]]
      }
      else{
        BTout <- build_tree(VWplus, rplus, gradplus, logu, v, j, epsilon, f, par, joint);
        VWplus <- BTout[[4]]
        rplus <- BTout[[5]]
        gradplus <- BTout[[6]]
        VWprime <- BTout[[7]]
        gradprime <- BTout[[8]]
        logpprime <- BTout[[9]]
        nprime <- BTout[[10]]
        sprime <- BTout[[11]]
        alpha <- BTout[[12]]
        nalpha <- BTout[[13]]
      }
      ## % Use Metropolis-Hastings to decide whether or not to move to a
      ## % point from the half-tree we just generated.
      if((sprime == 1) & (runif(1,0,1) < nprime/n)){
        VW <- VWprime
        logp <- logpprime
        grad <- gradprime
      }
      ##% Update number of valid points we've seen.
      n <- n + nprime;
      ## % Decide if it's time to stop.
      s <- sprime && stop_criterion(VWminus, VWplus, rminus, rplus)
      ## % Increment depth.
      j <- j + 1
      print(paste("j = ", j, sep=""))
      flush.console()
    }
    
    ##% Do adaptation of epsilon if we're still doing burn-in.'
    eta <- 1 / (m - 1 + t0)
    Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
    print(paste("Hbar = ", Hbar))
    if(m <= Madapt){
      epsilon <- exp(mu - sqrt(m-1)/gamma * Hbar)
      eta <- (m-1)^-kappa
      print(paste("eta = ", eta))
      epsilonbar <- exp((1 - eta) * log(epsilonbar) + eta * log(epsilon))
      print(paste("epsilonbar = ", epsilonbar))
    }
    else{
      epsilon <- epsilonbar
    }
    print(paste("m = ", m, sep=""))
    print(paste("epsilon = ", epsilon, sep=""))

    
    W <- exp(VW)
    theta <- thetagamtrans(gam, W)
    if(W==Wold){
      V <- Vold
    }
    else{
      Va <- av + T/2
      Vb <- bv + sum( (dat - theta[-1])^2 )/2
      V <- rinvgamma(1, Va, Vb)
    }
    print(c(V,W))
    
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    if(m > Madapt)
      samples[m-Madapt,] <- c(theta, V, W, time, epsilon)
  }
  return(samples)
}



nutsiter <- function(theta0, logp, grad, n, s, Hbar, epsilon, epsilonbar, m, Madapt, f, par, delta=0.6){



  ## Resample momenta.
  r0 <- rnorm(D, 0, 1);
  ## Joint log-probability of theta and momenta r.
  joint <- logp - 0.5 * sum(r0^2)    
  ## Resample u ~ uniform([0, exp(joint)]).
  ## Equivalent to (log(u) - joint) ~ exponential(1).
  logu <- joint - rexp(1)
  ## Initialize tree.
  thetaminus <- theta0
  thetaplus <- theta0
  rminus <- r0
  rplus <- r0
  gradminus <- grad
  gradplus <- grad
  ## Initial height j = 0.
  j <- 0
  ## If all else fails, the next sample is the previous sample.
  theta <- theta0
  ##% Initially the only valid point is the initial point.
  n <- 1
  
  ##% Main loop---keep going until the criterion s == 0.
  s <- 1
  while(s == 1){
    ## % Choose a direction. -1=backwards, 1=forwards.
    v <- 2*(runif(1,0,1) < 0.5)-1;
    ##% Double the size of the tree.
    if(v == -1){
      BTout <- build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f, par, joint)
      thetaminus <- BTout[[1]]
      rminus <- BTout[[2]]
      gradminus <- BTout[[3]]
      thetaprime <- BTout[[7]]
      gradprime <- BTout[[8]]
      logpprime <- BTout[[9]]
      nprime <- BTout[[10]]
      sprime <- BTout[[11]]
      alpha <- BTout[[12]]
      nalpha <- BTout[[13]]
    }
    else{
      BTout <- build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f, par, joint);
      thetaplus <- BTout[[4]]
      rplus <- BTout[[5]]
      gradplus <- BTout[[6]]
      thetaprime <- BTout[[7]]
      gradprime <- BTout[[8]]
      logpprime <- BTout[[9]]
      nprime <- BTout[[10]]
      sprime <- BTout[[11]]
      alpha <- BTout[[12]]
      nalpha <- BTout[[13]]
    }
    ## % Use Metropolis-Hastings to decide whether or not to move to a
    ## % point from the half-tree we just generated.
    if((sprime == 1) & (runif(1,0,1) < nprime/n)){
      theta <- thetaprime
      logp <- logpprime
      grad <- gradprime
    }
    ##% Update number of valid points we've seen.
    n <- n + nprime;
    ## % Decide if it's time to stop.
    s <- sprime && stop_criterion(thetaminus, thetaplus, rminus, rplus)
    ## % Increment depth.
    j <- j + 1
    print(paste("j = ", j, sep=""))
    flush.console()
  }
    
  ##% Do adaptation of epsilon if we're still doing burn-in.'
  eta <- 1 / (m - 1 + t0)
  Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
  if(m <= Madapt){
    epsilon <- exp(mu - sqrt(m-1)/gamma * Hbar)
    eta <- (m-1)^-kappa
    epsilonbar <- exp((1 - eta) * log(epsilonbar) + eta * log(epsilon))
  }
  else{
    epsilon <- epsilonbar
  }
  print(paste("m = ", m, sep=""))
  print(paste("epsilon = ", epsilon, sep=""))
  return(list(theta, epsilon, logp, grad, n, s, Hbar, epsilonbar))
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







leapfrog <- function(theta, r, grad, epsilon, f, par){
  rprime <- r + 0.5 * epsilon * grad
  thetaprime <- theta + epsilon * rprime
  fout <- f(thetaprime, par)
  logpprime <- fout$logp
  gradprime <- fout$grad
  rprime <- rprime + 0.5 * epsilon * gradprime
  return(list(thetaprime=thetaprime, rprime=rprime, gradprime=gradprime, logpprime=logpprime))
}


stop_criterion <- function(thetaminus, thetaplus, rminus, rplus){
  thetavec <- thetaplus - thetaminus
  criterion <- (sum(thetavec * rminus) >= 0) & (sum(thetavec * rplus) >= 0)
  return(criterion)
}


build_tree <- function(theta, r, grad, logu, v, j, epsilon, f, par, joint0){
  if (j == 0){
    ##% Base case: Take a single leapfrog step in the direction v.
    lfout <- leapfrog(theta, r, grad, v*epsilon, f, par)
    thetaprime <- lfout$thetaprime
    rprime <- lfout$rprime
    gradprime <- lfout$gradprime
    logpprime <- lfout$logpprime
    joint  <- logpprime - 0.5 * sum(rprime^2) 
    ##% Is the new point in the slice?
    nprime <- logu < joint
    ## Is the simulation wildly inaccurate?
    sprime <- logu - 1000 < joint
    ##% Set the return values---minus=plus for all things here, since the
    ##% "tree" is of depth 0.
    thetaminus <- thetaprime
    thetaplus <- thetaprime
    rminus <- rprime
    rplus <- rprime
    gradminus <- gradprime
    gradplus <- gradprime
    ##% Compute the acceptance probability.
    alphaprime <- min(1, exp(logpprime - 0.5 * sum(rprime^2) - joint0))
    nalphaprime <- 1
  }
  else{
    ## % Recursion: Implicitly build the height j-1 left and right subtrees.
    BTout <- build_tree(theta, r, grad, logu, v, j-1, epsilon, f, par, joint0)
    thetaminus <- BTout[[1]]
    rminus <- BTout[[2]]
    gradminus <- BTout[[3]]
    thetaplus <- BTout[[4]]
    rplus <- BTout[[5]]
    gradplus <- BTout[[6]]
    thetaprime <- BTout[[7]]
    gradprime <- BTout[[8]]
    logpprime <- BTout[[9]]
    nprime <- BTout[[10]]
    sprime <- BTout[[11]]
    alphaprime <- BTout[[12]]
    nalphaprime <- BTout[[13]]
    ##% No need to keep going if the stopping criteria were met in the first
    ##% subtree.
    if(sprime == 1){
      if(v == -1){
        BTout <- build_tree(thetaminus, rminus, gradminus, logu, v, j-1, epsilon, f, par, joint0)
        thetaminus <- BTout[[1]]
        rminus <- BTout[[2]]
        gradminus <- BTout[[3]]
        thetaprime2 <- BTout[[7]]
        gradprime2 <- BTout[[8]]
        logpprime2 <- BTout[[9]]
        nprime2 <- BTout[[10]]
        sprime2 <- BTout[[11]]
        alphaprime2 <- BTout[[12]]
        nalphaprime2 <- BTout[[13]]
      }
      else{
        BTout <- build_tree(thetaplus, rplus, gradplus, logu, v, j-1, epsilon, f, par, joint0);
        thetaplus <- BTout[[4]]
        rplus <- BTout[[5]]
        gradplus <- BTout[[6]]
        thetaprime2 <- BTout[[7]]
        gradprime2 <- BTout[[8]]
        logpprime2 <- BTout[[9]]
        nprime2 <- BTout[[10]]
        sprime2 <- BTout[[11]]
        alphaprime2 <- BTout[[12]]
        nalphaprime2 <- BTout[[13]]
      }
      ## % Choose which subtree to propagate a sample up from.
      ## (this npr bit needed because R can't compare NaN to anything)
      npr <- nprime2 / (nprime + nprime2 + (nprime + nprime2 == 0))
      if (runif(1,0,1) < npr){
        thetaprime <- thetaprime2
        gradprime <- gradprime2
        logpprime <- logpprime2
      }
      ## Update the number of valid points.
      nprime <- nprime + nprime2
      ##% Update the stopping criterion.
      sprime <- sprime & sprime2 & stop_criterion(thetaminus, thetaplus, rminus, rplus)
      ##% Update the acceptance probability statistics.
      alphaprime <- alphaprime + alphaprime2
      nalphaprime <- nalphaprime + nalphaprime2
    }
  }
  return(list(thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime))
}

find_reasonable_epsilon <- function(theta0, grad0, logp0, f, par){
  epsilon <- 1
  r0 <- rnorm(length(theta0), 0, 1)
  ##% Figure out what direction we should be moving epsilon.
  lfout <- leapfrog(theta0, r0, grad0, epsilon, f, par)
  rprime <- lfout[[2]]
  logpprime <- lfout[[4]]
  acceptprob <- exp(logpprime - logp0 - 0.5 * sum(rprime^2) - sum(r0^2))
  a = 2 * (acceptprob > 0.5) - 1
  ##% Keep moving epsilon in that direction until acceptprob crosses 0.5.
  while(acceptprob^a > 2^(-a)){
    epsilon <- c(epsilon * 2^a)
    lfout <- leapfrog(theta0, r0, grad0, epsilon, f, par)
    rprime <- lfout[[2]]
    logpprime <- lfout[[4]]
    acceptprob <- exp(logpprime - logp0 - 0.5 * sum(rprime^2) - sum(r0^2))
    ## if epsilon is already small, leave it alone
    if(epsilon < 10^(-4)){
      acceptprob <- a
      epsilon <- c(epsilon * 2^(-a))
    }
  }
  return(epsilon)
}







nuts_da <- function(f, M, Madapt, theta0, par, delta=0.6){
  ## assert(size(theta0, 1) == 1);
  D <- length(theta0)
  samples <- matrix(0, nrow=M+Madapt, ncol=D)
  fout <- f(theta0, par)
  logp <- fout$logp
  grad <- fout$grad
  samples[1,] <- theta0
  
  ## Choose a reasonable first epsilon by a simple heuristic.
  epsilon <- find_reasonable_epsilon(theta0, grad, logp, f, par)
  print("epsilon initialized")

  ## Initialize dual averaging algorithm.
  epsilonbar <- 1
  Hbar <- 0

  ## Parameters to the dual averaging algorithm.
  gamma <- 0.05
  t0 <- 10
  kappa <- 0.75
  mu <- log(10*epsilon)


  for(m in  2:(M+Madapt)){
    theta0 <- samples[m-1,]

    ## Resample momenta.
    r0 <- rnorm(D, 0, 1);
    ## Joint log-probability of theta and momenta r.
    joint <- logp - 0.5 * sum(r0^2)    
    ## Resample u ~ uniform([0, exp(joint)]).
    ## Equivalent to (log(u) - joint) ~ exponential(1).
    logu <- joint - rexp(1)
    ## Initialize tree.
    thetaminus <- theta0
    thetaplus <- theta0
    rminus <- r0
    rplus <- r0
    gradminus <- grad
    gradplus <- grad
    ## Initial height j = 0.
    j <- 0
    ## If all else fails, the next sample is the previous sample.
    theta <- theta0
    ##% Initially the only valid point is the initial point.
    n <- 1
    
    ##% Main loop---keep going until the criterion s == 0.
    s <- 1
    while(s == 1){
      ## % Choose a direction. -1=backwards, 1=forwards.
      v <- 2*(runif(1,0,1) < 0.5)-1;
      ##% Double the size of the tree.
      if(v == -1){
        BTout <- build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f, par, joint)
        thetaminus <- BTout[[1]]
        rminus <- BTout[[2]]
        gradminus <- BTout[[3]]
        thetaprime <- BTout[[7]]
        gradprime <- BTout[[8]]
        logpprime <- BTout[[9]]
        nprime <- BTout[[10]]
        sprime <- BTout[[11]]
        alpha <- BTout[[12]]
        nalpha <- BTout[[13]]
      }
      else{
        BTout <- build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f, par, joint);
        thetaplus <- BTout[[4]]
        rplus <- BTout[[5]]
        gradplus <- BTout[[6]]
        thetaprime <- BTout[[7]]
        gradprime <- BTout[[8]]
        logpprime <- BTout[[9]]
        nprime <- BTout[[10]]
        sprime <- BTout[[11]]
        alpha <- BTout[[12]]
        nalpha <- BTout[[13]]
      }
      ## % Use Metropolis-Hastings to decide whether or not to move to a
      ## % point from the half-tree we just generated.
      if((sprime == 1) & (runif(1,0,1) < nprime/n)){
        theta <- thetaprime
        logp <- logpprime
        grad <- gradprime
      }
      ##% Update number of valid points we've seen.
      n <- n + nprime;
      ## % Decide if it's time to stop.
      s <- sprime && stop_criterion(thetaminus, thetaplus, rminus, rplus)
      ## % Increment depth.
      j <- j + 1
      print(paste("j = ", j, sep=""))
      flush.console()
    }
    
    ##% Do adaptation of epsilon if we're still doing burn-in.'
    eta <- 1 / (m - 1 + t0)
    Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
    if(m <= Madapt){
      epsilon <- exp(mu - sqrt(m-1)/gamma * Hbar)
      eta <- (m-1)^-kappa
      epsilonbar <- exp((1 - eta) * log(epsilonbar) + eta * log(epsilon))
    }
    else{
      epsilon <- epsilonbar
    }
    print(paste("m = ", m, sep=""))
    print(paste("epsilon = ", epsilon, sep=""))
    
    
    samples[m,] <- theta
    
  }
  samples <- samples[(Madapt+1):(Madapt+M), ]
  return(list(samples=samples, epsilon=epsilon))
}
