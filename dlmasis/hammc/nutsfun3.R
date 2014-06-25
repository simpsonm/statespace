nuts_da <- function(f, M, Madapt, theta0, delta=0.6){
  ## assert(size(theta0, 1) == 1);
  D <- length(theta0)
  samples <- matrix(0, nrow=M+Madapt, ncol=D)
  fout <- f(theta0)
  logp <- fout$logp
  grad <- fout$grad
  samples[1,] <- theta0
  
  ## Choose a reasonable first epsilon by a simple heuristic.
  epsilon <- find_reasonable_epsilon(theta0, grad, logp, f, par)
  print("epsilon initialized")

  ## Parameters to the dual averaging algorithm.
  gamma <- 0.05
  t0 <- 10
  kappa <- 0.75
  mu <- log(10*epsilon)
  ## Initialize dual averaging algorithm.
  epsilonbar <- 1
  Hbar <- 0

  for(m in  2:(M+Madapt)){
    ## Resample momenta.
    r0 <- rnorm(D, 0, 1);
    ## Joint log-probability of theta and momenta r.
    joint <- logp - 0.5 * sum(r0^2)    
    ## Resample u ~ uniform([0, exp(joint)]).
    ## Equivalent to (log(u) - joint) ~ exponential(1).
    logu <- joint - rexp(1)
    ## Initialize tree.
    thetaminus <- samples[m-1, ]
    thetaplus <- samples[m-1, ]
    rminus <- r0
    rplus <- r0
    gradminus <- grad
    gradplus <- grad
    ## Initial height j = 0.
    j <- 0
    ## If all else fails, the next sample is the previous sample.
    samples[m, ] <- samples[m-1, ]
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
        samples[m, ] <- thetaprime
        logp <- logpprime
        grad <- gradprime
      }
      ##% Update number of valid points we've seen.
      n <- n + nprime;
      ## % Decide if it's time to stop.
      s <- sprime && stop_criterion(thetaminus, thetaplus, rminus, rplus)
      ## % Increment depth.
      j <- j + 1
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
  }
  samples <- samples[(Madapt+1):(Madapt+M), ]
  return(list(samples=samples, epsilon=epsilon))
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
      npr <- nprime2 / (nprime + nprime2)
      if(is.nan(npr))
        npr <- 0
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
    }
  }
  return(epsilon)
}


