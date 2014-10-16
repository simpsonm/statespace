## returns TRUE if target density of log-concave, FALSE otherwise
logcon <- function(a,b,c, eps=.1){
  RHS <- sqrt((a + 1)^3/b)*4/3*sqrt(2/3) + eps
  ## + eps to make sure ARS algorithm doesn't fail on
  ## near non-log-concave cases
  out <- (c > RHS)
  return(out)
}

## log of the conditional posterior density of W (V) given V (W), gamma (psi), data
logpiVW <- function(x, a, b, c, d){
  out <- -(a + 1)*log(x) - b/x + c*sqrt(x) - d*x
  return(out)
}

## first derivative of log of the conditional posterior of W (V)
logpiVWprime <- function(x, a, b, c, d){
  out <- - (a + 1)/x + b/(x^2) + c/2/sqrt(x) - d
  return(out)
}

## compute M in logtarget - logprop <= M
propM <- function(df, a, b, c, d, mn, propvar){
  M <- optimize(logpirej, c(-10^2,10^2), maximum=TRUE, a, b, c, d, mn, propvar, df)
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
logpilVW <- function(x, a, b, c, d){
  out <- -d*exp(x) + c*exp(x/2) - a*x - b*exp(-x)
  if(exp(x) == Inf)
      out <- -Inf
  return(out)
}

## first derivative of log of the conditional posterior of logW (logV)
logpilVWprime <- function(x, a, b, c, d){
  out <- -d*exp(x) + c*exp(x/2)/2 - a + b*exp(-x)
  return(out)
}

## difference between log conditional posterior of logW (logV) and the proposal density
logpirej <- function(x, a, b, c, d, mn, propvar, df){
  out <- logpilVW(x, a, b, c, d) - dtprop(x, mn, propvar, df)
  return(out)
}


## samples W conditional on V,gamma
Wgamiter <- function(a,b,c,d){
  lcon <- logcon(a, b, c)
  adrej <- lcon
  if(lcon){
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, c=c, d=d)$maximum
##    RHS <- sqrt((a + 1)^3/b)*4/3*sqrt(2/3)
##    pout <- matrix(c(mn,a,b,c,d,RHS),nrow=1)
##    colnames(pout) <- c("mn", "a", "b", "c", "d","RHS")
##    print(pout)
    try(W <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(mn/2, mn, mn*2),
                      lb=TRUE, xlb=0, a=a, b=b, c=c, d=d))
    if(W==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    W <- VWrejiter(a, b, c, d)
  }
  return(W)
}


## samples V conditional on W,psi 
Vpsiiter <- function(a,b,c,d){
  lcon <- logcon(a, b, c)
  adrej <- lcon
  if(lcon){
    mn <- optimize(logpiVW, c(0,10^10), maximum=TRUE, a=a, b=b, c=c, d=d)$maximum
##    RHS <- sqrt((a + 1)^3/b)*4/3*sqrt(2/3)
##    pout <- matrix(c(mn,a,b,c,d,RHS),nrow=1)
##    colnames(pout) <- c("mn", "a", "b", "c", "d", "RHS")
##    print(pout)
    try(V <- ars(n=1, logpiVW, logpiVWprime, ns=200, x=c(mn/2, mn, mn*2),
                 lb=TRUE, xlb=0, a=a, b=b, c=c, d=d))
    if(V==0){
      adrej <- FALSE
    }
  }
  if(!adrej){
    V <- VWrejiter(a, b, c, d)
  }
  return(V)
}

VWrejiter <- function(a, b, c, d){
  mn <- optimize(logpilVW, c(-10^2,10^2), maximum=TRUE, a=a, b=b, c=c, d=d)$maximum
  propvar <- - 1 /( -d*exp(mn) + c*exp(mn/2)/4 - b*exp(-mn) )
  df <- 1
  rej <- TRUE
  rejit <- 1
  M <- optimize(logpirej, c(-10^2, 10^2), a=a, b=b, c=c, d=d, mn=mn,
                propvar=propvar, df=df, maximum=TRUE)$objective
  while(rej){
    prop <- rtprop(1, mn, propvar, df)
    R <- logpirej(prop, a, b, c, d, mn, propvar, df) - M
    u <- runif(1,0,1)
    if(log(u)<R){
      W <- exp(prop)
      rej <- FALSE
    }
    rejit <- rejit + 1
  }
  return(W)
}

## transforms gamma to theta
thetagamtrans <- function(gam, W){
  N <- length(gam)-1
  theta <- rep(0,N+1)
  theta[1] <- gam[1]
  sW <- sqrt(W)
  for(i in 1:N){
    theta[i+1] <- sW*gam[i+1] + theta[i]
  }
  return(theta)
}



