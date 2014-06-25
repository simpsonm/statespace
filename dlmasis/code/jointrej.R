source("dlmasisfun.R")

logpijoint <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]
  if(lvw < 0){
    out <-  (ayou + T/2)*( log(a*exp(lvw) + 2*b*exp(lvw/2) + c)) + ame*lvw + bet*exp(-lvw)
  }
  else{
    out <-  (ayou + T/2)*( lvw + log(a + 2*b*exp(-lvw/2) + c*exp(-lvw))) + ame*lvw + bet*exp(-lvw)
  }
  return(-out)
}

gradlpjoint <- function(lvw, par){
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

logpirejjoint <- function(x, par, mn, std, df){
  out <- logpijoint(x, par) - dt((x-mn)/std, df, log=TRUE)-log(std)
  return(out)
}

VWgamiter <- function(dat, gam, theta, av, aw, bv, bw){
  T <- length(dat)
  cgam <- cumsum(gam[-1])
  gam0 <- gam[1]
  ag <- sum( cgam^2 )/2
  bg <- -sum( (dat - gam0) * cgam )/2
  cg <- bv + sum( (dat - gam0)^2 )/2
  par <- c(ag, bg, cg, aw, av, bw, T)
  Va <- av + T/2
  Vb <- bv + sum( (dat - theta[-1])^2 )/2
  V <- rinvgamma(1, Va, Vb)
  mn <- uniroot(gradlpjoint, c(-100,100), par=par)$root
  g <- ag*exp(mn) + 2*bg*exp(mn/2) + cg
  gp <- ag*exp(mn) + bg*exp(mn/2)
  gpp <- ag*exp(mn) + bg*exp(mn/2)/2
  prec <- abs((av + T/2)*(g*gpp - gp^2)/gpp^2 + bw*exp(-mn))
  std <- sqrt(1/prec)
  df <- 1
  M <- optimize(logpirejjoint, c(-100,100), maximum=TRUE, par=par, mn=mn, std=std, df=df)$objective
  rej <- TRUE
  rejit <- 1
  while(rej){
    prop <- rt(1, df)*std + mn
    R <- logpirejjoint(prop, par, mn, std, df) - M
    u <- runif(1,0,1)
    if(log(u)<R)
        rej <- FALSE
    rejit <- rejit + 1
  }
  W <- exp(prop)
  return(c(V,W))
}

VWpsiiter <- function(dat, psi, theta, av, aw, bv, bw){
  T <- length(dat)
  psi0 <- psi[1]
  psiLT <- c(0,psi[-1])
  Lpsi <- psiLT[-1] - psiLT[-(T+1)]
  ys <- c(psi0, dat)
  Ly <- ys[-1] - ys[-(T+1)]
  ap <- sum(Lpsi^2)/2
  bp <- -sum(Lpsi*Ly)/2
  cp <- bw + sum( Ly^2 )/2
  par <- c(ap, bp, cp, av, aw, bv, T)
  Wa <- aw + T/2
  Wb <- bw + sum( (theta[-1] - theta[-(T+1)])^2 )/2
  W <- rinvgamma(1, Wa, Wb)
  mn <- uniroot(gradlpjoint, c(-100,100), par=par)$root
  g <- ap*exp(mn) + 2*bp*exp(mn/2) + cp
  gp <- ap*exp(mn) + bp*exp(mn/2)
  gpp <- ap*exp(mn) + bp*exp(mn/2)/2
  prec <- abs((aw + T/2)*(g*gpp - gp^2)/gpp^2 + bv*exp(-mn))
  std <- sqrt(1/prec)
  df <- 1
  M <- optimize(logpirejjoint, c(-10,10), par=par, mn=mn, std=std, df=df, maximum=TRUE)$objective
  rej <- TRUE
  rejit <- 1
  while(rej){
    prop <- rt(1, df)*std + mn
    R <- logpirejjoint(prop, par, mn, std, df) - M
    u <- runif(1,0,1)
    if(log(u)<R)
        rej <- FALSE
    rejit <- rejit + 1
  }
  V <- exp(prop)

  return(c(V,W))
}

distjointsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    theta <- awolthsmooth(dat, V, W, m0, C0)
    gam <- gamtrans(theta, W)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VW <- VWgamiter(dat, gam, theta, av, aw, bv, bw)
    V <- VW[1]
    W <- VW[2]
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, theta)
  }
  return(out)
}


errorjointsam <- function(n, start, dat, av, aw, bv, bw, m0, C0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  out <- samoutsetup(n, T)
  for(i in 1:n){
    ptma <- proc.time()
    psi <- awolpssmooth(dat, V, W, m0, C0)
    theta <- thetapsitrans(dat, psi, V)
    ptmb <- proc.time()
    smoothtime <- ptmb[3]-ptma[3]
    VW <- VWpsiiter(dat, psi, theta, av, aw, bv, bw)
    V <- VW[1]
    W <- VW[2]
    out[i,] <- c(NA,NA,NA,NA,NA, smoothtime, V, W, theta)
  }
  return(out)
}


jointwrap <- function(par, n, samp){
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
  if(samp=="error")
      time <- system.time(out <- errorjointsam(n, start, dat, av, aw, bv, bw, m0, C0))
  if(samp=="dist")
      time <- system.time(out <- distjointsam(n, start, dat, av, aw, bv, bw, m0, C0))
  outdat <- data.frame(out)
  outdat$time <- time[1]
  cols <- ncol(outdat)
  outdat <- outdat[,c( cols, 1:(cols-1) )]
  print(paste(c(samp, " T=", T, " V=", start[1], " W=", start[2], " FINISHED"), collapse=""))
  return(outdat)

}

