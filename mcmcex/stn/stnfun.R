source("../mcmcexfun.R")

stnstatesam <- function(n, start, dat, av=0, aw=0, bv=0, bw=0){
  T <- length(dat)
  V <- start[1]
  W <- start[2]
  R <- W/V
  out <- mcmc(matrix(0, nrow=n, ncol=T+9))
  colnames(out) <- c(paste("theta",0:T,sep=""),"V","W", "time",
                     "logconV", "adrejV", "logconW", "adrejW", "kernel")
  for(i in 1:n){
    time <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    theta <- thetaiter(dat, V, W)
    VRiter <- VRthetaiter(dat, theta, av, aw, bv, bw)
    V <- VRiter[1]
    R <- VRiter[2]
    W <- R*V
    time2 <- sum(as.numeric(unlist(strsplit(format(Sys.time(), "%M:%OS3"), ":")))*c(60, 1))
    time <- time2-time
    out[i,] <- c(theta,V,W, time,c(NA,NA,NA,NA), NA)
  }
  return(out)
}


VRthetaiter <- function(dat, theta, av, aw, bv, bw){
  T <- length(dat)
  Va <- av + T/2
  Vb <- bv + sum((dat - theta[-1])^2)/2
  Wa <- aw + T/2
  Wb <- bw + sum((theta[-1]-theta[-(T+1)])^2)/2
  V <- rinvgamma(1, Va, Vb)
  R <- rinvgamma(1, Wa, Wb/V)
  return(c(V,R))
}



