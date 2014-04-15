source("../code/mcmcexfun.R")
source("../code/allatoncefun.R")
load("trialtSAM.RData")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
set.seed(234134530)
samshort$iter <- 1:length(samshort[,1])

itsim <- function(dat, y){
  V.T <- dat$V.T[1]
  W.T <- dat$W.T[1]
  T.T <- dat$T.T[1]
  V <- dat$V[1]
  W <- dat$W[1]
  yt <- y$y[y$V.T==V.T & y$W.T==W.T & y$T.T==T.T]
  theta <- awolthsmooth(yt, V, W, 0, 10^7)
  bv <- sum( (yt - theta[-1])^2 )
  bw <- sum( (theta[-1] - theta[-(T.T+1)])^2 )
  out <- dat
  out$bv <- bv
  out$bw <- bw
  gam <- gamtrans(theta, W)
  psi <- psitrans(yt, theta, V)
  gam0 <- gam[1]
  cgam <- cumsum(gam[-1])
  agam <- sum(cgam^2)/(2*V)
  bgam <- sum((yt-gam0)*cgam)/V
  ys <- c(psi[1], yt)
  Ly <- ys[-1] - ys[-(T.T+1)]
  psiLT <- c(0,psi[-1])
  Lpsi <- psiLT[-1] - psiLT[-(T.T+1)]
  apsi <- sum(Lpsi^2)/2/W
  bpsi <- sum(Lpsi*Ly)/W
  out$agam <- agam
  out$bgam <- bgam
  out$apsi <- apsi
  out$bpsi <- bpsi
  return(out)
}

parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(8)
}


samcor <- ddply(.data=samshort, .variables=.(iter), .fun=itsim, y=simdata, .parallel=parallel)
rm(samshort)

newcors <- function(samcor){
  Vbv <- cor(samcor$V, samcor$bv)
  Wbv <- cor(samcor$W, samcor$bv)
  Vbw <- cor(samcor$V, samcor$bw)
  Wbw <- cor(samcor$W, samcor$bw)
  Wagam <- cor(samcor$W, samcor$agam)
  Wbgam <- cor(samcor$W, samcor$bgam)
  Vapsi <- cor(samcor$V, samcor$apsi)
  Vbpsi <- cor(samcor$V, samcor$bpsi)

  out <- data.frame(V=samcor$V.T[1], W=samcor$W.T[1], T=samcor$T.T[1])
  out$Vbv <- Vbv
  out$Wbv <- Wbv
  out$Vbw <- Vbw
  out$Wbw <- Wbw
  out$Wagam <- Wagam
  out$Wbgam <- Wbgam
  out$Vapsi <- Vapsi
  out$Vbpsi <- Vbpsi
  return(out)
}

newpostcors <- ddply(.data=samcor, .variables=.(V.T, W.T, T.T), .fun=newcors, .parallel=parallel)
save(newpostcors, file="newpostcors.RData")
