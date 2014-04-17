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
samshort$bv <- 0
samshort$bw <- 0
samshort$agam <- 0
samshort$bgam <- 0
samshort$apsi <- 0
samshort$bpsi <- 0
n <- length(samshort[,1])
for(i in 1:n){
  print(c(i,n))
  V.T <- samshort$V.T[i]
  W.T <- samshort$W.T[i]
  T.T <- samshort$T.T[i]
  V <- samshort$V[i]
  W <- samshort$W[i]
  yt <- y$y[y$V.T==V.T & y$W.T==W.T & y$T.T==T.T]
  theta <- awolthsmooth(yt, V, W, 0, 10^7)
  bv <- sum( (yt - theta[-1])^2 )
  bw <- sum( (theta[-1] - theta[-(T.T+1)])^2 )
  samshort$bv[i] <- bv
  samshort$bw[i] <- bw
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
  samshort$agam[i] <- agam
  samshort$bgam[i] <- bgam
  samshort$apsi[i] <- apsi
  samshort$bpsi[i] <- bpsi
}

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

parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(10)
}

newpostcors <- ddply(.data=samshort, .variables=.(V.T, W.T, T.T), .fun=newcors, .parallel=parallel)
save(newpostcors, file="newpostcors.RData")
