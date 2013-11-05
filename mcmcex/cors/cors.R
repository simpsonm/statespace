source("../mcmcexfun.R")
load("../mixing/triintSAM.RData")
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
  mod <- dlmModPoly(order=1, dV=V, dW=W)
  filt <- dlmFilter(yt, mod)
  theta <- dlmBSample(filt)
  bv <- sum( (yt - theta[-1])^2 )
  bw <- sum( (theta[-1] - theta[-(T.T+1)])^2 )
  out <- dat
  out$bv <- bv
  out$bw <- bw
  return(out)
}

samcor <- ddply(.data=samshort2, .variables=.(iter), .fun=itsim, y=simdata)
rm(samshort)

newcors <- function(samcor){
  Vbv <- cor(samcor$V, samcor$bv)
  Wbv <- cor(samcor$W, samcor$bv)
  Vbw <- cor(samcor$V, samcor$bw)
  Wbw <- cor(samcor$W, samcor$bw)
  out <- data.frame(V=samcor$V.T[1], W=samcor$W.T[1], T=samcor$T.T[1])
  out$Vbv <- Vbv
  out$Wbv <- Wbv
  out$Vbw <- Vbw
  out$Wbw <- Wbw
  return(out)
}

parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(8)
}

newpostcors <- ddply(.data=samcor, .variables=.(V.T, W.T, T.T), .fun=newcors, .parallel=parallel)
save(newpostcors, file="newpostcors.RData")
