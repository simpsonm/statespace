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
y <- simdata

itsim <- function(dat, y){
  V <- dat$V.T[1]
  W <- dat$W.T[1]
  T <- dat$T.T[1]
  yt <- y$y[y$V.T==V & y$W.T==W & y$T.T==T]
  mod <- dlmModPoly(order=1, dV=V, dW=W)
  filt <- dlmFilter(yt, mod)
  theta <- dlmBSample(filt)
  v2 <- sum( (yt - theta[-1])^2 )
  w2 <- sum( (theta[-1] - theta[-(T+1)])^2 )
  out <- dat
  out$v2 <- v2
  out$w2 <- w2
  return(out)
}

samcor <- ddply(.data=samshort, .variables=.(iter), .fun=itsim, y=simdata)
rm(samshort)

newcors <- function(samcor){
  Vv2 <- cor(samcor$V, samcor$v2)
  Wv2 <- cor(samcor$W, samcor$v2)
  Vw2 <- cor(samcor$V, samcor$w2)
  Ww2 <- cor(samcor$W, samcor$w2)
  out <- data.frame(V=samcor$V.T[1], W=samcor$W.T[1], T=samcor$T.T[1])
  out$Vv2 <- Vv2
  out$Wv2 <- Wv2
  out$Vw2 <- Vw2
  out$Ww2 <- Ww2
  return(out)
}

parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(8)
}

newpostcors <- ddply(.data=samcor, .variables=.(V.T, W.T, T.T), .fun=newcors, .parallel=parallel)
save(newpostcors, file="newpostcors.RData")
