## code for running the simulations from the paper.
## WARNING: THIS WILL TAKE ON THE ORDER OF WEEKS TO COMPLETE ON A UNIVERSITY CLUSTER
source("dlmasisfun.R")
set.seed(152893627) ## needed to replicate my dataset
## set the time series lengths and the true values of V and W
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
## create all simulated datasets
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
## add hyperparameters to datasets
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
## list of samplers to use to sample each posterior
sams <- c("state", "dist", "error", "sdint", "seint", "deint", "triint",
          "sdalt", "sealt", "dealt", "trialt", "wdist", "werror", "fullcis")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
## If doMC package is installed, attempt to use 8 threads for parallel processing
parallel <- require(doMC, quietly=TRUE) 
if(parallel){
  registerDoMC(8)
}

## can randomize seed here for obtaining a different posterior sample
## from the same posterior distribution

## run the samplers
system.time(samout <- fullsim(samplers, simdata, n, burn, parallel))
save(samout, file="samout.RData")
