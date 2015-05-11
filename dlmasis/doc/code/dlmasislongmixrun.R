## code for running the simulations from the paper.
## WARNING: THIS WILL TAKE ON THE ORDER OF WEEKS TO COMPLETE ON A UNIVERSITY CLUSTER

source("dlmasislongfun.R")
set.seed(152893627) ## needed to replicate my dataset

## set the time series lengths and the true values of V and W
T <- c(10, 100, 1000)
V <- 10^(c(0:10) / 2 - 2)
W <- V

## create all simulated datasets
simgrid <- expand.grid(V.T = V, W.T = W, T.T = T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0 = 0, C0 = 1)

## add hyperparameters to datasets
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av - 1) * simdata$V.T
simdata$bw <- (simdata$aw - 1) * simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7

## list of samplers to use to sample each posterior
sams <- c("state", "dist", "error", "sdint", "seint", "deint", "triint",
          "sdalt", "sealt", "dealt", "trialt", "wdist", "werror", "fullcis")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 10500
burn <- 500

## If doParallel package is installed, attempt to use 8 threads for parallel processing
## Must use doParallel and 4 threads to replicate my posterior draws
parallel <- require(doParallel, quietly=TRUE) 
if(parallel){
  cl <- makeCluster(4, "FORK") ## Requires Unix system
  registerDoParallel(cl)
  clusterSetRNGStream(cl, iseed = 8910) ## Alter seed for posterior simulations here
  ## Sets up L'Ecuyer RNG on the clusters
}

## run the samplers
system.time(samout <- fullsim(samplers, simdata, n, burn, parallel))
## When parallel = TRUE this will throw a bunch of warnings of this form:
## "4: <anonymous>: ... may be used in an incorrect context: ‘.fun(piece, ...)’"
## This is due to a bug in plyr and should hopefully be fixed shortly

save(samout, file="samoutlongmix.RData")
stopCluster(cl)
