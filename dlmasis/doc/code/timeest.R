## code for estimating the amount of time to run the simulations in the paper
## runs a few small simulations and estimates a simple linear regression

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
ns <- c(10, 25, 50, 100)
outtime <- rep(0,length(ns))
burn <- 5 ## burn < n

## If doParallel package is installed, attempt to use 8 threads for parallel processing
## Must use doParallel and 8 threads to replicate my posterior draws
parallel <- require(doParallel, quietly=TRUE) 
if(parallel){
  cl <- makeCluster(8, "FORK") ## Requires Unix system
  registerDoParallel(cl)
  clusterSetRNGStream(cl, iseed = 32511) ## Alter seed here
  ## Sets up L'Ecuyer RNG on the clusters
}

## run the samplers
for(i in 1:length(ns)){
  n <- ns[i]
  outtime[i] <- system.time(samout <- fullsim(samplers, simdata, n, burn, parallel))[3]
}
stopCluster(cl)

## fit the regression
regdat <- data.frame(n=ns, time=outtime)
o <- lm(time~n, data=regdat)
regdat ## print the data
o ## prints the regression outputs

## predicts time, in hours, to run all samplers w/ sample sizes of
## 3500, 5500, 7500, and 10500 (including burn in)
## then prints these estimates.
predict(o, newdata=list(n=c(3500, 5500, 7500, 10500, 20500)))/3600

o2 <- lm(time~n + I(n^2), data=regdat)
o2
predict(o2, newdata=list(n=c(3500, 5500, 7500, 10500, 20500)))/3600


