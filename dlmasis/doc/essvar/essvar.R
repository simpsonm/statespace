## code for running the simulations from the paper.
## WARNING: THIS WILL TAKE ON THE ORDER OF WEEKS TO COMPLETE ON A UNIVERSITY CLUSTER

source("dlmasisfun.R")
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
sams <- c("state", "deint")
T <- 100
V <- c(.01, 100)
W <- V
library(plyr)
simdata <- subset(simdata, T.T==T & V.T %in% V & W.T %in% W)
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

K <- 50
samouts <- NULL

for(k in 1:K){
  system.time(samout <- fullsim(samplers, simdata, n, burn, parallel))
  samout$k <- k
  samouts <- rbind(samouts, samout[,c(1:4,12:14)])
}
stopCluster(cl)

save(samouts, file = "samouts.RData")

load("samouts.RData")

samouts$V.ESP <- samouts$V.ES/10000
samouts$W.ESP <- samouts$W.ES/10000

ESmns <- aggregate(cbind(V.ES, W.ES) ~ sampler + V.T + W.T + T.T, data=samouts, mean)
ESsds <- aggregate(cbind(V.ES, W.ES) ~ sampler + V.T + W.T + T.T, data=samouts, sd)
ESses <- ESsds
ESses$V.ES <- ESses$V.ES/sqrt(K)
ESses$W.ES <- ESses$W.ES/sqrt(K)

ESPmns <- aggregate(cbind(V.ESP, W.ESP) ~ sampler + V.T + W.T + T.T, data=samouts, mean)
ESPsds <- aggregate(cbind(V.ESP, W.ESP) ~ sampler + V.T + W.T + T.T, data=samouts, sd)
ESPses <- ESPsds
ESPses$V.ESP <- ESPses$V.ESP/sqrt(K)
ESPses$W.ESP <- ESPses$W.ESP/sqrt(K)

ESPs <- cbind(ESPmns, ESPsds[,c(5,6)])
ESPs$V.T <- factor(as.character(ESPs$V.T))
ESPs$W.T <- factor(as.character(ESPs$W.T))
colnames(ESPs)[7:8] <- c("V.SD", "W.SD")
colnames(ESPs)[2:4] <- c("V", "W", "T")
ESPs$sampler[ESPs$sampler == "state"] <- "State"
ESPs$sampler[ESPs$sampler == "deint"] <- "SD-SE GIS"
ESPs <- ESPs[order(ESPs$sampler),c(1,2,3,5,7,6,8)]

library(xtable)
xtable(ESPs, digits=4)

