source("../code/dlmasisfun.R")
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
sams <- c("fullcis")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(4)
}
system.time(fullcissamout <- fullsim(samplers, simdata, n, burn, parallel))
save(fullcissamout, file="cissamout.RData")

print("Mixing Simulations Complete")
