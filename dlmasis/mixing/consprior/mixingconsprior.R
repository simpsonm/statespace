source("../../mcmcexfun.R")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^c(-1, 0, 1)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1) ## sets bv = 4 for all sims, prior mean = 1
simdata$bw <- (simdata$aw-1) ## sets bw = 4 for all sims, prior mean = 1
sams <- c("state", "dist", "error")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(2)
}
system.time(samout <- fullsim(samplers, simdata, n, burn, parallel))
save(samout, file="samout.RData")

