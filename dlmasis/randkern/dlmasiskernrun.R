source("../code/dlmasisfun.R")
set.seed(152893627) ## same simulated data as before
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
initidx <- c(1,4,7,10)
V <- V[initidx]
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
sams <- c("sdkern", "sekern", "dekern", "trikern")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
a1 <- 5
a2 <- a1
set.seed(25212342) ## new seed for simulations
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(2)
}
system.time(samout <- fullsim(samplers, simdata, n, burn,
                              a1, a2, parallel))
save(samout, file="samoutinit.RData")

