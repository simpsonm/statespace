source("../mcmcexfun.r")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
M <- expand.grid(V.S=c(1/100, 100), W.S=c(1/100, 100))
M[5,] <- c(1,1)
simdatatemp <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata <- ddply(M, .(V.S, W.S), dfun, simgrid=simdatatemp)
K <- 100
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt")
samplers <- expand.grid(sampler=sams[1], iter=1:K)
samplers$sampler <- as.character(samplers$sampler)
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(4)
}
a1 <- 5
a2 <- a1
ns <- data.frame(n=c(50, 100, 500, 1000))
conv.diag <- ddply(ns, .(n), mean.diag, simdata=simdata, a1=a1, a2=a2,
                   parallel=parallel, .parallel=parallel)
save(conv.diag, "conv.diag.RData")




