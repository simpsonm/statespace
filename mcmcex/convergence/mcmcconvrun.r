source("../mcmcexfun.r")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V

simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
M <- expand.grid(ch=1, V.S=c(1/100, 100), W.S=c(1/100, 100))
M[5,] <- c(1,1,1)
simdatatemp <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata <- ddply(M, .(ch, V.S, W.S), dfun, simgrid=simdatatemp)
K <- 100
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt", "sdkern",
          "sekern", "dekern", "trikern")
samplers <- expand.grid(sampler=sams, iter=1:K)
samplers$sampler <- as.character(samplers$sampler)
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(4)
}
a1 <- 5
a2 <- a1
ns <- data.frame(n=c(50, 100, 500, 1000))
conv.diag <- ddply(ns, .(n), chains.diag, samplers=samplers, simdata=simdata, a1=a1, a2=a2,
                   parallel=parallel, .parallel=parallel)
save(conv.diag, "conv.diag.RData")
conv.diag$iter <- NULL
meanconv.diag <- ddply(conv.diag, .(sampler, V.T, W.T, T.T), function(x){
    data.frame(G.D.M=mean(x$G.D.M), G.D.V=mean(x$G.D.V), G.D.W=mean(x$G.D.W))
    }, .parallel=parallel)
save(meanconv.diag, "meanconv.diag.RData")
rm(conv.diag)
rm(meanconv.diag)

M <- data.frame(ch=1:5, V.S=c(1,1,1,1,1), W.S=c(1,1,1,1,1))
simdata <- ddply(M, .(ch, V.S, W.S), dfun, simgrid=simdatatemp)
K <- 10
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt", "sdkern",
          "sekern", "dekern", "trikern")
samplers <- expand.grid(sampler=sams[1], iter=1)
samplers$sampler <- as.character(samplers$sampler)
ns <- data.frame(n=3000)
mix.diag <- ddply(ns, .(n), chains.diag, samplers=samplers, simdata=simdata, a1=a1, a2=a2,
                   parallel=parallel, mix=TRUE, .parallel=parallel)
mix.diag$n <- mix.diag$n2
mix.diag$n2 <- NULL
save(mix.diag, "mix.diag.RData")

