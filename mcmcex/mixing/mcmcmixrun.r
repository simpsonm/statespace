source("../mcmcexfun.r")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
a1 <- 5
a2 <- a1
system.time(samout <- fullsim(samplers, simdata, n, burn, a1, a2))
save(samout, file="samout.RData")

