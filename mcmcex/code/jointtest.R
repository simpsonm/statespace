source("jointrej.R")
set.seed(152893627)
T <- 100
W <- 1000
V <- .0001
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=10^7)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
par <- simdata
n <- 1000

dtime <- system.time(dist <- jointwrap(par, n, "dist"))
etime <- system.time(error <- jointwrap(par, n, "error"))
stime <- system.time(state <- samwrap(par, n, "state"))

check <- rbind(apply(state, 2, mean, na.rm=TRUE)[c(1,7:9)],apply(dist, 2, mean, na.rm=TRUE)[c(1,7:9)],apply(error, 2, mean, na.rm=TRUE)[c(1,7:9)])
rownames(check) <- c("state", "dist", "error")
check

###basic problem - rej sampling extremely inefficient for extreme values
