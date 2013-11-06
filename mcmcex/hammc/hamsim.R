source("../mcmcexfun.R")
source("hamfun.R")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
sams <- c("dist", "error")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 3000
burn <- 500
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(4)
}
system.time(hamout <- hamfullsim(samplers, simdata, n, burn, parallel))
save(hamout, file="samout.RData")





T.T <- 1000
V.T <- 1000
W.T <- .01

dat <- simdata$y[simdata$V.T==V.T & simdata$W.T==W.T & simdata$T.T==T.T]
n <- 1000
av <- 5
aw <- 5
bv <- V.T*(av-1)
bw <- W.T*(aw-1)
eps <- 0.01349 * min(V.T/W.T, W.T/V.T)
L <- 100 * max(V.T/W.T, W.T/V.T)
start <- c(V.T,W.T)


dtest <- distjointsam(n, start, dat, av, aw, bv, bw, eps, L)
etest <- errorjointsam(n, start, dat, av, aw, bv, bw, eps, L)

Dtest <- data.frame(dtest)
Etest <- data.frame(etest)

dvw <- mcmc(cbind(Dtest$V, Dtest$W))
evw <- mcmc(cbind(Etest$V, Etest$W))

mean(Dtest$acc)
mean(Etest$acc)

summary(dvw)
summary(evw)


par(mfrow=c(2,2))
plot(ts(Dtest$V), ylab="VD")
plot(ts(Dtest$W), ylab="WD")
plot(ts(Etest$V), ylab="VE")
plot(ts(Etest$W), ylab="WE")
