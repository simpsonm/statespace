source("mcmcexfun.R")
source("hamfun.R")
set.seed(152893627)
T <- c(10, 100, 1000)
V <- 10^(c(0:10)/2-2)
W <- V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)

T.T <- 100
V.T <- 10
W.T <- 1

dat <- simdata$y[simdata$V.T==V.T & simdata$W.T==W.T & simdata$T.T==T.T]
n <- 1000
av <- 5
aw <- 5
bv <- V.T*(av-1)
bw <- W.T*(aw-1)
eps <- 0.0001
L <- 100
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
