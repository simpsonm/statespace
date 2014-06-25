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





T.T <- 10
W.T <- 1
V.T <- 1

dat <- simdata$y[simdata$V.T==V.T & simdata$W.T==W.T & simdata$T.T==T.T]
n <- 1000
av <- 5
aw <- 5
bv <- V.T*(av-1)
bw <- W.T*(aw-1)
eps <- 0.0001349 
L <- 1000
start <- c(V.T,W.T)


system.time(dtest <- distjointsam(n, start, dat, av, aw, bv, bw, eps, L))
system.time(etest <- errorjointsam(n, start, dat, av, aw, bv, bw, eps, L))

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





source("hamfun.R")
a <- 8
b <- -7
c <- 12
par <- c(8, -7, 12, 5,  5,  4, 10)
par <- c(8, -7, 12, 5,  5,  4, 100)

curve(exp(-lpr2(x, par)), xlim=c(-10,10))


lpr(1, par)
lpr2(1, par)
lpr(10000, par)



lvw <- 1000000

  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

(ayou + T/2)*( lvw +  log(a + 2*b*exp(-lvw/2) + c*exp(-lvw))) + ame*lvw + bet*exp(-lvw)
  
  return(out)
}
