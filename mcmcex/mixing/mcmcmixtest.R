source("../code/mcmcexfun.R")
set.seed(152893627)
T <- c(10)
V <- c(1/10)
W <- 1/V
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=1)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt",
          "sdkern", "sekern", "dekern", "trikern",
          "fullcis", "partialcis")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sampler <- sams
n <- 2000
burn <- 100

system.time(samout <- fullsim(samplers, simdata, n, burn, FALSE))


[1] "fullcis T=10 V=1 W=1"

Process R aborted (core dumped) at Wed Apr  2 18:08:42 2014
Trap in ARS: infinite while in update_ of ars.cpp near l. 810

terminate called after throwing an instance of 'returnR'



source("../code/mcmcexfun.R")
T <- 10
V <- 100
W <- 100
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=10^7)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
par <- simdata
n <- 3000


set.seed(12314124)
dtime <- system.time(dist <- samwrap(par, n, "dist"))

etime <- system.time(error <- samwrap(par, n, "error"))
dcon <- apply(dist[,c(6:9)], 2, mean)
econ <- apply(error[,c(2:5)], 2, mean)
V <- 1
W <- 1000
simgrid <- expand.grid(V.T=V, W.T=W, T.T=T)
simdata <- ddply(simgrid, .(V.T, W.T, T.T), lldsim, m0=0, C0=10^7)
simdata$av <- 5
simdata$aw <- 5
simdata$bv <- (simdata$av-1)*simdata$V.T
simdata$bw <- (simdata$aw-1)*simdata$W.T
simdata$m0 <- 0
simdata$C0 <- 10^7
par <- simdata
n <- 500
dtime2 <- system.time(dist <- samwrap(par, n, "dist"))
etime2 <- system.time(error <- samwrap(par, n, "error"))
dcon2 <- apply(dist[,c(6:9)], 2, mean)
econ2 <- apply(error[,c(2:5)], 2, mean)


etime
dtime
etime2
dtime2


disttime <- system.time(dist <- samwrap(par, n, "dist"))
errortime <- system.time(error <- samwrap(par, n, "error"))
deinttime <- system.time(deint <- samwrap(par, n, "deint"))
fullcistime <- system.time(fullcis <- samwrap(par, n, "fullcis"))

disttime2 <- system.time(dist <- samwrap(par, n, "dist"))
errortime2 <- system.time(error <- samwrap(par, n, "error"))
deinttime2 <- system.time(deint <- samwrap(par, n, "deint"))
fullcistime2 <- system.time(fullcis <- samwrap(par, n, "fullcis"))

disttime3 <- system.time(dist <- samwrap(par, n, "dist"))
errortime3 <- system.time(error <- samwrap(par, n, "error"))
deinttime3 <- system.time(deint <- samwrap(par, n, "deint"))
fullcistime3 <- system.time(fullcis <- samwrap(par, n, "fullcis"))

times1 <- c(disttime[3], errortime[3], deinttime[3], fullcistime[3])
times2 <- c(disttime2[3], errortime2[3], deinttime2[3], fullcistime2[3])
times3 <- c(disttime3[3], errortime3[3], deinttime3[3], fullcistime3[3])

normaltimes <- cbind(times1, times2, times3)
newtimes <- cbind(times1, times2, times3)

source("../code/mcmcexfun.R")
set.seed(152893627)
T <- 10
V <- 100
W <- 100
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

fullcis <- samwrap(par, n, "fullcis")
state <- samwrap(par, n, "state")
dist <- samwrap(par, n, "dist")
error <- samwrap(par, n, "error")
sdalt <- samwrap(par, n, "sdalt")
sdint <- samwrap(par, n, "sdint")
sealt <- samwrap(par, n, "sealt")
seint <- samwrap(par, n, "seint")
dealt <- samwrap(par, n, "dealt")
deint <- samwrap(par, n, "deint")
trialt <- samwrap(par, n, "trialt")
triint <- samwrap(par, n, "triint")
fullcis <- samwrap(par, n, "fullcis")
partialcis <- samwrap(par, n, "partialcis")
sdkern <- samwrap(par, n, "sdkern")
sekern <- samwrap(par, n, "sekern")
dekern <- samwrap(par, n, "dekern")
trikern <- samwrap(par, n, "trikern")

apply(state, 2, mean, na.rm=TRUE)[1:9]
apply(dist, 2, mean, na.rm=TRUE)[1:9]
apply(error, 2, mean, na.rm=TRUE)[1:9]
apply(sdalt, 2, mean, na.rm=TRUE)[1:9]
apply(sdint, 2, mean, na.rm=TRUE)[1:9]
apply(sealt, 2, mean, na.rm=TRUE)[1:9]
apply(seint, 2, mean, na.rm=TRUE)[1:9]
apply(dealt, 2, mean, na.rm=TRUE)[1:9]
apply(deint, 2, mean, na.rm=TRUE)[1:9]
apply(trialt, 2, mean, na.rm=TRUE)[1:9]
apply(triint, 2, mean, na.rm=TRUE)[1:9]
apply(fullcis, 2, mean, na.rm=TRUE)[1:9]
apply(partialcis, 2, mean, na.rm=TRUE)[1:9]
apply(sdkern, 2, mean, na.rm=TRUE)[1:9]
apply(sekern, 2, mean, na.rm=TRUE)[1:9]
apply(dekern, 2, mean, na.rm=TRUE)[1:9]
apply(trikern, 2, mean, na.rm=TRUE)[1:9]

par(mfrow=c(2,1))
plot(1:n, cumsum(state$V)/1:n, ylab="V", type="l", ylim=c(.5*V, 1.5*V))
lines(1:n, cumsum(dist$V)/1:n, col="red")
lines(1:n, cumsum(error$V)/1:n, col="blue")
lines(1:n, cumsum(sdalt$V)/1:n, col="green")
lines(1:n, cumsum(sdint$V)/1:n, col="orange")
lines(1:n, cumsum(sealt$V)/1:n, col="purple")
lines(1:n, cumsum(seint$V)/1:n, col="pink")
lines(1:n, cumsum(dealt$V)/1:n, col="yellow")
lines(1:n, cumsum(deint$V)/1:n, col="turquoise")
lines(1:n, cumsum(triint$V)/1:n, col="maroon")
lines(1:n, cumsum(trialt$V)/1:n, col="brown")
lines(1:n, cumsum(fullcis$V)/1:n, col="magenta")
lines(1:n, cumsum(partialcis$V)/1:n, col="darkblue")
lines(1:n, cumsum(sdkern$V)/1:n, col="darkorange")
lines(1:n, cumsum(sekern$V)/1:n, col="cyan")
lines(1:n, cumsum(dekern$V)/1:n, col="black")
lines(1:n, cumsum(trikern$V)/1:n, col="black")
plot(1:n, cumsum(state$W)/1:n, ylab="W", type="l", ylim=c(.5*W, 1.5*W))
lines(1:n, cumsum(dist$W)/1:n, col="red")
lines(1:n, cumsum(error$W)/1:n, col="blue")
lines(1:n, cumsum(sdalt$W)/1:n, col="green")
lines(1:n, cumsum(sdint$W)/1:n, col="orange")
lines(1:n, cumsum(sealt$W)/1:n, col="purple")
lines(1:n, cumsum(seint$W)/1:n, col="pink")
lines(1:n, cumsum(dealt$W)/1:n, col="yellow")
lines(1:n, cumsum(deint$W)/1:n, col="turquoise")
lines(1:n, cumsum(triint$W)/1:n, col="maroon")
lines(1:n, cumsum(trialt$W)/1:n, col="brown")
lines(1:n, cumsum(fullcis$W)/1:n, col="magenta")
lines(1:n, cumsum(partialcis$W)/1:n, col="darkblue")
lines(1:n, cumsum(sdkern$W)/1:n, col="darkorange")
lines(1:n, cumsum(sekern$W)/1:n, col="cyan")
lines(1:n, cumsum(dekern$W)/1:n, col="black")
lines(1:n, cumsum(trikern$W)/1:n, col="black")

par(mfrow=c(2,1))
plot(1:n, cumsum(state$V)/1:n, ylab="V", type="l", ylim=c(.5*V, 1.5*V))
lines(1:n, cumsum(dist$V)/1:n, col="red")
lines(1:n, cumsum(error$V)/1:n, col="blue")
lines(1:n, cumsum(dealt$V)/1:n, col="yellow", lty=2)
lines(1:n, cumsum(deint$V)/1:n, col="turquoise", lty=2)
plot(1:n, cumsum(state$W)/1:n, ylab="W", type="l", ylim=c(.5*W, 1.5*W))
lines(1:n, cumsum(dist$W)/1:n, col="red")
lines(1:n, cumsum(error$W)/1:n, col="blue")
lines(1:n, cumsum(dealt$W)/1:n, col="yellow", lty=2)
lines(1:n, cumsum(deint$W)/1:n, col="turquoise", lty=2)


Vmns <- c(mean(state$V), mean(dist$V), mean(error$V), mean(dealt$V), mean(deint$V), mean(sdalt$V), mean(sdint$V), mean(sealt$V), mean(seint$V), mean(triint$V), mean(trialt$V), mean(fullcis$V), mean(partialcis$V), mean(sdkern$V), mean(sekern$V), mean(dekern$V), mean(trikern$V))
Wmns <- c(mean(state$W), mean(dist$W), mean(error$W), mean(dealt$W), mean(deint$W), mean(sdalt$W), mean(sdint$W), mean(sealt$W), mean(seint$W), mean(triint$W), mean(trialt$W), mean(fullcis$W), mean(partialcis$W), mean(sdkern$W), mean(sekern$W), mean(dekern$W), mean(trikern$W))


mns <- cbind(Vmns, Wmns)
rownames(mns) <- c("state", "dist", "error", "dealt", "deint", "sdalt", "sdint", "sealt", "seint", "triint", "trialt", "fullcis", "partialcis", "sdkern", "sekern", "dekern", "trikern")

mns


