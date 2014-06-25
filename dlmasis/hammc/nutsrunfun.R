source("nutsfun2.R")
source("../mcmcexfun.R")
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




T.T <- 100
W.T <- 1
V.T <- 1
dat <- simdata$y[simdata$V.T==V.T & simdata$W.T==W.T & simdata$T.T==T.T]
av <- 5
aw <- 5
bv <- V.T*(av-1)
bw <- W.T*(aw-1)
start <- c(V.T,W.T)
M <- 100
Madapt <- 100

set.seed(2141423424)
system.time(test <- distjointsam(M, Madapt, start, dat, av, aw, bv, bw, 0.5))

sam <- data.frame(test)

par(mfrow=c(2,2))
plot(ts(sam$V))
plot(ts(sam$W))
hist(sam$V)
hist(sam$W)

summary(mcmc(data.frame(sam$V, sam$W)))



lpr2 <- function(theta){
  out <- dnorm(theta, mean=5, sd=2, log=TRUE)
  return(out)
}

gradlpr2 <- function(theta){
  out <- -(theta - 5)/4
  return(out)
}

f2 <- function(theta, par){
  logp <- lpr2(theta)
  grad <- gradlpr2(theta)
  out <- list(logp=logp, grad=grad)
  return(out)
}

M <- 10000
Madapt <- 1000


test <- nuts_da(f2, M, Madapt, 10, 0, 0.6)

sam <- test$samples
eps <- test$epsilon

sam2 <- rnorm(M, mean=5, sd=10)
summary(mcmc(data.frame(sam, sam2)))

mean(sam)
mean(sam2)
