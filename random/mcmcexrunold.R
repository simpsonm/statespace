source("mcmcexfunold.R")

n <- 5000
T <- 15
a1 <- 5
a2 <- 5
Vs <- c(.01, .1, 1, 10, 100, 1000)
Ws <- Vs
VWs <- length(Vs)
set.seed(232)
Vcors <- matrix(0, nrow=VWs, ncol=VWs)
Wcors <- Vcors
THcors <- array(0, dim=c(VWs, VWs, T+1))
sTH <- list()
sV <- list()
sW <- list()

for(i in 1:VWs){
  sTH[[i]] <- list()
  sV[[i]] <- list()
  sW[[i]] <- list()
  for(j in 1:VWs){
    V <- Vs[i]
    W <- Ws[j]
    b1 <- (a1-1)*V
    b2 <- (a2-1)*W
    start <- c(V,W)
    dat <- llsim(T, V, W, 0, 1)
    sam <- simsmooth(n, start, dat, a1, a2, b1, b2)
    sTH[[i]][[j]] <- sam[-c(1:500),1:(T+1)]
    sV[[i]][[j]] <- sam[-c(1:500),T+2]
    sW[[i]][[j]] <- sam[-c(1:500),T+3]
    THcors[i,j,] <- diag(cor(sTH[[i]][[j]][-1,], sTH[[i]][[j]][-(n-500),]))
    Vcors[i,j] <- cor(sV[[i]][[j]][-1], sV[[i]][[j]][-(n-500)])
    Wcors[i,j] <- cor(sW[[i]][[j]][-1], sW[[i]][[j]][-(n-500)])
  }
}

colnames(Vcors) <- paste("V=", Vs, sep="")
colnames(Wcors) <- paste("V=", Vs, sep="")
rownames(Vcors) <- paste("W=", Ws, sep="")
rownames(Wcors) <- paste("W=", Ws, sep="")

Vcors
Wcors

library(xtable)
xtable(Vcors, caption=c("First order autocorrlelation for samples of V"))
xtable(Wcors, caption=c("First order autocorrlelation for samples of W"))


n <- 5000
T <- 15
a1 <- 5
a2 <- 5
V <- 10000
W <- .01
b1 <- (a1-1)*V
b2 <- (a2-1)*W
start <- c(V,W)
set.seed(2342)
dat <- llsim(T, V, W, 0, 1)
sam1 <- distFFBS(n, start, dat, a1, a2, b1, b2, trans=FALSE)
sam2 <- simsmooth(n, start, dat, a1, a2, b1, b2)

round(cbind(summary(sam1)[[1]][,1],summary(sam2)[[1]][,1]),4)

par(mfrow=c(2,2))
plot(ts(sam1[,ncol(sam1)]))
plot(ts(sam2[,ncol(sam1)]))
plot(ts(sam1[,ncol(sam1)-1]))
plot(ts(sam2[,ncol(sam1)-1]))

par(mfrow=c(2,1))
plot(ts(cumsum(sam1[,ncol(sam1)])/1:n), ylim=c(0,10))
lines(ts(cumsum(sam2[,ncol(sam1)])/1:n), col="red")
plot(ts(cumsum(sam1[,ncol(sam1)-1])/1:n), ylim=c(0,15000))
lines(ts(cumsum(sam2[,ncol(sam1)-1])/1:n), col="red")


plot(ts(sam2[,ncol(sam1)]))
plot(ts(sam1[,ncol(sam1)-1]))
plot(ts(sam2[,ncol(sam1)-1]))


source("mcmcexfun.r")
n <- 5000
T <- 15
a1 <- 5
a2 <- 5
V <- 1000
W <- 0.1
b1 <- (a1-1)*V
b2 <- (a2-1)*W
start <- c(V,W)
set.seed(234223123)
dat <- llsim(T, V, W, 0, 1)

sam1 <- sdsmooth(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2, jacob=TRUE)
sam2 <- sdsmooth(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2, jacob=FALSE)
##sam3 <- dissmoothAR(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
sam3 <- simsmooth(n, start, dat, a1, a2, b1, b2)


round(cbind(summary(sam1)[[1]][,1],summary(sam2)[[1]][,1],summary(sam3)[[1]][,1]),4)


round(cbind(summary(sam1)[[1]][,1],summary(sam3)[[1]][,1]),4)


plot(ts(dat))
lines(summary(sam1)[[1]][2:(T+1),1], col="red")
lines(summary(sam3)[[1]][2:(T+1),1], col="blue")


k1 <- T+2
k2 <- T+3

dV <- sam1[,k1]
dW <- sam1[,k2]
ddV <- sam2[,k1]
ddW <- sam2[,k2]
sV <- sam3[,k1]
sW <- sam3[,k2]


par(mfrow=c(3,2))
acf(dV)
acf(dW)
acf(ddV)
acf(ddW)
acf(sV)
acf(sW)

cor(dV[-1],dV[-n])
cor(dW[-1],dW[-n])
cor(ddV[-1],ddV[-n])
cor(ddW[-1],ddW[-n])
cor(sV[-1],sV[-n])
cor(sW[-1],sW[-n])

par(mfrow=c(3,2))
plot(ts(dV))
plot(ts(dW))
plot(ts(ddV))
plot(ts(ddW))
plot(ts(sV))
plot(ts(sW))




source("mcmcexfun.r")
set.seed(123)
n <- 5000
chains <- 3
T <- 10
a1 <- 5
a2 <- 5
V <- 1000
W <- 1
starts <- c(.01, 1, 100)
startV <- V*starts
startW <- W*starts
b1 <- (a1-1)*V
b2 <- (a2-1)*W
dat <- llsim(T, V, W, 0, 1)
simlist <- list()
dislist <- list()

for(j in 1:chains){
  print(j)
  start <- c(startV[j], startW[j])
  dis <- distsmooth(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
  sim <- simsmooth(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
  dislist[[j]] <- dis
  simlist[[j]] <- sim
}

##plot ergodic means from multiple chains
par(mfrow=c(2,1))
plot(ts(cumsum(dislist[[1]][,T+2])/1:n), main="V")#, ylim=c(.5,1.5))
for(i in 2:chains){
  lines(ts(cumsum(dislist[[i]][,T+2])/1:n))
}
for(i in 1:chains){
  lines(ts(cumsum(simlist[[i]][,T+2])/1:n), col="red")
}
plot(ts(cumsum(dislist[[1]][,T+3])/1:n), main="W")#, ylim=c(8,12))
for(i in 2:chains){
  lines(ts(cumsum(dislist[[i]][,T+3])/1:n))
}
for(i in 1:chains){
  lines(ts(cumsum(simlist[[i]][,T+3])/1:n), col="red")
}

par(mfrow=c(4,3))
for(i in 1:chains){
  plot(ts(dislist[[i]][,T+2]), main="dis V", ylim=c(0,4))
  abline(h=mean(dislist[[i]][,T+2]), col="red")
}
for(i in 1:chains){
  plot(ts(simlist[[i]][,T+2]), main="sim V", ylim=c(0,4))
  abline(h=mean(simlist[[i]][,T+2]), col="red")
}
for(i in 1:chains){
  plot(ts(dislist[[i]][,T+3]), main="dis W", ylim=c(0,4))
  abline(h=mean(dislist[[i]][,T+3]), col="red")
}
for(i in 1:chains){
  plot(ts(simlist[[i]][,T+3]), main="sim W", ylim=c(0,4))
  abline(h=mean(simlist[[i]][,T+3]), col="red")
}



### Gelman Rubin
gd1 <- list()
gd2 <- list()
for(i in 1:I){
    a1 <- disRlist[[i]][[1]]
    a2 <- disRlist[[i]][[2]]
    a3 <- disRlist[[i]][[3]]
    gd1[[i]] <- gelman.diag(mcmc.list(a1,a2,a3))
    a1 <- simlist[[i]][[1]]
    a2 <- simlist[[i]][[2]]
    a3 <- simlist[[i]][[3]]
    gd2[[i]] <- gelman.diag(mcmc.list(a1,a2,a3))
  }


