source("varfun.r")

dat2010 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2010.csv")
dat2011 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2011.csv")

## finds the 3 heated 25cm temperature 2011 time series and averages them
idx <- grep("H.25cm.T", colnames(dat2011))
datheat2011 <- dat2011[,idx]
datMH11 <- apply(datheat2011,1,mean)[100:200]

## Fit model 3b using MLE
## note: parameter estimates are on log scale
MLEest11 <- dlmMLE(datMH11, rep(0,4), mymodMLE, other=c(2,24,3))
MLEstart <- exp(MLEest11$par) # my estimates don't replicate

priorIG <- list(rep(.1,8), rep(.01,8), rep(.001,8))
priorU <- list(rep(10000,4), rep(10000000,4), rep(10000000000,4))
priorHT <- list(c(rep(25,4), rep(1,4)), c(rep(25,4), rep(10,4)), c(rep(10,4), rep(1,4)))
chains <- 3
n <- 100
IGsam <- data.frame(NULL)
URsam <- data.frame(NULL)
UMsam <- data.frame(NULL)
HTsam <- data.frame(NULL)
starts <- list()
starts[[1]] <- MLEstart*1/100
starts[[2]] <- MLEstart
starts[[3]] <- MLEstart*100

## inverse cdf is probably better for UR and UM sam
system.time(
  for(j in 1:3){
    for(i in 1:chains){
      print(c(i,j))
      IGsamtemp <- postsamIG(n, datMH11, priorIG[[j]], starts[[i]])
      IGsam <- rbind(IGsam, data.frame(IGsamtemp, prior="IG", chain=i, hyper=j))
      print("IG finished")
      URsamtemp <- postsamUR(n, datMH11, priorU[[j]], MLEstart)
      URsam <- rbind(URsam, data.frame(URsamtemp, prior="UR", chain=i, hyper=j))
      print("UR finished")
      UMsamtemp <- postsamUM(n, datMH11, priorU[[j]], starts[[i]])
      UMsam <- rbind(UMsam, data.frame(UMsamtemp, prior="UM", chain=i, hyper=j))
      print("UM finished")
      HTsamtemp <- postsamHT(n, datMH11, priorHT[[j]], starts[[i]])
      HTsam <- rbind(HTsam, data.frame(HTsamtemp, prior="HT", chain=i, hyper=j))
      print("HT finished")
    }
  }
  )

##plot ergodic means from multiple chains
k <- 1
par(mfrow=c(2,2))
nams <- c("V","Wm", "Wb", "Ws")
for(j in 1:3){
  top <- max(c(cumsum(IGsam[[j]][[1]][,k])/1:n,cumsum(IGsam[[j]][[2]][,k])/1:n,cumsum(IGsam[[j]][[3]][,k])/1:n))
  plot(ts(cumsum(IGsam[[j]][[1]][,k])/1:n), main=nams[k], ylim=c(0,top))
  lines(ts(cumsum(IGsam[[j]][[2]][,k])/1:n), col="red")
  lines(ts(cumsum(IGsam[[j]][[3]][,k])/1:n), col="blue")
}

k <- 1
par(mfrow=c(2,2))
nams <- c("V","Wm", "Wb", "Ws")
for(j in 1:3){
  top <- max(c(cumsum(UMsam[[j]][[1]][,k])/1:n,cumsum(UMsam[[j]][[2]][,k])/1:n,cumsum(UMsam[[j]][[3]][,k])/1:n))
  plot(ts(cumsum(UMsam[[j]][[1]][,k])/1:n), main=nams[k], ylim=c(0,top))
  lines(ts(cumsum(UMsam[[j]][[2]][,k])/1:n), col="red")
  lines(ts(cumsum(UMsam[[j]][[3]][,k])/1:n), col="blue")
}

k <- 1
par(mfrow=c(2,2))
nams <- c("V","Wm", "Wb", "Ws")
for(j in 1:3){
  top <- max(c(cumsum(HTsam[[j]][[1]][,k])/1:n,cumsum(HTsam[[j]][[2]][,k])/1:n,cumsum(HTsam[[j]][[3]][,k])/1:n))
  plot(ts(cumsum(HTsam[[j]][[1]][,k])/1:n), main=nams[k], ylim=c(0,top))
  lines(ts(cumsum(HTsam[[j]][[2]][,k])/1:n), col="red")
  lines(ts(cumsum(HTsam[[j]][[3]][,k])/1:n), col="blue")
}

lims <- c(.005, 3e-11, .005, .0000000001)
k <- 4

par(mfrow=c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(ts(HTsam[[i]][[j]][,k]), main=paste(i,j), ylim=c(0,lims[k]))
  }
}


par(mfrow=c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(ts(URsam[[i]][[j]][,k]), main=paste(i,j), ylim=c(0,lims[k]))
  }
}


par(mfrow=c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(ts(UMsam[[i]][[j]][,k]), main=paste(i,j), ylim=c(0,lims[k]))
  }
}

lims <- c(.005, 3e-11, .005, .001)
par(mfrow=c(3,3))
for(i in 1:3){
  for(j in 1:3){
    plot(ts(IGsam[[i]][[j]][,k]), main=paste(i,j), ylim=c(0,lims[k]))
  }
}



plot(ts(HTsam[[1]][,2]))

X11.options(type="Xlib")
source("varfun.r")
V <- 1
W <- 1
T <- 100
dat <- llsim(T, V, W, 0, 1)
start <- c(V,W)
n <- 1000
prior <- c(1,1,25,25)
test <- mcmc(statesam(1000, start, dat, prior))
plot(test)

