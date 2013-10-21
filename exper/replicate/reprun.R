## code snippets for dealing with the data and fitting model 3b

source("repfun.R")

dat2010 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2010.csv")
dat2011 <- read.csv("../data/GYE_SoilMoistureTemperatureData_2011.csv")

## finds the 3 heated 25cm temperature 2011 time series and averages them
idx <- grep("H.25cm.T", colnames(dat2011))
datheat2011 <- dat2011[,idx]
datMH11 <- apply(datheat2011,1,mean)

## plotting each of the time series
plot(ts(datheat2011[,1]), ylim=c(-5,30))
lines(datheat2011[,2], col="red")
lines(datheat2011[,3], col="blue")
lines(datMH11, col="green")


## Fit model 3b using MLE
## note: parameter estimates are on log scale
system.time(MLEest11 <- dlmMLE(datMH11, rep(0,4), mymodMLE, other=c(2,24,3)))
##   user  system elapsed
##  64.68    0.00   64.81
MLEstart <- exp(MLEest11$par) # my estimates don't replicate

## Fit model 3b using the 3 different sets of priors originally used
## - my trace plots don't replicate
prior1 <- c(10, 100, 10, 100, 10, 100, 10, 100)
prior2 <- c(10, 1000, 10, 1000, 10, 1000, 10, 1000)
prior3 <- c(1, 1000, 1, 1000, 1, 1000, 1, 1000)
priors <- list(prior1, prior2, prior3)
n <- 1000
chains <- list()
times <- list()
for(i in 1:3){
  print(i)
  prior <- priors[[i]]
  tim <- system.time(sam <- postsam(n, datMH11, prior, MLEstart))
  times[[i]] <- tim
  chains[[i]] <- sam
}
save(chains, file="chains.RData")

## Fitting model 3b again with different priors - checking to see
## if there was a confusion in the parameterization of the gamma dist
## again my trace plots don't replicate
prior1 <- c(10, 1/100, 10, 1/100, 10, 1/100, 10, 1/100)
prior2 <- c(10, 1/1000, 10, 1/1000, 10, 1/1000, 10, 1/1000)
prior3 <- c(1, 1/1000, 1, 1/1000, 1, 1/1000, 1, 1/1000)
priors <- list(prior1, prior2, prior3)
n <- 1000
chains2 <- list()
times2 <- list()
for(i in 1:3){
  print(i)
  prior <- priors[[i]]
  tim <- system.time(sam <- postsam(n, datMH11, prior, MLEstart))
  times2[[i]] <- tim
  chains2[[i]] <- sam
}
save(chains2, file="chains2.RData")

## Fitting model 3b with the IG(e,e) prior, letting e shrink to check for
## sensitivity to e. Very sensitive.
prior1 <- c(1, 1, 1, 1, 1, 1, 1, 1)
prior2 <- rep(.1,8)
prior3 <- rep(.01, 8)
priors <- list(prior1, prior2, prior3)
n <- 1000
chains3 <- list()
times3 <- list()
for(i in 1:3){
  print(i)
  prior <- priors[[i]]
  tim <- system.time(sam <- postsam(n, datMH11, prior, MLEstart))
  times3[[i]] <- tim
  chains3[[i]] <- sam
}
save(chains3, file="chains3.RData")

## Plotting trace plots for all three chains for all three sets of priors

load("chains.Rdata")
load("chains2.Rdata")
load("chains3.Rdata")

yls <- c(3,5,5,2)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(ts(chains[[1]][,i]), ylab=colnames(chains[[1]])[i], ylim=c(0,yls[i]))
  lines(chains[[2]][,i], col="red")
  lines(chains[[3]][,i], col="blue")
  abline(h=MLEstart[i], lty=2)
}

yls <- c(.001,.0022,.0015,.00015)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(ts(chains2[[1]][,i]), ylab=colnames(chains[[1]])[i], ylim=c(0,yls[i]))
  lines(chains2[[2]][,i], col="red")
  lines(chains2[[3]][,i], col="blue")
  abline(h=MLEstart[i], lty=2)
}

yls <- rep(1/2,4)
par(mfrow=c(2,2))
for(i in 1:4){
  plot(ts(chains3[[1]][,i]), ylab=colnames(chains3[[1]])[i])#, ylim=c(0,yls[i]))
  lines(chains3[[2]][,i], col="red")
  lines(chains3[[3]][,i], col="blue")
  abline(h=MLEstart[i], lty=2)
}


