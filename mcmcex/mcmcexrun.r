## Simulates data from 36 local level models at 3 different lengths of the time series
## then fits a local level model to each one using both the state sampler and the scaled
## disturbance sampler. The samples and the autocorrelations from each chain are stored.

## This takes about 18-24 hours to run

## setup and data simulation
source("mcmcexfun.r")
n <- 2000
burn <- 500
Ts <- c(10,100,1000)
LT <- length(Ts)
a1 <- 5
a2 <- 5
Vs <- c(.01, .1, 1, 10, 100, 1000)
LVW <- length(Vs)
Ws <- Vs

set.seed(152893627)
simdata <- list()
for(k in 1:LT){
  T <- Ts[k]
  simdata[[k]] <- list()
  for(i in 1:LVW){
    V <- Vs[i]
    simdata[[k]][[i]] <- list()
    for(j in 1:LVW){
      W <- Ws[j]
      simdata[[k]][[i]][[j]] <- llsim(T, V, W, 0, 1)
    }
  }
}

sVcors <- list()
sWcors <- list()
sTHcors <- list()
ssam <- list()
stimes <- list()


## first do the state sampler
for(k in 1:LT){
  T <- Ts[k]
  sVcors[[k]] <- matrix(0, LVW, LVW)
  sWcors[[k]] <- matrix(0, LVW, LVW)
  sTHcors[[k]] <- array(0, dim=c(LVW, LVW, T+1))
  ssam[[k]] <- list()
  stimes[[k]] <- list()
  colnames(sVcors[[k]]) <- paste("V=", Vs, sep="")
  colnames(sWcors[[k]]) <- paste("V=", Vs, sep="")
  rownames(sVcors[[k]]) <- paste("W=", Ws, sep="")
  rownames(sWcors[[k]]) <- paste("W=", Ws, sep="")
  for(i in 1:LVW){
    ssam[[k]][[i]] <- list()
    stimes[[k]][[i]] <- list()
    V <- Vs[i]
    for(j in 1:LVW){
      print(c(k,i,j))
      W <- Ws[j]
      b1 <- (a1-1)*V
      b2 <- (a2-1)*W
      start <- c(V,W)
      dat <- simdata[[k]][[i]][[j]]
      stateout <- statesam(n, start, dat, a1, a2, b1, b2, time=TRUE)
      stimes[[k]][[i]][[j]] <- stateout[[2]]
      state <- stateout[[1]][-c(1:burn),]
      ssam[[k]][[i]][[j]] <- state
      sTHcors[[k]][i,j,] <- diag(acf(state[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
      sVcors[[k]][i,j] <- acf(state[,T+2], plot=FALSE, lag.max=1)$acf[2]
      sWcors[[k]][i,j] <- acf(state[,T+3], plot=FALSE, lag.max=1)$acf[2]
    }
  }
}
scors <- list(sTHcors, sVcors, sWcors)
save(stimes, file="stimes.RData")
save(scors, file="scors.RData")
save(ssam, file="statsam.RData")
rm(ssam)



## next do the scaled disturbance sampler
set.seed(413013)
dVcors <- list()
dWcors <- list()
dTHcors <- list()
dsam <- list()
dtimes <- list()

for(k in 1:LT){
  T <- Ts[k]
  dVcors[[k]] <- matrix(0, LVW, LVW)
  dWcors[[k]] <- matrix(0, LVW, LVW)
  dTHcors[[k]] <- array(0, dim=c(LVW, LVW, T+1))
  dsam[[k]] <- list()
  dtimes[[k]] <- list()
  colnames(dVcors[[k]]) <- paste("V=", Vs, sep="")
  colnames(dWcors[[k]]) <- paste("V=", Vs, sep="")
  rownames(dVcors[[k]]) <- paste("W=", Ws, sep="")
  rownames(dWcors[[k]]) <- paste("W=", Ws, sep="")
  for(i in 1:LVW){
    dsam[[k]][[i]] <- list()
    dtimes[[k]][[i]] <- list()
    V <- Vs[i]
    for(j in 1:LVW){
      print(c(k,i,j))
      W <- Ws[j]
      b1 <- (a1-1)*V
      b2 <- (a2-1)*W
      start <- c(V,W)
      dat <- simdata[[k]][[i]][[j]]
      distout <- distsam(n, start, dat, a1, a2, b1, b2, time=TRUE)
      dtimes[[k]][[i]][[j]] <- distout[[2]]
      dist <- distout[[1]][-c(1:burn),]
      dsam[[k]][[i]][[j]] <- dist
      dTHcors[[k]][i,j,] <- diag(acf(dist[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
      dVcors[[k]][i,j] <- acf(dist[,T+2], plot=FALSE, lag.max=1)$acf[2]
      dWcors[[k]][i,j] <- acf(dist[,T+3], plot=FALSE, lag.max=1)$acf[2]
    }
  }
}

dcors <- list(dTHcors, dVcors, dWcors)
save(dcors, file="dcors.RData")
save(dsam, file="distsam.RData")
save(dtimes, file="dtimes.RData")
rm(dsam)

set.seed(413412)
eVcors <- list()
eWcors <- list()
eTHcors <- list()
esam <- list()
etimes <- list()

## finally do the scaled error sampler
for(k in 1:LT){
  T <- Ts[k]
  eVcors[[k]] <- matrix(0, LVW, LVW)
  eWcors[[k]] <- matrix(0, LVW, LVW)
  eTHcors[[k]] <- array(0, dim=c(LVW, LVW, T+1))
  esam[[k]] <- list()
  etimes[[k]] <- list()
  colnames(eVcors[[k]]) <- paste("V=", Vs, sep="")
  colnames(eWcors[[k]]) <- paste("V=", Vs, sep="")
  rownames(eVcors[[k]]) <- paste("W=", Ws, sep="")
  rownames(eWcors[[k]]) <- paste("W=", Ws, sep="")
  for(i in 1:LVW){
    esam[[k]][[i]] <- list()
    etimes[[k]][[i]] <- list()
    V <- Vs[i]
    for(j in 1:LVW){
      print(c(k,i,j))
      W <- Ws[j]
      b1 <- (a1-1)*V
      b2 <- (a2-1)*W
      start <- c(V,W)
      dat <- simdata[[k]][[i]][[j]]
      errout <- errorsam(n, start, dat, a1, a2, b1, b2, time=TRUE)
      etimes[[k]][[i]][[j]] <- errout[[2]]
      err <- errout[[1]][-c(1:burn),]
      esam[[k]][[i]][[j]] <- err
      eTHcors[[k]][i,j,] <- diag(acf(err[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
      eVcors[[k]][i,j] <- acf(err[,T+2], plot=FALSE, lag.max=1)$acf[2]
      eWcors[[k]][i,j] <- acf(err[,T+3], plot=FALSE, lag.max=1)$acf[2]
    }
  }
}

ecors <- list(eTHcors, eVcors, eWcors)
save(ecors, file="ecors.RData")
save(esam, file="errorsam.RData")
save(etimes, file="etimes.Rdata")
rm(esam)

