## Simulates data from 9 local level models at 3 different lengths of the time series
## then fits a local level model to each one using both the state sampler and the scaled
## disturbance sampler with 3 different sets of starting values. The samples and the
## autocorrelations from each chain are stored.

## takes around 40 hours to run. Ick.

## setup and data simulation
source("../mcmcexfun.r")
n <- 2000
Ts <- c(10,100,1000)
LT <- length(Ts)
a1 <- 5
a2 <- 5
Vs <- c(1, 10, 100)
LVW <- length(Vs)
Ws <- Vs
starts <- c(1/100, 1, 100)
chains <- length(starts)

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

## first do the state sampler
set.seed(1890)
for(ch in 1:chains){
  sVcors[[ch]] <- list()
  sWcors[[ch]] <- list()
  sTHcors[[ch]] <- list()
  ssam[[ch]] <- list()
  for(k in 1:LT){
    T <- Ts[k]
    sVcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    sWcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    sTHcors[[ch]][[k]] <- array(0, dim=c(LVW, LVW, T+1))
    ssam[[ch]][[k]] <- list()
    colnames(sVcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    colnames(sWcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    rownames(sVcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    rownames(sWcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    for(i in 1:LVW){
      ssam[[ch]][[k]][[i]] <- list()
      V <- Vs[i]
      for(j in 1:LVW){
        print(c(ch,k,i,j))
        W <- Ws[j]
        b1 <- (a1-1)*V
        b2 <- (a2-1)*W
        start <- c(V,W)*starts[ch]
        dat <- simdata[[k]][[i]][[j]]
        state <- statesam(n, start, dat, a1, a2, b1, b2)
        ssam[[ch]][[k]][[i]][[j]] <- state
        sTHcors[[ch]][[k]][i,j,] <- diag(acf(state[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
        sVcors[[ch]][[k]][i,j] <- acf(state[,T+2], plot=FALSE, lag.max=1)$acf[2]
        sWcors[[ch]][[k]][i,j] <- acf(state[,T+3], plot=FALSE, lag.max=1)$acf[2]
      }
    }
  }
}
scors <- list(sTHcors, sVcors, sWcors)
save(scors, file="scors.RData")
save(ssam, file="statsam.RData")
rm(ssam)



## next do the scaled disturbance sampler
dVcors <- list()
dWcors <- list()
dTHcors <- list()
dsam <- list()

set.seed(09142)
for(ch in 1:chains){
  dVcors[[ch]] <- list()
  dWcors[[ch]] <- list()
  dTHcors[[ch]] <- list()
  dsam[[ch]] <- list()
  for(k in 1:LT){
    T <- Ts[k]
    dVcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    dWcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    dTHcors[[ch]][[k]] <- array(0, dim=c(LVW, LVW, T+1))
    dsam[[ch]][[k]] <- list()
    colnames(dVcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    colnames(dWcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    rownames(dVcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    rownames(dWcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    for(i in 1:LVW){
      dsam[[ch]][[k]][[i]] <- list()
      V <- Vs[i]
      for(j in 1:LVW){
        print(c(ch,k,i,j))
        W <- Ws[j]
        b1 <- (a1-1)*V
        b2 <- (a2-1)*W
        start <- c(V,W)*starts[ch]
        dat <- simdata[[k]][[i]][[j]]
        dist <- distsam(n, start, dat, a1, a2, b1, b2)
        dsam[[ch]][[k]][[i]][[j]] <- dist
        dTHcors[[ch]][[k]][i,j,] <- diag(acf(dist[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
        dVcors[[ch]][[k]][i,j] <- acf(dist[,T+2], plot=FALSE, lag.max=1)$acf[2]
        dWcors[[ch]][[k]][i,j] <- acf(dist[,T+3], plot=FALSE, lag.max=1)$acf[2]
      }
    }
  }
}
dcors <- list(dTHcors, dVcors, dWcors)
save(dcors, file="dcors.RData")
save(dsam, file="distsam.RData")
rm(dsam)


eVcors <- list()
eWcors <- list()
eTHcors <- list()
esam <- list()


## finally do the scaled error sampler
set.seed(423142)
for(ch in 1:chains){
  eVcors[[ch]] <- list()
  eWcors[[ch]] <- list()
  eTHcors[[ch]] <- list()
  esam[[ch]] <- list()
  for(k in 1:LT){
    T <- Ts[k]
    eVcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    eWcors[[ch]][[k]] <- matrix(0, LVW, LVW)
    eTHcors[[ch]][[k]] <- array(0, dim=c(LVW, LVW, T+1))
    esam[[ch]][[k]] <- list()
    colnames(eVcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    colnames(eWcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
    rownames(eVcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    rownames(eWcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
    for(i in 1:LVW){
      esam[[ch]][[k]][[i]] <- list()
      V <- Vs[i]
      for(j in 1:LVW){
        print(c(ch,k,i,j))
        W <- Ws[j]
        b1 <- (a1-1)*V
        b2 <- (a2-1)*W
        start <- c(V,W)*starts[ch]
        dat <- simdata[[k]][[i]][[j]]
        err <- errorsam(n, start, dat, a1, a2, b1, b2)
        esam[[ch]][[k]][[i]][[j]] <- err
        eTHcors[[ch]][[k]][i,j,] <- diag(acf(err[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
        eVcors[[ch]][[k]][i,j] <- acf(err[,T+2], plot=FALSE, lag.max=1)$acf[2]
        eWcors[[ch]][[k]][i,j] <- acf(err[,T+3], plot=FALSE, lag.max=1)$acf[2]
      }
    }
  }
}

ecors <- list(eTHcors, eVcors, eWcors)
save(ecors, file="ecors.RData")
save(esam, file="errorsam.RData")
rm(esam)

