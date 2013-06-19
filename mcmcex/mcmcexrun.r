## Simulates data from 36 local level models at 3 different lengths of the time series
## then fits a local level model to each one using both the state sampler and the scaled
## disturbance sampler. The samples and the autocorrelations from each chain are stored.

## This takes about 18-24 hours to run

source("mcmcexfun.r")
set.seed(15236289)
n <- 5000
burn <- 500
Ts <- c(10,100,500)
LT <- length(Ts)
a1 <- 5
a2 <- 5
Vs <- c(.01, .1, 1, 10, 100, 1000)
LVW <- length(Vs)
Ws <- Vs
dVcors <- list()
dWcors <- list()
dTHcors <- list()
sVcors <- list()
sWcors <- list()
sTHcors <- list()
dsam <- list()
ssam <- list()

for(k in 1:LT){
  T <- Ts[k]
  dVcors[[k]] <- matrix(0, LVW, LVW)
  dWcors[[k]] <- matrix(0, LVW, LVW)
  sVcors[[k]] <- matrix(0, LVW, LVW)
  sWcors[[k]] <- matrix(0, LVW, LVW)
  dTHcors[[k]] <- array(0, dim=c(LVW, LVW, T+1))
  sTHcors[[k]] <- array(0, dim=c(LVW, LVW, T+1))
  dsam[[k]] <- list()
  ssam[[k]] <- list()
  colnames(sVcors[[k]]) <- paste("V=", Vs, sep="")
  colnames(sWcors[[k]]) <- paste("V=", Vs, sep="")
  rownames(sVcors[[k]]) <- paste("W=", Ws, sep="")
  rownames(sWcors[[k]]) <- paste("W=", Ws, sep="")
  colnames(dVcors[[k]]) <- paste("V=", Vs, sep="")
  colnames(dWcors[[k]]) <- paste("V=", Vs, sep="")
  rownames(dVcors[[k]]) <- paste("W=", Ws, sep="")
  rownames(dWcors[[k]]) <- paste("W=", Ws, sep="")
  for(i in 1:LVW){
    dsam[[k]][[i]] <- list()
    ssam[[k]][[i]] <- list()
    for(j in 1:LVW){
      print(c(k,i,j))
      V <- Vs[i]
      W <- Ws[j]
      b1 <- (a1-1)*V
      b2 <- (a2-1)*W
      start <- c(V,W)
      dat <- llsim(T, V, W, 0, 1)
      sim <- statesam(n, start, dat, a1, a2, b1, b2)
      sim <- sim[-c(1:burn),]
      dist <- distsam(n, start, dat, a1, a2, b1, b2)
      dist <- dist[-c(1:burn),]
      dsam[[k]][[i]][[j]] <- dist
      ssam[[k]][[i]][[j]] <- sim
      sTHcors[[k]][i,j,] <- diag(acf(sim[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
      sVcors[[k]][i,j] <- acf(sim[,T+2], plot=FALSE, lag.max=1)$acf[2]
      sWcors[[k]][i,j] <- acf(sim[,T+3], plot=FALSE, lag.max=1)$acf[2]
      dTHcors[[k]][i,j,] <- diag(acf(dist[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
      dVcors[[k]][i,j] <- acf(dist[,T+2], plot=FALSE, lag.max=1)$acf[2]
      dWcors[[k]][i,j] <- acf(dist[,T+3], plot=FALSE, lag.max=1)$acf[2]
    }
  }
}

sims <- list(ssam, dsam)
cors <- list(sTHcors, sVcors, sWcors, dTHcors, dVcors, dWcors)
save(cors, file="cors.RData")
memory.limit(4095)  # needed to have enough space to save the draws
save(ssam, file="simsam.RData")
save(dsam, file="distsam.RData")

