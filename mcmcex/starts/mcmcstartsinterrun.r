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

truth <- c(TRUE, FALSE)

cornam <- matrix("", nrow=4, ncol=2)
samnam <- matrix("", nrow=4, ncol=2)
samps <- c("sd", "se", "de", "tr")
ints <- c("i", "a")
cend <- "cors.RData"
send <- "sam.RData"

for(samp in 1:4){
  for(int in 1:2){
    tnam <- paste(samps[samp], ints[int], sep="")
    cornam[samp,int] <- paste(tnam, cend, sep="")
    samnam[samp,int] <- paste(tnam, send, sep="")
  }
}

## first do the state sampler
set.seed(1890)
for(samp in 1:4){
  for(int in 1:2){
    Vcors <- list()
    Wcors <- list()
    THcors <- list()
    sam <- list()
    for(ch in 1:chains){
      Vcors[[ch]] <- list()
      Wcors[[ch]] <- list()
      THcors[[ch]] <- list()
      sam[[ch]] <- list()
      for(k in 1:LT){
        T <- Ts[k]
        Vcors[[ch]][[k]] <- matrix(0, LVW, LVW)
        Wcors[[ch]][[k]] <- matrix(0, LVW, LVW)
        THcors[[ch]][[k]] <- array(0, dim=c(LVW, LVW, T+1))
        sam[[ch]][[k]] <- list()
        colnames(Vcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
        colnames(Wcors[[ch]][[k]]) <- paste("V=", Vs, sep="")
        rownames(Vcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
        rownames(Wcors[[ch]][[k]]) <- paste("W=", Ws, sep="")
        for(i in 1:LVW){
          sam[[ch]][[k]][[i]] <- list()
          V <- Vs[i]
          for(j in 1:LVW){
            print(c(ch,k,i,j))
            W <- Ws[j]
            b1 <- (a1-1)*V
            b2 <- (a2-1)*W
            start <- c(V,W)*starts[ch]
            dat <- simdata[[k]][[i]][[j]]
            tempsam <- samwrapper(n, start, dat, a1, a2, b1, b2, truth, samp)
            sam[[ch]][[k]][[i]][[j]] <- tempsam
            THcors[[ch]][[k]][i,j,] <- diag(acf(tempsam[,1:(T+1)], plot=FALSE, lag.max=1)$acf[2,,])
            Vcors[[ch]][[k]][i,j] <- acf(tempsam[,T+2], plot=FALSE, lag.max=1)$acf[2]
            Wcors[[ch]][[k]][i,j] <- acf(tempsam[,T+3], plot=FALSE, lag.max=1)$acf[2]
          }
        }
      }
    }
    cors <- list(THcors, Vcors, Wcors)
    save(cors, file=cornam[samp, int])
    save(sam, file=samnam[samp, int])
    rm(sam)
  }
}



