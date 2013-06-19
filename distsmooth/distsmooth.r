## Code snippets for implementing an algorithm based on Koopman (1993) for simulating
## the system disturbances - an alternative to using FFBS to sample the states,
## then transforming. Only implemented for local level model. Not interesting at the
## moment because it can't help convergence, but ultimately might be interesting for
## computational speed.

library(dlm)
library(coda)

## simulate from a local level model
V <- 2
W <- 3
T <- 10
theta0 <- 0
dat <- rep(0,T)
theta <- c(theta0, rep(0,T))
for(t in 1:T){
  theta[t+1] <- theta[t] + rnorm(1, 0, sd=sqrt(W))
  dat[t] <- theta[t+1] + rnorm(1, 0, sd=sqrt(V))
}

## run the kalman filter
mod <- dlmModPoly(1, V, W, 0, 1)
filt <- dlmFilter(dat, mod)

##############################################################################
## From here we're basically recursively finding the means, variances, and
## covariances of the system disturbances, conditional on the data
##############################################################################
v <- dat - filt$f
F <- rep(0,T)
P <- F
K <- F
for(t in 1:T){
  U.R <- filt$U.R[[t]]
  D.R <- filt$D.R[t,]
  P[t] <-  U.R^2*D.R^2
  F[t] <- P[t] + V
  K[t] <- P[t]/F[t]
}

r <- rep(0,T+1)
e <- rep(0,T)
D <- e
N <- r
VHu <- e
for(t in T:1){
  e[t] <- v[t]/F[t] - K[t]*r[t+1]
  r[t] <- e[t] + r[t+1]
  D[t] <-  1/F[t] + N[t+1]*K[t]^2
  N[t] <- D[t]-2*N[t+1]*K[t] + N[t+1]
}

L <- 1 - K
Ls <- matrix(0,T,T)
for(t in 1:T){
  for(s in 1:t){
    if(t==s){
      Ls[t,s] <- 1
    }
    else{
      Ls[t,s] <- prod(L[(t-1):(s)])
    }
  }
}

VHu <- W*(1-W*N[1:T])


VSig <- diag(VHu)
for(t in 2:T){
  for(s in 1:(t-1)){
    VSig[t,s] <- -W^2*N[t]*(1-K[t-1])*Ls[t-1,s]
    VSig[s,t] <- VSig[t,s]
  }
}

###############################################################
## This finishes finding all of that information

###############################################################
## This is me toying round with the results
mns <- W*r
CK <- chol(VSig)
esam <- mns[1:T] + t(CK)%*%rnorm(T)

thetasmooth <- dlmSmooth(dat, mod)$s

thetasmooth[-1]-thetasmooth[1:T]
mns[1:T]

testsam <- dlmBSample(filt)
###############################################################

## An extremely large sample test to see if the sampler works by
## checking if it agrees with FFBS
I <- 1000000
sam1 <- matrix(0,ncol=T+1,nrow=I)
sam2 <- sam1
for(i in 1:I){
  print(i)

  ## first sample the disturbances
  tempa <- mns[1:T] + t(CK)%*%rnorm(T)

  ## next sample the initial state (not in Koopman, but easy)
  ## then put it together
  mn <- T/V*mean(dat-cumsum(tempa)) + mod$m0/mod$C0
  vr <- 1/(1/mod$C0 + T/V)
  mn <- vr*mn
  sam1[i,] <- c(rnorm(1,mn,sqrt(vr)),tempa)

  ## now sample the states using FFBS and transform into
  ## the disturbances
  temp <- dlmBSample(filt)
  sam2[i,] <- c(temp[1],temp[-1]-temp[1:T])
}

## some diagnostic stuff to look at the difference between the two sampling methods
## agree
sum1 <- summary(mcmc(sam1))
sum2 <- summary(mcmc(sam2))

sum1[[1]]
sum2[[1]]

par(mfrow=c(1,2))
plot(density(sam1[,1]), xlim=c(-5,5), ylim=c(0,0.5))
plot(density(sam2[,1]), xlim=c(-5,5), ylim=c(0,0.5))

par(mfrow=c(1,2))
hist(sam1[,1], freq=FALSE)
plot(sam2[,1]), xlim=c(-5,5), ylim=c(0,0.5))

ks.test(sam1[,1],sam2[,1])

cbind(sum1[[1]][,1:2],sum2[[1]][,1:2])[,c(1,3,2,4)]
(sum1[[1]]-sum2[[1]])/((sum1[[1]] + sum2[[1]])/2)
##WORKING!!!!
