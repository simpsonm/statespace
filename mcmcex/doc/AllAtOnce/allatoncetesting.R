source("../../mcmcexfun.R")

awolthsmooth <- function(y, V, W, m0, C0){
  n <- length(y)
  O0 <- 1/C0 + 1/W
  Ot <- 1/V + 2/W
  OT <- 1/V + 1/W
  Ott1 <- -1/W
  a <- rep(0,n+1)
  E <- a
  m <- a
  a[1] <- m0/C0
  E[1] <- 1/O0
  m[1] <- E[1]*a[1]
  for(t in 1:(n-1)){
    a[t+1] <- y[t]/V
    E[t+1] <- 1/(Ot - Ott1^2*E[t])
    m[t+1] <- E[t+1]*(a[t+1]-Ott1*m[t])
  }
  a[n+1] <- y[n]/V
  E[n+1] <- 1/(OT - Ott1^2*E[n])
  m[n+1] <- E[n+1]*(a[n+1]-Ott1*m[n])
  theta <- rnorm(n+1)
  theta[n+1] <- m[n+1] + sqrt(E[t+1])*theta[n+1]
  for(t in n:1){
    theta[t] <- m[t] - E[t]*Ott1*theta[t+1] + sqrt(E[t])*theta[t]
  }
  return(theta)
}

awolpssmooth <- function(y, V, W, m0, C0){
  n <- length(y)
  O0 <- 1/W + 1/C0
  Ot <- 2*V/W + 1
  OT <- V/W + 1
  sV <- sqrt(V)
  O01 <- sV/W
  Ost <- -V/W
  a <- rep(0,n+1)
  E <- a
  m <- a
  a[1] <- y[1]/W + m0/C0
  a[2] <- sV/W*(2*y[1] - y[2])
  E[1] <- 1/O0
  m[1] <- E[1]*a[1]
  E[2] <- 1/(Ot - O01^2*E[1])
  m[2] <- E[2]*(a[2] - O01*m[1])
  for(t in 2:n){
    a[t+1] <- sV/W*(2*y[t] - y[t-1] - y[t+1])
    E[t+1] <- 1/(Ot - Ost^2*E[t])
    m[t+1] <- E[t+1]*(a[t+1] - Ost*m[t])
  }
  a[n+1] <- sV/W*(y[n] - y[n-1])
  E[n+1] <- 1/(OT - Ost^2*E[n])
  m[n+1] <- E[n+1]*(a[n+1] - Ost*m[n])
  psi <- rnorm(n+1)
  psi[n+1] <- m[n+1] + sqrt(E[t+1])*psi[n+1]
  for(t in n:2){
    psi[t] <- m[t] - E[t]*Ost*psi[t+1] + sqrt(E[t])*psi[t]
  }
  psi[1] <- m[1] - E[1]*O01*psi[2] + sqrt(E[1])*psi[1]
  return(psi)
}

awolpssmooth2 <- function(y, V, W, m0, C0){
  n <- length(y)
  O0 <- 1/W + 1/C0
  Ot <- 2*V/W + 1
  OT <- V/W + 1
  sV <- sqrt(V)
  O01 <- sV/W
  Ost <- -V/W
  a <- rep(0,n+1)
  a[1] <- y[1]/W + m0/C0
  a[2] <- sV/W*(2*y[1] - y[2])
   for(t in 2:n){
    a[t+1] <- sV/W*(2*y[t] - y[t-1] - y[t+1])
   }
  a[n+1] <- sV/W*(y[n] - y[n-1])
  psi <- rnorm(n+1)
  O <- diag(c(O0, rep(Ot, n-1), OT))
  diag(O[-1,-(n+1)]) <- Ost
  diag(O[-(n+1),-1]) <- Ost
  O[1,2] <- O01
  O[2,1] <- O01
  E <- solve(O)
  K <- chol(E)
  psi <- c(E%*%matrix(a, ncol=1) + t(K)%*%matrix(psi, ncol=1))
  return(psi)
}

V <- 1
W <- 100
n <- 20
dat <- rnorm(n, 0, sqrt(V)) + cumsum(rnorm(n, 0, sqrt(W)))
m0 <- 0
C0 <- 10000

k <- 10000
psitest <- matrix(0,nrow=k, ncol=n+1)
thetest <- psitest
control <- psitest
psitest2 <- psitest

for(i in 1:k){
  psi <- awolpssmooth(dat, V, W, m0, C0)
  theta <- awolthsmooth(dat, V, W, m0, C0)
  mod <- dlmModPoly(order=1, dV=V, dW=W)
  filt <- dlmFilter(dat, mod)
  thetatrue <- dlmBSample(filt)
  control[i,] <- thetatrue
  thetest[i,] <- theta
  thetapsi <- thetapsitrans(dat, psi, V)
  psitest[i,] <- thetapsi
  psi2 <- awolpssmooth2(dat, V, W, m0, C0)
  thetapsi2 <- thetapsitrans(dat, psi2, V)
  psitest2[i,] <- thetapsi2
}

mns <- cbind(summary(mcmc(control))[[1]][,1],summary(mcmc(thetest))[[1]][,1],summary(mcmc(psitest))[[1]][,1],summary(mcmc(psitest2))[[1]][,1])
colnames(mns) <- c("control", "thetaAWOL", "psiAWOL", "psiAWOL2")
mns

s <- 1
par(mfrow=c(2,2))
plot(1:k, cumsum(control[,s])/c(1:k), type="l")
plot(1:k, cumsum(thetest[,s])/c(1:k), type="l", col="red")
plot(1:k, cumsum(psitest[,s])/c(1:k), type="l", col="blue")
plot(1:k, cumsum(psitest2[,s])/c(1:k), type="l", col="green")

Ss <- 1:9
par(mfrow=c(3,3))
for(s in Ss){
  plot(1:k, cumsum(control[,s])/c(1:k), type="l")
  lines(1:k, cumsum(thetest[,s])/c(1:k), type="l", col="red")
  lines(1:k, cumsum(psitest[,s])/c(1:k), type="l", col="blue")
  lines(1:k, cumsum(psitest2[,s])/c(1:k), type="l", col="green")
}



    
