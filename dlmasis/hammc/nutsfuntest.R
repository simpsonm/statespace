source("nutsfun.R")

lpr <- function(theta){
  Sigma <- matrix(0,2,2)
  sig <- c(1,1)
  rho <- 0
  Sigma[1,1] <- sig[1]
  Sigma[2,2] <- sig[2]
  Sigma[1,2] <- sqrt(sig[1]*sig[2]*rho)
  Sigma[2,1] <- Sigma[1,2]
  mu <- c(0,0)
  diff <- t(c(theta-mu))
  Siginv <- solve(Sigma)
  detSig <- det(Sigma)
  out <- - log(detSig)/2 - diff%*%Siginv%*%t(diff)/2
  return(out)
}

gradlpr <- function(theta){
  Sigma <- matrix(0,2,2)
  sig <- c(1,1)
  rho <- 0
  Sigma[1,1] <- sig[1]
  Sigma[2,2] <- sig[2]
  Sigma[1,2] <- sqrt(sig[1]*sig[2]*rho)
  Sigma[2,1] <- Sigma[1,2]
  mu <- c(0,0)
  diff <- t(c(theta-mu))
  Siginv <- solve(Sigma)
  out <- -diff%*%Siginv
  return(out)
}


M <- 10
Madapt <- 5
theta0 <- c(0,0)
delta <- 0.6
delmax <- 1000

test <- NUTSda(M, Madapt, theta0, delta, lpr, gradlpr, delmax)
