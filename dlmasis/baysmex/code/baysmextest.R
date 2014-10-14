source("baysmexfun.R")
hllmsim <- function(N,J,U,V,W,beta){
  theta <- matrix(0,nrow=N+1,ncol=J)
  dat <- matrix(0,nrow=N,ncol=J)
  u <- rnorm(N,0,sqrt(U))
  mu <- cumsum(c(0,u))
  for(t in 1:N){
    theta[t+1,] <- mu[t+1] + beta*theta[t,] + rnorm(J)*sqrt(W)
    dat[t,] <- theta[t+1,] + rnorm(J)*sqrt(V)
  }
  return(dat)
}

N <- 40
J <- 4
U <- 1
V <- rep(1,J)
W <- rep(10,J)
beta <- rep(.2,J)
testdat <- hllmsim(N,J,U,V,W,beta)


start <- c(U,V,W,beta,rep(0,N+1),rep(0,(N+1)*J))
au <- 5
av <- rep(5,J)
aw <- rep(5,J)
bu <- (au-1)*U
bv <- (av-1)*V
bw <- (aw-1)*W
V0 <- rep(100,J)
W0 <- rep(100,J)
m0 <- 0
U0 <- 100
beta0 <- 0
B <- 100
n <- 10000

testout <- naivegibbs(n, start, testdat, au, bu, av, bv, aw, bw, V0, W0, m0, U0, beta0, B)

Uout <- testout[,1]
Vout <- testout[,2:(J+1)]
Wout <- testout[,(J+2):(2*J+1)]
betaout <- testout[,(2*J+2):(3*J+1)]
muout <- testout[,(3*J+2):(3*J+2+N)]
thetaout <- matrix(testout[,(3*J+2+N+1):(3*J + (J+1)*(N+1)+1)], byrow=FALSE, ncol=J)





plot(ts(Uout))

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(Vout[,j]))
}

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(Wout[,j]))
}

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(betaout[,j]),ylim=c(-2,2))
}







plot(ts(cumsum(Uout)/1:n))

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(cumsum(Vout[,j])/(1:n)))
}

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(cumsum(Wout[,j])/(1:n)))
}

par(mfrow=c(2,2))
for(j in 1:J){
  plot(ts(cumsum(betaout[,j])/(1:n)), ylim=c(-3,3))
}


