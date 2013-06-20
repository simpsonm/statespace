## Code snippets for testing the convergence of the samplers
## mainly for plotting ergodic means

source("mcmcexfun.r")
n <- 1000
chains <- 3
T <- 10
a1 <- 5
a2 <- 5
V <- 1
W <- 10
starts <- c(.01, 1, 100)
startV <- V*starts
startW <- W*starts
b1 <- (a1-1)*V
b2 <- (a2-1)*W
dat <- llsim(T, V, W, 0, 1)
statelist <- list()
distlist <- list()
errorlist <- list()

for(j in 1:chains){
  print(j)
  start <- c(startV[j], startW[j])
  dist <- distsam(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
  state <- statesam(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
  error <- errorsam(n, start, dat, a1=a1, a2=a2, b1=b1, b2=b2)
  distlist[[j]] <- dist
  statelist[[j]] <- state
  errorlist[[j]] <- error
}

##plot ergodic means from multiple chains
par(mfrow=c(2,1))
plot(ts(cumsum(distlist[[1]][,T+2])/1:n), main="V", ylim=c(V*.5,V/.5))
for(i in 2:chains){
  lines(ts(cumsum(distlist[[i]][,T+2])/1:n))
}
for(i in 1:chains){
  lines(ts(cumsum(statelist[[i]][,T+2])/1:n), col="red")
}
for(i in 1:chains){
  lines(ts(cumsum(errorlist[[i]][,T+2])/1:n), col="blue")
}
legend("topright", c("dist", "state", "error"), lty=c(1,1,1), col=c("black", "red", "blue"))
plot(ts(cumsum(distlist[[1]][,T+3])/1:n), main="W", ylim=c(.5*W,W/.5))
for(i in 2:chains){
  lines(ts(cumsum(distlist[[i]][,T+3])/1:n))
}
for(i in 1:chains){
  lines(ts(cumsum(statelist[[i]][,T+3])/1:n), col="red")
}
for(i in 1:chains){
  lines(ts(cumsum(errorlist[[i]][,T+3])/1:n), col="blue")
}
legend("topright", c("dist", "state", "error"), lty=c(1,1,1), col=c("black", "red", "blue"))
