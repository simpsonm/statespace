source("baysmexfun.R")
hllmsim <- function(N,J,U,V,W){
  dat <- matrix(0,nrow=J, ncol=N)
  u <- rnorm(N,0,sqrt(U))
  mu <- cumsum(c(0,u))
  for(j in 1:J){
    theta <- cumsum(c(0,rnorm(N, 0, sqrt(W[j]))))
    dat[j,] <- theta[-1] + mu[-1] + rnorm(N,0,sqrt(V[j]))
  }
  return(dat)
}


N <- 20
J <- 1
U <- 1
V <- rep(1,J)
W <- rep(1,J)
testdat <- hllmsim(N,J,U,V,W)

start <- c(V,U,W)
Qu <- 100
Qv <- rep(100,J)
Qw <- rep(100,J)
m0 <- rep(0,J+1)
C0 <- diag(100,J+1)
n <- 2000

par(mfrow=c(J,1))
for(j in 1:J){
  plot(ts(testdat[j,]))
}



set.seed(235)
gibbsout <- naivegibbs(n, start, testdat, Qv, Qu, Qw, m0, C0)

interout <- disterrorinter(n, start, testdat, Qv, Qu, Qw, m0, C0)


R <- 5
gl <- list()
il <- list()
system.time(
    for(r in 1:R){
      print(r)
      ##testdat <- hllmsim(N,J,U,V,W)
      newstart <- start*exp(rnorm(length(start)))
      gibbsout <- naivegibbs(n, newstart, testdat, Qv, Qu, Qw, m0, C0)
      gl[[r]] <- gibbsout[,1:(2*J+2)]
##      interout <- disterrorinter(n, newstart, testdat, Qv, Qu, Qw, m0, C0)
##      il[[r]] <- interout[,1:(2*J+2)]
      print("finished")
    }
    )



m <- 1
ll <- c(rep(0,2*J+1),-5)
ul <- c(c(V,U,W)*3,5)

par(mfrow=c(J+1,2))
for(j in 1:(2*J + 2)){
  plot(ts(cumsum(gl[[1]][-c(1:m),j])/1:(n-m)), ylab=colnames(gl[[1]])[j], ylim=c(ll[j],ul[j]))
##  lines(1:(n-m), cumsum(il[[1]][-c(1:m),j])/1:(n-m), col="red")
  for(r in 2:R){
    lines(1:(n-m), cumsum(gl[[r]][-c(1:m),j])/1:(n-m))
##    lines(1:(n-m), cumsum(il[[r]][-c(1:m),j])/1:(n-m), col="red")
  }
}

r <- 3
par(mfrow=c(2*J+1,1))
for(j in 1:(2*J + 1)){
##  plot(ts(il[[r]][,j]),ylab=colnames(il[[r]])[j])
  plot(ts(gl[[r]][,j]),ylab=colnames(gl[[r]])[j],col="red")
}

r <- 1
par(mfrow=c(2*J+1,2))
for(j in 1:(2*J + 1)){
  plot(ts(cumsum(il[[r]][,j])/1:n),ylab=colnames(il[[r]])[j])
  plot(ts(cumsum(gl[[r]][,j])/1:n),ylab=colnames(gl[[r]])[j],col="red")
}



r <- 3
par(mfrow=c(4,2))
j <- 1
plot(ts(il[[1]][,j]),ylab=colnames(il[[r]])[j])
plot(ts(gl[[1]][,j]),ylab=colnames(gl[[r]])[j],col="red")
j <- j+J+1
plot(ts(il[[1]][,j]),ylab=colnames(il[[r]])[j])
plot(ts(gl[[1]][,j]),ylab=colnames(gl[[r]])[j],col="red")
j <- 2
plot(ts(il[[1]][,j]),ylab=colnames(il[[r]])[j])
plot(ts(gl[[1]][,j]),ylab=colnames(gl[[r]])[j],col="red")
j <- j+J+1
plot(ts(il[[1]][,j]),ylab=colnames(il[[r]])[j])
plot(ts(gl[[1]][,j]),ylab=colnames(gl[[r]])[j],col="red")



j <- 2
par(mfrow=c(R,2))
for(r in 1:R){
  plot(ts(cumsum(il[[r]][,j])/1:n),ylab=colnames(il[[r]])[j])
  plot(ts(cumsum(gl[[r]][,j])/1:n),ylab=colnames(gl[[r]])[j],col="red")
}



library(plyr)
library(reshape2)
effdat <- read.csv("Efficiency.csv")[,-1]
effmelt <- melt(effdat, id=c("period", "rep", "trt", "id"))
KStealdat <- dcast(effmelt, formula=period~rep, subset=.(trt=="KSteal"))[,-1]



###################################
## For looking at the chains

par(mfrow=c(2,1))
j <- 2
for(k in 1:1){
  traceplot(gibbs[[j]][,k],ylim=c(0,1), ylab="gibbs")
  traceplot(inter[[j]][,k],ylim=c(0,1), ylab="inter")
}

par(mfrow=c(2,1))
j <- 2
for(k in 1:1){
  traceplot(gibbs[,k,drop=FALSE],ylim=c(0,1), ylab="gibbs")
  traceplot(inter[,k,drop=FALSE],ylim=c(0,1), ylab="inter")
}


effectiveSize(gibbs[,1:13,drop=FALSE])/sum(gibbstime[,3])
effectiveSize(inter[,1:13,drop=FALSE])/sum(intertime[,3])

