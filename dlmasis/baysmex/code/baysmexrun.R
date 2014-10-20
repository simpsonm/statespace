source("baysmexfun.R")

library(plyr)
library(reshape2)
library(coda)
effdat <- read.csv("Efficiency.csv")[,-1]
effmelt <- melt(effdat, id=c("period", "rep", "trt", "id"))
logit <- function(x){
  return(log(x/(1-x)))
}
KStealdat <- logit(t(dcast(effmelt, formula=period~rep, subset=.(trt=="KSteal"))[,-1]))
J <- nrow(KStealdat)
N <- ncol(KStealdat)

C <- 5 ## number of chains
n <- 20000 ## number of sims per chain

start <- matrix(0, nrow=C, ncol=2*J+1)
start[1,] <- .5
for(i in 2:5){
  start[i,] <- start[1,]*exp(rnorm(2*J+1, sd=5))
}

au <- 1.5
av <- rep(1.5,J)
aw <- rep(1.5,J)
bu <- (au-1)*.5
bv <- (av-1)*.5
bw <- (aw-1)*.5
m0 <- rep(0,J+1)
C0 <- diag(100,J+1)

gibbs <- mcmc.list()
inter <- mcmc.list()
gibbstime <- matrix(0,nrow=C,ncol=5)
intertime <- gibbstime

for(i in 1:C){
  print(i)
  gibbstime[i,] <- system.time(gibbs[[i]] <- mcmc(naivegibbs(n, start[i,], KStealdat, au, bu, av, bv, aw, bw, m0, C0)))
  print("gibbs done")
  print(gibbstime)
  intertime[i,] <- system.time(inter[[i]] <- mcmc(disterrorinter(n, start[i,], KStealdat, au, bu, av, bv, aw, bw, m0, C0)))
  print("inter done")
  print(intertime)
}

save(gibbs, file="gibbs.R")
save(inter, file="inter.R")
save(gibbstime, file="gibbstime.R")
save(intertime, file="intertime.R")
