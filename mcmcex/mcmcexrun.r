
## a wrapper for simulating from any of the interweaving/alternating
## functions quickly in a loop

source("mcmcexfun.r")



set.seed(152893627)
simdata <- data.frame(y=NULL, t=NULL, V=NULL, W=NULL, T=NULL)
for(k in 1:LT){
  T <- Ts[k]
  for(i in 1:LVW){
    V <- Vs[i]
    for(j in 1:LVW){
      W <- Ws[j]
      y <- llsim(T, V, W, 0, 1)
      newsim <- data.frame(y=y, t=1:length(y), V.T=V, W.T=W, T.T=T)
      simdata <- merge(simdata, newsim, all=TRUE)
    }
  }
}

simdata1 <- data.frame(simdata, ch=1)
simdata2 <- data.frame(simdata, ch=2)
simdata3 <- data.frame(simdata, ch=3)
simdat <- merge(merge(simdata1, simdata2, all=TRUE), simdata3, all=TRUE)




sams <- c("state", "dist", "error", "sdint", "seint", "deint", "triint", "sdalt", "sealt", "dealt", "trialt")
sam <- ddply(simdata, .(V.T, W.T, T.T), samwrap, n=5, a1=5, a2=5,
b1=5, b2=5, samp="state", ch=1)






