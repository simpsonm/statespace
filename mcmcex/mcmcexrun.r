source("mcmcexfun.r")
set.seed(152893627)
simdata <- data.frame(y=NULL, t=NULL, V.T=NULL, W.T=NULL, T.T=NULL)
Ts <- c(10, 100, 1000)
Vs <- 10^(c(0:10)/2-2)
Ws <- Vs
LT <- length(Ts)
LVW <- length(Vs)
for(k in 1:LT){
  T <- Ts[k]
  for(i in 1:LVW){
    V <- Vs[i]
    for(j in 1:LVW){
      W <- Ws[j]
      y <- llsim(T, V, W, 0, 1)
      newsim <- data.frame(y=y, t=1:length(y), V.T=V, W.T=W, T.T=T)
      simdata <- rbind(simdata, newsim)
    }
  }
}
simdata1 <- data.frame(simdata, ch=1)
simdata2 <- data.frame(simdata, ch=2)
simdata3 <- data.frame(simdata, ch=3)
simdat <- rbind(simdata1, simdata2, simdata3)
sams <- c("state", "dist", "error", "sdint", "seint", "deint",
          "triint", "sdalt", "sealt", "dealt", "trialt")
samplers <- data.frame(sams=rep(1,length(sams)))
samplers$sams <- sams

system.time(test <- fullsim(samplers, simdata1, 4, 1, 5, 5))
save(test, file="summary.RData")

