###################################
## For looking at the chains

library(plyr)
library(reshape2)
library(coda)
library(ggplot2)
effdat <- read.csv("Efficiency.csv")[,-1]
effmelt <- melt(effdat, id=c("period", "rep", "trt", "id"))
logit <- function(x){
  return(log(x/(1-x)))
}
KStealdat <- logit(t(dcast(effmelt, formula=period~rep, subset=.(trt=="KSteal"))[,-1]))
J <- nrow(KStealdat)
N <- ncol(KStealdat)


load("gibbs.R")
load("gibbstime.R")
load("inter.R")
load("intertime.R")


gelman.plot(gibbs[,1:13,drop=FALSE], ask=TRUE)
gelman.plot(inter[,1:13,drop=FALSE], ask=TRUE)


allchainsburn <- mcmc.list()
for(i in 1:10){
  if(i < 6){
    all <- gibbs[[i]][-c(1:5000),]
  }
  if(i > 5){
    all <- inter[[i-5]][-c(1:5000),]
  }
  mus <- all[,14:(13+36)]
  thets <- all[,-c(1:(13+36))]
  for(j in 1:J){
    for(t in 1:36){
      thets[,(t + (j-1)*36)] <- thets[,(t + (j-1)*36)] + mus[,t]
    }
  }
  all[,-c(1:(13+36))] <- thets
  allchainsburn[[i]] <- mcmc(all) 
}

allsum <- summary(allchainsburn)

allmns <- allsum[[1]]
allquants <- allsum[[2]]
parmns <- allmns[1:13,]
parquants <- allquants[1:13,]
statequants <- allquants[-c(1:13),]
statemeds <- matrix(statequants[,3], byrow=FALSE, nrow=36)
colnames(statemeds) <- c("mu_t", paste("theta_", 1:6, "t", sep=""))
statemns <- matrix(allmns[-c(1:13),1], byrow=FALSE, nrow=36)
colnames(statemns) <- c("mu_t", paste("theta_", 1:6, "t", sep=""))

thetmeds <- data.frame(period=rep(0:35,6), phi=c(statemeds[,-1]), rep=factor(paste("Replication", rep(1:6,each=36), sep=" ")))
mumeds <- data.frame(period=0:35, mu=statemeds[,1])
kstealdatf <- data.frame(period=rep(1:35,6), leff=c(t(KStealdat)), rep=factor(paste("Replication", rep(1:6,each=35), sep=" ")))


statell <- matrix(statequants[,1], byrow=FALSE, nrow=36)
stateul <- matrix(statequants[,5], byrow=FALSE, nrow=36)
mucis <- data.frame(period=0:35, ll=statell[,1], ul=stateul[,1], mu=mumeds$mu)
thetcis <- data.frame(period=rep(0:35,6), rep=rep(1:6,each=36), ll=c(statell[,-1]), ul=c(stateul[,-1]), phi=thetmeds$phi)



phiplot <- ggplot(data=thetmeds, aes(x=period, y=phi)) + geom_line(aes(size="phi")) + facet_wrap(~rep, ncol=2) + geom_line(data=kstealdatf, aes(x=period, y=leff, size="data")) + geom_line(data=mumeds, aes(x=period, y=mu, size="mu"), lty=2) + scale_size_manual(values=c(0.25, 0.5, 0.5), labels=c(expression(y[j][,][t]), expression(mu[t]), expression(phi[j][,][t])), guide=guide_legend(title=NULL, override.aes=list(size=c(0.1,0.5,0.5), linetype=c(1,2,1)))) + ylab("logit efficiency")

mudatplot <- qplot(period, mu, data=mumeds, geom="line", lwd=I(1)) + geom_line(data=kstealdatf, aes(x=period, y=leff, lty=rep),lwd=0.25)
muphiplot <- qplot(period, mu, data=mumeds, geom="line", lwd=I(1)) + geom_line(data=thetmeds, aes(x=period, y=phi, lty=rep),lwd=0.25)

ggsave("phiplot.png", phiplot, width=10, height=5)

ggsave("phiplot.pdf", phiplot, width=10, height=5)
##ggsave("mudatplot.pdf", mudatplot)
##ggsave("muphiplot.pdf", muphiplot)


gibbsburn <- mcmc.list()
interburn <- mcmc.list()
for(i in 1:5){
  gibbsburn[[i]] <- mcmc(allchainsburn[[i]])
  interburn[[i]] <- mcmc(allchainsburn[[i+5]])
}

gibbseff <- effectiveSize(gibbsburn[,1:13, drop=FALSE])
intereff <- effectiveSize(interburn[,1:13, drop=FALSE])

gibbstpd <- sum(gibbstime[,3])/5/20000*15000/gibbseff*1000
intertpd <- sum(intertime[,3])/5/20000*15000/intereff*1000

library(xtable)

timetab <- rbind(gibbseff, intereff, gibbstpd, intertpd)
timetab <- timetab[,c(1:6,8:13,7)]
colnames(timetab) <- c(paste("$V_", 1:6, "$", sep=""), "$U$", paste("$W_", 1:6, "$", sep=""))
rownames(timetab) <- NULL

xtable(timetab)

parests <- cbind(parmns[,1], parquants[,c(3,1,5)])
colnames(parests) <- c("Mean", "Median", "2.5\\%", "97.5\\%")
rownames(parests) <- c(paste("$V_",1:6, "$",sep=""),"$U$",paste("$W_",1:6,"$",sep=""))

xtable(parests, digits=3)

## effsizes, effsize/time or time/effsize, plots of rep and trt
## gelman-rubin
