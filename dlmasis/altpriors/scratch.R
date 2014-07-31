rm(samout)
load("../../altpriors/priorsims/stateOUT.RData")
stateout <- out
load("../../altpriors/priorsims/distOUT.RData")
distout <- out
load("../../altpriors/priorsims/errorOUT.RData")
errorout <- out
rm(out)
stateout$sampler <- "state"
distout$sampler <- "dist"
errorout$sampler <- "error"
samout <- rbind(stateout, distout, errorout)
base <- c("error", "state", "dist")
samout$Wefftime <- samout$time*3000/samout$W.ES
samout$Vefftime <- samout$time*3000/samout$V.ES
samout$type <- "Base" 
samout$samplers <- "Base"
samout$samplers[samout$sampler=="error"] <- "SE"
samout$samplers[samout$sampler=="dist"] <- "SD"
samout$samplers[samout$sampler=="state"] <- "State"
samlevels <- c("State", "SD", "SE", "State-SD", "State-SE", "SD-SE", 
               "Triple", "CIS", "PartialCIS", "WSD", "WSE")
samout$samplers <- factor(samout$samplers, levels=samlevels)
samout$V.time <- samout$time/samout$V.ES
samout$W.time <- samout$time/samout$W.ES
meltedsam <- melt(samout, id=c("type", "samplers", "sampler", "V.T", "W.T", 
                            "T.T"))



## base plots for alternative priors
vars <- c("V.ES", "W.ES")
title <- "ESP for V and W in the base algorithms, T="
p1 <- plotfun(meltedsam, vars, c(base), 10, title)
p2 <- plotfun(meltedsam, vars, c(base), 100, title)
p3 <- plotfun(meltedsam, vars, c(base), 1000, title)
ggsave(filename="basealtpriorESplot10.pdf", plot=p1, width=8, height=3.75)
ggsave(filename="basealtpriorESplot100.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="basealtpriorESplot1000.pdf", plot=p3, width=8, height=3.75)



plotfun(meltedsam, vars, "state", 10, title)
