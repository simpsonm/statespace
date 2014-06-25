## Plot Setup

library(ggplot2)
library(scales)
library(plyr)
library(xtable)
library(reshape2)
library(gridExtra)
load("../../mixing/samout.RData")
samout2 <- samout
load("../../wrongscale/OldDAs/samout.RData")
samout$stime <- 0
samout <- rbind(samout2, samout)
rm(samout2)
base <- c("error", "state", "dist")
alts <- c("sdalt", "sealt", "dealt", "trialt")
ints <- c("sdint", "seint", "deint", "triint")
kerns <- c("sdkern", "sekern", "dekern", "trikern")
cis <- c("fullcis", "partialcis")
wrongs <- c("errorda", "distda")
samout$V.ES[samout$sampler %in% kerns] <- samout$V.ES[samout$sampler %in% kerns]*2
samout$W.ES[samout$sampler %in% kerns] <- samout$W.ES[samout$sampler %in% kerns]*2
samout$V.ES[samout$sampler == "trikern"] <- samout$V.ES[samout$sampler == "trikern"]*(3/2)
samout$W.ES[samout$sampler == "trikern"] <- samout$W.ES[samout$sampler == "trikern"]*(3/2)
samout$type <- "Base" #$
samout$type[samout$sampler %in% alts] <- "Alt" 
samout$type[samout$sampler %in% ints] <- "GIS" 
samout$type[samout$sampler %in% kerns] <- "RKern" 
samout$type[samout$sampler %in% cis] <- "CIS" 
samout$type[samout$sampler %in% wrongs] <- "W-Base" 
samout$samplers <- "Base"
samout$samplers[substr(samout$sampler, 1, 2)=="sd"] <- "State-Dist" 
samout$samplers[substr(samout$sampler, 1, 2)=="se"] <- "State-Error" 
samout$samplers[substr(samout$sampler, 1, 2)=="de"] <- "Dist-Error" 
samout$samplers[substr(samout$sampler, 1, 3)=="tri"] <- "Triple" 
samout$samplers[samout$sampler=="fullcis"] <- "FullCIS"
samout$samplers[samout$sampler=="partialcis"] <- "PartialCIS"
samout$samplers[samout$sampler=="error"] <- "Error"
samout$samplers[samout$sampler=="dist"] <- "Dist"
samout$samplers[samout$sampler=="state"] <- "State"
samout$samplers[samout$sampler=="errorda"] <- "W-Error"
samout$samplers[samout$sampler=="distda"] <- "W-Dist"
samlevels <- c("State", "Dist", "Error", "State-Dist", "State-Error", "Dist-Error", 
               "Triple", "FullCIS", "PartialCIS", "W-Dist", "W-Error")
samout$samplers <- factor(samout$samplers, levels=samlevels)
samout$V.time <- samout$time/samout$V.ES
samout$W.time <- samout$time/samout$W.ES
meltedsam <- melt(samout, id=c("type", "samplers", "sampler", "V.T", "W.T", 
                            "T.T"))
Vs <- unique(meltedsam$V.T)[1:9] #$
Ws <- Vs
breaks <- Vs[seq(1,9,2)]
label_both_parsed <- function(variable, value){
  llply(as.character(paste(variable, value, sep = ": ")), function(x) parse(text = x))
}
label_both_parsed_split <- function(variable, value){
  llply(as.character(paste(variable, value, sep = ": ")), 
        function(x) parse(text = strsplit(x, "\\.")[[1]][1]))
}
label_parsed_split <- function(variable, value){
  llply(as.character(value), function(x) parse(text = strsplit(x, "\\.")[[1]][1]))
}
plotfun <- function(meltedsam, vars, sams, T, title){
  castedsam <- dcast(meltedsam, formula=sampler + V.T + W.T + variable + samplers ~ ., 
                     subset=.(variable %in% vars  & T.T==T & sampler %in% sams &
                       V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T, fill=value/2500)) + #$
         geom_tile() +
         scale_fill_gradient("ESP", low=muted("red"), high="white",
           guide=guide_colorbar(barheight=10),
           limits=c(0,1), na.value="white", trans="sqrt") +
         facet_grid(variable~samplers, scales="free", labeller=label_parsed_split) +
         scale_x_log10("V = noise", breaks=breaks) + scale_y_log10("W = signal", breaks=breaks) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  return(out)
}
load("../../cors/newpostcors.RData")
newpostcors <- newpostcors[newpostcors$V.T <= 10^2 & newpostcors$W.T <= 10^2,]
plotfuntime <- function(meltedsam, vars, sams, T, title, lims){
  castedsam <- dcast(meltedsam, 
                     formula=sampler + V.T + W.T + variable + 
                     samplers ~ ., 
                     subset=.(variable %in% vars  & 
                         T.T==T & sampler %in% sams &
                         V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T, 
                    fill=log(value*1000/60))) + #$
         geom_tile() +
         scale_fill_gradient("Log min", 
                             high=muted("red"), low="white",
                             guide=guide_colorbar(barheight=10), 
                             limits=lims, na.value="red") +
         facet_grid(variable~samplers, scales="free", 
                    labeller=label_parsed_split) +
         scale_x_log10("V = noise", breaks=breaks) + 
             scale_y_log10("W = signal", breaks=breaks) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, 
                   hjust = 1, vjust=0.5))
  return(out)
  ##, limits=c(0,top)
}
plotfuncor <- function(newpostcors, var, title){
  dat <- newpostcors[newpostcors$T.T!=10,]
  id <- which(colnames(newpostcors)==var)
  colnames(dat)[id] <- "value"
  out <- ggplot(data=dat, aes(x=V.T, y=W.T, fill=value)) +
      geom_tile() +
      scale_fill_gradient2("Corr", low=muted("blue"), high=muted("red"),
         limits=c(-1,1), mid="white") +
      facet_grid(.~T, scales="free", labeller=label_both) +
      scale_x_log10("V = noise", breaks=breaks) +
      scale_y_log10("W = signal", breaks=breaks) +
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  return(out)
}


## baseESplot, fig.cap=cap, echo=FALSE, fig.height=3.75, fig.width=8, out.width=".7\\textwidth"
vars <- c("V.ES", "W.ES")
title <- "ESP for V and W in the base algorithms, T="
##p1 <- plotfun(meltedsam, vars, c(base,wrongs), 10, title)
p2 <- plotfun(meltedsam, vars, c(base,wrongs), 100, title)
p3 <- plotfun(meltedsam, vars, c(base,wrongs), 1000, title)
##p1
ggsave(filename="baseESplot1.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="baseESplot2.pdf", plot=p3, width=8, height=3.75)



## corplot, echo=FALSE, fig.height=3, fig.width=5, out.width=".245\\textwidth", fig.cap=cap
title <- expression(paste("Posterior Correlation Between V and ",b[V], sep=""))
pvv <- plotfuncor(newpostcors, "Vbv", title)
title <- expression(paste("Posterior Correlation Between W and ",b[W], sep=""))
pww <- plotfuncor(newpostcors, "Wbw", title)
title <- expression(paste("Posterior Correlation Between V and ",b[W], sep=""))
pvw <- plotfuncor(newpostcors, "Vbw", title)
title <- expression(paste("Posterior Correlation Between W and ",b[V], sep=""))
pwv <- plotfuncor(newpostcors, "Wbv", title)
title <- expression(paste("Posterior Correlation Between V and ",a[psi], sep=""))
pva <- plotfuncor(newpostcors, "Vapsi", title)
title <- expression(paste("Posterior Correlation Between W and ",a[gamma], sep=""))
pwa <- plotfuncor(newpostcors, "Wagam", title)
title <- expression(paste("Posterior Correlation Between V and ",b[psi], sep=""))
pvb <- plotfuncor(newpostcors, "Vbpsi", title)
title <- expression(paste("Posterior Correlation Between W and ",b[gamma], sep=""))
pwb <- plotfuncor(newpostcors, "Wbgam", title)
ggsave(filename="corplot1.pdf", plot=pvv, width=5, height=3)
ggsave(filename="corplot2.pdf", plot=pvw, width=5, height=3)
ggsave(filename="corplot3.pdf", plot=pva, width=5, height=3)
ggsave(filename="corplot4.pdf", plot=pvb, width=5, height=3)
ggsave(filename="corplot5.pdf", plot=pww, width=5, height=3)
ggsave(filename="corplot6.pdf", plot=pwv, width=5, height=3)
ggsave(filename="corplot7.pdf", plot=pwa, width=5, height=3)
ggsave(filename="corplot8.pdf", plot=pwb, width=5, height=3)


## intESplot, fig.cap=cap, echo=FALSE, fig.height=3.75, fig.width=8, out.width=".7\\textwidth"
vars <- c("V.ES", "W.ES")
sams <- c("deint", "seint", "sdint", "triint", "fullcis")
title <- "ESP for V and W in the GIS and CIS algorithms, T="
##p1 <- plotfun(meltedsam, vars, sams, 10, title)
p2 <- plotfun(meltedsam, vars, sams, 100, title)
p3 <- plotfun(meltedsam, vars, sams, 1000, title)
##p1
ggsave(filename="intESplot1.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="intESplot2.pdf", plot=p3, width=8, height=3.75)



## altESplot, fig.cap=cap, echo=FALSE, fig.height=3.75, fig.width=7, out.width=".7\\textwidth"
vars <- c("V.ES", "W.ES")
sams <- c(ints)
title <- "ESP for V and W in the alternating algorithms, T="
##p1 <- plotfun(meltedsam, vars, sams, 10, title)
p2 <- plotfun(meltedsam, vars, sams, 100, title)
p3 <- plotfun(meltedsam, vars, sams, 1000, title)
##p1
ggsave(filename="altESplot1.pdf", plot=p2, width=7, height=3.75)
ggsave(filename="altESplot2.pdf", plot=p3, width=7, height=3.75)



## baseinttimeplot, fig.cap=cap, echo=FALSE, fig.width=10, fig.height=3.25, out.width=".8\\textwidth"
vars <- c("V.time", "W.time")
sams <- c("dist", "error", "deint", "state", "seint", "sdint", "triint", "fullcis")
title <- "Log of time (minutes) per 1000 effective draws for base and interweaving samplers, T="
##p1 <- plotfuntime(meltedsam, vars, sams, 10, title, 25*60/1000)
p2 <- plotfuntime(meltedsam, vars, sams, 100, title, c(-3.5,5))
p3 <- plotfuntime(meltedsam, vars, sams, 1000, title, c(-1,8))
##p1
ggsave(filename="baseinttimeplot1.pdf", plot=p2, width=10, height=3.25)
ggsave(filename="baseinttimeplot2.pdf", plot=p3, width=10, height=3.25)



## altinttimeplot, fig.cap=cap, echo=FALSE, fig.width=8, fig.height=3.75, out.width=".49\\textwidth"
sams <- c(alts,ints)
title <- "Log of time (minutes) per 1000 effective draws for alternating samplers, T="
##p1 <- plotfuntime(meltedsam, vars, sams, 10, title, 25*60/1000)
p2 <- plotfuntime(meltedsam, vars, sams, 100, title, c(-3.5,5))
p3 <- plotfuntime(meltedsam, vars, sams, 1000, title, c(-1,8))
ggsave(filename="altinttimeplot1.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="altinttimeplot2.pdf", plot=p3, width=8, height=3.75)

