## Plot Setup

library(ggplot2)
library(scales)
library(plyr)
library(xtable)
library(reshape2)
library(gridExtra)
load("../../cors/newpostcors.RData")
newpostcors <- newpostcors[newpostcors$V.T <= 10^2 & newpostcors$W.T <= 10^2,]
load("../../mixing/samout.RData")
samout2 <- samout
load("../../wrongscale/OldDAs/samout.RData")
samout$stime <- 0
samout2 <- rbind(samout2, samout)
load("../../cis/fullcissamout.RData")
samout <- rbind(samout2[samout2$sampler != "fullcis",], fullcissamout)
rm(samout2)

## Correctly name things
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
samout$samplers[substr(samout$sampler, 1, 2)=="sd"] <- "State-SD" 
samout$samplers[substr(samout$sampler, 1, 2)=="se"] <- "State-SE" 
samout$samplers[substr(samout$sampler, 1, 2)=="de"] <- "SD-SE" 
samout$samplers[substr(samout$sampler, 1, 3)=="tri"] <- "Triple"
samout$samplers[samout$sampler=="fullcis"] <- "CIS"
samout$samplers[samout$sampler=="partialcis"] <- "PartialCIS"
samout$samplers[samout$sampler=="error"] <- "SE"
samout$samplers[samout$sampler=="dist"] <- "SD"
samout$samplers[samout$sampler=="state"] <- "State"
samout$samplers[samout$sampler=="errorda"] <- "WSE"
samout$samplers[samout$sampler=="distda"] <- "WSD"
samlevels <- c("State", "SD", "SE", "State-SD", "State-SE", "SD-SE", 
               "Triple", "CIS", "PartialCIS", "WSD", "WSE")
samout$samplers <- factor(samout$samplers, levels=samlevels)

## time per effective draw of V and W
samout$V.time <- samout$time/samout$V.ES
samout$W.time <- samout$time/samout$W.ES

## Used to create the plots - puts the dataframe in a convenient form
meltedsam <- melt(samout, id=c("type", "samplers", "sampler", "V.T", "W.T", 
                            "T.T"))

## Breaks for both axes
Vs <- unique(meltedsam$V.T)[1:9] 
breaks <- Vs[seq(1,9,2)]
## labels for both axes
labs <- c("0.01", "0.1", "1", "10", "100")

## Function for taking the variable names, splitting off irrelevant text, and
## parsing the result as Latex.
label_parsed_split <- function(variable, value){
  llply(as.character(value), function(x) parse(text = strsplit(x, "\\.")[[1]][1]))
}

## Function for creating each of the ESP plots.
## meltedsam = dataframe that has been melted appropriately
## vars = variables to include in the plot (vector of strings)
## sams = samplers used in the plots (vector of strings)
## T = length of the time series (10, 100 or 1000)
## title = title of the plot
## guide = TRUE to include a legend
## type = TRUE to facet by sampler type (for GIS vs Alt plots)
##   else facets by variable (V or W)
## Returns a plot created by ggplot
plotfunES <- function(meltedsam, vars, sams, T, title, guide, type){
  if(guide){
    guide <- guide_colorbar(barheight=10)
  }
  if(type){
    facgrid <- facet_grid(type~samplers, scales="free", labeller=label_parsed_split)
  }
  else{
    facgrid <- facet_grid(variable~samplers, scales="free", labeller=label_parsed_split)
  }
  castedsam <- dcast(meltedsam, formula=samplers + V.T + W.T + variable + type ~ ., 
                     subset=.(variable %in% vars  & T.T==T & sampler %in% sams &
                         V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T, fill=value/2500)) + geom_tile() +
         scale_fill_gradient("ESP", low=muted("red"), high="white", guide=guide,
                             limits=c(0,1), na.value="white") +
         facgrid +
         scale_x_log10("V = noise", breaks=breaks, labels=labs) + 
         scale_y_log10("W = signal", breaks=breaks, labels=labs) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.8),
              axis.text.y = element_text(angle = 0, hjust = 1.1, vjust=0.8))
  return(out)
}

## Function for creating the computation time plots.
## meltedsam = dataframe that has been melted appropriately
## vars = variables to include in the plot (vector of strings)
## sams = samplers used in the plots (vector of strings)
## T = length of the time series (10, 100 or 1000)
## title = title of the plot
## guide = TRUE to include a legend
## type = TRUE to facet by sampler type (for GIS vs Alt plots)
##   else facets by variable (V or W)
## lims = lower and upper limits of the scale
## Returns a plot created by ggplot
plotfuntime <- function(meltedsam, vars, sams, T, title, lims, guide, type){
  if(guide){
    guide <- guide_colorbar(barheight=10)
  }
  if(type){
    facgrid <- facet_grid(type~samplers, scales="free", labeller=label_parsed_split)
  }
  else{
    facgrid <- facet_grid(variable~samplers, scales="free", labeller=label_parsed_split)
  }
  castedsam <- dcast(meltedsam, formula=samplers + V.T + W.T + variable + type ~ ., 
                     subset=.(variable %in% vars  & T.T==T & sampler %in% sams &
                         V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T,
                    fill=log(value*1000/60))) + 
         geom_tile() +
         scale_fill_gradient("Log min", high=muted("red"), low="white", guide=guide,
                             limits=lims, na.value=muted("red")) +
         facgrid +
         scale_x_log10("V = noise", breaks=breaks, labels=labs) + 
         scale_y_log10("W = signal", breaks=breaks, labels=labs) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.8),
               axis.text.y = element_text(angle = 0, hjust = 1.1, vjust=0.8))
  return(out)
}

## Function for creating the posterior correlation plots for T=100
## postcors = data frame of posterior correlation information
## var = variable to include in the plot (string) 
## title = title of the plot
## Returns a plot created by ggplot
plotfuncor <- function(newpostcors, var, title){
  dat <- newpostcors[newpostcors$T.T==100,]
  id <- which(colnames(newpostcors)==var)
  colnames(dat)[id] <- "value"
  out <- ggplot(data=dat, aes(x=V.T, y=W.T, fill=value)) +
      geom_tile() +
      scale_fill_gradient2("Corr", low=muted("blue"), high=muted("red"),
         limits=c(-1,1), mid="white") +
      facet_grid(.~T, scales="free", labeller=label_both) +
      scale_x_log10("V = noise", breaks=breaks, labels=labs) +
      scale_y_log10("W = signal", breaks=breaks, labels=labs) +
      ggtitle(title) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.8),
            axis.text.y = element_text(angle = 0, hjust = 1.1, vjust=0.8))
  return(out)
}


## ESP plots, for creating Figures 1, H.1, H.2 and H.3.
titlea <- "ESP for V in Alt and GIS samplers, T="
titleb <- "ESP for W in Alt and GIS samplers, T="
titlec <- "ESP for V and W in base and CIS samplers, T="
p1a <- plotfunES(meltedsam, "V.ES", c(alts,ints), 10, titlea, FALSE, TRUE)
p1b <- plotfunES(meltedsam, "W.ES", c(alts,ints), 10, titleb, FALSE, TRUE)
p1c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs,"fullcis"), 10, titlec, TRUE, FALSE)
p2a <- plotfunES(meltedsam, "V.ES", c(alts,ints), 100, titlea, FALSE, TRUE)
p2b <- plotfunES(meltedsam, "W.ES", c(alts,ints), 100, titleb, FALSE, TRUE)
p2c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs,"fullcis"), 100, titlec, TRUE, FALSE)
p3a <- plotfunES(meltedsam, "V.ES", c(alts,ints), 1000, titlea, FALSE, TRUE)
p3b <- plotfunES(meltedsam, "W.ES", c(alts,ints), 1000, titleb, FALSE, TRUE)
p3c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs,"fullcis"), 1000, titlec, TRUE, FALSE)
ggsave(filename="altintESplotV10.pdf", plot=p1a, width=6, height=3.75)
ggsave(filename="altintESplotW10.pdf", plot=p1b, width=6, height=3.75)
ggsave(filename="basecisESplot10.pdf", plot=p1c, width=9.4, height=3.75)
ggsave(filename="altintESplotV100.pdf", plot=p2a, width=6, height=3.75)
ggsave(filename="altintESplotW100.pdf", plot=p2b, width=6, height=3.75)
ggsave(filename="basecisESplot100.pdf", plot=p2c, width=9.4, height=3.75)
ggsave(filename="altintESplotV1000.pdf", plot=p3a, width=6, height=3.75)
ggsave(filename="altintESplotW1000.pdf", plot=p3b, width=6, height=3.75)
ggsave(filename="basecisESplot1000.pdf", plot=p3c, width=9.4, height=3.75)

## log time plots, for creating Figures 2, H.4, H.5 and H.6
titlea <- "Time per 1000 eff. draws in base and CIS samplers, T="
titleb <- "Time per 1000 eff. draws for V in Alt and GIS samplers, T="
titlec <- "Time per 1000 eff. draws for W in Alt and GIS samplers, T="
p1a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs,"fullcis"), 10, titlea, c(-5,2), TRUE, FALSE)
p1b <- plotfuntime(meltedsam, "V.time", c(alts,ints), 10, titleb, c(-5,2), FALSE, TRUE)
p1c <- plotfuntime(meltedsam, "W.time", c(alts,ints), 10, titlec, c(-5,2), FALSE, TRUE)
p2a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs,"fullcis"), 100, titlea, c(-3.5,3), TRUE, FALSE)
p2b <- plotfuntime(meltedsam, "V.time", c(alts,ints), 100, titleb, c(-3.5,3), FALSE, TRUE)
p2c <- plotfuntime(meltedsam, "W.time", c(alts,ints), 100, titlec, c(-3.5,3), FALSE, TRUE)
p3a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs,"fullcis"), 1000, titlea, c(-1,6), TRUE, FALSE)
p3b <- plotfuntime(meltedsam, "V.time", c(alts,ints), 1000, titleb, c(-1,6), FALSE, TRUE)
p3c <- plotfuntime(meltedsam, "W.time", c(alts,ints), 1000, titlec, c(-1,6), FALSE, TRUE)
ggsave(filename="basecistimeplot10.pdf", plot=p1a, width=9.4, height=3.75)
ggsave(filename="altgisVtimeplot10.pdf", plot=p1b, width=6, height=3.75)
ggsave(filename="altgisWtimeplot10.pdf", plot=p1c, width=6, height=3.75)
ggsave(filename="basecistimeplot100.pdf", plot=p2a, width=9.4, height=3.75)
ggsave(filename="altgisVtimeplot100.pdf", plot=p2b, width=6, height=3.75)
ggsave(filename="altgisWtimeplot100.pdf", plot=p2c, width=6, height=3.75)
ggsave(filename="basecistimeplot1000.pdf", plot=p3a, width=9.4, height=3.75)
ggsave(filename="altgisVtimeplot1000.pdf", plot=p3b, width=6, height=3.75)
ggsave(filename="altgisWtimeplot1000.pdf", plot=p3c, width=6, height=3.75)

## corplot, for creating Figure G.1
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
ggsave(filename="corplot1.pdf", plot=pvv, width=4, height=3)
ggsave(filename="corplot2.pdf", plot=pvw, width=4, height=3)
ggsave(filename="corplot3.pdf", plot=pva, width=4, height=3)
ggsave(filename="corplot4.pdf", plot=pvb, width=4, height=3)
ggsave(filename="corplot5.pdf", plot=pww, width=4, height=3)
ggsave(filename="corplot6.pdf", plot=pwv, width=4, height=3)
ggsave(filename="corplot7.pdf", plot=pwa, width=4, height=3)
ggsave(filename="corplot8.pdf", plot=pwb, width=4, height=3)
