## Plot Setup
library(ggplot2)
library(scales)
library(plyr)
library(xtable)
library(reshape2)
library(gridExtra)
###########################################################################
###### make sure that postcors.RData and samout.RData are in this directory
###########################################################################
load("../../cors/newpostcors.RData")
postcors <- newpostcors[newpostcors$V.T <= 10^2 & newpostcors$W.T <= 10^2,]
load("samoutlongmix.RData")

## Correctly name things
base <- c("error", "state", "dist")
alts <- c("sdalt", "sealt", "dealt", "trialt")
ints <- c("sdint", "seint", "deint", "triint")
kerns <- c("sdkern", "sekern", "dekern", "trikern")
cis <- c("fullcis")
wrongs <- c("wdist", "werror")
samout$type <- "Base" #$
samout$type[samout$sampler %in% alts] <- "Alt" 
samout$type[samout$sampler %in% ints] <- "GIS" 
samout$type[samout$sampler %in% cis] <- "CIS" 
samout$type[samout$sampler %in% wrongs] <- "W-Base" 
samout$samplers <- "Base"
samout$samplers[substr(samout$sampler, 1, 2)=="sd"] <- "State-SD" 
samout$samplers[substr(samout$sampler, 1, 2)=="se"] <- "State-SE" 
samout$samplers[substr(samout$sampler, 1, 2)=="de"] <- "SD-SE" 
samout$samplers[substr(samout$sampler, 1, 3)=="tri"] <- "Triple"
samout$samplers[samout$sampler=="fullcis"] <- "CIS"
samout$samplers[samout$sampler=="error"] <- "SE"
samout$samplers[samout$sampler=="dist"] <- "SD"
samout$samplers[samout$sampler=="state"] <- "State"
samout$samplers[samout$sampler=="werror"] <- "WSE"
samout$samplers[samout$sampler=="wdist"] <- "WSD"
samlevels <- c("State", "SD", "SE", "WSD", "WSE", "State-SD", "State-SE", "SD-SE", "Triple", "CIS")
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
## parsing the result as Latex code.
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
## Returns a plot created by ggplot
plotfunES <- function(meltedsam, vars, sams, T, title, guide, type){
  if(guide){
    guide <- guide_colorbar(barheight=10)
  }
  facgrid <- facet_grid(variable~samplers, scales="free", labeller=label_parsed_split)
  castedsam <- dcast(meltedsam, formula=samplers + V.T + W.T + variable + type ~ ., 
                     subset=.(variable %in% vars  & T.T==T & sampler %in% sams &
                         V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T, fill=value/10000)) + geom_tile() +
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
## lims = lower and upper limits of the scale
## Returns a plot created by ggplot
plotfuntime <- function(meltedsam, vars, sams, T, title, lims, guide){
  if(guide){
    guide <- guide_colorbar(barheight=10)
  }
  facgrid <- facet_grid(variable~samplers, scales="free", labeller=label_parsed_split)
  castedsam <- dcast(meltedsam, formula=samplers + V.T + W.T +
                         variable + type ~ ., 
                     subset=.(variable %in% vars  & T.T==T &
                                  sampler %in% sams &
                                      V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  castedsam$value[log(castedsam$value*1000/60) > lims[2]] <-
    exp(lims[2])*60/1000 - .0001
  castedsam$value[log(castedsam$value*1000/60) < lims[1]] <-
    exp(lims[1])*60/1000 + .0001
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T,
                    fill=log(value*1000/60))) + 
      geom_tile() +
      scale_fill_gradient("Log min", high=muted("red"),
                              low="white", guide=guide,
                              limits=lims, na.value="grey") +
      facgrid +
      scale_x_log10("V = noise", breaks=breaks, labels=labs) + 
      scale_y_log10("W = signal", breaks=breaks, labels=labs) +
      ggtitle(paste(title, T, sep="")) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,
                vjust=0.8),
            axis.text.y = element_text(angle = 0, hjust = 1.1,
                vjust=0.8))
  return(out)
}

## Function for creating the posterior correlation plots for T=100
## postcors = data frame of posterior correlation information
## var = variable to include in the plot (string) 
## title = title of the plot
## Returns a plot created by ggplot
plotfuncor <- function(postcors, var, title){
  dat <- postcors[postcors$T.T==100,]
  id <- which(colnames(postcors)==var)
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
titlea <- "ESP for V and W GIS and CIS samplers, T="
titleb <- "ESP for V and W in Alt samplers, T="
titlec <- "ESP for V and W in base samplers, T="
p1a <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(ints,"fullcis"), 10, titlea, FALSE)
p1b <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(alts), 10, titleb, FALSE)
p1c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs), 10, titlec, TRUE)
p2a <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(ints,"fullcis"), 100, titlea, FALSE)
p2b <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(alts), 100, titleb, FALSE)
p2c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs), 100, titlec, TRUE)
p3a <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(ints,"fullcis"), 1000, titlea, FALSE)
p3b <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(alts), 1000, titleb, FALSE)
p3c <- plotfunES(meltedsam, c("V.ES", "W.ES"), c(base,wrongs), 1000, titlec, TRUE)

ggsave(filename="altintESplotV10.pdf", plot=p1a, width=6.7, height=3.5)
ggsave(filename="altintESplotW10.pdf", plot=p1b, width=5.7, height=3.5)
ggsave(filename="basecisESplot10.pdf", plot=p1c, width=8, height=3.5)
ggsave(filename="altintESplotV100.pdf", plot=p2a, width=6.7, height=3.5)
ggsave(filename="altintESplotW100.pdf", plot=p2b, width=5.7, height=3.5)
ggsave(filename="basecisESplot100.pdf", plot=p2c, width=8, height=3.5)
ggsave(filename="altintESplotV1000.pdf", plot=p3a, width=6.7, height=3.5)
ggsave(filename="altintESplotW1000.pdf", plot=p3b, width=5.7, height=3.5)
ggsave(filename="basecisESplot1000.pdf", plot=p3c, width=8, height=3.5)

## log time plots, for creating Figures 2, H.4, H.5 and H.6
titlea <- "Time per 1000 eff. draws in base samplers, T="
titleb <- "Time per 1000 eff. draws GIS and CIS samplers, T="
titlec <- "Time per 1000 eff. draws in Alt, T="
p1a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs), 10, titlea, c(-5,0), TRUE)
p1b <- plotfuntime(meltedsam, c("V.time", "W.time"), c(ints, "fullcis"), 10, titleb, c(-5,0), FALSE)
p1c <- plotfuntime(meltedsam, c("V.time", "W.time"), alts, 10, titlec, c(-5,0), FALSE)
p2a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs), 100, titlea, c(-4,3), TRUE)
p2b <- plotfuntime(meltedsam, c("V.time", "W.time"), c(ints,"fullcis"), 100, titleb, c(-4,3), FALSE)
p2c <- plotfuntime(meltedsam, c("V.time", "W.time"), alts, 100, titlec, c(-4,3), FALSE)
p3a <- plotfuntime(meltedsam, c("V.time", "W.time"), c(base,wrongs), 1000, titlea, c(-2,6), TRUE)
p3b <- plotfuntime(meltedsam, c("V.time", "W.time"), c(ints,"fullcis"), 1000, titleb, c(-2,6), FALSE)
p3c <- plotfuntime(meltedsam, c("V.time", "W.time"), alts, 1000, titlec, c(-2,6), FALSE)

ggsave(filename="basecistimeplot10.pdf", plot=p1a, width=8, height=3.5)
ggsave(filename="altgisVtimeplot10.pdf", plot=p1b, width=6.7, height=3.5)
ggsave(filename="altgisWtimeplot10.pdf", plot=p1c, width=5.7, height=3.5)
ggsave(filename="basecistimeplot100.pdf", plot=p2a, width=8, height=3.5)
ggsave(filename="altgisVtimeplot100.pdf", plot=p2b, width=6.7, height=3.5)
ggsave(filename="altgisWtimeplot100.pdf", plot=p2c, width=5.7, height=3.5)
ggsave(filename="basecistimeplot1000.pdf", plot=p3a, width=8, height=3.5)
ggsave(filename="altgisVtimeplot1000.pdf", plot=p3b, width=6.7, height=3.5)
ggsave(filename="altgisWtimeplot1000.pdf", plot=p3c, width=5.7, height=3.5)

## corplot, for creating Figure G.1
title <- expression(paste("Posterior Correlation Between V and ",b[V], sep=""))
pvv <- plotfuncor(postcors, "Vbv", title)
title <- expression(paste("Posterior Correlation Between W and ",b[W], sep=""))
pww <- plotfuncor(postcors, "Wbw", title)
title <- expression(paste("Posterior Correlation Between V and ",b[W], sep=""))
pvw <- plotfuncor(postcors, "Vbw", title)
title <- expression(paste("Posterior Correlation Between W and ",b[V], sep=""))
pwv <- plotfuncor(postcors, "Wbv", title)
title <- expression(paste("Posterior Correlation Between V and ",a[psi], sep=""))
pva <- plotfuncor(postcors, "Vapsi", title)
title <- expression(paste("Posterior Correlation Between W and ",a[gamma], sep=""))
pwa <- plotfuncor(postcors, "Wagam", title)
title <- expression(paste("Posterior Correlation Between V and ",b[psi], sep=""))
pvb <- plotfuncor(postcors, "Vbpsi", title)
title <- expression(paste("Posterior Correlation Between W and ",b[gamma], sep=""))
pwb <- plotfuncor(postcors, "Wbgam", title)

ggsave(filename="corplot1.pdf", plot=pvv, width=4, height=3)
ggsave(filename="corplot2.pdf", plot=pvw, width=4, height=3)
ggsave(filename="corplot3.pdf", plot=pva, width=4, height=3)
ggsave(filename="corplot4.pdf", plot=pvb, width=4, height=3)
ggsave(filename="corplot5.pdf", plot=pww, width=4, height=3)
ggsave(filename="corplot6.pdf", plot=pwv, width=4, height=3)
ggsave(filename="corplot7.pdf", plot=pwa, width=4, height=3)
ggsave(filename="corplot8.pdf", plot=pwb, width=4, height=3)



###############################################################
### Now make sure that samoutlong.RData is in the directory
###############################################################

load("samoutlong.RData")
## log time in minutes per 1000 effective draws
samout$V.ET <- samout$time/samout$V.ES/60*1000
samout$W.ET <- samout$time/samout$W.ES/60*1000

samout <- samout[,c(1,4,14,15)]
samout$sampler[samout$sampler=="dealt"] <- "SD-SE Alt"
samout$sampler[samout$sampler=="deint"] <- "SD-SE GIS"
colnames(samout)[2:4] <- c("T", "V", "W")
meltedsam <- melt(samout, id=c("sampler", "T"))
pgisalt <- qplot(T, value, data=meltedsam, facets=variable~., linetype=sampler, geom="line", size=I(2)) + ylab("time (minutes) per 1000 effective draws")

ggsave(filename="gisaltplot.pdf", plot=pgisalt, width=10, height=4)



