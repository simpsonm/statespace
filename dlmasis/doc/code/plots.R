## Plot Setup

library(ggplot2)
library(scales)
library(plyr)
library(xtable)
library(reshape2)
library(gridExtra)
## make sure that postcors.RData and samout.RData are in this directory
load("postcors.RData")
postcors <- postcors[postcors$V.T <= 10^2 & postcors$W.T <= 10^2,]
load("samout.RData")

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
## Returns a plot created by ggplot
plotfun <- function(meltedsam, vars, sams, T, title){
  castedsam <- dcast(meltedsam, formula=sampler + V.T + W.T + variable + samplers ~ ., 
                     subset=.(variable %in% vars  & T.T==T & sampler %in% sams &
                       V.T<=10^2 & W.T<=10^2))
  colnames(castedsam)[6] <- "value"
  ## value = elements of variable. They should all be effective sample sizes,
  ## so value/2500 computes effective sample proportion
  out <- ggplot(data=castedsam, aes(x=V.T, y=W.T, fill=value/2500)) + 
         geom_tile() +
         scale_fill_gradient("ESP", low=muted("red"), high="white",
           guide=guide_colorbar(barheight=10),
           limits=c(0,1), na.value="white", trans="sqrt") +
         facet_grid(variable~samplers, scales="free", labeller=label_parsed_split) +
         scale_x_log10("V = noise", breaks=breaks, labels=labs) + scale_y_log10("W = signal", breaks=breaks, labels=labs) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  return(out)
}

## Function for creating the computation time plots.
## meltedsam = dataframe that has been melted appropriately
## vars = variables to include in the plot (vector of strings)
## sams = samplers used in the plots (vector of strings)
## T = length of the time series (10, 100 or 1000)
## title = title of the plot
## lims = lower and upper limits of the scale
## Returns a plot created by ggplot
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
         scale_x_log10("V = noise", breaks=breaks, labels=labs) + 
         scale_y_log10("W = signal", breaks=breaks, labels=labs) +
         ggtitle(paste(title, T, sep="")) +
         theme(axis.text.x = element_text(angle = 90, 
                   hjust = 1, vjust=0.5))
  return(out)
}

## Function for creating the posterior correlation plots for T=100
## postcors = data frame of posterior correlation information
## var = variable to include in the plot (string) 
## title = title of the plot
## Returns a plot created by ggplot
plotfuncor <- function(postcors, var, title){
  dat <- postcors[postcors$T.T==100,] ## change this if you want to see T=10 or T=1000
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
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
  return(out)
}


## create plot of base sampler ESPs, Figures 1 and G.1
vars <- c("V.ES", "W.ES") ## effective sample size for V and W
title <- "ESP for V and W in the base algorithms, T="
## plots include base samplers (state, SD, and SE) and
## wrongly scaled samplers (W-SD, W-SE)
p1 <- plotfun(meltedsam, vars, c(base,wrongs), 10, title)
p2 <- plotfun(meltedsam, vars, c(base,wrongs), 100, title)
p3 <- plotfun(meltedsam, vars, c(base,wrongs), 1000, title)
ggsave(filename="baseESplot10.pdf", plot=p1, width=8, height=3.75)
ggsave(filename="baseESplot100.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="baseESplot1000.pdf", plot=p3, width=8, height=3.75)



# create plot of posterior correlations, Figure F.1
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


## create plot of interweaving ESPs, Figures 2 and G.2
vars <- c("V.ES", "W.ES") ## effective sample size for V and W
## Includes SD-SE Int, SE-State Int, SD-State Int, Triple Int, and CIS samplers
sams <- c("deint", "seint", "sdint", "triint", "fullcis") 
title <- "ESP for V and W in the GIS and CIS algorithms, T="
p1 <- plotfun(meltedsam, vars, sams, 10, title)
p2 <- plotfun(meltedsam, vars, sams, 100, title)
p3 <- plotfun(meltedsam, vars, sams, 1000, title)
ggsave(filename="intESplot10.pdf", plot=p1, width=8, height=3.75)
ggsave(filename="intESplot100.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="intESplot1000.pdf", plot=p3, width=8, height=3.75)



## create plot of alternating ESPs, Figures 3 and G.3
vars <- c("V.ES", "W.ES") ## effective sample size for V and W
## Includes SD-SE Alt, SE-State Alt, SD-State Alt, and Triple Alt
sams <- c(alts)
title <- "ESP for V and W in the alternating algorithms, T="
p1 <- plotfun(meltedsam, vars, sams, 10, title)
p2 <- plotfun(meltedsam, vars, sams, 100, title)
p3 <- plotfun(meltedsam, vars, sams, 1000, title)
ggsave(filename="altESplot10.pdf", plot=p1, width=7, height=3.75)
ggsave(filename="altESplot100.pdf", plot=p2, width=7, height=3.75)
ggsave(filename="altESplot1000.pdf", plot=p3, width=7, height=3.75)



## create plot of log time per 1000 effective draws for
## base and interweaving samplers, Figures 4 and G.4
vars <- c("V.time", "W.time") ## time per effective draw for V and W
## Includes SD, SE, DE Int, State, SE Int, SD Int, Triple Int and CIS samplers
sams <- c("dist", "error", "deint", "state", "seint", "sdint", "triint", "fullcis")
title <- "Log minutes per 1000 effective draws for base and interweaving samplers, T="
p1 <- plotfuntime(meltedsam, vars, sams, 10, title, c(-5,3.5))
p2 <- plotfuntime(meltedsam, vars, sams, 100, title, c(-3.5,5))
p3 <- plotfuntime(meltedsam, vars, sams, 1000, title, c(-1,8))
ggsave(filename="baseinttimeplot10.pdf", plot=p1, width=10, height=3.25)
ggsave(filename="baseinttimeplot100.pdf", plot=p2, width=10, height=3.25)
ggsave(filename="baseinttimeplot1000.pdf", plot=p3, width=10, height=3.25)



## create plot of log time per 1000 effective draws for
## alternating samplers, Figure G.5
vars <- c("V.time", "W.time") ## time per effective draw for V and W
## Includes DE Alt, SE Alt, SD Alt, and Triple Alt samplers
sams <- c(alts)
title <- "Log minutes per 1000 effective draws for alternating samplers, T="
p1 <- plotfuntime(meltedsam, vars, sams, 10, title, c(-5,1))
p2 <- plotfuntime(meltedsam, vars, sams, 100, title, c(-3.5,5))
p3 <- plotfuntime(meltedsam, vars, sams, 1000, title, c(-1,8))
ggsave(filename="altinttimeplot10.pdf", plot=p1, width=8, height=3.75)
ggsave(filename="altinttimeplot100.pdf", plot=p2, width=8, height=3.75)
ggsave(filename="altinttimeplot1000.pdf", plot=p3, width=8, height=3.75)

