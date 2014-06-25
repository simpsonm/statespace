## A set of functions for simulating from and fitting local level models
library(dlm)
library(coda)
library(MCMCpack)
library(ars)
library(plyr)
source("dlmasisfun.R")

chains.diag <- function(ns, samplers, simdata, av, aw, parallel, mix=FALSE){
  n <- ns$n[1]
  G.D.full <- ddply(samplers, .(sampler, iter), fullsam.diag, simdata=simdata,
                    n=n, av=av, aw=aw, parallel=parallel, .parallel=parallel, mix=mix)
  return(G.D.full)
}


fullsam.diag <- function(samplers, simdata, n, av, aw, parallel=parallel, mix=FALSE){
  sampler <- samplers$sampler[1]
  sam <- ddply(simdata, .(V.T, W.T, T.T, V.S, W.S, ch), samwrapstart,
               n=n, av=av, aw=aw, samp=sampler)
  samdiag <- ddply(sam, .(V.T, W.T, T.T), sam.diag, .parallel=parallel,
                   parallel=parallel, mix=mix)
  return(samdiag)
}

sam.diag <- function(sam, parallel, mix=FALSE){
  T <- sam$T.T[1]
  namid <- grep("(V.T|W.T|T.T|time)", colnames(sam))
  sam.par <- sam[,-namid]
  sam.par <- sam.par[,1:(T+1+2+2+1)]
  if(!mix){
    sam.list <- dlply(sam.par, .(ch, V.S, W.S), function(x){
      namid <- grep("(V.S|W.S|ch)", colnames(x))
      x.par <- x[,-namid]
      return(mcmc(x.par))
    }, .parallel=parallel)
    GD <- gelman.diag(mcmc.list(sam.list))
    G.D <- data.frame(G.D.M=GD[[2]], G.D.V=GD[[1]][1,1], G.D.W=GD[[1]][2,1])
  }
  else{
    ns <- seq(from=100, to=dim(sam.par)[1]/5, by= 100)
    G.D <- NULL
    for(n in ns){
      sam.list <- dlply(sam.par, .(ch, V.S, W.S), function(x,n){
        namid <- grep("(V.S|W.S|ch)", colnames(x))
        x.par <- x[1:n,-namid]
        return(mcmc(x.par))
      }, .parallel=parallel, n=n)
      GD <- gelman.diag(mcmc.list(sam.list))
      G.D <- rbind(G.D, data.frame(n2=n, G.D.M=GD[[2]], G.D.V=GD[[1]][1,1], G.D.W=GD[[1]][2,1]))
    }
  }
  return(G.D)
}

## Wrapper for quickly simulating from all samplers w/ multiple
## chains at diff starting values
samwrapstart <- function(par, n, av, aw, samp){
  dat <- par$y[order(par$t)]
  T <- length(dat)
  start <- c(par$V.S[1]*par$V.T[1], par$W.S[1]*par$W.T[1])
  bv <- (av-1)*par$V.T[1]
  bw <- (aw-1)*par$W.T[1]
  if(samp=="state")
    time <- system.time(out <- statesam(n, start, dat, av, aw, bv, bw))
  if(samp=="dist")
    time <- system.time(out <- distsam(n, start, dat, av, aw, bv, bw))
  if(samp=="error")
    time <- system.time(out <- errorsam(n, start, dat, av, aw, bv, bw))
  if(samp=="sdint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 1))
  if(samp=="seint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 2))
  if(samp=="deint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, TRUE, 3))
  if(samp=="triint")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(TRUE, TRUE), 4))
  if(samp=="sdalt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 1))
  if(samp=="sealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 2))
  if(samp=="dealt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, FALSE, 3))
  if(samp=="trialt")
    time <- system.time(out <- samwrapper(n, start, dat, av, aw, bv, bw, c(FALSE, FALSE), 4))
  if(samp=="sdkern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 1/2, 0)))
  if(samp=="sekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(1/2, 0, 1/2)))
  if(samp=="dekern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw, c(0, 1/2, 1/2)))
  if(samp=="trikern")
    time <- system.time(out <- randkernsam(n, start, dat, av, aw, bv, bw))
  if(samp=="partialcis")
    time <- system.time(out <- partialcissam(n, start, dat, av, aw, bv, bw))
  if(samp=="fullcis")
    time <- system.time(out <- fullcissam(n, start, dat, av, aw, bv, bw))
  outdat <- data.frame(out)
  outdat$time <- time[1]
  return(outdat[,c(T+9, (T+2):(T+8), 1:(T+1))])
}
