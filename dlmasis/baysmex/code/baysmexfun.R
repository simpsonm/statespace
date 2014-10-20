library(ars)
source("mcfa.R")
source("rej.R")

naivegibbs <- function(n, start, dat, au, bu, av, bv, aw, bw, m0, C0){
  J <- nrow(dat)
  N <- ncol(dat)
  F <- cbind(rep(1,J),diag(1,J))
  Vs <- start[1:J]
  Ws <- start[(J+1):(2*J+1)]
  out <- matrix(0, nrow=n, ncol=2*J + (J+1)*(N+1) + 1)
  thetanames <- NULL
  for(j in 1:J){
    thetanames <- c(thetanames, paste("theta", j, 0:N, sep="."))
  }
  colnames(out) <- c(paste("V", 1:J, sep="."), "U", paste("W", 1:J, sep="."),
                     paste("mu", 0:N, sep="."), thetanames)
  aut <- au + N/2
  avt <- av + N/2
  awt <- aw + N/2
  but <- bu
  bvt <- bv
  bwt <- bw
  for(i in 1:n){
    if(J>1){
      V <- diag(Vs)
    }
    if(J==1){
      V <- Vs
    }
    W <- diag(Ws)
    phi <- mcfaphi(dat, V, W, m0, C0, F, N, J)
    mu <- phi[1,]
    for(j in 1:J){
      thetaj <- phi[j+1,]
      yj <- dat[j,]
      bvt[j] <- bv[j] + sum((yj - mu[-1] - thetaj[-1])^2)/2
      bwt[j] <- bw[j] + sum(diff(thetaj)^2)/2
      Vs[j] <- 1/rgamma(1,shape=avt[j], rate=bvt[j])
      Ws[j+1] <- 1/rgamma(1,shape=awt[j], rate=bwt[j])
    }
    but <- bu + sum(diff(mu)^2)/2
    Ws[1] <- 1/rgamma(1, shape=aut, rate=but)
    out[i,] <- c(Vs,Ws,t(phi))
  }
  return(out)
}


disterrorinter <- function(n, start, dat, au, bu, av, bv, aw, bw, m0, C0){
  J <- nrow(dat)
  N <- ncol(dat)
  F <- cbind(rep(1,J),diag(1,J))
  Vs <- start[1:J]
  Ws <- start[(J+1):(2*J+1)]
  out <- matrix(0, nrow=n, ncol=2*J + (J+1)*(N+1) + 1)
  thetanames <- NULL
  for(j in 1:J){
    thetanames <- c(thetanames, paste("theta", j, 0:N, sep="."))
  }
  colnames(out) <- c(paste("V", 1:J, sep="."), "U", paste("W", 1:J, sep="."),
                     paste("mu", 0:N, sep="."), thetanames)
  aut <- au + N/2
  avt <- av + N/2
  awt <- aw + N/2
  but <- bu
  bvt <- bv
  bwt <- bw
  for(i in 1:n){
    V <- diag(Vs)
    W <- diag(Ws)
    phi <- mcfaphi(dat, V, W, m0, C0, F, N, J)
    mu <- phi[1,]
    for(j in 1:J){
      thetaj <- phi[j+1,]
      yj <- dat[j,]
      bvt[j] <- bv[j] + sum((yj - mu[-1] - thetaj[-1])^2)/2
      Vs[j] <- 1/rgamma(1,shape=avt[j], rate=bvt[j])
      gamj <- c(thetaj[1], diff(thetaj)/sqrt(Ws[j+1]))
      cwj <- sum((yj - mu[-1] - gamj[1])*cumsum(gamj[-1]))/Vs[j]
      dwj <- sum(cumsum(gamj[-1])^2)/(2*Vs[j])
      Ws[j+1] <- Wgamiter(aw[j], bw[j], cwj, dwj)
      thetaj <- thetagamtrans(gamj, Ws[j+1])
      psij <- c(thetaj[1], (yj - mu[-1] - thetaj[-1])/sqrt(Vs[j]))
      cvj1 <-  sum(diff(psij[-1])*(diff(yj) - diff(mu[-1])))
      cvj2 <- psij[2]*(yj[1]-mu[2]-psij[1])
      cvj <- (cvj1 + cvj2)/Ws[j+1]
      dvj <- (psij[2]^2 + sum(diff(psij[-1])^2))/(2*Ws[j+1])
      Vs[j] <- Vpsiiter(av[j], bv[j], cvj, dvj)
      thetaj <- c(psij[1], yj - mu[-1] - sqrt(Vs[j])*psij[-1])
      phi[j+1,] <- thetaj
      bwt[j] <- bw[j] + sum(diff(thetaj)^2)/2
      Ws[j+1] <- 1/rgamma(1,shape=awt[j], rate=bwt[j])
    }
    but <- bu + sum(diff(mu)^2)/2
    Ws[1] <- 1/rgamma(1, shape=aut, rate=but)
    out[i,] <- c(Vs,Ws,t(phi))
  }
  return(out)
}
