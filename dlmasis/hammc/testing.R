source("nutsfun2.R")

lpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]

  out <- (ayou + T/2)*( lvw +  log(a + 2*b*exp(-lvw/2) + c*exp(-lvw))) + ame*lvw + bet*exp(-lvw)
  return(-out)
}

pr <- function(lvw, par){
  logp <- lpr(lvw,par)
  out <- exp(logp)
  return(out)
}

gradlpr <- function(lvw, par){
  a <- par[1]
  b <- par[2]
  c <- par[3]
  ame <- par[4]
  ayou <- par[5]
  bet <- par[6]
  T <- par[7]
  g <- a*exp(lvw) + 2*b*exp(lvw/2) + c
  gp <- a*exp(lvw) + b*exp(lvw/2)
  out <- - (ayou + T/2) * gp / g - ame + bet*exp(-lvw)
  return(out)
}

f <- function(lvw, par){
  logprob <- lpr(lvw, par)
  grad <- gradlpr(lvw, par)
  out <- list(logprob=logprob, grad=grad)
  return(out)
}



V.T <- 100
W.T <- 1
T.T <- 100
w <- rnorm(T.T + 1, 0, sqrt(W.T))
theta <- cumsum(w)
v <- rnorm(T.T, 0, sqrt(V.T))
y <- theta[-1] + v
av <- 5
aw <- 5
bw <- (aw-1)*W.T
bv <- (av-1)*V.T
gam <- gamtrans(theta, W.T)
cgam <- cumsum(gam[-1])
gam0 <- gam[1]
a <- sum( cgam^2 )/2
b <- -sum( (y - gam0) * cgam )/2
c <- bv + sum( (y - gam0)^2 )/2
par <- c(a, b, c, aw, av, bw, T)

delta <- 0.6

M <- 10000
Madapt <- 1000
nutsout <- nuts_da(f, M, Madapt, W.T, par, delta=0.6)

gradroot <- uniroot(gradlpr, c(-100, 100), par=par)
mn <- gradroot$root

g <- a*exp(mn) + 2*b*exp(mn/2) + c
gp <- a*exp(mn) + b*exp(mn/2)
gpp <- a*exp(mn) + b*exp(mn/2)/2

prec <- (av + T/2)*(g*gpp - gp^2)/gpp^2 + bw*exp(-mn)
std <- sqrt(1/prec)

metout <- rep(0, M)
Wold <- log(W.T)
for(m in 1:(Madapt + M)){
  prop <- rnorm(1, mn, std)
  num1 <- lpr(prop, par)
  denom1 <- lpr(Wold, par)
  num2 <- dnorm(Wold, mn, std, log=TRUE)
  denom2 <- dnorm(prop, mn, std, log=TRUE)
  u <- runif(1,0,1)
  if(log(u)< num1 - denom1 + num2 - denom2)
    Wold <- prop
  if(m>Madapt)
    metout[m-Madapt] <- Wold
}

sams <- mcmc(data.frame(nuts=nutsout$samples, met=metout))
rejectionRate(sams)
summary(sams)
samsD <- data.frame(sams)

par(mfrow=c(2,1))
plot(ts(sams[,1]))
plot(ts(sams[,2]))

library(ggplot2)
library(gridExtra)

lims <- c(min(sams), max(sams))
normpr <- integrate(pr, -Inf, Inf, par=par, rel.tol=10^(-50))
C <- normpr$value
x <- seq(lims[1],lims[2],length.out=1000)
curvedat <- data.frame(x=x, y=pr(x,par)/C, prop=dnorm(x,mn,std))

p1 <- ggplot(samsD, aes(x=nuts)) +
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_line(aes(x=x, y=y), data=curvedat, color="red")
p2 <- ggplot(samsD, aes(x=met)) +
  geom_histogram(aes(y=..density..), binwidth=0.01) +
  geom_line(aes(x=x, y=y), data=curvedat, color="red")
grid.arrange(p1, p2)






p1 <- qplot(nuts, geom="histogram", binwidth=0.01, xlim=lims, data=samsD, margins=FALSE) + theme(axis.text.y = element_blank())
p2 <- qplot(met, geom="histogram", binwidth=0.01, xlim=lims, data=samsD, margins=FALSE) + theme(axis.text.y = element_blank())
x <- seq(lims[1],lims[2],length.out=1000)
curvedat <- data.frame(x=x, y=exp(lpr(x, par)), prop=dnorm(x,mn,std))
p3 <- qplot(x, y, data=curvedat, geom="line", margins=FALSE) + theme(axis.text.y = element_blank())
p4 <- qplot(x, prop, data=curvedat, geom="line", margins=FALSE) + theme(axis.text.y = element_blank())


grid.arrange(p1, p3, p2, p4, ncol=2)
