source("scratch.R")

lpr2 <- function(theta){
  out <- dnorm(theta, mean=5, sd=10, log=TRUE)
  return(out)
}

gradlpr2 <- function(theta){
  out <- -(theta - 5)/100
  return(out)
}

f2 <- function(theta, par){
  logp <- lpr(theta)
  grad <- gradlpr(theta)
  out <- list(logp=logp, grad=grad)
  return(out)
}


lpr2 <- function(theta){
  out <- dnorm(theta, mean=5, sd=2, log=TRUE)
  return(out)
}

gradlpr2 <- function(theta){
  out <- -(theta - 5)/4
  return(out)
}

f2 <- function(theta){
  logp <- lpr2(theta)
  grad <- gradlpr2(theta)
  out <- list(logp=logp, grad=grad)
  return(out)
}

M <- 10000
Madapt <- 1000
##theta0 <- c(0.1,0.2)
theta0 <- 10
delta <- 0.6

test <- nuts_da(f2, M, Madapt, theta0, delta)

sam <- test$samples
eps <- test$epsilon
sam2 <- rnorm(10000,5,10)

par(mfrow=c(1,2))
hist(sam)
hist(sam2)

mean(sam)
mean(sam2)

theta <- c(-0.2058817, 0.2003900)
r <- c(-0.5612439, 0.5462734)
u <- 0.8468824
v <- 1
j <- 1
em <- .345
theta0 <- c(0, 0)
r0 <- c(-0.5967585,    0.5808407 )
delmax <- 1000.0000000


test <- BuildTree(theta, r, u, v, j, e, theta0, r0, lpr, gradlpr, delmax)
