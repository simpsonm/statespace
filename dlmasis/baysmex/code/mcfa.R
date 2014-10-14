mcfatheta <- function(y, V, W, mu, W0, beta){
  n <- length(y)
  O0 <- 1/W0 + beta^2/W
  Ot <- 1/V + (1 + beta^2)/W
  OT <- 1/V + 1/W
  Ott1 <- -beta/W
  o <- rep(0,n+1)
  E <- o
  m <- o
  o[1] <- mu[1]/W0 + mu[2]*beta/W
  E[1] <- 1/O0
  m[1] <- E[1]*o[1]
  for(t in 1:(n-1)){
    o[t+1] <- y[t]/V + (mu[t+1] - mu[t+2]*beta)/W
    E[t+1] <- 1/(Ot - Ott1^2*E[t])
    m[t+1] <- E[t+1]*(o[t+1]-Ott1*m[t])
  }
  o[n+1] <- y[n]/V + mu[n+1]/W
  E[n+1] <- 1/(OT - Ott1^2*E[n])
  m[n+1] <- E[n+1]*(o[n+1]-Ott1*m[n])
  theta <- rnorm(n+1)
  theta[n+1] <- m[n+1] + sqrt(E[t+1])*theta[n+1]
  for(t in n:1){
    theta[t] <- m[t] - E[t]*Ott1*theta[t+1] + sqrt(E[t])*theta[t]
  }
  return(theta)
}


mcfapsi <- function(y, V, W, mu, W0, beta){
  n <- length(y)
  O0 <- 1/W0 + beta^2/W
  Ot <- 1 + V/W*(1 + beta^2)
  OT <- 1 + V/W
  O10 <- -sqrt(V)*beta/W
  Ott1 <- -V*beta/W
  o <- rep(0,n+1)
  E <- o
  m <- o
  o[1] <- mu[1]/W0 + (y[1] - mu[2])*beta/W
  E[1] <- 1/O0
  m[1] <- E[1]*o[1]
  o[2] <- sqrt(V)/W*(mu[2] - y[1] + beta*(y[2]-beta*y[1]-mu[3]))
  E[2] <- 1/(Ot - O10^2*E[1])
  m[2] <- E[2](o[2] - O10*m[1])
  for(t in 2:(n-1)){
    o[t+1] <- sqrt(V)/W*(mu[t+1] + beta*y[t-1] - y[t] + beta*(y[t+1] - beta*y[t] - mu[t+2]))
    E[t+1] <- 1/(Ot - Ott1^2*E[t])
    m[t+1] <- E[t+1]*(o[t+1]-Ott1*m[t])
  }
  o[n+1] <- sqrt(V)/W*(mu[n+1] + beta*y[n-1] - y[n])
  E[n+1] <- 1/(OT - Ott1^2*E[n])
  m[n+1] <- E[n+1]*(o[n+1]-Ott1*m[n])
  psi <- rnorm(n+1)
  psi[n+1] <- m[n+1] + sqrt(E[t+1])*psi[n+1]
  for(t in n:2){
    psi[t] <- m[t] - E[t]*Ott1*psi[t+1] + sqrt(E[t])*psi[t]
  }
  psi[1] <-  m[1] - E[t]*O10*psi[t+1] + sqrt(E[t])*psi[t]
  return(psi)
}


mcfamu <- function(W, W0, U, U0, m0, theta, beta){
  n <- nrow(theta) - 1
  J <- ncol(theta)
  O0 <- 1/U0 + sum(1/W0) + 1/U
  Ot <- sum(1/W) + 2/U
  OT <- sum(1/W) + 1/U
  Ott1 <- -1/U
  o <- rep(0,n+1)
  E <- o
  m <- o
  o[1] <- m0/U0 + sum(theta[1,]/W0)
  E[1] <- 1/O0
  m[1] <- E[1]*o[1]
  for(t in 1:(n-1)){
    o[t+1] <- sum((theta[t+1,] - theta[t,]*beta)/W)
    E[t+1] <- 1/(Ot - Ott1^2*E[t])
    m[t+1] <- E[t+1]*(o[t+1]-Ott1*m[t])
  }
  o[n+1] <- sum((theta[n+1,] - theta[n,]*beta)/W)
  E[n+1] <- 1/(OT - Ott1^2*E[n])
  m[n+1] <- E[n+1]*(o[n+1]-Ott1*m[n])
  mu <- rnorm(n+1)
  mu[n+1] <- m[n+1] + sqrt(E[t+1])*mu[n+1]
  for(t in n:1){
    mu[t] <- m[t] - E[t]*Ott1*mu[t+1] + sqrt(E[t])*mu[t]
  }
  return(mu)
}
