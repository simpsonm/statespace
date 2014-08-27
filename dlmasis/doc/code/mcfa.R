## Function for applying the MCFA algorithm to simulate
## the latent states (theta's) in the local level model
## y = (y_1, y_2, ... , y_T) is a vector of the observed time series
## V is the observation level variance
## W is the system level variance
## m0 is the prior mean for theta_0
## C0 is the prior variance for theta_0
## returns a draw from p(theta|V,W,y)
mcfathsmooth <- function(y, V, W, m0, C0){
  n <- length(y)
  O0 <- 1/C0 + 1/W
  Ot <- 1/V + 2/W
  OT <- 1/V + 1/W
  Ott1 <- -1/W
  a <- rep(0,n+1)
  E <- a
  m <- a
  ##create the parameters of the distribution of theta
  ## E is Sigma in the definition from the MCFA
  a[1] <- m0/C0
  E[1] <- 1/O0
  m[1] <- E[1]*a[1]
  for(t in 1:(n-1)){
    a[t+1] <- y[t]/V
    E[t+1] <- 1/(Ot - Ott1^2*E[t])
    m[t+1] <- E[t+1]*(a[t+1]-Ott1*m[t])
  }
  a[n+1] <- y[n]/V
  E[n+1] <- 1/(OT - Ott1^2*E[n])
  m[n+1] <- E[n+1]*(a[n+1]-Ott1*m[n])
  ##initialize output vector with N(0,1) draws
  theta <- rnorm(n+1) 
  ## Recursively obtain each theta_t from the initialization and the later theta_t's
  theta[n+1] <- m[n+1] + sqrt(E[t+1])*theta[n+1]
  for(t in n:1){
    theta[t] <- m[t] - E[t]*Ott1*theta[t+1] + sqrt(E[t])*theta[t]
  }
  return(theta)
}

## Function for applying the MCFA algorithm to simulate
## the scaled errors (psi's) in the local level model
## y = (y_1, y_2, ... , y_T) is a vector of the observed time series
## V is the observation level variance
## W is the system level variance
## m0 is the prior mean for theta_0
## C0 is the prior variance for theta_0
## returns a draw from p(psi|V,W,y)
mcfapssmooth <- function(y, V, W, m0, C0){
  n <- length(y)
  O0 <- 1/W + 1/C0
  Ot <- 2*V/W + 1
  OT <- V/W + 1
  sV <- sqrt(V)
  O01 <- sV/W
  Ost <- -V/W
  a <- rep(0,n+1)
  E <- a
  m <- a
  ##create the parameters of the distribution of theta
  ## E is Sigma in the definition from the MCFA
  a[1] <- y[1]/W + m0/C0
  a[2] <- sV/W*(2*y[1] - y[2])
  E[1] <- 1/O0
  m[1] <- E[1]*a[1]
  E[2] <- 1/(Ot - O01^2*E[1])
  m[2] <- E[2]*(a[2] - O01*m[1])
  for(t in 2:n){
    a[t+1] <- sV/W*(2*y[t] - y[t-1] - y[t+1])
    E[t+1] <- 1/(Ot - Ost^2*E[t])
    m[t+1] <- E[t+1]*(a[t+1] - Ost*m[t])
  }
  a[n+1] <- sV/W*(y[n] - y[n-1])
  E[n+1] <- 1/(OT - Ost^2*E[n])
  m[n+1] <- E[n+1]*(a[n+1] - Ost*m[n])
  ##initialize output vector with N(0,1) draws
  psi <- rnorm(n+1)
  ## Recursively obtain each psi_t from the initialization and the later psi_t's
  psi[n+1] <- m[n+1] + sqrt(E[t+1])*psi[n+1]
  for(t in n:2){
    psi[t] <- m[t] - E[t]*Ost*psi[t+1] + sqrt(E[t])*psi[t]
  }
  psi[1] <- m[1] - E[1]*O01*psi[2] + sqrt(E[1])*psi[1]
  return(psi)
}
