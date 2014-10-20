mcfaphi <- function(y, V, W, m0, C0, F, N, J){
  C0inv <- solve(C0)
  Vinv <- solve(V)
  Winv <- solve(W)
  tF <- t(F)
  O0 <- C0inv + Winv
  OT <- tF%*%Vinv%*%F + Winv
  Ot <- OT + Winv
  Ott1 <- -Winv
  o <- matrix(0, nrow=(J+1), ncol=(N+1))
  m <- o
  E <- array(0,c(J+1,J+1,N+1))
  o[,1] <- C0inv%*%m0
  E[,,1] <- solve(O0)
  m[,1] <- E[,,1]%*%o[,1]
  for(t in 1:N){
    o[,t+1] <- tF%*%Vinv%*%y[,t]
    E[,,t+1] <- solve(Ot - Ott1%*%E[,,t]%*%Ott1)
    m[,t+1] <- E[,,t+1]%*%(o[,t+1] - Ott1%*%m[,t])
  }
  phi <- matrix(rnorm((J+1)*(N+1)),ncol=(N+1))
  phi[,N+1] <- m[,N+1] + t(chol(E[,,N+1]))%*%phi[,N+1]
  for(t in N:1){
    phi[,t] <- m[,t] - E[,,t]%*%Ott1%*%phi[,t+1] + t(chol(E[,,t]))%*%phi[,t]
  }
  return(phi)
}
