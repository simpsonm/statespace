library(ggplot2)

set.seed(1235341)
n <- 1000
sig2s <- c(.001, 1000)
dat <- data.frame(mu=NULL, t=NULL, sig2=NULL, DA=NULL)
for(j in 1:2){
  sig2 <- sig2s[j]
  theta1 <- sqrt(sig2)*rnorm(1)
  y <- theta1 + rnorm(1)
  out1 <- rep(0,n)
  out2 <- out1
  mu <- 0
  for(i in 1:n){
    thmn <- (mu + sig2*y)/(1 + sig2)
    thsd <- sqrt(sig2 / (1 + sig2))
    theta1 <- rnorm(1, thmn, thsd)
    mu <- rnorm(1, theta1, sqrt(sig2))
    out1[i] <- mu
  }
  mu <- 0
  for(i in 1:n){
    thmn <- sig2*(y - mu)/(1 + sig2)
    thsd <- sqrt(sig2 / (1 + sig2))
    theta2 <- rnorm(1, thmn, thsd)
    mu <- rnorm(1,y-theta2, 1)
    out2[i] <- mu
  }
  tempdat <- data.frame(mu=c(out1,out2), t=rep(1:n,2),sig2=sig2,DA=c(rep("theta[1]",n),rep("theta[2]",n)))
  dat <- data.frame(rbind(dat,tempdat))
}


p1 <- ggplot(data=dat[dat$sig2==sig2s[1],], aes(x=t, y=mu)) +
    geom_line() +
    facet_grid(facets=.~DA, labeller=label_parsed) +
    ylab(expression(mu))

p2 <- ggplot(data=dat[dat$sig2==sig2s[2],], aes(x=t, y=mu)) +
    geom_line() +
    facet_grid(facets=.~DA, labeller=label_parsed) +
    ylab(expression(mu))

ggsave(filename="trace1.pdf", plot=p1, width=6, height=2)
ggsave(filename="trace2.pdf", plot=p2, width=6, height=2)

