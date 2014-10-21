## my first jags file ##

model {
  
  for (i in 1:N) {
  K[i] ~ dbin(theta, Nc[i])    #syntax was wrong on this line.  dbin(p,N) not dbinom(N,p)
  theta <- ilogit(beta[0]+(beta[1]*X[i])+(beta[2]*Y[i]))
  }

  beta ~ dnorm(mu, tau)
  #mu ~ dnorm(0,1) should this be a prior? Should this be mu <- mean(X)?
  tau <- 1 / pow(sigma, 2)
  sigma ~ dunif(0,100)
  #tau ~ dgamma(0.001, 0.001)
  #sigma <- 1.0/sqrt(tau)
}
