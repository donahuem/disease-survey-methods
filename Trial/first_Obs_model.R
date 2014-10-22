## my first jags file ##

model {
  for (i in 1:N) { 
  K[i] ~ dbin(theta, Nc[i]) # Number of diseased colonies is distributed binomially with prevalence and total number of colonies counted
  theta <- ilogit(b0+b[1]*X[i]+b[2]*Y[i]) # transforms prevalence from 0-1
  }
  
  #priors
  #theta ~ dbeta(1,1)
  b0 ~ dnorm(mu, tau)
  mu ~ dnorm(0,1)
  tau <- 1/pow(sigma, 2)
  sigma ~ dunif(0,100)
  #b0 ~ dnorm(0, 1.0E-6)
  #b[1] ~ dnorm(0, 1.0E-6)
  #b[2] ~ dnorm(0, 1.0E-6)
  #tau ~ dgamma(0.001, 0.001)
  #sigma <- 1.0/sqrt(tau)
}
