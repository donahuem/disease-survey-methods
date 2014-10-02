## my first jags file ##

model 
{
  theta ~ dbeta(1,1)
  
  #alpha ~ unif(0,1)
  #beta ~ unif(0,1)

  for (i in 1:N) {
    #Krep[i] <- dbinom(Nc[i], theta)
    K[i] ~ dbinom(Nc[i], theta)    
  }
}
