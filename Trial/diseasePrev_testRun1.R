rm(list=ls())
setwd("C:/Users/Iain/Desktop/Jamie")

library(rjags)
library(coda)
library(xtable)

Num_Colonies <- sample(seq(30,60,1), 50, replace=TRUE)
Diseased <- rbinom(Num_Colonies, Num_Colonies, .2)
Depth <- sample(seq(5,40,1), 50, replace=TRUE)
SizeFreq <- sample(seq(5,160,1), 50, replace=TRUE)
simData <- data.frame(Num_Colonies, Diseased, Depth, SizeFreq)

ds.data <- list(N=nrow(simData), K=simData$Diseased, Nc=simData$Num_Colonies, X=simData$Depth, Y=simData$SizeFreq)
ds.init <- list(list(theta = .5, b = .1), list(theta=.2, b=.4))

n.adapt=1000
n.iter=20000
n.chains=length(ds.init)

prevalence = jags.model("first_Obs_model.R",
                        data=ds.data,
                        inits=ds.init,
                        n.chains=n.chains,
                        n.adapt=n.adapt)

load.module("dic")
prev.jags <- jags.samples(prevalence, variable.names=c("theta"), n.iter=n.iter, n.thin=1)
prev.coda <- coda.samples(prevalence, variable.names=c("deviance", "theta"), n.iter=n.iter)
prev.df=as.data.frame(rbind(prev.coda))

#some plots (I thinned these just so they would plot faster)
xyplot(window(prev.coda, thin=10))
densityplot(window(prev.coda, thin=10))
autocorr.plot(prev.coda) 
  hist(prev.jags$theta, freq=TRUE, breaks=500)
