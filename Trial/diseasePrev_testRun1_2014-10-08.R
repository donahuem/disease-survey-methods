rm(list=ls()) # clear workspace

# load libraries
library(rjags)
library(coda)
library(xtable)

belts <- read.csv('model1data.csv') # read document

#Create some simulated data where we know the prevalence
#Num_Colonies <- sample(seq(30,60,1),50,replace=TRUE)
#Diseased <- rbinom(Num_Colonies,Num_Colonies,.2) #prevalence of 0.2
#simData <- data.frame(Num_Colonies,Diseased)
# rm(Num_Colonies,Diseased) #remove these local variables to avoid confusion

######
ds.data <- list(N=nrow(belts), K=belts$Diseased, Nc=belts$NumColonies)
#ds.data <- list(N=nrow(simData), K=simData$Diseased, Nc=simData$Num_Colonies)
ds.init <- list(list(theta = .5),list(theta=0.2))
n.adapt=1000
n.iter=20000
n.chains=length(ds.init)
prevalence = jags.model("first_Obs_model.R",
                        data=ds.data,
                        inits=ds.init,
                        n.chains=n.chains,
                        n.adapt=n.adapt)
load.module("dic") # monitor to calculate summary statistics
prev.jags <- jags.samples(prevalence,variable.names=c("theta"), n.iter=n.iter, n.thin=1) # thin defines the monitor interval recording every nth value
prev.coda <- coda.samples(prevalence,variable.names=c("deviance","theta"),n.iter=n.iter) # stores monitored values in coda file
prev.df=as.data.frame(rbind(prev.coda)) 
#some plots (I thinned these just so they would plot faster)
xyplot(window(prev.coda, thin=10))
densityplot(window(prev.coda, thin=10))
autocorr.plot(prev.coda) #(no autocorrelation issues in the chain - no need to thin) ??
hist(prev.jags$theta,freq=TRUE, breaks=500)
