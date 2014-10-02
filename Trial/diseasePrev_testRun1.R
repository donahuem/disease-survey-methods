rm(list=ls())

library(rjags)
library(coda)
library(xtable)

disprev <- read.csv('COUCH_HICORDIS_EXPORT.csv')
attach(disprev)

library(plyr)
library(tcltk)

#summarize <- ddply(disprev, .(Date, Site, Transect_Number),
 #                  summarise,
  #                 NumColonies = length(Disease_Type),
   #                Diseased = sum(Disease_Type == "Healthy"))

summarize <- ddply(disprev, .(Transect_Number),
                   summarise,
                   NumColonies = length(Disease_Type),
                   Diseased = sum(Disease_Type == "Healthy"))

newData <- summarize[c(1,2,4:6), ]

ds.data <- list(N=nrow(newData), K=newData$Diseased, Nc=newData$NumColonies)

ds.init <- list(list(theta = .5),list(theta=0.2))

n.adapt=1000
n.iter=20000
n.chains=length(ds.init)

prevalence = jags.model("first_Obs_model.R",
                        data=ds.data,
                        inits=ds.init,
                        n.chains=n.chains,
                        n.adapt=n.adapt)