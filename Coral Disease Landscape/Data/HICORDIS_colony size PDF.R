rm(list=ls()) #cleared workspace
setwd("~/Dropbox/R Files/DZ_methods")#set directory
data<-read.csv("hicordis_baysianData.csv", header=T)
attach(data)

area<-3.14*(Colony_Length/2)^2

a<-mean(Colony_Length)
b<-sd(Colony_Length)

###Histograms and PDFs of colony length for ALL SPECIES

##Plot log normal PDF
h<-hist(Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="All Colonies", ylim=c(0, 300000),) 
xfit<-seq(min(Colony_Length),max(Colony_Length),length=40) 
yfit<-dlnorm(xfit,meanlog=mean(log(Colony_Length)),sdlog=sd(log(Colony_Length))) 
yfit <- yfit*diff(h$mids[1:2])*length(Colony_Length) 
lines(xfit, yfit, col="black", lwd=2)

###Normal PDF of log transformed colony size to confirm that log normal distribution is ok
ll<-log(Colony_Length)
hist(ll)
h<-hist(ll, breaks=10, density=10, col="blue", xlab="Log Colony Size", main="All Colonies", ylim=c(0, 70000),) 
xfit<-seq(min(ll),max(ll),length=40) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)

####Gamma
rate<-mean(Colony_Length)/var(Colony_Length)
shape<-rate*mean(Colony_Length)
scale<-1/rate

h<-hist(Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="All Colonies", ylim=c(0, 1000000),) 
xfit<-seq(min(Colony_Length),max(Colony_Length),length=40) 
Sid_GPDF_yfit<-dgamma(xfit,scale=scale,shape=shape) #PDF given shape and rate
yfit <- Sid_GPDF_yfit*diff(h$mids[1:2])*length(Colony_Length) #transform pdf to fit over histogram
lines(xfit, yfit, col="black", lwd=2)


#### Gamma Area
rate<-mean(area)/var(area)
shape<-rate*mean(area)
scale<-1/rate

h<-hist(area, breaks=10, density=10, col="blue", xlab="Area", main="All Colonies", ylim=c(0, 200000),) 
xfit<-seq(min(area),max(area),length=40) 
area_GPDF_yfit<-dgamma(xfit,scale=scale,shape=shape) #PDF given shape and rate
yfit <- area_GPDF_yfit*diff(h$mids[1:2])*length(area) #transform pdf to fit over histogram
lines(xfit, yfit, col="black", lwd=2)



##############
##########
#Look at histograms and PDFs for 4 most abundant coral species
plob<-subset(data,Species=="Porites_lobata")
pcom<-subset(data,Species=="Porites_compressa")
mcap<-subset(data,Species=="Montipora_capitata")
pmean<-subset(data,Species=="Pocillopora_meandrina")

par(mfrow=c(2,2))

##############
#Porites lobata
#Log normal
h<-hist(plob$Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="Porites lobata", ylim=c(0, 200000),) 
xfit<-seq(min(plob$Colony_Length),max(plob$Colony_Length),length=40) 
yfit<-dlnorm(xfit,meanlog=mean(log(plob$Colony_Length)),sdlog=sd(log(plob$Colony_Length))) 
yfit <- yfit*diff(h$mids[1:2])*length(plob$Colony_Length) 
lines(xfit, yfit, col="black", lwd=2)

#Normal 
ll<-log(plob$Colony_Length)
h<-hist(ll, breaks=10, density=10, col="blue", xlab="Log Colony Size", main="All Colonies", ylim=c(0, 25000),) 
xfit<-seq(min(ll),max(ll),length=40) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)

#Porites compressa
h<-hist(pcom$Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="Porites compressa", ylim=c(0, 30000),) 
xfit<-seq(min(pcom$Colony_Length),max(pcom$Colony_Length),length=40) 
yfit<-dlnorm(xfit,meanlog=mean(log(pcom$Colony_Length)),sdlog=sd(log(pcom$Colony_Length))) 
yfit <- yfit*diff(h$mids[1:2])*length(pcom$Colony_Length) 
lines(xfit, yfit, col="black", lwd=2)

#Normal 
ll<-log(pcom$Colony_Length)
h<-hist(ll, breaks=10, density=10, col="blue", xlab="Log Colony Size", main="All Colonies", ylim=c(0, 10000),) 
xfit<-seq(min(ll),max(ll),length=40) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)


#Pocillopora meandrina
h<-hist(pmean$Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="Pocillopora meandrina", ylim=c(0, 30000),) 
xfit<-seq(min(pmean$Colony_Length),max(pmean$Colony_Length),length=40) 
yfit<-dlnorm(xfit,meanlog=mean(log(pmean$Colony_Length)),sdlog=sd(log(pmean$Colony_Length))) 
yfit <- yfit*diff(h$mids[1:2])*length(pmean$Colony_Length) 
lines(xfit, yfit, col="black", lwd=2)

#Normal 
ll<-log(pmean$Colony_Length)
h<-hist(ll, breaks=10, density=10, col="blue", xlab="Log Colony Size", main="Pocillopora meandrina", ylim=c(0, 10000),) 
xfit<-seq(min(ll),max(ll),length=40) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)

#Montipora capitata
h<-hist(mcap$Colony_Length, breaks=10, density=10, col="blue", xlab="Colony Size", main="Montipora capitata", ylim=c(0, 60000),) 
xfit<-seq(min(mcap$Colony_Length),max(mcap$Colony_Length),length=40) 
yfit<-dlnorm(xfit,meanlog=mean(log(mcap$Colony_Length)),sdlog=sd(log(mcap$Colony_Length))) 
yfit <- yfit*diff(h$mids[1:2])*length(mcap$Colony_Length) 
lines(xfit, yfit, col="black", lwd=2)

#Normal 
ll<-log(mcap$Colony_Length)
h<-hist(ll, breaks=10, density=10, col="blue", xlab="Log Colony Size", main="Montipora capitata", ylim=c(0, 10000),) 
xfit<-seq(min(ll),max(ll),length=40) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)
