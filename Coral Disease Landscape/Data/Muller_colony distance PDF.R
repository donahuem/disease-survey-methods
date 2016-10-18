rm(list=ls()) #cleared workspace
setwd("~/Dropbox/R Files/DZ_methods")#set directory
dis<-read.csv("Mullercolonydistance.csv", header=T)
attach(dis)

#Subset Montastrea colonies
new<-subset(Data,Species=="Montastrea")
hist(new$distance,breaks=14)


##Plot log normal PDF for Montastrea
h<-hist(new$distance, breaks=10, density=10, col="blue", xlab="Distance", main="Log Normal PDF_Montastrea", ylim=c(0, 40),) 
xfit<-seq(min(new$distance),max(new$distance),length=100) 
Mon_PDF_yfit<-dlnorm(xfit,meanlog=mean(log(new$distance)),sdlog=sd(log(new$distance))) 
yfit <- Mon_PDF_yfit*diff(h$mids[1:2])*length(new$distance) 
lines(xfit, yfit, col="black", lwd=2)

##Log transform and plot normal PDF confirms log normal distribution is appropriate
ll<-log(new$distance)
hist(ll)

h<-hist(ll, breaks=6, col="lightgray", xlab="Log Distance", main="Montastrea", ylim=c(0,30),) 
xfit<-seq(min(ll),max(ll),length=100) 
yfit<-dnorm(xfit,mean=mean(ll),sd=sd(ll)) 
yfit <- yfit*diff(h$mids[1:2])*length(ll) 
lines(xfit, yfit, col="black", lwd=2)

####
#Subset Siderastrea colonies
new<-subset(Data,Species=="Siderastrea")
hist(new$distance,breaks=14)


##Plot normal PDF for Siderastrea (It's not a great fit, but definiatley not log normal)
h<-hist(new$distance, breaks=10, density=10, col="blue", xlab="Distance", main="Normal PDF_Siderastrea", ylim=c(0, 25),) 
xfit<-seq(min(new$distance),max(new$distance),length=100) 
Sid_NPDF_yfit<-dnorm(xfit,mean=mean(new$distance),sd=sd(new$distance)) #PDF given mean and sd
yfit <- Sid_NPDF_yfit*diff(h$mids[1:2])*length(new$distance) #transform pdf to fit over histogram
lines(xfit, yfit, col="black", lwd=2)


##Plot gamma PDF for Siderastrea (I think this fits the data better than a normal distribution)
rate<-mean(new$distance)/var(new$distance)
shape<-rate*mean(new$distance)

h<-hist(new$distance, breaks=10, density=10, col="blue", xlab="Distance", main="Gamma PDF_Siderastrea", ylim=c(0, 25),) 
xfit<-seq(min(new$distance),max(new$distance),length=100) 
Sid_GPDF_yfit<-dgamma(xfit,rate=rate,shape=shape) #PDF given shape and rate
yfit <- Sid_GPDF_yfit*diff(h$mids[1:2])*length(new$distance) #transform pdf to fit over histogram
lines(xfit, yfit, col="black", lwd=2)

#extract PDFs
pdfcomp<-cbind(xfit,Mon_PDF_yfit,Sid_NPDF_yfit,Sid_GPDF_yfit)
write.csv(pdfcomp,"PDF_Mullercolonydist.csv")
