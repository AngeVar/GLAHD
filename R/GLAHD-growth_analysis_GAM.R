#-------------------------------------------------------------------------------------
#- This script provides an alternative curve fit for growth analysis via a GAM.
#- 
#-------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")
source("R/fitGAM/derivSimulCI.R")
source("R/fitGAM/plotCIdate.R")
source("R/fitGAM/smoothplot.R")

# for GAM
library(mgcv)








#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1


#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------








#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- fit gam model to each plant one at a time, predict the derivatives
growth.l <- split(dat3,dat3$Code)
log3fits <- output <- data.out <- list()
kgam=5
gamfits <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  
  #- fit the gam
  g <- gam(lnTotMass ~ s(Time, k=kgam), data=tofit)

  #- plot fit
  #smoothplot(Time, lnTotMass, data=tofit, kgam=kgam)

  #- create a vector of "dates" on which to estimate the derivative 
  dates <- seq(min(tofit$Time), max(tofit$Time), by=1)
  
  #- extract the derivative 
  fd <- derivSimulCI(g, samples = 10000, n=length(dates))
  dydt <- fd[[1]]$deriv[,1]

  #- put derivatives into a list of dataframes
  gamfits[[i]] <- data.frame(Code=tofit$Code[1],Time=dates,dydt=dydt)
  
}

#- merge dataframes, combine with treatment key
gamfits.df <- do.call(rbind,gamfits)
key <- unique(subset(dat3,select=c("Code","Species","Treatment","Location","Taxa","Range")))
gamfits2 <- merge(gamfits.df,key,by=c("Code"),all.x=T)
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- plot some results

g.trt <- summaryBy(dydt~Time+Treatment+Location+Range,data=gamfits2,FUN=mean,keep.names=T)
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(dydt~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)

#maximum of RGR
maxdydt<-gamfits2[ gamfits2$dydt %in% tapply(gamfits2$dydt, gamfits2$Code, max), ]#max RGR per Code
maxdydt$combotrt <- as.factor(paste(maxdydt$Location,maxdydt$Range,maxdydt$Treatment,sep="_"))

ylims=c(0,0.2)
xlims=c(0,35)
colors <- c("blue","red")
windows(10,5);par(mfrow=c(1,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)
boxplot(dydt~Treatment+Taxa,data=maxdydt,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(MaxRGR),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(maxdydt$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)

#- compare statistically 
#Warming increases MAXRGR in temperate but not tropical taxa, slightly more in wide than narrow (P=0.124).
#Warming causes an earlier peak in Tropical Taxa but not Temperate taxa #poly4 did not find this

library(nlme)
maxdydt$Sp_RS_EN <- as.factor(with(maxdydt,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
maxdydt$Prov_Sp_EN <- as.factor(with(maxdydt,paste(Taxa,Species)))
fm.maxdydt<- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt,log(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt,log(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt$residuals[,1])
anova(fm.maxdydt)  
plot(allEffects(fm.maxdydt))

fm.tmaxdydt<- lme(sqrt(Time)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.tmaxdydt,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.tmaxdydt,sqrt(Time)+1~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.tmaxdydt,sqrt(Time)+1~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.tmaxdydt, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.tmaxdydt$residuals[,1])
anova(fm.tmaxdydt) 
plot(allEffects(fm.tmaxdydt))
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- Calculate interval-based RGR for each plant
dat3 <- dat3[with(dat3,order(Code,Date)),] #- make sure the observations for each plant are in order
growth.l <- split(dat3,dat3$Code)          #- split into a list for each plant
for (i in 1:length(growth.l)){
  #- use diff() to calculate RGR. Note that I've added a trailing "NA", or else the vector would be too short to fit
  #-   the data frame.
  growth.l[[i]]$RGR <- c(diff(growth.l[[i]]$lnTotMass),NA)/c(unname(diff(growth.l[[i]]$Date)),NA) 
  growth.l[[i]]$dDiameter <- c(diff(growth.l[[i]]$Diameter),NA) 
  
}
dat4 <- do.call(rbind,growth.l)
windows(30,40)
plotBy(RGR~jitter(as.numeric(Date))|Code,type="b",data=dat4,legend=F);abline(h=0) 
smoothScatter(x=dat4$Time,y=dat4$RGR,ylab="RGR",xlab="Time since treatment began (days)",cex.lab=1.5,
              xlim=c(0,60));abline(h=0)

#- plot RGR via the segment method for each taxa-combination. Overlay model outputs.
dat4$combotrt <- as.factor(paste(dat4$Location,dat4$Range,dat4$Treatment,sep="_"))
dat4.l <- split(dat4,dat4$combotrt)
names(g.trt)[1] <- "Time"
rates.trt.l <- split(g.trt,g.trt$combotrt)
windows(40,40);par(mfrow=c(4,2),mar=c(0.2,0.2,0.2,0.2),oma=c(7,7,3,3))
ylims=c(0,0.3)
for (i in 1:length(dat4.l)){
  smoothScatter(x=dat4.l[[i]]$Time,y=dat4.l[[i]]$RGR,ylab="",xlab="",ylim=ylims,axes=F,
                xlim=c(0,60));abline(h=0)
  lines(rates.trt.l[[i]]$dydt~rates.trt.l[[i]]$Time,lwd=3)
  legend("top",legend=dat4.l[[i]]$combotrt[1],bty="n")
  if(i%%2==1) axis(2,las=1)
  if(i%%2==0) axis(4,las=1)
  if(i>6) axis(1)
}
mtext(side=1,"Time (days after treatment began)",outer=T,line=4,cex=1.5,xpd=NA)
mtext(side=2,"RGR",outer=T,line=3,cex=1.5)

#-----------------------------------------------------------------------------------------
#Relative enhancement of RGR
ber<- summaryBy(dydt~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]/bio[,5]
BER<- bio[,c(1:5,11)]
bersum<- summaryBy(ber~Range+Location+Time, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(1,1.5), xlim=c(0,67))
mtext(text=expression(RGR[W]:RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(1,1.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)

#maximum enhancement of RGR
sber<- subset(BER, Location =="S")
nber<-subset(BER, Location =="N")
sber[ sber$ber %in% tapply(sber$ber, sber$Range, max), ]
nber[ nber$ber %in% tapply(nber$ber, nber$Range, max), ]

#Absolute enhancement of RGR
ber<- summaryBy(dydt~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]-bio[,5]
BER<- bio[,c(1:5,11)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-0.02,0.04), xlim=c(0,67))
mtext(text=expression(RGR[W]-RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-0.02,0.04),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

#--------------------------------------------------------------------------------------------------


