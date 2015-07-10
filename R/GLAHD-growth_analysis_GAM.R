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
  
  #- get the predicted mass
  newDF <- data.frame(Time=dates) ## needs to be a data frame for predict
  X0 <- exp(predict(g, newDF))    ## exp() needed to convert from ln(mass) to mass
  
  #- put mass into the dataframe
  gamfits[[i]]$predMass <- X0
}

#- merge dataframes, combine with treatment key
gamfits.df <- do.call(rbind,gamfits)
key <- unique(subset(dat3,select=c("Code","Species","Treatment","Location","Taxa","Range")))
gamfits2 <- merge(gamfits.df,key,by=c("Code"),all.x=T)
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

gamfits2$AGR <- gamfits2$dydt*gamfits2$predMass



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- plot some results

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=mean,keep.names=T)
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(dydt~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)
plotBy(AGR~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
plotBy(predMass~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)

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

#-------------------------------------------------------------------------------------
#- plot mass and RGR over time for each taxa-treatment combination, overlay data and model output
dat4$RGR_time <- dat4$Time+5
dat4.l.taxa <- split(dat4,dat4$Taxa)
rates.taxa <- summaryBy(predMass+dydt+AGR~Time+Treatment+Location+Range+Taxa,data=gamfits2,FUN=mean,keep.names=T)
rates.taxa.l <- split(rates.taxa,rates.taxa$Taxa)
rates.l <- split(gamfits2,gamfits2$Taxa)

# loop over each taxa, plot a set of four figures, export to a pdf
pdf(file="C:/Repos/GLAHD/Output/RGR_M_Taxa_filtered_GAM.pdf")

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
       widths=rep(1,4), heights=rep(1,4))
for (i in 1:length(dat4.l.taxa)){
  ylims=c(0,max(dat4.l.taxa[[i]]$TotMass+5))
  #- plot home mass over time
  plot(TotMass~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
  plotBy(predMass~Time|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(predMass~Time,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
  
  #- plot warmed mass over time
  plot(TotMass~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Warmed",sep='_'))
  plotBy(predMass~Time|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(predMass~Time,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
  
  #- plot home RGR over time
  ylims=c(0,max(dat4.l.taxa[[i]]$RGR+0.1,na.rm=T))
  
  plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",xlim=c(0,60),
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
  plotBy(dydt~Time|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(dydt~Time,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
  
  #- plot warmed RGR over time
  plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",xlim=c(0,60),
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Home",sep='_'))
  plotBy(dydt~Time|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(dydt~Time,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
  
}
dev.off()

#-------------------------------------------------------------------------------------------------------
#- Do analysis suggested by Mark, where interval-based RGR estimates are compared to model estimates
#-  The basic question: is a 4th order polynomial model sufficient to describe the data,
#-      particularly the peaked RGR response to warming during the second growth interval?


#- the interval-based RGR estimate for the segment in question
rgr.seg2 <- subset(dat4,Date==as.Date("2014-11-17"))
names(rgr.seg2)[18] <- "RGR.seg"

#- get the modeled RGR for each day of the growth interval in question for each plant
#remember that time "1" is 2014-11-07
table(dat2$Date) # second growth interval is from 2014-11-17 to 2014-11-26, which is "Time" 11 to 20
rgr.mod.seg2 <- subset(gamfits2,Time>=11 & Time <=20)
rgr.mod.seg2.plant <- summaryBy(predMass+AGR+dydt~Code+Species+Treatment+Location+Taxa+Range,data=rgr.mod.seg2,keep.names=T)

#- put them together
rgr2 <- merge(rgr.seg2,rgr.mod.seg2.plant,by=c("Code","Species","Treatment","Location","Taxa","Range"))

#- plot both estimates of RGR for each taxa

#plot RGR seg, then RGR modeled
ylims=c(0,0.2)
xlims=c(0,35)
colors <- c("blue","red")
windows(10,10);par(mfrow=c(2,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)

#RGR segmented
boxplot(RGR.seg~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[seg]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=F,las=2,cex.axis=1.5)
abline(v=16.4)

#RGR modeled
boxplot(dydt~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[mod]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(rgr2$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)


#- by-plot
windows(20,20)
plotBy(dydt~RGR.seg|Location,data=rgr2,xlab="RGR-seg",ylab="RGR-mod",xlim=c(0,0.22),ylim=c(0,0.22),cex.lab=2)
abline(0,1)

#- compare statistical analyes
library(nlme)
rgr2$Sp_RS_EN <- as.factor(with(rgr2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rgr2$Prov_Sp_EN <- as.factor(with(rgr2,paste(Taxa,Species)))
fm.rgr.mod <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
plot(fm.rgr.mod,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.rgr.mod,log(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.rgr.mod,log(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.rgr.mod, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.rgr.mod$residuals[,1])
anova(fm.rgr.mod)  

fm.rgr.seg <- lme(RGR.seg~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
plot(fm.rgr.seg,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.rgr.seg,RGR.seg~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.rgr.seg,RGR.seg~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.rgr.seg, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.rgr.seg$residuals[,1])
anova(fm.rgr.seg)  


#-----------------------------------------------------------------------------------------
#Relative enhancement of biomass
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]/bio[,5]
BER<- bio[,c(1:5,11)]
bersum<- summaryBy(ber~Range+Location+Time, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,2), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,2),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)

#Absolute enhancement of biomass
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]-bio[,5]
BER<- bio[,c(1:5,11)]
bersum<- summaryBy(ber~Range+Location+Time, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-5,20), xlim=c(0,67))
mtext(text=expression(M[W]-M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-5,20),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

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
       legend=F,type="o", main="South", ylim=c(0.5,1.5), xlim=c(0,67))
mtext(text=expression(RGR[W]:RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,1.5),yaxt='n', xlim=c(0,67))
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
#LAR

rate <- merge(gamfits2,dat2,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,19)]# merge with dat2 to get LA
rate$LAR <- rate$leafArea/rate$predMass
rate.trt <- summaryBy(LAR+predMass+dydt+AGR~Time+Treatment+Location+Range,data=rate,FUN=mean,keep.names=T)
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(LAR~Time|combotrt,type='o',data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)
plotBy(LAR~dydt|combotrt, type='o', data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="bottomright",pch=15, ylim=c(60,144), xlab="RGR")


