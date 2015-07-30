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
gamfits2$AGR <- gamfits2$dydt*gamfits2$predMass
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#Compares model outputs against Interval-based data



# #- Calculate interval-based RGR for each plant
# dat3 <- dat3[with(dat3,order(Code,Date)),] #- make sure the observations for each plant are in order
# growth.l <- split(dat3,dat3$Code)          #- split into a list for each plant
# for (i in 1:length(growth.l)){
#   #- use diff() to calculate RGR. Note that I've added a trailing "NA", or else the vector would be too short to fit
#   #-   the data frame.
#   growth.l[[i]]$RGR <- c(diff(growth.l[[i]]$lnTotMass),NA)/c(unname(diff(growth.l[[i]]$Date)),NA) 
#   growth.l[[i]]$dDiameter <- c(diff(growth.l[[i]]$Diameter),NA) 
#   
# }
# dat4 <- do.call(rbind,growth.l)
# windows(30,40)
# plotBy(RGR~jitter(as.numeric(Date))|Code,type="b",data=dat4,legend=F);abline(h=0) 
# smoothScatter(x=dat4$Time,y=dat4$RGR,ylab="RGR",xlab="Time since treatment began (days)",cex.lab=1.5,
#               xlim=c(0,60));abline(h=0)
# 
# #- plot RGR via the segment method for each taxa-combination. Overlay model outputs.
# dat4$combotrt <- as.factor(paste(dat4$Location,dat4$Range,dat4$Treatment,sep="_"))
# dat4.l <- split(dat4,dat4$combotrt)
# names(g.trt)[1] <- "Time"
# rates.trt.l <- split(g.trt,g.trt$combotrt)
# windows(40,40);par(mfrow=c(4,2),mar=c(0.2,0.2,0.2,0.2),oma=c(7,7,3,3))
# ylims=c(0,0.3)
# for (i in 1:length(dat4.l)){
#   smoothScatter(x=dat4.l[[i]]$Time,y=dat4.l[[i]]$RGR,ylab="",xlab="",ylim=ylims,axes=F,
#                 xlim=c(0,60));abline(h=0)
#   lines(rates.trt.l[[i]]$dydt~rates.trt.l[[i]]$Time,lwd=3)
#   legend("top",legend=dat4.l[[i]]$combotrt[1],bty="n")
#   if(i%%2==1) axis(2,las=1)
#   if(i%%2==0) axis(4,las=1)
#   if(i>6) axis(1)
# }
# mtext(side=1,"Time (days after treatment began)",outer=T,line=4,cex=1.5,xpd=NA)
# mtext(side=2,"RGR",outer=T,line=3,cex=1.5)
# 
# #- plot AGR via the segment method for each taxa-combination. Overlay model outputs.
# dat4$AGR <- with(dat4, RGR*TotMass)
# dat4$combotrt <- as.factor(paste(dat4$Location,dat4$Range,dat4$Treatment,sep="_"))
# dat4.l <- split(dat4,dat4$combotrt)
# names(g.trt)[1] <- "Time"
# rates.trt.l <- split(g.trt,g.trt$combotrt)
# windows(40,40);par(mfrow=c(4,2),mar=c(0.2,0.2,0.2,0.2),oma=c(7,7,3,3))
# ylims=c(0,5)
# for (i in 1:length(dat4.l)){
#   smoothScatter(x=dat4.l[[i]]$Time,y=dat4.l[[i]]$AGR,ylab="",xlab="",ylim=ylims,axes=F,
#                 xlim=c(0,60));abline(h=0)
#   lines(rates.trt.l[[i]]$AGR~rates.trt.l[[i]]$Time,lwd=3)
#   legend("top",legend=dat4.l[[i]]$combotrt[1],bty="n")
#   if(i%%2==1) axis(2,las=1)
#   if(i%%2==0) axis(4,las=1)
#   if(i>6) axis(1)
# }
# mtext(side=1,"Time (days after treatment began)",outer=T,line=4,cex=1.5,xpd=NA)
# mtext(side=2,"AGR",outer=T,line=3,cex=1.5)
# 
# #-------------------------------------------------------------------------------------
# #- plot mass and RGR over time for each taxa-treatment combination, overlay data and model output
# dat4$RGR_time <- dat4$Time+5
# dat4.l.taxa <- split(dat4,dat4$Taxa)
# rates.taxa <- summaryBy(predMass+dydt+AGR~Time+Treatment+Location+Range+Taxa,data=gamfits2,FUN=mean,keep.names=T)
# rates.taxa.l <- split(rates.taxa,rates.taxa$Taxa)
# rates.l <- split(gamfits2,gamfits2$Taxa)
# 
# # loop over each taxa, plot a set of four figures, export to a pdf
# pdf(file="C:/Repos/GLAHD/Output/RGR_M_Taxa_filtered_GAM_logM.pdf")
# 
# layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
#        widths=rep(1,4), heights=rep(1,4))
# for (i in 1:length(dat4.l.taxa)){
#   ylims=c(log(min(dat4.l.taxa[[i]]$TotMass)),log(max(dat4.l.taxa[[i]]$TotMass+5)))
#   #- plot home mass over time
#   plot(log(TotMass)~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",
#        main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
#   plotBy(log(predMass)~Time|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
#   lines(log(predMass)~Time,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
#   
#   #- plot warmed mass over time
#   plot(log(TotMass)~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",
#        main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Warmed",sep='_'))
#   plotBy(log(predMass)~Time|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
#   lines(log(predMass)~Time,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
#   
#   #- plot home RGR over time
#   ylims=c(0,max(dat4.l.taxa[[i]]$RGR+0.1,na.rm=T))
#   
#   plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",xlim=c(0,60),
#        main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
#   plotBy(dydt~Time|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
#   lines(dydt~Time,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
#   
#   #- plot warmed RGR over time
#   plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",xlim=c(0,60),
#        main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Home",sep='_'))
#   plotBy(dydt~Time|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
#   lines(dydt~Time,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
#   
# }
# dev.off()
# 
# #-------------------------------------------------------------------------------------------------------
# #- Do analysis suggested by Mark, where interval-based RGR estimates are compared to model estimates
# #-  The basic question: is a 4th order polynomial model sufficient to describe the data,
# #-      particularly the peaked RGR response to warming during the second growth interval?
# 
# 
# #- the interval-based RGR estimate for the segment in question
# rgr.seg2 <- subset(dat4,Date==as.Date("2014-11-17"))
# names(rgr.seg2)[18] <- "RGR.seg"
# 
# #- get the modeled RGR for each day of the growth interval in question for each plant
# #remember that time "1" is 2014-11-07
# table(dat2$Date) # second growth interval is from 2014-11-17 to 2014-11-26, which is "Time" 11 to 20
# rgr.mod.seg2 <- subset(gamfits2,Time>=11 & Time <=20)
# rgr.mod.seg2.plant <- summaryBy(predMass+AGR+dydt~Code+Species+Treatment+Location+Taxa+Range,data=rgr.mod.seg2,keep.names=T)
# 
# #- put them together
# rgr2 <- merge(rgr.seg2,rgr.mod.seg2.plant,by=c("Code","Species","Treatment","Location","Taxa","Range"))
# 
# #- plot both estimates of RGR for each taxa
# 
# #plot RGR seg, then RGR modeled
# ylims=c(0,0.2)
# xlims=c(0,35)
# colors <- c("blue","red")
# windows(10,10);par(mfrow=c(2,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)
# 
# #RGR segmented
# boxplot(RGR.seg~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
# title(ylab=expression(RGR[seg]),cex.lab=2,line=2.5)
# axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=F,las=2,cex.axis=1.5)
# abline(v=16.4)
# 
# #RGR modeled
# boxplot(dydt~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
# title(ylab=expression(RGR[mod]),cex.lab=2,line=2.5)
# axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(rgr2$Taxa),las=2,cex.axis=1.5)
# abline(v=16.4)
# 
# 
# #- by-plot
# windows(20,20)
# plotBy(dydt~RGR.seg|Location,data=rgr2,xlab="RGR-seg",ylab="RGR-mod",xlim=c(0,0.22),ylim=c(0,0.22),cex.lab=2)
# abline(0,1)
# 
# #- compare statistical analyes
# library(nlme)
# rgr2$Sp_RS_EN <- as.factor(with(rgr2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
# rgr2$Prov_Sp_EN <- as.factor(with(rgr2,paste(Taxa,Species)))
# fm.rgr.mod <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
# plot(fm.rgr.mod,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm.rgr.mod,log(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
# plot(fm.rgr.mod,log(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
# qqnorm(fm.rgr.mod, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
# hist(fm.rgr.mod$residuals[,1])
# anova(fm.rgr.mod)  
# 
# fm.rgr.seg <- lme(RGR.seg~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
# plot(fm.rgr.seg,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm.rgr.seg,RGR.seg~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
# plot(fm.rgr.seg,RGR.seg~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
# qqnorm(fm.rgr.seg, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
# hist(fm.rgr.seg$residuals[,1])
# anova(fm.rgr.seg)  
#
#
##Seems to perform very well
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------




#- plot Mass, RGR (dydt) and AGR against Time

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=mean,keep.names=T)
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(30,40);par(mfrow=c(3,1), mar=c(0.5,6,0.5,3), oma=c(5,1,1,1))
plotBy(predMass~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",legend.cex=2, pch=15, xaxt='n', ylab="", type="o")
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, ylim=c(0.04,0.12),ylab= "",xaxt='n', type="o")
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(AGR~Time|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15,ylab= "",xaxt='n', type="o")
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)



#Find and plot maximum of RGR
maxdydt<-gamfits2[ gamfits2$dydt %in% tapply(gamfits2$dydt, gamfits2$Code, max), ]#max RGR per Code

ylims=c(0,0.2)
xlims=c(0,35)
colors <- c("blue","red")
windows(10,5);par(mfrow=c(1,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)
boxplot(dydt~Treatment+Taxa,data=maxdydt,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[max]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(maxdydt$Taxa),las=2,cex.axis=1.5)
abline(v=16.5)

#- Did maxRGR increase with warming?
#Warming increases maxRGR in temperate but not tropical taxa (Treatment:Location).

library(nlme)
library(effects)
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

#Did warming make the peak occur earlier?
#Model assumptions are difficult to meet to test for differences in timing of maximum RGR
#since maxRGR often occur at Time 1

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------



# Plot RGR and AGR over mass
g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=mean,keep.names=T)
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(dydt~predMass|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15, ylim=c(0.04,0.12), ylab="RGR", xlab="Mass")
plotBy(dydt~predMass|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15, ylim=c(0.04,0.12),ylab="RGR", xlab="Mass (g))", log="x")

plotBy(AGR~predMass|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15, ylim=c(0,4))
plotBy(AGR~predMass|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15, ylim=c(0,2), xlim=c(0,30),xlab="Mass")
plotBy(dydt~AGR|combotrt,data=g.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15, ylim=c(0.04,0.12), ylab="RGR", xlab="AGR", log="x")



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Relative enhancement of biomass
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:3,5,8)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="South", ylim=c(0.5,2), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "North", ylim=c(0.5,2),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(55,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=0.8)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~predMass|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.8,2), xlim=c(0,36))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~predMass|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.8,2),yaxt='n', xlim=c(0,80))
mtext(text=expression(M[H]~(g)),side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)


#Absolute enhancement of biomass
bio$ber2<- with(bio, predMassWarm-predMass)
BER<- bio[,c(1:3,5,9)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber2~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-5,20), xlim=c(0,67))
mtext(text=expression(M[W]-M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber2~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-5,20),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber2~predMass|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-5,20), xlim=c(0,36))
mtext(text=expression(M[W]-M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber2~predMass|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-5,20),yaxt='n', xlim=c(0,80))
mtext(text=expression(M[H]~(g)),side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Relative enhancement of RGR
rer<- summaryBy(dydt~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
rerh <- subset(rer, Treatment == "Home")
rerw <- subset (rer, Treatment == "Warmed")
names(rerw)<- c("Time","Range","Location","Treatment","dydtWarm" )
rel <- merge(rerh,rerw, by=c("Time","Range","Location"))
rel$rer<- with(rel, dydtWarm/dydt)
RER<- rel[,c(1:3,8)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(rer~Time|Range,data=subset(RER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,1.5), xlim=c(0,67))
mtext(text=expression(RGR[W]:RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(rer~Time|Range,data=subset(RER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,1.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)

#maximum enhancement of RGR
srer<- subset(RER, Location =="S")
nrer<-subset(RER, Location =="N")
srer[ srer$rer %in% tapply(srer$rer, srer$Range, max), ]
nrer[ nrer$rer %in% tapply(nrer$rer, nrer$Range, max), ]

#Absolute enhancement of RGR
rel$rer2<- with(rel, dydtWarm-dydt)
RER<- rel[,c(1:3,9)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(rer2~Time|Range,data=subset(RER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-0.02,0.04), xlim=c(0,67))
mtext(text=expression(RGR[W]-RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(rer2~Time|Range,data=subset(RER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-0.02,0.04),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)



#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#RGR per taxa
sumrate<-summaryBy(dydt~Time+Taxa+Treatment+Location, data=gamfits2, FUN=c(mean,standard.error))
rate.l <- split(sumrate,sumrate$Taxa)
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,0.15)
xlims <- c(0,60)
palette(c("black","red"))
for (i in 1:length(rate.l)){
  toplot <- rate.l[[i]]
  plotBy(dydt.mean~Time|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Time,y=toplot$dydt.mean,
                                  SE=toplot$dydt.standard.error,direction="updown",col=c("black","red")))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(rate.l)){
    mtext(text="RGR",side=2,outer=T,line=2,cex=2)
    mtext(text="Time",side=1,outer=T,line=2,cex=1)
    
  }      
}

legend(x=80,y=0.10,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#--------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Extract data from specific times and test for significant differences

dat5<- subset(gamfits2,Time==60)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.agr <- lme(sqrt(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 

plot(allEffects(fm1.agr))      

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)                 

plot(allEffects(fm1.rgr))     

fm1.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)                 

fm1.mass <- lme(predMass~Treatment+Location+Range+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass))     




#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#RGR = LAR * ULR   (Leaf area ratio and Unit leaf rate)
#Is the change in RGR due to leaf areaor photosynthesis per unit leaf area?


#Get LAR and ULR
rate <- merge(gamfits2,dat2,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,18:19)]# merge with dat2 to get LA
rate$LAR <- rate$leafArea/rate$TotMass
rate$ULR <- with(rate,dydt/LAR)

#LAR per taxa
sumrate<-summaryBy(LAR~Time+Taxa+Treatment+Location, data=rate, FUN=c(mean,standard.error))
rate.l <- split(sumrate,sumrate$Taxa)
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,180)
xlims <- c(0,60)
palette(c("black","red"))
for (i in 1:length(rate.l)){
  toplot <- rate.l[[i]]
  plotBy(LAR.mean~Time|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Time,y=toplot$LAR.mean,
                                  SE=toplot$LAR.standard.error,direction="updown", col=c("black","red")))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(rate.l)){
    mtext(text="LAR",side=2,outer=T,line=2,cex=2)
    mtext(text="Time",side=1,outer=T,line=2,cex=1)
    
  }      
}
legend(x=80,y=180,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#LAR per taxa over mass
sumrate<-summaryBy(LAR~predMass+Taxa+Treatment+Location, data=rate, FUN=c(mean,standard.error))
rate.l <- split(sumrate,sumrate$Taxa)
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,180)
xlims <- c(0,80)
palette(c("black","red"))
for (i in 1:length(rate.l)){
  toplot <- rate.l[[i]]
  plotBy(LAR.mean~predMass|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$predMass,y=toplot$LAR.mean,
                                  SE=toplot$LAR.standard.error,direction="updown", col=c("black","red")))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(rate.l)){
    mtext(text="LAR",side=2,outer=T,line=2,cex=2)
    mtext(text="Mass",side=1,outer=T,line=2,cex=1)
    
  }      
}
legend(x=100,y=180,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#Absolute enhancement of LAR
ler<- summaryBy(LAR~Time+Range+Location+Treatment, data=rate, FUN=mean, keep.names=T) 
lerh <- subset(ler, Treatment == "Home")
lerw <- subset (ler, Treatment == "Warmed")
names(lerw)<- c("Time","Range","Location","Treatment","LARWarm" )
lea <- merge(lerh,lerw, by=c("Time","Range","Location"))
lea$ler<- with(lea, LARWarm-LAR)
LER<- lea[,c(1:3,8)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ler~Time|Range,data=subset(LER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-30,20), xlim=c(0,67))
mtext(text=expression(LAR[W]-LAR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ler~Time|Range,data=subset(LER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topleft",type="o", main = "North", ylim=c(-30,20),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of LAR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

#Relative enhancement of LAR

lea$ler2<- with(lea, LARWarm/LAR)
LER<- lea[,c(1:3,9)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ler2~Time|Range,data=subset(LER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,1.5), xlim=c(0,67))
mtext(text=expression(LAR[W]-LAR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ler2~Time|Range,data=subset(LER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topleft",type="o", main = "North",ylim=c(0.5,1.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of LAR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)


#--------------------------------------------------------------------------------------------------

#Extract data from specific times and test for significant differences
#  1 11 20 25 32 39 46 53 60
ratelar<- subset(rate,Time==60)

ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme(log(LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR)                 


fm1.LAR <- lme(log(LAR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(allEffects(fm1.LAR))     



#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


#Co-plot LAR
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Time+Treatment+Location+Range,data=rate,FUN=mean,keep.names=T)
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(LAR~Time|combotrt,type='o',data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15, xlim=c(0,60), ylim=c(50,145))

plotBy(LAR~predMass|combotrt,type='o',data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15, ylim=c(70,145), xlab="Mass (g)")


rate$Sp_RS_EN <- as.factor(with(rate,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rate$Prov_Sp_EN <- as.factor(with(rate,paste(Taxa,Species)))


fm.lar.rgr <- lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rate)
plot(fm.lar.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.lar.rgr,log(LAR)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.lar.rgr,log(LAR)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.lar.rgr, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.lar.rgr$residuals[,1])
anova(fm.lar.rgr) 
plot(allEffects(fm.lar.rgr))

fm.lar.rgr <- lme(log(LAR)~dydt+Treatment+dydt:Location+
                    Treatment:Location+dydt:Range+Treatment:Range+Location:Range+
                    dydt:Treatment:Location+dydt:Treatment:Range+
                    dydt:Treatment:Location:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rate)



fm.lar.rgr <- lm(log(LAR)~dydt*combotrt,data=rate.trt)
anova(fm.lar.rgr) 
fm.lar.rgr2 <- lm(log(LAR)~dydt+combotrt,data=rate.trt)
anova(fm.lar.rgr2) 
plot(allEffects(fm.lar.rgr2))


#------------------------------------------------------------------------------

#ULR


#ULR per taxa
sumrate3<-summaryBy(ULR~Time+Taxa+Treatment+Location, data=rate, FUN=c(mean,standard.error))
rate.l <- split(sumrate3,sumrate3$Taxa)
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,0.002)
xlims <- c(0,60)
palette(c("black","red"))
for (i in 1:length(rate.l)){
  toplot <- rate.l[[i]]
  plotBy(ULR.mean~Time|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Time,y=toplot$ULR.mean,
                                  SE=toplot$ULR.standard.error,direction="updown", col=c("black","red")))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(rate.l)){
    mtext(text="ULR",side=2,outer=T,line=2,cex=2)
    mtext(text="Time",side=1,outer=T,line=2,cex=1)
    
  }      
}

legend(x=80,y=180,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#ULR over mass per taxa
sumrate3<-summaryBy(ULR~predMass+Taxa+Treatment+Location, data=rate, FUN=c(mean,standard.error))
rate.l <- split(sumrate3,sumrate3$Taxa)
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,0.002)
xlims <- c(0,80)
palette(c("black","red"))
for (i in 1:length(rate.l)){
  toplot <- rate.l[[i]]
  plotBy(ULR.mean~predMass|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$predMass,y=toplot$ULR.mean,
                                  SE=toplot$ULR.standard.error,direction="updown", col=c("black","red")))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(rate.l)){
    mtext(text="ULR",side=2,outer=T,line=2,cex=2)
    mtext(text="Mass",side=1,outer=T,line=2,cex=1)
    
  }      
}

legend(x=100,y=180,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#Absolute enhancement of ULR
ler<- summaryBy(ULR~Time+Range+Location+Treatment, data=rate, FUN=mean, keep.names=T) 
lerh <- subset(ler, Treatment == "Home")
lerw <- subset (ler, Treatment == "Warmed")
names(lerw)<- c("Time","Range","Location","Treatment","ULRWarm" )
lea <- merge(lerh,lerw, by=c("Time","Range","Location"))
lea$ler<- with(lea, ULRWarm-ULR)
LER<- lea[,c(1:3,8)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ler~Time|Range,data=subset(LER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-0.0003,0.0003), xlim=c(0,67))
mtext(text=expression(ULR[W]-ULR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ler~Time|Range,data=subset(LER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-0.0003,0.0003),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of ULR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

#Relative enhancement of ULR

lea$ler2<- with(lea, ULRWarm/ULR)
LER<- lea[,c(1:3,9)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ler2~Time|Range,data=subset(LER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,1.5), xlim=c(0,67))
mtext(text=expression(ULR[W]-ULR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ler2~Time|Range,data=subset(LER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North",ylim=c(0.5,1.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of ULR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)

#--------------------------------------------------------------------------------------------------

#Extract data from specific times and test for significant differences
#  1 11 20 25 32 39 46 53 60
rateulr<- subset(rate,Time==1)

rateulr$Location <- factor(rateulr$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
rateulr$Sp_RS_EN <- as.factor(with(rateulr,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rateulr$Prov_Sp_EN <- as.factor(with(rateulr,paste(Taxa,Species)))
rateulr$Sp_Loc_EN <- as.factor(with(rateulr,paste(Species,Location)))

fm1.ULR <- lme(ULR~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rateulr)
plot(fm1.ULR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.ULR,ULR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.ULR,ULR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.ULR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.ULR$residuals[,1])
anova(fm1.ULR)                 


fm1.LAR <- lme(log(LAR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rateulr)
plot(allEffects(fm1.ULR))     



#--------------------------------------------------------------------------------------------------


#Co-plot RGR and ULR
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Time+Treatment+Location+Range,data=rate,FUN=mean,keep.names=T)
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(ULR~dydt|combotrt, type='o', data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
         legendwhere="topleft",pch=15, xlim=c(0,0.12), ylim=c(0,0.0015), xlab="RGR")

fm.ulr.rgr <- lm(ULR~dydt*combotrt,data=rate.trt)
anova(fm.lar.rgr) 
fm.lar.rgr2 <- lm(ULR~dydt+combotrt,data=rate.trt)
anova(fm.lar.rgr2) 
plot(allEffects(fm.lar.rgr2))
  
rate$Sp_RS_EN <- as.factor(with(rate,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rate$Prov_Sp_EN <- as.factor(with(rate,paste(Taxa,Species)))
  
  
fm.ulr.rgr <- lme(ULR~dydt*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rate)
plot(fm.ulr.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.ulr.rgr,ULR~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.ulr.rgr,ULR~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.ulr.rgr, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.ulr.rgr$residuals[,1])
anova(fm.ulr.rgr) 
plot(allEffects(fm.ulr.rgr))


windows(40,30);par(mfrow=c(1,1))
plotBy(ULR~predMass|combotrt,type='o',data=rate.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15,  xlab="Mass (g)")


#------------------------------------------------------------------------------

#growth response coefficients

grc<-summaryBy(dydt+LAR+ULR~Treatment+Range+Location+Time, data=rate, FUN=mean, keep.names=T)  
grcwarm<- subset(grc, Treatment=="Warmed")
grchome<- subset(grc, Treatment=="Home")
names(grcwarm)<- c("Treatment","Range","Location","Time","dydtwarm","LARwarm","ULRwarm"  )
grc2<- merge(grchome,grcwarm, by=c("Range","Location","Time"))
grc2$dLAR<- with(grc2, ((LARwarm-LAR)/LAR))  
grc2$dRGR<- with(grc2, ((dydtwarm-dydt)/dydt)) 
grc2$dULR<- with(grc2, ((ULRwarm-ULR)/ULR))
grc2$LARbrc<- with(grc2, dLAR/dRGR)  
grc2$ULRbrc<- with(grc2, dULR/dRGR)
grc2$sum <- with(grc2, LARbrc+ULRbrc)

grc.trt <- summaryBy(LARbrc+ULRbrc~Time+Treatment+Location+Range,data=grc2,FUN=mean,keep.names=T)
grc.trt$combotrt <- as.factor(paste(grc.trt$Location,grc.trt$Range,sep="_"))

windows(40,30);par(mfrow=c(2,1),mar=c(1,3,1,1), oma=c(3,2,1,1))
plotBy(LARbrc~Time|combotrt, type='o', data=grc.trt,col=c("red","blue","orange","grey"),
       legendwhere="topleft",pch=15,  xlab="", ylab="")
abline(h=1)
mtext(text="dLAR/dRGR", side=2, line=3)
plotBy(ULRbrc~Time|combotrt, type='o', data=grc.trt,col=c("red","blue","orange","grey"),
       legend=F,pch=15, xlim=c(0,60), ylim=c(0.5,1.5), xlab="Time", ylab="")
abline(h=1)
mtext(text="dULR/dRGR", side=2, line=2.5)
mtext(text="Time", side=1, line=2.5)
#------------------------------------------------------------------------------

#Make bar graph
grctrt2<- grc.trt[,c(1:5,7)]; grctrt2$Category<-as.factor("LAR");names(grctrt2)<-c("Time","Treatment","Location","Range","value","combotrt","Category" )
grctrt3<- grc.trt[,c(1:4,6:7)]; grctrt3$Category<-as.factor("ULR");names(grctrt3)<-c("Time","Treatment","Location","Range","value","combotrt","Category") 
grctrt4<- rbind(grctrt2,grctrt3)
Nnarrow<-subset(grctrt4,combotrt=="N_narrow_Home"|combotrt=="N_narrow_Warmed")
Nwide<-subset(grctrt4,combotrt=="N_wide_Home"|combotrt=="N_wide_Warmed")
Snarrow<-subset(grctrt4,combotrt=="S_narrow_Home"|combotrt=="S_narrow_Warmed")
Swide<-subset(grctrt4,combotrt=="S_wide_Home"|combotrt=="S_wide_Warmed")

barplot(Nnarrow$value, beside =T, col=c("darkblue","red"))
abline(h=c(4,8,12,16,20,24)


#Big figure on how mean LAR and RGR change with warming per time point
windows(40,30);par(mfrow=c(3,4), mar=c(1,2,1,0.8), oma=c(5,6,2,2.5))
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==1),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),xlab="")
segments(subset(rate.trt, Time==1)[1,8],subset(rate.trt, Time==1)[1,5],subset(rate.trt, Time==1)[5,8],subset(rate.trt, Time==1)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==1)[2,8],subset(rate.trt, Time==1)[2,5],subset(rate.trt, Time==1)[6,8],subset(rate.trt, Time==1)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==1)[3,8],subset(rate.trt, Time==1)[3,5],subset(rate.trt, Time==1)[7,8],subset(rate.trt, Time==1)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==1)[4,8],subset(rate.trt, Time==1)[4,5],subset(rate.trt, Time==1)[8,8],subset(rate.trt, Time==1)[8,5], col="grey")#S-Wide
mtext(text="1", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==11),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="",xlab="")
segments(subset(rate.trt, Time==11)[1,8],subset(rate.trt, Time==11)[1,5],subset(rate.trt, Time==11)[5,8],subset(rate.trt, Time==11)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==11)[2,8],subset(rate.trt, Time==11)[2,5],subset(rate.trt, Time==11)[6,8],subset(rate.trt, Time==11)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==11)[3,8],subset(rate.trt, Time==11)[3,5],subset(rate.trt, Time==11)[7,8],subset(rate.trt, Time==11)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==11)[4,8],subset(rate.trt, Time==11)[4,5],subset(rate.trt, Time==11)[8,8],subset(rate.trt, Time==11)[8,5], col="grey")#S-Wide
mtext(text="11", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==20),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="",xlab="")
segments(subset(rate.trt, Time==20)[1,8],subset(rate.trt, Time==20)[1,5],subset(rate.trt, Time==20)[5,8],subset(rate.trt, Time==20)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==20)[2,8],subset(rate.trt, Time==20)[2,5],subset(rate.trt, Time==20)[6,8],subset(rate.trt, Time==20)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==20)[3,8],subset(rate.trt, Time==20)[3,5],subset(rate.trt, Time==20)[7,8],subset(rate.trt, Time==20)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==20)[4,8],subset(rate.trt, Time==20)[4,5],subset(rate.trt, Time==20)[8,8],subset(rate.trt, Time==20)[8,5], col="grey")#S-Wide
mtext(text="20", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==25),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),xlab="")
segments(subset(rate.trt, Time==25)[1,8],subset(rate.trt, Time==25)[1,5],subset(rate.trt, Time==25)[5,8],subset(rate.trt, Time==25)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==25)[2,8],subset(rate.trt, Time==25)[2,5],subset(rate.trt, Time==25)[6,8],subset(rate.trt, Time==25)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==25)[3,8],subset(rate.trt, Time==25)[3,5],subset(rate.trt, Time==25)[7,8],subset(rate.trt, Time==25)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==25)[4,8],subset(rate.trt, Time==25)[4,5],subset(rate.trt, Time==25)[8,8],subset(rate.trt, Time==25)[8,5], col="grey")#S-Wide
mtext(text="25", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==32),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="",xlab="")
segments(subset(rate.trt, Time==32)[1,8],subset(rate.trt, Time==32)[1,5],subset(rate.trt, Time==32)[5,8],subset(rate.trt, Time==32)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==32)[2,8],subset(rate.trt, Time==32)[2,5],subset(rate.trt, Time==32)[6,8],subset(rate.trt, Time==32)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==32)[3,8],subset(rate.trt, Time==32)[3,5],subset(rate.trt, Time==32)[7,8],subset(rate.trt, Time==32)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==32)[4,8],subset(rate.trt, Time==32)[4,5],subset(rate.trt, Time==32)[8,8],subset(rate.trt, Time==32)[8,5], col="grey")#S-Wide
mtext(text="32", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==39),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="",xlab="")
segments(subset(rate.trt, Time==39)[1,8],subset(rate.trt, Time==39)[1,5],subset(rate.trt, Time==39)[5,8],subset(rate.trt, Time==39)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==39)[2,8],subset(rate.trt, Time==39)[2,5],subset(rate.trt, Time==39)[6,8],subset(rate.trt, Time==39)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==39)[3,8],subset(rate.trt, Time==39)[3,5],subset(rate.trt, Time==39)[7,8],subset(rate.trt, Time==39)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==39)[4,8],subset(rate.trt, Time==39)[4,5],subset(rate.trt, Time==39)[8,8],subset(rate.trt, Time==39)[8,5], col="grey")#S-Wide
mtext(text="39", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==46),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150), xlab="RGR")
segments(subset(rate.trt, Time==46)[1,8],subset(rate.trt, Time==46)[1,5],subset(rate.trt, Time==46)[5,8],subset(rate.trt, Time==46)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==46)[2,8],subset(rate.trt, Time==46)[2,5],subset(rate.trt, Time==46)[6,8],subset(rate.trt, Time==46)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==46)[3,8],subset(rate.trt, Time==46)[3,5],subset(rate.trt, Time==46)[7,8],subset(rate.trt, Time==46)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==46)[4,8],subset(rate.trt, Time==46)[4,5],subset(rate.trt, Time==46)[8,8],subset(rate.trt, Time==46)[8,5], col="grey")#S-Wide
mtext(text="46", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==53),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="", xlab="RGR")
segments(subset(rate.trt, Time==53)[1,8],subset(rate.trt, Time==53)[1,5],subset(rate.trt, Time==53)[5,8],subset(rate.trt, Time==53)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==53)[2,8],subset(rate.trt, Time==53)[2,5],subset(rate.trt, Time==53)[6,8],subset(rate.trt, Time==53)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==53)[3,8],subset(rate.trt, Time==53)[3,5],subset(rate.trt, Time==53)[7,8],subset(rate.trt, Time==53)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==53)[4,8],subset(rate.trt, Time==53)[4,5],subset(rate.trt, Time==53)[8,8],subset(rate.trt, Time==53)[8,5], col="grey")#S-Wide
mtext(text="53", side=3, line=-1.5)
plotBy(LAR~dydt|combotrt, type='o', data=subset(rate.trt, Time==60),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=c(50,150),ylab="", xlab="RGR")
segments(subset(rate.trt, Time==60)[1,8],subset(rate.trt, Time==60)[1,5],subset(rate.trt, Time==60)[5,8],subset(rate.trt, Time==60)[5,5], col="red")#N-Narrow
segments(subset(rate.trt, Time==60)[2,8],subset(rate.trt, Time==60)[2,5],subset(rate.trt, Time==60)[6,8],subset(rate.trt, Time==60)[6,5], col="blue")#N-Wide
segments(subset(rate.trt, Time==60)[3,8],subset(rate.trt, Time==60)[3,5],subset(rate.trt, Time==60)[7,8],subset(rate.trt, Time==60)[7,5], col="orange")#S-Narrow
segments(subset(rate.trt, Time==60)[4,8],subset(rate.trt, Time==60)[4,5],subset(rate.trt, Time==60)[8,8],subset(rate.trt, Time==60)[8,5], col="grey")#S-Wide
mtext(text="60", side=3, line=-1.5)
mtext(text="LAR",side=2,outer=T,line=2,cex=2)
mtext(text="RGR",side=1,outer=T,line=2,cex=1)
legend(x=0.15,y=145,legend=c(levels(rate.trt$combotrt)),pch=15,cex=1.5,xpd=NA,col=c("red","black","blue","green","orange","cyan","grey","yellow"))


#Big figure on how mean ULR and RGR change with warming per time point
windows(40,30);par(mfrow=c(3,4), mar=c(1,2,1,0.8), oma=c(5,6,2,2.5))
ylims<-c(0,0.0015)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==1),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,xlab="")
segments(subset(rate.trt, Time==1)[1,8],subset(rate.trt, Time==1)[1,6],subset(rate.trt, Time==1)[5,8],subset(rate.trt, Time==1)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==1)[2,8],subset(rate.trt, Time==1)[2,6],subset(rate.trt, Time==1)[6,8],subset(rate.trt, Time==1)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==1)[3,8],subset(rate.trt, Time==1)[3,6],subset(rate.trt, Time==1)[7,8],subset(rate.trt, Time==1)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==1)[4,8],subset(rate.trt, Time==1)[4,6],subset(rate.trt, Time==1)[8,8],subset(rate.trt, Time==1)[8,6], col="grey")#S-Wide
mtext(text="1", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==11),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="",xlab="")
segments(subset(rate.trt, Time==11)[1,8],subset(rate.trt, Time==11)[1,6],subset(rate.trt, Time==11)[5,8],subset(rate.trt, Time==11)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==11)[2,8],subset(rate.trt, Time==11)[2,6],subset(rate.trt, Time==11)[6,8],subset(rate.trt, Time==11)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==11)[3,8],subset(rate.trt, Time==11)[3,6],subset(rate.trt, Time==11)[7,8],subset(rate.trt, Time==11)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==11)[4,8],subset(rate.trt, Time==11)[4,6],subset(rate.trt, Time==11)[8,8],subset(rate.trt, Time==11)[8,6], col="grey")#S-Wide
mtext(text="11", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==20),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="",xlab="")
segments(subset(rate.trt, Time==20)[1,8],subset(rate.trt, Time==20)[1,6],subset(rate.trt, Time==20)[5,8],subset(rate.trt, Time==20)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==20)[2,8],subset(rate.trt, Time==20)[2,6],subset(rate.trt, Time==20)[6,8],subset(rate.trt, Time==20)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==20)[3,8],subset(rate.trt, Time==20)[3,6],subset(rate.trt, Time==20)[7,8],subset(rate.trt, Time==20)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==20)[4,8],subset(rate.trt, Time==20)[4,6],subset(rate.trt, Time==20)[8,8],subset(rate.trt, Time==20)[8,6], col="grey")#S-Wide
mtext(text="20", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==25),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,xlab="")
segments(subset(rate.trt, Time==25)[1,8],subset(rate.trt, Time==25)[1,6],subset(rate.trt, Time==25)[5,8],subset(rate.trt, Time==25)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==25)[2,8],subset(rate.trt, Time==25)[2,6],subset(rate.trt, Time==25)[6,8],subset(rate.trt, Time==25)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==25)[3,8],subset(rate.trt, Time==25)[3,6],subset(rate.trt, Time==25)[7,8],subset(rate.trt, Time==25)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==25)[4,8],subset(rate.trt, Time==25)[4,6],subset(rate.trt, Time==25)[8,8],subset(rate.trt, Time==25)[8,6], col="grey")#S-Wide
mtext(text="25", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==32),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="",xlab="")
segments(subset(rate.trt, Time==32)[1,8],subset(rate.trt, Time==32)[1,6],subset(rate.trt, Time==32)[5,8],subset(rate.trt, Time==32)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==32)[2,8],subset(rate.trt, Time==32)[2,6],subset(rate.trt, Time==32)[6,8],subset(rate.trt, Time==32)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==32)[3,8],subset(rate.trt, Time==32)[3,6],subset(rate.trt, Time==32)[7,8],subset(rate.trt, Time==32)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==32)[4,8],subset(rate.trt, Time==32)[4,6],subset(rate.trt, Time==32)[8,8],subset(rate.trt, Time==32)[8,6], col="grey")#S-Wide
mtext(text="32", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==39),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="",xlab="")
segments(subset(rate.trt, Time==39)[1,8],subset(rate.trt, Time==39)[1,6],subset(rate.trt, Time==39)[5,8],subset(rate.trt, Time==39)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==39)[2,8],subset(rate.trt, Time==39)[2,6],subset(rate.trt, Time==39)[6,8],subset(rate.trt, Time==39)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==39)[3,8],subset(rate.trt, Time==39)[3,6],subset(rate.trt, Time==39)[7,8],subset(rate.trt, Time==39)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==39)[4,8],subset(rate.trt, Time==39)[4,6],subset(rate.trt, Time==39)[8,8],subset(rate.trt, Time==39)[8,6], col="grey")#S-Wide
mtext(text="39", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==46),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims, xlab="RGR")
segments(subset(rate.trt, Time==46)[1,8],subset(rate.trt, Time==46)[1,6],subset(rate.trt, Time==46)[5,8],subset(rate.trt, Time==46)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==46)[2,8],subset(rate.trt, Time==46)[2,6],subset(rate.trt, Time==46)[6,8],subset(rate.trt, Time==46)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==46)[3,8],subset(rate.trt, Time==46)[3,6],subset(rate.trt, Time==46)[7,8],subset(rate.trt, Time==46)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==46)[4,8],subset(rate.trt, Time==46)[4,6],subset(rate.trt, Time==46)[8,8],subset(rate.trt, Time==46)[8,6], col="grey")#S-Wide
mtext(text="46", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==53),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="", xlab="RGR")
segments(subset(rate.trt, Time==53)[1,8],subset(rate.trt, Time==53)[1,6],subset(rate.trt, Time==53)[5,8],subset(rate.trt, Time==53)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==53)[2,8],subset(rate.trt, Time==53)[2,6],subset(rate.trt, Time==53)[6,8],subset(rate.trt, Time==53)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==53)[3,8],subset(rate.trt, Time==53)[3,6],subset(rate.trt, Time==53)[7,8],subset(rate.trt, Time==53)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==53)[4,8],subset(rate.trt, Time==53)[4,6],subset(rate.trt, Time==53)[8,8],subset(rate.trt, Time==53)[8,6], col="grey")#S-Wide
mtext(text="53", side=3, line=-1.5)
plotBy(ULR~dydt|combotrt, type='o', data=subset(rate.trt, Time==60),col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15, xlim=c(0.04,0.12), ylim=ylims,ylab="", xlab="RGR")
segments(subset(rate.trt, Time==60)[1,8],subset(rate.trt, Time==60)[1,6],subset(rate.trt, Time==60)[5,8],subset(rate.trt, Time==60)[5,6], col="red")#N-Narrow
segments(subset(rate.trt, Time==60)[2,8],subset(rate.trt, Time==60)[2,6],subset(rate.trt, Time==60)[6,8],subset(rate.trt, Time==60)[6,6], col="blue")#N-Wide
segments(subset(rate.trt, Time==60)[3,8],subset(rate.trt, Time==60)[3,6],subset(rate.trt, Time==60)[7,8],subset(rate.trt, Time==60)[7,6], col="orange")#S-Narrow
segments(subset(rate.trt, Time==60)[4,8],subset(rate.trt, Time==60)[4,6],subset(rate.trt, Time==60)[8,8],subset(rate.trt, Time==60)[8,6], col="grey")#S-Wide
mtext(text="60", side=3, line=-1.5)
mtext(text="ULR",side=2,outer=T,line=2,cex=2)
mtext(text="RGR",side=1,outer=T,line=2,cex=1)
legend(x=0.15,y=0.0015,legend=c(levels(rate.trt$combotrt)),pch=15,cex=1.5,xpd=NA,col=c("red","black","blue","green","orange","cyan","grey","yellow"))


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#growth response coefficients V2 (seperately per treatment)


#- Calculate incremental increase in RGR, LAR and ULR for each plant
rate2 <- rate[with(rate,order(Code,Time)),] #- make sure the observations for each plant are in order

#Function to divide T2 by T1
for (i in 1:nrow(rate2)) {
rate2[i, "dLAR"] <- (rate2[i+1,13]-rate2[i,13])/rate2[i,13]
rate2[i, "dULR"] <- (rate2[i+1,14]-rate2[i,14])/rate2[i,14]
}
for (i in 1:nrow(rate2)) {
  rate2[i, "dRGR"] <- (rate2[i+1,8]-rate2[i,8])/rate2[i,8]
}
rate2$LARbrc <- with(rate2, dLAR/dRGR)
rate2$ULRbrc <- with(rate2, dULR/dRGR)

#remove last timepoint as it is not relevant
rate3<-rate2[!(rate2$Time == 60),]

#summarise per combotreatment
grc.trt <- summaryBy(LARbrc+ULRbrc~Time+Treatment+Location+Range,data=rate3,FUN=mean,keep.names=T)
grc.trt$combotrt <- as.factor(paste(grc.trt$Location,grc.trt$Range,grc.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(2,1),mar=c(1,3,1,1), oma=c(3,2,1,1))
plotBy(LARbrc~Time|combotrt, type='p', data=grc.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15,  xlab="", ylab="")
abline(h=1)
mtext(text="dLAR/dRGR", side=2, line=3)
plotBy(ULRbrc~Time|combotrt, type='p', data=grc.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15,  xlab="Time", ylab="")
abline(h=1)
mtext(text="dULR/dRGR", side=2, line=2.5)
mtext(text="Time", side=1, line=2.5)

#- Calculate incremental increase in RGR, LAR and ULR for each taxa
rate4<- summaryBy(dydt+LAR+ULR~Time+Taxa, data=rate, FUN=mean, keep.names=T)
rate4 <- rate[with(rate,order(Taxa,Time)),] #- make sure the observations for each plant are in order


#Function to divide T2 by T1
for (i in 1:nrow(rate4)) {
  rate4[i, "dLAR"] <- (rate4[i+1,13]-rate4[i,13])/rate4[i,13]
  rate4[i, "dULR"] <- (rate4[i+1,14]-rate4[i,14])/rate4[i,14]
}
for (i in 1:nrow(rate4)) {
  rate4[i, "dRGR"] <- (rate4[i+1,8]-rate4[i,8])/rate4[i,8]
}
rate4$LARbrc <- with(rate4, dLAR/dRGR)
rate4$ULRbrc <- with(rate4, dULR/dRGR)

#remove last timepoint as it is not relevant
rate5<-rate4[!(rate4$Time == 60),]

#summarise per combotreatment
grc.trt <- summaryBy(LARbrc+ULRbrc~Time+Treatment+Location+Range,data=rate5,FUN=mean,keep.names=T)
grc.trt$combotrt <- as.factor(paste(grc.trt$Location,grc.trt$Range,grc.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(2,1),mar=c(1,3,1,1), oma=c(3,2,1,1))
plotBy(LARbrc~Time|combotrt, type='p', data=grc.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15,  xlab="", ylab="")
abline(h=1)
mtext(text="dLAR/dRGR", side=2, line=3)
plotBy(ULRbrc~Time|combotrt, type='p', data=grc.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legend=F,pch=15,  xlab="Time", ylab="")
abline(h=1)
mtext(text="dULR/dRGR", side=2, line=2.5)
mtext(text="Time", side=1, line=2.5)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#Following Renton and Poorter 2011
# log(LAR) and log(ULR) for each combination (8) at each time point (9) - slopes are covariance factor
windows(40,30);par(mfrow=c(4,2),mar=c(1,3,1,1), oma=c(3,2,1,1))
ylims<-c(3,6)
xlims<-c(-9,-6)
plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"),
       col=rainbow(9),legendwhere="topleft",pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"&Time==60)),col="#FF00AAFF")

plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Warmed"&Time==60)),col="#FF00AAFF")

plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Home"&Time==60)),col="#FF00AAFF")

plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="N"& Range=="wide"&Treatment=="Warmed"&Time==60)),col="#FF00AAFF")

plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Home"&Time==60)),col="#FF00AAFF")
plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="narrow"&Treatment=="Warmed"&Time==60)),col="#FF00AAFF")
plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Home"&Time==60)),col="#FF00AAFF")
plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"),
       col=rainbow(9),pch=15,  xlab="", ylab="", xlim=xlims, ylim=ylims, legend=F)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Location=="S"& Range=="wide"&Treatment=="Warmed"&Time==60)),col="#FF00AAFF")


# log(LAR) and log(dydt) for each combination (8) at each time point (9) - slope+half of covariance factor is LAR contribution

# log(ULR) and log(dydt) for each combination (8) at each time point (9) - slope+half of covariance factor is ULR contribution


# Too complicated start simple
fm.lar.rgr <- lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))
summary(fm.lar.rgr)

fm.ulr.rgr <- lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))
summary(fm.ulr.rgr)

#Simple linear models

Time <-rep(unique(rate$Time), each=2)
GRC<- as.data.frame(Time)
GRC$Type <- as.factor(c("LAR","ULR"))
GRC[1,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==1)))$coefficients[2,1]
GRC[2,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==1)))$coefficients[2,1]
GRC[3,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==11)))$coefficients[2,1]
GRC[4,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==11)))$coefficients[2,1]
GRC[5,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==20)))$coefficients[2,1]
GRC[6,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==20)))$coefficients[2,1]
GRC[7,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==25)))$coefficients[2,1]
GRC[8,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==25)))$coefficients[2,1]
GRC[9,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==32)))$coefficients[2,1]
GRC[10,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==32)))$coefficients[2,1]
GRC[11,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==39)))$coefficients[2,1]
GRC[12,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==39)))$coefficients[2,1]
GRC[13,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==46)))$coefficients[2,1]
GRC[14,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==46)))$coefficients[2,1]
GRC[15,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==53)))$coefficients[2,1]
GRC[16,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==53)))$coefficients[2,1]
GRC[17,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(rate, Time==60)))$coefficients[2,1]
GRC[18,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(rate, Time==60)))$coefficients[2,1]

library(reshape2)
library(zoo) 
df1 <- transform(GRC[order(GRC$Time),], 
                 Time=factor(Time, levels=unique(Time)))
m1 <- acast(df1, Type~Time, value.var='V', fill=0)
windows(10,10); par(mar=c(4,4,4,4))
barplot(m1,beside=T, ylim=c(0,1), col=c("green","orange"), ylab="GRC",xlab="Time")
legend(1, 1, c("LAR","ULR"), cex=0.8, fill=c("green","orange"))
box()

#More complex models

Time <-rep(unique(rate$Time), each=2)
GRC<- as.data.frame(Time)
GRC$Type <- as.factor(c("LAR","ULR"))
GRC[1,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==1)))$coefficients$fixed[[2]]
GRC[2,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==1)))$coefficients$fixed[[2]]
GRC[3,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11)))$coefficients$fixed[[2]]
GRC[4,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11)))$coefficients$fixed[[2]]
GRC[5,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==20)))$coefficients$fixed[[2]]
GRC[6,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==20)))$coefficients$fixed[[2]]
GRC[7,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==25)))$coefficients$fixed[[2]]
GRC[8,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==25)))$coefficients$fixed[[2]]
GRC[9,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==32)))$coefficients$fixed[[2]]
GRC[10,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==32)))$coefficients$fixed[[2]]
GRC[11,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==39)))$coefficients$fixed[[2]]
GRC[12,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==39)))$coefficients$fixed[[2]]
GRC[13,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==46)))$coefficients$fixed[[2]]
GRC[14,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==46)))$coefficients$fixed[[2]]
GRC[15,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==53)))$coefficients$fixed[[2]]
GRC[16,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==53)))$coefficients$fixed[[2]]
GRC[17,"V"] <- summary(lme(log(LAR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==60)))$coefficients$fixed[[2]]
GRC[18,"V"] <- summary(lme(log(ULR)~log(dydt)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==60)))$coefficients$fixed[[2]]

library(reshape2)
library(zoo) 
df1 <- transform(GRC[order(GRC$Time),], 
                 Time=factor(Time, levels=unique(Time)))
m1 <- acast(df1, Type~Time, value.var='V', fill=0)
windows(10,10); par(mar=c(4,4,4,4))
barplot(m1,beside=T, ylim=c(-.5,1.5), col=c("green","orange"), ylab="GRC",xlab="Time")
legend(1, 1.5, c("LAR","ULR"), cex=0.8, fill=c("green","orange"))
box()

windows(40,30);par(mfrow=c(1,1))
ylims<-c(3,6)
xlims<-c(-9,-6)
plotBy(log(LAR)~log(ULR)|Time, data=subset(rate, Location=="N"& Range=="narrow"&Treatment=="Home"),
       col=rainbow(9),legendwhere="topright",pch=15,  xlab="log(ULR)", ylab="log(LAR)", xlim=xlims, ylim=ylims)
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==1)), col="#FF0000FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==11)),col="#FFAA00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==20)), col="#AAFF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==25)),col="#00FF00FF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==32)),col="#00FFAAFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==39)),col="#00AAFFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==46)),col="#0000FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==53)),col="#AA00FFFF")
abline(lm(log(LAR)~log(ULR),data=subset(rate, Time==60)),col="#FF00AAFF")


#
#Simple models for home and warmed seperately
#Simple linear models

Time <-rep(unique(rate$Time), each=2)
GRC<- as.data.frame(Time)
GRC$Type <- as.factor(c("LAR","ULR"))

ratewarm<- subset(rate, Location=="N"& Range == "narrow" &Treatment == "Home" )
GRC[1,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==1)))$coefficients[2,1]
GRC[2,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==1)))$coefficients[2,1]
GRC[3,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==11)))$coefficients[2,1]
GRC[4,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==11)))$coefficients[2,1]
GRC[5,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==20)))$coefficients[2,1]
GRC[6,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==20)))$coefficients[2,1]
GRC[7,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==25)))$coefficients[2,1]
GRC[8,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==25)))$coefficients[2,1]
GRC[9,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==32)))$coefficients[2,1]
GRC[10,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==32)))$coefficients[2,1]
GRC[11,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==39)))$coefficients[2,1]
GRC[12,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==39)))$coefficients[2,1]
GRC[13,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==46)))$coefficients[2,1]
GRC[14,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==46)))$coefficients[2,1]
GRC[15,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==53)))$coefficients[2,1]
GRC[16,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==53)))$coefficients[2,1]
GRC[17,"V"] <- summary(lm(log(LAR)~log(dydt),data=subset(ratewarm, Time==60)))$coefficients[2,1]
GRC[18,"V"] <- summary(lm(log(ULR)~log(dydt),data=subset(ratewarm, Time==60)))$coefficients[2,1]

library(reshape2)
library(zoo) 
df1 <- transform(GRC[order(GRC$Time),], 
                 Time=factor(Time, levels=unique(Time)))
m1 <- acast(df1, Type~Time, value.var='V', fill=0)
windows(10,10); par(mar=c(4,4,4,4))
barplot(m1,beside=T, ylim=c(-0.5,1.5), col=c("green","orange"), ylab="GRC",xlab="Time", main="S_Wide_Home")
legend(24, 1.5, c("LAR","ULR"), cex=0.8, fill=c("green","orange"))
box()

