#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plot the climate data for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")


#------------------------------------------------------------------------------------------------------------------------------s
#- read in the "fast" datasets (airT, RH, and PAR)
files <- list.files(path="W:/WorkingData/GHS39/GLAHD/Share/Data/Climate/heatwave climate data/",pattern="fast",full.names=T)

dat <- list()
for (i in 1:length(files)){
  dat[[i]] <- readTOA5(files[i])
}
fast.dat.raw <- do.call(rbind,dat)
fast.dat.raw$Room <- as.factor(substr(fast.dat.raw$Source,start=16,stop=16))


#- subset to data just for our experiment. (there is a lot of data preceding our work, from the last experment in S39)
startdate <- as.Date("2014-11-8") # this is the date of planting
fastdat <- subset(fast.dat.raw,Date >=startdate)

#- exclude room 1?
fastdat.GLAHD <- subset(fastdat,Room!=1)
fastdat.GLAHD$Room <- factor(fastdat.GLAHD$Room)
#------------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------
#- what proportion of daylight hours is PAR saturating (>1500)? It's a small fraction
length(which(fastdat.GLAHD$PAR_Avg>1500))/length(which(fastdat.GLAHD$PAR_Avg>4))
hist(fastdat.GLAHD$PAR_Avg,xlab="PAR",main="PAR histogram")
#------------------------------------------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------------------------------------------
#- read in the "slow" datasets, which include soil temperature
files.slow <- list.files(path="W:/WorkingData/GHS39/GLAHD/Share/Data/Climate/heatwave climate data/",pattern="Tsoil",full.names=T)
dat.slow <- list()
for (i in 1:length(files.slow)){
  dat.slow[[i]] <- readTOA5(files.slow[i])
}
slow.dat.raw <- do.call(rbind,dat.slow)
slow.dat.raw$Room <- as.factor(substr(slow.dat.raw$Source,start=16,stop=16))
slow.dat.raw$Tsoil_Avg.5. <- slow.dat.raw$Tsoil_Avg.6.<- slow.dat.raw$Tsoil_Avg.7. <- slow.dat.raw$Tsoil_Avg.8. <- NULL
slow.dat.raw$Date <- as.Date(slow.dat.raw$DateTime)

startdate <- as.Date("2014-11-5") # this is the date of planting
slowdat <- subset(slow.dat.raw,Date >=startdate)
#------------------------------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------
#- get daily means for air temperature and RH, daily sums for PAR
fast.day <- dplyr::summarize(group_by(fastdat.GLAHD,Date,Room),
                             Tair = mean(Tair_Avg,na.rm=T),
                             RH = mean(RH_Avg,na.rm=T),
                             PAR = sum(PAR_Avg,na.rm=T))

means <- summaryBy(Tair+RH+PAR~Room,data=subset(as.data.frame(fast.day),Date>=as.Date("2014-11-12") & Date <=as.Date("2014-11-22")))

#- get daily means for soil T
slow.day <- dplyr::summarize(group_by(slowdat,Date,Room),
                             Tsoil = mean(c(Tsoil_Avg.1.,Tsoil_Avg.2.,Tsoil_Avg.3.,Tsoil_Avg.4.),na.rm=T))

means <- summaryBy(Tair+RH+PAR~Room,data=subset(as.data.frame(fast.day),Date>=as.Date("2014-11-12") & Date <=as.Date("2014-11-22")))
#------------------------------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------------------------------
#- plot daily data
#colors <- colorRampPalette(colors=c("blue","yellow","orange","red"))(4)
colors <- c("black","grey","red","forestgreen","purple")
windows(40,30)
par(mfrow=c(3,1), mar=c(0.3,2,0.3,0.8), oma=c(5,6,6,2.5))
size=1.5
# plot airT
plotBy(Tair~Date|Room,data=fast.day,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fast.day$Date),max(fast.day$Date),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean air T (deg C)",side=2,outer=F,line=5,cex=1.2)
legend(x=as.Date("2014-11-20"),y=42,xpd=NA,legend=1:5,col=colors,pch=15,ncol=5,cex=1.5,title="Room")
# plot RH
plotBy(RH~Date|Room,data=fast.day,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fast.day$Date),max(fast.day$Date),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean RH (%)",side=2,outer=F,line=5,cex=1.2)

# plot PAR
plotBy(PAR~Date|Room,data=fast.day,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fast.day$Date),max(fast.day$Date),by="day"),labels=T,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Sum PAR (umol m-2)",side=2,outer=F,line=5,cex=1.2)
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Climate_daily_S39.pdf")

#------------------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------------------
#- get hourly means for air temperature, RH, and PAR
fastdat.GLAHD$DateTime_hr <- round.POSIXct(fastdat.GLAHD$DateTime,"hours")
fast.hour <- dplyr::summarize(group_by(fastdat.GLAHD,DateTime_hr,Room),
                             Tair = mean(Tair_Avg,na.rm=T),
                             RH = mean(RH_Avg,na.rm=T),
                             PAR = mean(PAR_Avg,na.rm=T))

#- plot hourly data
colors <- c("black","grey","red","forestgreen","purple")
windows(40,30)
par(mfrow=c(3,1), mar=c(2,2,0.3,0.8), oma=c(5,6,6,2.5))
size=1.5
# plot airT
plotBy(Tair~DateTime_hr|Room,data=fast.hour,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast.hour$DateTime_hr),max(fast.hour$DateTime_hr),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean air T (deg C)",side=2,outer=F,line=5,cex=1.2)
legend(x=as.POSIXct("2014-11-12 00:00:00"),y=62,xpd=NA,legend=1:5,col=colors,lty=1,lwd=2,ncol=5,cex=1.5,title="Room")

# plot RH
plotBy(RH~DateTime_hr|Room,data=fast.hour,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast.hour$DateTime_hr),max(fast.hour$DateTime_hr),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean RH (%)",side=2,outer=F,line=5,cex=1.2)

# plot PAR
plotBy(PAR~DateTime_hr|Room,data=fast.hour,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast.hour$DateTime_hr),max(fast.hour$DateTime_hr),by="day"),format="%d/%b",labels=T,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="PAR (umol m-2 s-1)",side=2,outer=F,line=5,cex=1.2)
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Climate_hourly_S39.pdf")
#------------------------------------------------------------------------------------------------------------------------------

# extract just the bay 1 data for Chris and Danielle
write.csv(subset(fast.hour,Room==1),file="W:/WorkingData/GHS39/GLAHD/Share/Output/Climate_hourly_room1.csv",row.names=F)







#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
# What were the climate condititons for each treatment throughout the experiment? This analysis must account for the rotation
# of treatments across glasshouse bays
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

#- read in the "key" for when each bay was in each treatment
roomkey <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/Climate/climate_bay_key.csv")
roomkey$Date <- as.Date(roomkey$Date,format="%d/%m/%Y")
roomkey$Treat <- factor(roomkey$Treat,levels=c("SouthHome","SouthWarmed","NorthHome","NorthWarmed"))
roomkey$Room <- as.factor(roomkey$Room)
#add extra dates to graph
Date<- c("2014-12-01","2014-12-01","2014-12-01","2014-12-01","2014-12-23","2014-12-23","2014-12-23","2014-12-23")
Room <- c("2","3","4","5","2","3","4","5")
Treat<- c("NorthHome","NorthWarmed","SouthHome","SouthWarmed","SouthWarmed","NorthHome","NorthWarmed","SouthHome")
extra<-data.frame(Date,Room,Treat)
extra$Date <- as.Date(extra$Date, format="%Y-%m-%d")
key<-rbind(roomkey,extra)


#- put a date object in the hourly air temperature and RH dataset
fast.hour$Date <- as.Date(fast.hour$DateTime_hr)
fast.hour$VPD <- getVPD(Ta=fast.hour$Tair,RH=fast.hour$RH)
fast.hour$hour <- hour(fast.hour$DateTime_hr)

#- merge
fast <- merge(fast.hour,key,by=c("Date","Room"))

boxplot(Tair~Treat,data=fast)
boxplot(VPD~Treat,data=fast)

# midday means (noon to 4pm)
md.means <- summaryBy(Tair+RH+VPD~Treat,data=subset(fast,hour>=12 & hour <=16),FUN=c(mean,sd))
# warming treatment in S
md.means$Tair.mean[2]-md.means$Tair.mean[1]
# warming treatment in N
md.means$Tair.mean[4]-md.means$Tair.mean[3]

# nightime means (midnight to 5am)
nt.means <- summaryBy(Tair+RH+VPD~Treat,data=subset(fast, hour <=5),FUN=c(mean,sd))
# warming treatment in S
nt.means$Tair.mean[2]-nt.means$Tair.mean[1]
# warming treatment in N
nt.means$Tair.mean[4]-nt.means$Tair.mean[3]

#Daily means
d.means <- summaryBy(Tair+RH+VPD~Treat,data=exp,FUN=c(mean,standard.error))
#------------------------------------------------------------------------------------------------------------------------------
d.means$Tair.mean-d.means$Tair.standard.error
d.means$Tair.mean+d.means$Tair.standard.error
plot(Tair.mean~Date, data=d.means)
lines((Tair.mean-Tair.sd)~Date,data=subset(d.means, d.means$Treat == "NorthWarmed" ))


#----------
colors <- c("black","red","forestgreen","purple")
windows(60,30)
par(mfrow=c(3,1), mar=c(2,2,0.3,0.8), oma=c(5,6,6,2.5))
size=1.5
# plot airT
plotBy(Tair~DateTime_hr|Treat,data=fast,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast$DateTime_hr),max(fast$DateTime_hr),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean air T (deg C)",side=2,outer=F,line=5,cex=1.2)
legend(x=as.POSIXct("2014-11-12 00:00:00"),y=52,xpd=NA,legend=levels(fast$Treat),col=colors,lty=1,lwd=2,ncol=4,cex=1.5,title="Treatment")

# plot RH
plotBy(RH~DateTime_hr|Treat,data=fast,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast$DateTime_hr),max(fast$DateTime_hr),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean RH (%)",side=2,outer=F,line=5,cex=1.2)

# plot PAR
plotBy(PAR~DateTime_hr|Treat,data=fast,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(fast$DateTime_hr),max(fast$DateTime_hr),by="day"),format="%d/%b",labels=T,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="PAR (umol m-2 s-1)",side=2,outer=F,line=5,cex=1.2)


#- put a date object in the hourly air temperature and RH dataset
fast.day$Date <- as.Date(fast.day$DateTime_hr)
fast.day$VPD <- getVPD(Ta=fast.day$Tair,RH=fast.day$RH)
fast.day$hour <- hour(fast.day$DateTime_hr)

#- merge
fastday <- merge(fast.day,key,by=c("Date","Room"))

#- plot daily data
#colors <- colorRampPalette(colors=c("blue","yellow","orange","red"))(4)
colors <- c("black","grey","red","forestgreen","purple")
windows(40,30)
par(mfrow=c(3,1), mar=c(0.3,2,0.3,0.8), oma=c(5,6,6,2.5))
size=1.5
# plot airT
plotBy(Tair~Date|Treat,data=fastday,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fastday$Date),max(fastday$Date),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean air T (deg C)",side=2,outer=F,line=5,cex=1.2)
legend(x=as.Date("2014-11-20"),y=38,xpd=NA,legend=levels(fastday$Treat),col=colors,pch=15,ncol=4,cex=1.5,title="Treatment")
# plot RH
plotBy(RH~Date|Treat,data=fastday,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fastday$Date),max(fastday$Date),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Mean RH (%)",side=2,outer=F,line=5,cex=1.2)

# plot PAR
plotBy(PAR~Date|Treat,data=fastday,pch=15,type="b",legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,0),box=T,las=1,cex.axis=1.5)
axis.Date(side=1,at=seq(from=min(fastday$Date),max(fastday$Date),by="day"),labels=T,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="Sum PAR (umol m-2)",side=2,outer=F,line=5,cex=1.2)

