#Supplemental climate figure

metdata<- read.csv("Data/GHS39_GLAHD_MAIN_MET-AIR_20141108-20150118_L1.csv")
metdata$DateTime<- as.POSIXct(metdata$DateTime)
metdata$Date<- as.Date(metdata$Date)

#- get hourly means for air temperature, VPD, and PAR
metdata$DateTime_hr <- round.POSIXct(metdata$DateTime,"hours")
met.hour <- dplyr::summarize(group_by(metdata,DateTime_hr,T_regime),
                              Tair = mean(Tair,na.rm=T),
                              PAR = mean(PAR,na.rm=T),
                             VPD = mean(VPD,na.rm=T))

startdate<- as.Date("2014-11-15")
enddate<- as.Date("2014-11-22")
met.hour_sub<- droplevels(subset(met.hour, DateTime_hr > startdate & DateTime_hr < enddate ))

# #- plot hourly data
colors <- c("black","grey","forestgreen","red")
windows(20,15)
par(mfrow=c(3,1), mar=c(2,2,0.3,0.8), oma=c(8,8,2,5))
size=1.5
# plot airT
plotBy(Tair~DateTime_hr|T_regime,data=met.hour_sub,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size,
       ylim=c(15,43))
magaxis(side=c(2,4),labels=c(1,1),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(met.hour_sub$DateTime_hr),max(met.hour_sub$DateTime_hr+1),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text=expression(paste("Tair (",degree,"C)")),side=2,outer=F,line=5,cex=1.2)
box()
legend("top",legend=c("Tropical-Warmed","Tropical-Home","Temperate-Warmed","Temperate_Home"),col=c("grey","black","red","forestgreen"),
       lty=1,lwd=2,ncol=4,bty='n',cex=1.5)
legend("topright", "a",bty='n',cex=1.8)
# plot VPD
plotBy(VPD~DateTime_hr|T_regime,data=met.hour_sub,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,1),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(met.hour_sub$DateTime_hr),max(met.hour_sub$DateTime_hr),by="day"),labels=F,las=2,cex.axis=1.5,tcl=0.5)
mtext(text="VPD (kPa)",side=2,outer=F,line=5,cex=1.2)
box()
legend("topright","b",bty='n',cex=1.8)
# plot PAR
plotBy(PAR~DateTime_hr|T_regime,data=met.hour_sub,pch=15,type="l",lwd=2,legend=F,axes=F,col=colors,cex=size)
magaxis(side=c(2,4),labels=c(1,1),box=T,las=1,cex.axis=1.5)
axis.POSIXct(side=1,at=seq(from=min(met.hour_sub$DateTime_hr),max(met.hour_sub$DateTime_hr),by="day"),format="%Y-%m-%d",labels=T,las=2,cex.axis=1.5,tcl=0.5)
mtext(text=expression(paste("PPFD (",mu,"mol",m^-2~s^-1,")")),side=2,outer=F,line=5,cex=1.2)
box()
legend("topright","c",bty='n',cex=1.8)

