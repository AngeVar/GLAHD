

#Figure 2: Mass over Time
g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

windows(9.352,6.616)
#par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(4,6,2,1),cex.axis=1)
par(fig=c(0,0.5,0.4,1), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))

plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", yaxt='n', type="l",ylim=c(0.1,70),xlim=c(0,64),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.S$Time,y=g.trt.S$predMass.mean,
                                SE=g.trt.S$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
#mtext(text="Time (Days)", side=1, line=3, cex=1.2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
mtext("Temperate",3,line=-2,at=14,cex=1.5, outer=F)
legend(2,60, legend=c("Narrow Home",expression(paste("Narrow +3.5",degree,"C")),"Wide Home",expression(paste("Wide +3.5",degree,"C"))),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext(text="Total biomass (g)", outer=F, side=2, line=3, cex=1.2)

par(fig=c(0.5,1,0.4,1), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))
plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', yaxt='n',ylab="", type="l",ylim=c(0.1,70),xlim=c(0,64),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.N$Time,y=g.trt.N$predMass.mean,
                                SE=g.trt.N$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
#mtext(text="Time (Days)", side=1, line=3, cex=1.2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext("Tropical",3,line=-2,at=12.5, cex=1.5, outer=F)

############################################
#Enhancement ratios
#get estimate of error using bootstrapping?
#library(boot)


par(fig=c(0,0.5,0,0.4), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, (predMassWarm-predMass)/predMass)
BER<- bio[,c(1:3,5,8)];SBERn<- subset(BER, Location =="S"&Range=="narrow"); SBERn<-SBERn[order(SBERn$Time),]
SBERw<- subset(BER, Location =="S"&Range=="wide"); SBERw<-SBERw[order(SBERw$Time),]
NBERn<- subset(BER,Location =="N"&Range=="narrow"); NBERn<-NBERn[order(NBERn$Time),]
NBERw<- subset(BER,Location =="N"&Range=="wide"); NBERw<-NBERw[order(NBERw$Time),]
plotBy(ber~Time,data=SBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main="", ylim=c(-0.25,1.1), xaxt='n',xlim=c(0,64),yaxp  = c(-0.5, 1, 3),yaxs="i",xaxs="i",)
lines(ber~Time,data=SBERw,col="black",lty=1,lwd=2)
#mtext(side=3,outer=T,text="Temperate", adj=0.2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.015)
abline(h=0)
mtext(expression(Delta~"Biomass (%)"), side=2, line=3, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2)

axis(side = 2, at = seq(from=-1,to=2,by=0.5), labels = T, tck = 0.015)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.015)
axis(side = 4, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.015)

par(fig=c(0.5,1,0,0.4), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))
plotBy(ber~Time,data=NBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main = "", ylim=c(-0.25,1.1),yaxt='n',xaxt='n', xlim=c(0,64),yaxp  = c(-0.5, 1, 3),yaxs="i",xaxs="i",)
lines(ber~Time,data=NBERw,col="black",lty=1,lwd=2)
#mtext(side=3,outer=T,text="Tropical", adj=0.8)
abline(h=0)
#legend(35,2, legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1,lwd=2, bty="n")
#text(2,1.95,labels="B",cex=2,adj=0.5)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.015)
mtext(text="Time (Days)", side=1, line=3, cex=1.2)
axis(side = 2, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.015)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.015)
axis(side = 4, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.015)

################################################################################
#Final biomass per taxa
g.trt.60<-subset(g.trt, Time=="60")
g.trt.60$combo <- as.factor(paste(g.trt.60$Location,g.trt.60$Range,sep="_"))
windows(5.513333,4.135);par(mar=c(6,6,6,4))
target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
g.trt.60<-g.trt.60[match(target, g.trt.60$combotrt),]

bar<-barplot(g.trt.60$predMass.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0,110), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.8)),4)), 
             axis.lty=0)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar(bar,g.trt.60$predMass.mean,g.trt.60$predMass.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.25,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.25,at=8.5,cex=1.15, outer=F)
legend(-0.1,95, legend=c("H",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),"red"),bty="n")
mtext(text="Final Biomass (g)", outer=T, side=2, line=-3.5)

#add dots with errorbars for individual species

g.sp <- summaryBy(dydt+predMass+AGR~Time+Taxa+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.sp.60<-subset(g.sp, Time=="60")
g.sp.Sn<- subset(g.sp.60,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw<- subset(g.sp.60,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw<- subset(g.sp.60,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww<- subset(g.sp.60,Location=="S"&Range=="wide"&Treatment=="Warmed")
g.sp.Nn<- subset(g.sp.60,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw<- subset(g.sp.60,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw<- subset(g.sp.60,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww<- subset(g.sp.60,Location=="N"&Range=="wide"&Treatment=="Warmed")

points(rep(0.5,length(g.sp.Sn$predMass.mean)),g.sp.Sn$predMass.mean, pch=16) #Sn home
points(rep(1.5,length(g.sp.Snw$predMass.mean)),g.sp.Snw$predMass.mean, pch=16) #Sn warmed
points(rep(3.5,length(g.sp.Sw$predMass.mean)),g.sp.Sw$predMass.mean, pch=16) #Sw home
points(rep(4.5,length(g.sp.Sww$predMass.mean)),g.sp.Sww$predMass.mean, pch=16) #Sw warmed
points(rep(6.5,length(g.sp.Nn$predMass.mean)),g.sp.Nn$predMass.mean, pch=16) #Nn home
points(rep(7.5,length(g.sp.Nnw$predMass.mean)),g.sp.Nnw$predMass.mean, pch=16) #Nn warmed 
points(rep(9.5,length(g.sp.Nw$predMass.mean)),g.sp.Nw$predMass.mean, pch=16) #Nw home
points(rep(10.5,length(g.sp.Nww$predMass.mean)),g.sp.Nww$predMass.mean, pch=16) #Nw warmed

means<-c(g.sp.Sn$predMass.mean,g.sp.Snw$predMass.mean,g.sp.Sw$predMass.mean,g.sp.Sww$predMass.mean,g.sp.Nn$predMass.mean,
         g.sp.Nnw$predMass.mean,g.sp.Nw$predMass.mean,g.sp.Nww$predMass.mean)
ses<-c(g.sp.Sn$predMass.standard.error,g.sp.Snw$predMass.standard.error,g.sp.Sw$predMass.standard.error,
       g.sp.Sww$predMass.standard.error,g.sp.Nn$predMass.standard.error,
       g.sp.Nnw$predMass.standard.error,g.sp.Nw$predMass.standard.error,g.sp.Nww$predMass.standard.error)
error.bar(c(rep(0.5,3),rep(1.5,3),rep(3.5,5),rep(4.5,5),rep(6.5,3),rep(7.5,3),rep(9.5,6),rep(10.5,6)),means,ses, length=0)

################################################################################
#RGR day 15
g.trt.15<-subset(g.trt, Time=="15")
g.trt.15$combo <- as.factor(paste(g.trt.15$Location,g.trt.15$Range,sep="_"))
windows(5.513333,4.135);par(mar=c(6,6,6,4))
target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
g.trt.15<-g.trt.15[match(target, g.trt.15$combotrt),]

bar<-barplot(g.trt.15$dydt.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0.05,0.18), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.8)),4)), 
             axis.lty=0)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar(bar,g.trt.15$dydt.mean,g.trt.15$dydt.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
legend(-0.1,0.16, legend=c("H",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),"red"),bty="n")
mtext(text="RGR Day 15", outer=T, side=2, line=-2)
mtext(text=expression((g~g^-1~day^-1)), outer=T, side=2, line=-3.5)

#add dots with errorbars for individual species

g.sp2 <- summaryBy(dydt+predMass+AGR~Time+Taxa+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.sp2.15<-subset(g.sp2, Time=="15")
g.sp.Sn2<- subset(g.sp2.15,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw2<- subset(g.sp2.15,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw2<- subset(g.sp2.15,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww2<- subset(g.sp2.15,Location=="S"&Range=="wide"&Treatment=="Warmed")
g.sp.Nn2<- subset(g.sp2.15,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw2<- subset(g.sp2.15,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw2<- subset(g.sp2.15,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww2<- subset(g.sp2.15,Location=="N"&Range=="wide"&Treatment=="Warmed")

points(rep(0.5,length(g.sp.Sn2$dydt.mean)),g.sp.Sn2$dydt.mean, pch=16) #Sn home
points(rep(1.5,length(g.sp.Snw2$dydt.mean)),g.sp.Snw2$dydt.mean, pch=16) #Sn warmed
points(rep(3.5,length(g.sp.Sw2$dydt.mean)),g.sp.Sw2$dydt.mean, pch=16) #Sw home
points(rep(4.5,length(g.sp.Sww2$dydt.mean)),g.sp.Sww2$dydt.mean, pch=16) #Sw warmed
points(rep(6.5,length(g.sp.Nn2$dydt.mean)),g.sp.Nn2$dydt.mean, pch=16) #Nn home
points(rep(7.5,length(g.sp.Nnw2$dydt.mean)),g.sp.Nnw2$dydt.mean, pch=16) #Nn warmed 
points(rep(9.5,length(g.sp.Nw2$dydt.mean)),g.sp.Nw2$dydt.mean, pch=16) #Nw home
points(rep(10.5,length(g.sp.Nww2$dydt.mean)),g.sp.Nww2$dydt.mean, pch=16) #Nw warmed

means<-c(g.sp.Sn2$dydt.mean,g.sp.Snw2$dydt.mean,g.sp.Sw2$dydt.mean,g.sp.Sww2$dydt.mean,g.sp.Nn2$dydt.mean,
         g.sp.Nnw2$dydt.mean,g.sp.Nw2$dydt.mean,g.sp.Nww2$dydt.mean)
ses<-c(g.sp.Sn2$dydt.standard.error,g.sp.Snw2$dydt.standard.error,g.sp.Sw2$dydt.standard.error,
       g.sp.Sww2$dydt.standard.error,g.sp.Nn2$dydt.standard.error,
       g.sp.Nnw2$dydt.standard.error,g.sp.Nw2$dydt.standard.error,g.sp.Nww2$dydt.standard.error)
error.bar(c(rep(0.5,3),rep(1.5,3),rep(3.5,5),rep(4.5,5),rep(6.5,3),rep(7.5,3),rep(9.5,6),rep(10.5,6)),means,ses, length=0)


################################################################################
#Respiration and Asat on mass basis with taxa dots
Asatm$Photom<-with(Asatm,Photo*SLA/10000)
asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=mean,keep.names=T)
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)

colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)

Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=c(mean,standard.error),keep.names=T)
#Rmass
ylims=c(5,20)
boxplot(Rmass.mean~Treatment*Range,data=subset(Rmass.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)

mtext(text=expression(R["mass"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=2.5)
legend(-0.2,21.2,"Temperate", bty="n", cex=1.3)

g.sp.Sn3<- subset(Rmass.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw3<- subset(Rmass.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw3<- subset(Rmass.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww3<- subset(Rmass.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn3$Rmass.mean)),g.sp.Sn3$Rmass.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw3$Rmass.mean)),g.sp.Snw3$Rmass.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw3$Rmass.mean)),g.sp.Sw3$Rmass.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww3$Rmass.mean)),g.sp.Sww3$Rmass.mean, pch=16) #Sw warmed
meansR1<-c(g.sp.Sn3$Rmass.mean,g.sp.Snw3$Rmass.mean,g.sp.Sw3$Rmass.mean,g.sp.Sww3$Rmass.mean)
sesR1<-c(g.sp.Sn3$Rmass.standard.error,g.sp.Snw3$Rmass.standard.error,g.sp.Sw3$Rmass.standard.error,
       g.sp.Sww3$Rmass.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansR1,sesR1, length=0)

boxplot(Rmass.mean~Treatment*Range,data=subset(Rmass.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
legend(-0.2,21.2,"Tropical", bty="n", cex=1.3)
legend("topright",legend=c("H","W"),ncol=1, fill=colors,cex=1, bty="n")

g.sp.Nn3<- subset(Rmass.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw3<- subset(Rmass.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw3<- subset(Rmass.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww3<- subset(Rmass.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn3$Rmass.mean)),g.sp.Nn3$Rmass.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw3$Rmass.mean)),g.sp.Nnw3$Rmass.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw3$Rmass.mean)),g.sp.Nw3$Rmass.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww3$Rmass.mean)),g.sp.Nww3$Rmass.mean, pch=16) #Nw warmed

meansR2<-c(g.sp.Nn3$Rmass.mean,g.sp.Nnw3$Rmass.mean,g.sp.Nw3$Rmass.mean,g.sp.Nww3$Rmass.mean)
sesR2<-c(g.sp.Nn3$Rmass.standard.error,g.sp.Nnw3$Rmass.standard.error,g.sp.Nw3$Rmass.standard.error,g.sp.Nww3$Rmass.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansR2,sesR2, length=0)

asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=c(mean,standard.error),keep.names=T)

#Asat mass
ylims=c(0.25,0.75)
boxplot(Photom.mean~Treatment*Range,data=subset(asatm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
axis(side=1,at=c(1.5,3.5),labels=levels(asatm.tm$Range),las=1,cex.axis=1.5)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=2.5)

g.sp.Sn4<- subset(asatm.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw4<- subset(asatm.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw4<- subset(asatm.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww4<- subset(asatm.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn4$Photom.mean)),g.sp.Sn4$Photom.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw4$Photom.mean)),g.sp.Snw4$Photom.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw4$Photom.mean)),g.sp.Sw4$Photom.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww4$Photom.mean)),g.sp.Sww4$Photom.mean, pch=16) #Sw warmed

meansA1<-c(g.sp.Sn4$Photom.mean,g.sp.Snw4$Photom.mean,g.sp.Sw4$Photom.mean,g.sp.Sww4$Photom.mean)
sesA1<-c(g.sp.Sn4$Photom.standard.error,g.sp.Snw4$Photom.standard.error,g.sp.Sw4$Photom.standard.error,
       g.sp.Sww4$Photom.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)

boxplot(Photom.mean~Treatment*Range,data=subset(asatm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(asatm.tm$Range),las=1,cex.axis=1.5)
g.sp.Nn4<- subset(asatm.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw4<- subset(asatm.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw4<- subset(asatm.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww4<- subset(asatm.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn4$Photom.mean)),g.sp.Nn4$Photom.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw4$Photom.mean)),g.sp.Nnw4$Photom.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw4$Photom.mean)),g.sp.Nw4$Photom.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww4$Photom.mean)),g.sp.Nww4$Photom.mean, pch=16) #Nw warmed

meansA2<-c(g.sp.Nn4$Photom.mean,g.sp.Nnw4$Photom.mean,g.sp.Nw4$Photom.mean,g.sp.Nww4$Photom.mean)
sesA2<-c(g.sp.Nn4$Photom.standard.error,g.sp.Nnw4$Photom.standard.error,g.sp.Nw4$Photom.standard.error,g.sp.Nww4$Photom.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)
