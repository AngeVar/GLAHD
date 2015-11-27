

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
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.02)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.02)
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
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = F, tck = 0.02)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.02)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
mtext("Tropical",3,line=-2,at=12.5, cex=1.5, outer=F)

############################################
#Enhancement ratios
#get estimate of error using taxa means?


par(fig=c(0,0.5,0,0.4), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))
ber<- summaryBy(predMass~Taxa+Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Taxa","Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Taxa","Time","Range","Location"))
bio$ber<- with(bio, (predMassWarm-predMass)/predMass)
bio$combo<-as.factor(paste(bio$Location,bio$Range,sep="_"))
BERsm<- summaryBy(ber~Time+combo, data=bio, FUN=c(mean,standard.error))

SBERn<- subset(BERsm, combo =="S_narrow"); SBERn<-SBERn[order(SBERn$Time),]
SBERw<- subset(BERsm, combo =="S_wide"); SBERw<-SBERw[order(SBERw$Time),]
NBERn<- subset(BERsm, combo =="N_narrow"); NBERn<-NBERn[order(NBERn$Time),]
NBERw<- subset(BERsm, combo =="N_wide"); NBERw<-NBERw[order(NBERw$Time),]

plotBy(ber.mean~Time,data=SBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main="", ylim=c(-0.25,1.1), xaxt='n',xlim=c(0,64),yaxp  = c(-0.5, 1, 3),yaxs="i",xaxs="i",ylab="", yaxt='n')
lines(ber.mean~Time,data=SBERw,col="black",lty=1,lwd=2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.02)
abline(h=0)
mtext(expression(Delta~"Biomass (%)"), side=2, line=3, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2)
#legend("topleft", legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1,lwd=2, bty="n")
  
axis(side = 2, at = seq(from=-1,to=2,by=0.5), labels = T, tck = 0.02)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
axis(side = 4, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.02)
SBERn$upper <- with(SBERn, ber.mean+ber.standard.error)
SBERn$lower <- with(SBERn, ber.mean-ber.standard.error)
xcord<-c(SBERn$Time,rev(SBERn$Time))
ycord<-c(SBERn$lower,rev(SBERn$upper))
polygon(x = xcord, y = ycord, col = alpha(grey,0.3), border = NULL,xpd=F)

par(fig=c(0.5,1,0,0.4), new=TRUE,mar=c(0,0,0,0),oma=c(4,6,3,3))
plotBy(ber.mean~Time,data=NBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main = "", ylim=c(-0.25,1.1),yaxt='n',xaxt='n', xlim=c(0,64),yaxp  = c(-0.5, 1, 3),yaxs="i",xaxs="i")
lines(ber.mean~Time,data=NBERw,col="black",lty=1,lwd=2)
#mtext(side=3,outer=T,text="Tropical", adj=0.8)
abline(h=0)
legend("topright", legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1,lwd=2, bty="n")
#text(2,1.95,labels="B",cex=2,adj=0.5)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.02)
mtext(text="Time (Days)", side=1, line=3, cex=1.2)
axis(side = 2, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.02)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.02)
axis(side = 4, at = seq(from=-1,to=2,by=0.5), labels = F, tck = 0.02)

################################################################################
#Final biomass per taxa
g.trt.60<-subset(g.trt, Time=="60")
g.trt.60$combo <- as.factor(paste(g.trt.60$Location,g.trt.60$Range,sep="_"))
windows(5.513333,4.135);par(mar=c(6,6,6,4))
target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
g.trt.60<-g.trt.60[match(target, g.trt.60$combotrt),]

bar<-barplot(g.trt.60$predMass.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0,88), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
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
# mtext("Temperate ",3,line=-1.25,at=2.5,cex=1.15, outer=F)
# mtext("Tropical ",3,line=-1.25,at=8.5,cex=1.15, outer=F)
#legend(-0.1,80, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="Final biomass", outer=T, side=2, line=-1.5, cex=1.3)
mtext(text="(g)", outer=T, side=2, line=-3.5)
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

bar<-barplot(g.trt.15$dydt.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0.04,0.16), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
             axis.lty=0)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar(bar,g.trt.15$dydt.mean,g.trt.15$dydt.standard.error )
#mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
legend(-0.1,0.145, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="RGR Day 15", outer=T, side=2, line=-1.5, cex=1.3)
mtext(text=expression((g~g^-1~day^-1)), outer=T, side=2, line=-3.5)

#% increase
(0.10174500-0.08872640)/0.08872640 #14.7% increase Sn
(0.10372688-0.08076066)/0.08076066 #28.4% increase Sw
(0.10204003-0.10442129)/0.10442129 #-2.2% Nn
(0.11170755-0.11088574)/0.11088574 #+0.7% Nw
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
#LMF,SLA and LAR on Day 15 from allometric data 
#- does not capture the same results as the regressions across the entire experiment

rate.trt <- summaryBy(LMF+SLA+LAR+ULR~Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))
rate11<-subset(rate.trt, Time == 11)
target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
rate11<-rate11[match(target, rate11$combotrt),]

windows(5.513333,4.135);par(mar=c(6,6,6,4))
bar<-barplot(rate11$LMF.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0.42,0.65), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
             axis.lty=0)
error.bar(bar,rate11$LMF.mean,rate11$LMF.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
legend(8.2,0.6, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="LMF Day 15", outer=T, side=2, line=-3)

windows(5.513333,4.135);par(mar=c(6,6,6,4))
bar<-barplot(rate11$SLA.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(150,260), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
             axis.lty=0)
error.bar(bar,rate11$SLA.mean,rate11$SLA.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
# legend(-0.1,240, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="SLA Day 15", outer=T, side=2, line=-2)

windows(5.513333,4.135);par(mar=c(6,6,6,4))
bar<-barplot(rate11$LAR.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(82,148), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
             axis.lty=0)
error.bar(bar,rate11$LAR.mean,rate11$LAR.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
#legend(-0.1,140, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="LAR Day 15", outer=T, side=2, line=-2)

windows(5.513333,4.135);par(mar=c(6,6,6,4))
bar<-barplot(rate11$ULR.mean, space=c(0,0,1,0,1,0,1,0),ylim=c(0.0005,0.0013), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.6)),4)), 
             axis.lty=0)
error.bar(bar,rate11$ULR.mean,rate11$ULR.standard.error )
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate ",3,line=-1.5,at=2.5,cex=1.15, outer=F)
mtext("Tropical ",3,line=-1.5,at=8.5,cex=1.15, outer=F)
#legend(-0.1,140, legend=c("Home",expression(paste("+3.5",degree,"C"))),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),alpha("red",0.6)),bty="n")
mtext(text="NAR Day 15", outer=T, side=2, line=-2)


fm1m1 <- lme(1/(LMF)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))#, method="ML")
plot(fm1m1,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m1,1/(LMF)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m1,1/(LMF)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m1, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m1$residuals[,1])
anova(fm1m1)
plot(effect("Treatment",fm1m1),multiline=T)

fm1m2 <- lme(SLA~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))#, method="ML")
plot(fm1m2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m2,SLA~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m2,SLA~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m2$residuals[,1])
anova(fm1m2)
plot(effect("Treatment:Location",fm1m2),multiline=T) 

fm1m3 <- lme(LAR~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))#, method="ML")
plot(fm1m3,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m3,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m3,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m3, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m3$residuals[,1])
anova(fm1m3)
plot(effect("Treatment:Location",fm1m3),multiline=T)  #<0.0001
plot(effect("Treatment:Range",fm1m3),multiline=T)     #0.0545

fm1m4 <- lme(ULR~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=subset(rate, Time==11))#, method="ML")
plot(fm1m4,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m4,ULR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m4,ULR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m4, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m4$residuals[,1])
anova(fm1m4)
plot(effect("Treatment",fm1m4),multiline=T)   #<0.0001
plot(effect("Location",fm1m4),multiline=T)     #0.0185

################################################################################
#Respiration and Asat on mass basis with taxa dots
Asatm$Photom<-with(Asatm,Photo*SLA/10000)
asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=c(mean,standard.error))
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=c(mean,standard.error))

colors <- c(alpha("black",0.6),alpha("red",0.6))
windows(5.513333,4.135);par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(6,6,6,4),cex.axis=1.2)
asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=c(mean,standard.error),keep.names=T)

#Asat mass
ylims=c(0.35,0.7)
boxplot(Photom.mean~Treatment*Range,data=subset(asatm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=0,majorn=4, labels=c(1,0),frame.plot=T,las=1)
#mtext(text=expression(A["sat"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text="Photosynthesis",side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,line=2.5)
#legend("top","Temperate", bty="n", cex=1.15)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1)
# g.sp.Sn4<- subset(asatm.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
# g.sp.Snw4<- subset(asatm.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
# g.sp.Sw4<- subset(asatm.tm,Location=="S"&Range=="wide"&Treatment=="Home")
# g.sp.Sww4<- subset(asatm.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
# points(rep(1,length(g.sp.Sn4$Photom.mean)),g.sp.Sn4$Photom.mean, pch=16) #Sn home
# points(rep(2,length(g.sp.Snw4$Photom.mean)),g.sp.Snw4$Photom.mean, pch=16) #Sn warmed
# points(rep(3,length(g.sp.Sw4$Photom.mean)),g.sp.Sw4$Photom.mean, pch=16) #Sw home
# points(rep(4,length(g.sp.Sww4$Photom.mean)),g.sp.Sww4$Photom.mean, pch=16) #Sw warmed
# 
# meansA1<-c(g.sp.Sn4$Photom.mean,g.sp.Snw4$Photom.mean,g.sp.Sw4$Photom.mean,g.sp.Sww4$Photom.mean)
# sesA1<-c(g.sp.Sn4$Photom.standard.error,g.sp.Snw4$Photom.standard.error,g.sp.Sw4$Photom.standard.error,
#          g.sp.Sww4$Photom.standard.error)
# error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)
# 
boxplot(Photom.mean~Treatment*Range,data=subset(asatm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=0,majorn=4,labels=c(0,0),frame.plot=T,las=1)
#legend("top","Tropical", bty="n", cex=1.15)
#legend(2,0.8,legend=c("Home",expression(paste("+3.5",degree,"C"))),cex=1, fill=colors,bty="n")
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1)
# g.sp.Nn4<- subset(asatm.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
# g.sp.Nnw4<- subset(asatm.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
# g.sp.Nw4<- subset(asatm.tm,Location=="N"&Range=="wide"&Treatment=="Home")
# g.sp.Nww4<- subset(asatm.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
# points(rep(1,length(g.sp.Nn4$Photom.mean)),g.sp.Nn4$Photom.mean, pch=16) #Nn home
# points(rep(2,length(g.sp.Nnw4$Photom.mean)),g.sp.Nnw4$Photom.mean, pch=16) #Nn warmed 
# points(rep(3,length(g.sp.Nw4$Photom.mean)),g.sp.Nw4$Photom.mean, pch=16) #Nw home
# points(rep(4,length(g.sp.Nww4$Photom.mean)),g.sp.Nww4$Photom.mean, pch=16) #Nw warmed
# 
# meansA2<-c(g.sp.Nn4$Photom.mean,g.sp.Nnw4$Photom.mean,g.sp.Nw4$Photom.mean,g.sp.Nww4$Photom.mean)
# sesA2<-c(g.sp.Nn4$Photom.standard.error,g.sp.Nnw4$Photom.standard.error,g.sp.Nw4$Photom.standard.error,g.sp.Nww4$Photom.standard.error)
# error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)
windows(5.513333,4.135);par(mfrow=c(1,2),mar=c(0,0,0,0),oma=c(6,6,6,4),cex.axis=1.2)
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=c(mean,standard.error),keep.names=T)
#Rmass
ylims=c(5,18)
boxplot(Rmass.mean~Treatment*Range,data=subset(Rmass.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=c(5,5),majorn=3,labels=c(1,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1)
#mtext(text=expression(R["mass"]),side=2,outer=F,cex=1.3,line=4)
mtext(text="Respiration",side=2,outer=F,cex=1.2,line=4)
mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=F,cex=1.1,line=2.5)
mtext(text="Temperate", side=3, cex=1.15)
# g.sp.Sn3<- subset(Rmass.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
# g.sp.Snw3<- subset(Rmass.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
# g.sp.Sw3<- subset(Rmass.tm,Location=="S"&Range=="wide"&Treatment=="Home")
# g.sp.Sww3<- subset(Rmass.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
# points(rep(1,length(g.sp.Sn3$Rmass.mean)),g.sp.Sn3$Rmass.mean, pch=16) #Sn home
# points(rep(2,length(g.sp.Snw3$Rmass.mean)),g.sp.Snw3$Rmass.mean, pch=16) #Sn warmed
# points(rep(3,length(g.sp.Sw3$Rmass.mean)),g.sp.Sw3$Rmass.mean, pch=16) #Sw home
# points(rep(4,length(g.sp.Sww3$Rmass.mean)),g.sp.Sww3$Rmass.mean, pch=16) #Sw warmed
# meansR1<-c(g.sp.Sn3$Rmass.mean,g.sp.Snw3$Rmass.mean,g.sp.Sw3$Rmass.mean,g.sp.Sww3$Rmass.mean)
# sesR1<-c(g.sp.Sn3$Rmass.standard.error,g.sp.Snw3$Rmass.standard.error,g.sp.Sw3$Rmass.standard.error,
#        g.sp.Sww3$Rmass.standard.error)
# error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansR1,sesR1, length=0)

boxplot(Rmass.mean~Treatment*Range,data=subset(Rmass.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=c(5,5),majorn=3,labels=c(0,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1)
mtext(text="Tropical", side=3, cex=1.15)
legend("topright",legend=c("Home",expression(paste("+3.5",degree,"C"))),cex=1, fill=colors,bty="n")
# g.sp.Nn3<- subset(Rmass.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
# g.sp.Nnw3<- subset(Rmass.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
# g.sp.Nw3<- subset(Rmass.tm,Location=="N"&Range=="wide"&Treatment=="Home")
# g.sp.Nww3<- subset(Rmass.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
# points(rep(1,length(g.sp.Nn3$Rmass.mean)),g.sp.Nn3$Rmass.mean, pch=16) #Nn home
# points(rep(2,length(g.sp.Nnw3$Rmass.mean)),g.sp.Nnw3$Rmass.mean, pch=16) #Nn warmed 
# points(rep(3,length(g.sp.Nw3$Rmass.mean)),g.sp.Nw3$Rmass.mean, pch=16) #Nw home
# points(rep(4,length(g.sp.Nww3$Rmass.mean)),g.sp.Nww3$Rmass.mean, pch=16) #Nw warmed
# 
# meansR2<-c(g.sp.Nn3$Rmass.mean,g.sp.Nnw3$Rmass.mean,g.sp.Nw3$Rmass.mean,g.sp.Nww3$Rmass.mean)
# sesR2<-c(g.sp.Nn3$Rmass.standard.error,g.sp.Nnw3$Rmass.standard.error,g.sp.Nw3$Rmass.standard.error,g.sp.Nww3$Rmass.standard.error)
# error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansR2,sesR2, length=0)


################################################################################
#Jmax and Vcmax with dots
jmax.tm <- summaryBy(Jmax~Taxa+Treatment+Location+Range,data=acifit,FUN=c(mean,standard.error))
vcmax.tm <- summaryBy(Vcmax~Taxa+Treatment+Location+Range,data=acifit,FUN=c(mean,standard.error))
jmaxm.tm <- summaryBy(Jmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=c(mean,standard.error))
vcmaxm.tm <- summaryBy(Vcmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=c(mean,standard.error))

colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)
#Jmax
ylims=c(90,250)
boxplot(Jmax.mean~Treatment*Range,data=subset(jmax.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(J["max"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.9,line=2.5)
legend(-0,260,"Temperate", bty="n", cex=1.3)

g.sp.Sn4<- subset(jmax.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw4<- subset(jmax.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw4<- subset(jmax.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww4<- subset(jmax.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn4$Jmax.mean)),g.sp.Sn4$Jmax.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw4$Jmax.mean)),g.sp.Snw4$Jmax.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw4$Jmax.mean)),g.sp.Sw4$Jmax.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww4$Jmax.mean)),g.sp.Sww4$Jmax.mean, pch=16) #Sw warmed

meansA1<-c(g.sp.Sn4$Jmax.mean,g.sp.Snw4$Jmax.mean,g.sp.Sw4$Jmax.mean,g.sp.Sww4$Jmax.mean)
sesA1<-c(g.sp.Sn4$Jmax.standard.error,g.sp.Snw4$Jmax.standard.error,g.sp.Sw4$Jmax.standard.error,
         g.sp.Sww4$Jmax.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)

boxplot(Jmax.mean~Treatment*Range,data=subset(jmax.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
legend(-0,260,"Tropical", bty="n", cex=1.3)
g.sp.Nn4<- subset(jmax.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw4<- subset(jmax.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw4<- subset(jmax.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww4<- subset(jmax.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn4$Jmax.mean)),g.sp.Nn4$Jmax.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw4$Jmax.mean)),g.sp.Nnw4$Jmax.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw4$Jmax.mean)),g.sp.Nw4$Jmax.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww4$Jmax.mean)),g.sp.Nww4$Jmax.mean, pch=16) #Nw warmed

meansA2<-c(g.sp.Nn4$Jmax.mean,g.sp.Nnw4$Jmax.mean,g.sp.Nw4$Jmax.mean,g.sp.Nww4$Jmax.mean)
sesA2<-c(g.sp.Nn4$Jmax.standard.error,g.sp.Nnw4$Jmax.standard.error,g.sp.Nw4$Jmax.standard.error,g.sp.Nww4$Jmax.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)


#Jmax mass
ylims=c(2,4.3)
boxplot(Jmaxm.mean~Treatment*Range,data=subset(jmaxm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
axis(side=1,at=c(1.5,3.5),labels=levels(jmaxm.tm$Range),las=1,cex.axis=1.5)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(J["max"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=2.5)
g.sp.Sn4<- subset(jmaxm.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw4<- subset(jmaxm.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw4<- subset(jmaxm.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww4<- subset(jmaxm.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn4$Jmaxm.mean)),g.sp.Sn4$Jmaxm.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw4$Jmaxm.mean)),g.sp.Snw4$Jmaxm.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw4$Jmaxm.mean)),g.sp.Sw4$Jmaxm.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww4$Jmaxm.mean)),g.sp.Sww4$Jmaxm.mean, pch=16) #Sw warmed

meansA1<-c(g.sp.Sn4$Jmaxm.mean,g.sp.Snw4$Jmaxm.mean,g.sp.Sw4$Jmaxm.mean,g.sp.Sww4$Jmaxm.mean)
sesA1<-c(g.sp.Sn4$Jmaxm.standard.error,g.sp.Snw4$Jmaxm.standard.error,g.sp.Sw4$Jmaxm.standard.error,
         g.sp.Sww4$Jmaxm.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)

boxplot(Jmaxm.mean~Treatment*Range,data=subset(jmaxm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(jmaxm.tm$Range),las=1,cex.axis=1.5)
g.sp.Nn4<- subset(jmaxm.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw4<- subset(jmaxm.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw4<- subset(jmaxm.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww4<- subset(jmaxm.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn4$Jmaxm.mean)),g.sp.Nn4$Jmaxm.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw4$Jmaxm.mean)),g.sp.Nnw4$Jmaxm.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw4$Jmaxm.mean)),g.sp.Nw4$Jmaxm.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww4$Jmaxm.mean)),g.sp.Nww4$Jmaxm.mean, pch=16) #Nw warmed

meansA2<-c(g.sp.Nn4$Jmaxm.mean,g.sp.Nnw4$Jmaxm.mean,g.sp.Nw4$Jmaxm.mean,g.sp.Nww4$Jmaxm.mean)
sesA2<-c(g.sp.Nn4$Jmaxm.standard.error,g.sp.Nnw4$Jmaxm.standard.error,g.sp.Nw4$Jmaxm.standard.error,g.sp.Nww4$Jmaxm.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)

###########################################################################################
colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)
#Vcmax
ylims=c(60,240)
boxplot(Vcmax.mean~Treatment*Range,data=subset(vcmax.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.9,line=2.5)
legend(-0,250,"Temperate", bty="n", cex=1.3)

g.sp.Sn4<- subset(vcmax.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw4<- subset(vcmax.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw4<- subset(vcmax.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww4<- subset(vcmax.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn4$Vcmax.mean)),g.sp.Sn4$Vcmax.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw4$Vcmax.mean)),g.sp.Snw4$Vcmax.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw4$Vcmax.mean)),g.sp.Sw4$Vcmax.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww4$Vcmax.mean)),g.sp.Sww4$Vcmax.mean, pch=16) #Sw warmed

meansA1<-c(g.sp.Sn4$Vcmax.mean,g.sp.Snw4$Vcmax.mean,g.sp.Sw4$Vcmax.mean,g.sp.Sww4$Vcmax.mean)
sesA1<-c(g.sp.Sn4$Vcmax.standard.error,g.sp.Snw4$Vcmax.standard.error,g.sp.Sw4$Vcmax.standard.error,
         g.sp.Sww4$Vcmax.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)

boxplot(Vcmax.mean~Treatment*Range,data=subset(vcmax.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
legend(-0,250,"Tropical", bty="n", cex=1.3)
g.sp.Nn4<- subset(vcmax.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw4<- subset(vcmax.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw4<- subset(vcmax.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww4<- subset(vcmax.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn4$Vcmax.mean)),g.sp.Nn4$Vcmax.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw4$Vcmax.mean)),g.sp.Nnw4$Vcmax.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw4$Vcmax.mean)),g.sp.Nw4$Vcmax.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww4$Vcmax.mean)),g.sp.Nww4$Vcmax.mean, pch=16) #Nw warmed

meansA2<-c(g.sp.Nn4$Vcmax.mean,g.sp.Nnw4$Vcmax.mean,g.sp.Nw4$Vcmax.mean,g.sp.Nww4$Vcmax.mean)
sesA2<-c(g.sp.Nn4$Vcmax.standard.error,g.sp.Nnw4$Vcmax.standard.error,g.sp.Nw4$Vcmax.standard.error,g.sp.Nww4$Vcmax.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)


#Vcmax mass
ylims=c(1,4.8)
boxplot(Vcmaxm.mean~Treatment*Range,data=subset(vcmaxm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
axis(side=1,at=c(1.5,3.5),labels=levels(vcmaxm.tm$Range),las=1,cex.axis=1.5)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=2.5)
g.sp.Sn4<- subset(vcmaxm.tm,Location=="S"&Range=="narrow"&Treatment=="Home")
g.sp.Snw4<- subset(vcmaxm.tm,Location=="S"&Range=="narrow"&Treatment=="Warmed")
g.sp.Sw4<- subset(vcmaxm.tm,Location=="S"&Range=="wide"&Treatment=="Home")
g.sp.Sww4<- subset(vcmaxm.tm,Location=="S"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Sn4$Vcmaxm.mean)),g.sp.Sn4$Vcmaxm.mean, pch=16) #Sn home
points(rep(2,length(g.sp.Snw4$Vcmaxm.mean)),g.sp.Snw4$Vcmaxm.mean, pch=16) #Sn warmed
points(rep(3,length(g.sp.Sw4$Vcmaxm.mean)),g.sp.Sw4$Vcmaxm.mean, pch=16) #Sw home
points(rep(4,length(g.sp.Sww4$Vcmaxm.mean)),g.sp.Sww4$Vcmaxm.mean, pch=16) #Sw warmed

meansA1<-c(g.sp.Sn4$Vcmaxm.mean,g.sp.Snw4$Vcmaxm.mean,g.sp.Sw4$Vcmaxm.mean,g.sp.Sww4$Vcmaxm.mean)
sesA1<-c(g.sp.Sn4$Vcmaxm.standard.error,g.sp.Snw4$Vcmaxm.standard.error,g.sp.Sw4$Vcmaxm.standard.error,
         g.sp.Sww4$Vcmaxm.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,5),rep(4,5)),meansA1,sesA1, length=0)

boxplot(Vcmaxm.mean~Treatment*Range,data=subset(vcmaxm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(jmaxm.tm$Range),las=1,cex.axis=1.5)
g.sp.Nn4<- subset(vcmaxm.tm,Location=="N"&Range=="narrow"&Treatment=="Home")
g.sp.Nnw4<- subset(vcmaxm.tm,Location=="N"&Range=="narrow"&Treatment=="Warmed")
g.sp.Nw4<- subset(vcmaxm.tm,Location=="N"&Range=="wide"&Treatment=="Home")
g.sp.Nww4<- subset(vcmaxm.tm,Location=="N"&Range=="wide"&Treatment=="Warmed")
points(rep(1,length(g.sp.Nn4$Vcmaxm.mean)),g.sp.Nn4$Vcmaxm.mean, pch=16) #Nn home
points(rep(2,length(g.sp.Nnw4$Vcmaxm.mean)),g.sp.Nnw4$Vcmaxm.mean, pch=16) #Nn warmed 
points(rep(3,length(g.sp.Nw4$Vcmaxm.mean)),g.sp.Nw4$Vcmaxm.mean, pch=16) #Nw home
points(rep(4,length(g.sp.Nww4$Vcmaxm.mean)),g.sp.Nww4$Vcmaxm.mean, pch=16) #Nw warmed

meansA2<-c(g.sp.Nn4$Vcmaxm.mean,g.sp.Nnw4$Vcmaxm.mean,g.sp.Nw4$Vcmaxm.mean,g.sp.Nww4$Vcmaxm.mean)
sesA2<-c(g.sp.Nn4$Vcmaxm.standard.error,g.sp.Nnw4$Vcmaxm.standard.error,g.sp.Nw4$Vcmaxm.standard.error,g.sp.Nww4$Vcmaxm.standard.error)
error.bar(c(rep(1,3),rep(2,3),rep(3,6),rep(4,6)),meansA2,sesA2, length=0)
