#Figures for GAM output
source("R/GLAHD_gamfits.R") #gives you the gamfits2 output
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
###############################################################################################################

#Figure 2: Mass over Time
g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(1,70),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.N$Time,y=g.trt.N$predMass.mean,
                                SE=g.trt.N$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)

axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)

legend(0,60, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=12.5, outer=F)

plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.S$Time,y=g.trt.S$predMass.mean,
                                SE=g.trt.S$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
mtext(text="Time (Days)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=14, outer=F)
mtext(text="Total biomass (g)", outer=T, side=2, line=1)

###############################################################################################################
#Figure 3:RGR over Time


windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(dydt.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0.01,0.13),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.N$Time,y=g.trt.N$dydt.mean,
                                SE=g.trt.N$dydt.standard.error,direction="updown",
                                col="black",0))
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
abline(v=15, lty=2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.12,by=0.02), labels = T, tck = 0.01)
axis(side = 4, at = seq(from=0,to=0.12,by=0.02), labels = F, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)

legend(0,60, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=46.5, outer=F)

plotBy(dydt.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0.01,0.13),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.S$Time,y=g.trt.S$dydt.mean,
                                SE=g.trt.S$dydt.standard.error,direction="updown",
                                col="black",0))
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
abline(v=15, lty=2)
mtext(text="Time (Days)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.12,by=0.02), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=0.12,by=0.02), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=45, outer=F)
mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=0.5)

###############################################################################################################
#Figure 4: Bargraph of RGR at 15 days

g.trt.15<-subset(g.trt, Time=="15")
g.trt.15$combo <- as.factor(paste(g.trt.15$Location,g.trt.15$Range,sep="_"))

target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
g.trt.15<-g.trt.15[match(target, g.trt.15$combotrt),]

bar<-barplot(g.trt.15$dydt.mean, space=c(0,0,0.5,0,0.5,0,0.5),ylim=c(0.05,0.12), xpd = F, col=c(rep(c("black","red"),4)), axis.lty=1, 
        names.arg=c("H","W","H","W","H","W","H","W"))
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar(bar,g.trt.15$dydt.mean,g.trt.15$dydt.standard.error )
mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=-1.5)
mtext("  Narrow                Wide                 Narrow              Wide    ", side=1, line = 2.5)
abline(v=4.75)
box()
mtext("Temperate Eucalyptus",3,line=-1.5,at=1.7, outer=F)
mtext("Tropical Eucalyptus",3,line=-1.5,at=6.5, outer=F)

#LA over Total mass
#- load libraries from script
source("R/loadLibraries.R");source("R/gamplotfunctions.R");library(scales)
#- read in the data, do a few conversions
dat <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Totmass <- base::rowSums(dat[,11:13]) #total mass is the sum of leaf, stem, and root mass
dat$Treat <- as.factor(ifelse(dat$Pot < 20, "Home",
                              ifelse(dat$Pot>=40,"Pre","Warmed")))
dat$Taxa <- factor(dat$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
dat[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
dat2 <- subset(dat,Treat!="Pre");dat2$Treat <- factor(dat2$Treat)
dat2$Range<- ifelse(dat2$Species == "CAM"|dat2$Species == "TER", "wide","narrow")


la.trt <- summaryBy(Leafarea~Totmass+Treatment+Location+Range,data=dat2,FUN=c(mean,standard.error))
la.trt$combotrt <- as.factor(paste(la.trt$Location,la.trt$Range,la.trt$Treatment,sep="_"))
la.trt.S<- subset(la.trt, Location == "S")
la.trt.N<- subset(la.trt, Location == "N")

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0,10000),xlim=c(0,110),pch=1,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=la.trt.N$Totmass,y=la.trt.N$Leafarea.mean,
                                SE=la.trt.N$Leafarea.standard.error,direction="updown",
                                col="black",0))
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=1,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)

axis(side = 1, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=10000,by=1000), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=10000,by=1000), labels = F, tck = 0.01)


legend(0,60, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=22, outer=F)

plotBy(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0,10000),xlim=c(0,110),pch=1,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=la.trt.S$Totmass,y=la.trt.S$Leafarea.mean,
                                SE=la.trt.S$Leafarea.standard.error,direction="updown",
                                col="black",0))
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=1,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
mtext(text="Total biomass (g)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=110,by=10), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=10000,by=1000), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=10000,by=1000), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=25, outer=F)
mtext(text="Leaf area", outer=T, side=2, line=1)

###############################################################################################################

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

rate <- merge(gamfits2,dat2,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,18:19)]# merge with dat2 to get LA
rate$LAR <- rate$leafArea/rate$TotMass
rate$ULR <- with(rate,dydt/LAR)
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))

windows(60,40);par(mfrow=c(4,2), mar=c(0.5,6,0.5,3), oma=c(5,1,1,1))

palettS <- c("black","grey","blue","cyan"); palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|combotrt,data=sg.trt,col=palettS,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|combotrt,data=sg.trt,col=palettS,,
       legend=F,pch=3, ylab= "",xaxt='n', type="l",ylim=c(0.04,0.12),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(dydt.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3, ylab= "",xaxt='n', type="l",ylim=c(0.04,0.12),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(AGR.mean~Time|combotrt,data=sg.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",,log="y",
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))
  

plotBy(AGR.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",log="y",
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|combotrt,data=srate.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|combotrt,data=nrate.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)
 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#RGR and AGR over Mass
palettS <- c("black","grey","blue","cyan"); palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

windows(60,30);par(mfrow=c(3,2), mar=c(0.5,6,0.5,3), oma=c(5,1,1,1))

palettS <- c("black","grey","blue","cyan")
palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"))
ng.trt<- droplevels(subset(g.trt, Location == "N"))

plotBy(dydt.mean~predMass.mean|combotrt,data=sg.trt,col=palettS,,
       legendwhere="topright",pch=3, ylim=c(0.04,0.12),ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=sg.trt$predMass.mean,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(dydt.mean~predMass.mean|combotrt,data=ng.trt,col=palettN,
       legendwhere="topright",pch=3, ylim=c(0.04,0.12),ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=ng.trt$predMass.mean,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(AGR.mean~predMass.mean|combotrt,data=sg.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=sg.trt$predMass.mean,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(AGR.mean~predMass.mean|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=ng.trt$predMass.mean,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~predMass.mean|combotrt,data=srate.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=srate.trt$predMass.mean,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Mass (g)", side=1, line=3)

plotBy(LAR.mean~predMass.mean|combotrt,data=nrate.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=nrate.trt$predMass.mean,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Mass (g)", side=1, line=3)



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

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.9,2), xlim=c(0,64),yaxp  = c(0.9, 2, 11))
text(2,1.95,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.9,2),yaxt='n', xlim=c(0,64),yaxp  = c(0.9, 2, 11))
abline(h=1)
legend(45,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
text(2,1.95,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~predMass,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.9,2), xlim=c(0,37),xaxp  = c(0, 35, 7),yaxp  = c(0.9, 2, 11))
mtext(text="Biomass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~predMass,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(2,1.95,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~predMass,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
legend=F,type="p", main = "", ylim=c(0.9,2),yaxt='n', xlim=c(0,75),xaxp  = c(0, 80, 8),yaxp  = c(0.9, 2, 11))
points(ber~predMass,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.95,labels="D",cex=2,adj=0.5)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]


#Relative enhancement of biomass using rate data
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=rate, FUN=mean, keep.names=T) 
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
legend(53,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=0.8)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

#Relative enhancement of biomass per species
ber<- summaryBy(predMass~Time+Species+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Species", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Species","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))
SBER$Species2<-factor(SBER$Species,levels(SBER$Species)[c(2,5,1,3,4)]);NBER$Species2<-factor(NBER$Species,levels(NBER$Species)[c(2,5,1,3,4)])

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Species2,data=SBER,col=c("black","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,2.5), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(54,2.5, legend=c(levels(SBER$Species2)), pch=16, col=c("black","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Species2,data=NBER,col=c("black","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,2.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(54,2.5, legend=c(levels(NBER$Species2)), pch=16, col=c("black","orange","grey","red","magenta"), cex=0.8)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Relative enhancement of biomass per taxa
ber<- summaryBy(predMass~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,3.5), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(54,3.7, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,3.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(54,3.7, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Relative enhancement of biomass per taxa
ber<- summaryBy(predMass~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~predMass|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),log="x",
       legend=F,type="p", main="South", ylim=c(0.5,3.5), xlim=c(0,45))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(35,3.7, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~predMass|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),log="x",
       legend=F,type="p", main = "North", ylim=c(0.5,3.5),yaxt='n', xlim=c(0,90))
mtext(text="Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(70,3.7, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Relative enhancement of RGR per taxa
ber<- summaryBy(dydt~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","dydtWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, dydtWarm/dydt)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,1.8), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(50,1.8, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,1.8), xlim=c(0,67),yaxt='n')
mtext(text="Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(50,1.8, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(40,40);
split.screen( figs = c( 2, 1 ) )
split.screen( figs = c( 1, 2 ), screen = 2 )
palett <- c("red","black","blue","green","orange","cyan","grey","yellow")

screen(1)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(predMass.mean~Time|combotrt,data=g.trt,col=palett,
       legendwhere="topleft", cex=0.5, pch=3, xaxt='n', ylab="", type="o",xlab="",
       panel.first=adderrorbars(x=g.trt$Time,y=g.trt$predMass.mean,
                                SE=g.trt$predMass.standard.error,direction="updown",
                                col=c("red","blue","orange","grey","black","green","cyan","yellow")))
mtext(side=2, text="Total biomass (g)",  line=3)
mtext(side=1, text="Time (Days)",  line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)

screen(3)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p",  ylim=c(0.8,2), xlim=c(0,67), ylab="",)
mtext(text=expression(M[W]:M[H]),side=2,line=3)
mtext(side=3, text="South")
abline(h=1)

screen(4)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p",  ylim=c(0.8,2), xlim=c(0,67), yaxt='n')
mtext(side=3, text="North")
abline(h=1)
legend(45,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#Mass, RGR and AGR over Time per species
g.trt <- summaryBy(dydt+predMass+AGR~Species+Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Species+Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))

palettS <- c("black","blue","green","cyan","purple"); palettN <- c("black","orange","tan4","red","magenta")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|Species,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|Species,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|Species,data=sg.trt,col=palettS,,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(dydt.mean~Time|Species,data=ng.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))
plotBy(AGR.mean~Time|Species,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))


plotBy(AGR.mean~Time|Species,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|Species,data=srate.trt,col=palettS,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|Species,data=nrate.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#Mass, RGR and AGR over Time per taxa
g.trt <- summaryBy(dydt+predMass+AGR~Taxa+Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Taxa+Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))

palettS <- c("black","black","blue","blue","blue","green","cyan","purple"); palettN <- c("black","black","black","orange","orange","orange","tan4","red","magenta")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|Taxa,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|Taxa,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|Taxa,data=sg.trt,col=palettS,,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(dydt.mean~Time|Taxa,data=ng.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))
plotBy(AGR.mean~Time|Taxa,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))


plotBy(AGR.mean~Time|Taxa,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|Taxa,data=srate.trt,col=palettS,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|Taxa,data=nrate.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)


# ltys = c("22", "44", "13", "1343", "73", "2262",
#          "12223242", "F282", "F4448444", "224282F2", "F1")

windows(8,8);par(mfrow=c(1,1), mar=c(0.5,3,0.5,0.5), oma=c(3,1,1,1))
plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0,70),lty="22",lwd=3,
       panel.first=adderrorbars(x=g.trt$Time,y=g.trt$predMass.mean,
                                SE=g.trt$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="22",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="44",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="44",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="13",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="13",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="solid",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="solid",lwd=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)
mtext(text="Time (Days)", side=1, line=2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
legend("topleft", legend=c("Tropical Narrow","Tropical Narrow Warmed","Tropical Wide","Tropical Wide Warmed","Temperate Narrow","Temperate Narrow Warmed",
                           "Temperate Wide","Temperate Wide Warmed"), 
       col=c("black","red","black","red",alpha("black",0.6),alpha("red",0.6),alpha("black",0.6),alpha("red",0.6)),
       lty=c("22","22","44","44","13","13","solid","solid"), lwd=3, cex=0.8)

#Size over time

g.trt <- summaryBy(d2h+TotMass~Time+Treatment+Location+Range,data=dat3,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(18,12);par(mfrow=c(1,2), mar=c(2,0,2,0),oma=c(2,4,2,2))
  plotBy(d2h.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
         legend=F, xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3,
         panel.first=adderrorbars(x=g.trt$Time,y=g.trt$d2h.mean,
                                  SE=g.trt$d2h.standard.error,direction="updown",
                                  col="black",0))
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = FALSE, tck = 0.01)
  mtext(side=2, text=expression(paste('d'^'2','h',sep=' ')), line=2.5)
  mtext(text="Time (Days)", side=1, line=2)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = T, tck = 0.01)
  legend("topleft", ncol=2, legend=c("Tropical Narrow","Tropical Narrow Warmed","Tropical Wide","Tropical Wide Warmed","Temperate Narrow","Temperate Narrow Warmed",
                                     "Temperate Wide","Temperate Wide Warmed"), 
         col=c("black","red","black","red",alpha("black",0.6),alpha("red",0.6),alpha("black",0.6),alpha("red",0.6)),
         lty=c("22","22","44","44","13","13","solid","solid"), lwd=3, cex=0.8)
  
  #Size over mass
  plotBy(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",yaxt='n',
         legend=F, xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3,
         panel.first=adderrorbars(x=g.trt$TotMass.mean,y=g.trt$d2h.mean,
                                  SE=g.trt$d2h.standard.error,direction="updown",
                                  col="black",0))
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = FALSE, tck = 0.01)
  mtext(text="TotMass (g)", side=1, line=2)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = T, tck = 0.01)
#par(xpd=T)

  