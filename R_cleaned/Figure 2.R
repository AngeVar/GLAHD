#Figure 2
source("R/create_datasets.R")

#Figure 2: Mass over time
g.trt <- summaryBy(predMass~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

CI<-1.645 #90% CI #1.96 95% CI

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))
NnH<-subset(g.trt, combotrt=="N_narrow_Home")
NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
NwH<-subset(g.trt, combotrt=="N_wide_Home")
NwW<- subset(g.trt, combotrt=="N_wide_Warmed")

plotBy(predMass.mean~Time,data=NnH,legend=F,type="l",las=1,#yaxs="i",xaxs="i",
       ylim=c(1,70),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,
       ylab=expression(Total~mass~(g)),
       xlab="")
lines(predMass.mean~Time, data=NnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=NwH,col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=NwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)

NnH$high<- with(NnH,predMass.mean+predMass.standard.error*CI )
NnH$low<- with(NnH,predMass.mean-predMass.standard.error*CI )
NnW$high<- with(NnW,predMass.mean+predMass.standard.error*CI )
NnW$low<- with(NnW,predMass.mean-predMass.standard.error*CI )
NwH$high<- with(NwH,predMass.mean+predMass.standard.error*CI )
NwH$low<- with(NwH,predMass.mean-predMass.standard.error*CI )
NwW$high<- with(NwW,predMass.mean+predMass.standard.error*CI )
NwW$low<- with(NwW,predMass.mean-predMass.standard.error*CI )

polygon(x = c(NnH$Time, rev(NnH$Time)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NnW$Time, rev(NnW$Time)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
polygon(x = c(NwH$Time, rev(NwH$Time)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NwW$Time, rev(NwW$Time)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)

legend(0,60, legend=c("Narrow Home","Narrow +3.5","Wide Home","Wide +3.5"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical",3,line=-2,at=12.5, cex=1.5, outer=F)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
        
SnH<-subset(g.trt, combotrt=="S_narrow_Home")
SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
SwH<-subset(g.trt, combotrt=="S_wide_Home")
SwW<- subset(g.trt, combotrt=="S_wide_Warmed")

plotBy(predMass.mean~Time,data=SnH,legend=F,type="l",las=1,#yaxs="i",xaxs="i",
       ylim=c(1,70),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,
       ylab=expression(Total~mass~(g)),
       xlab="")
lines(predMass.mean~Time, data=SnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=SwH,col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=SwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)

SnH$high<- with(SnH,predMass.mean+predMass.standard.error*CI )
SnH$low<- with(SnH,predMass.mean-predMass.standard.error*CI )
SnW$high<- with(SnW,predMass.mean+predMass.standard.error*CI )
SnW$low<- with(SnW,predMass.mean-predMass.standard.error*CI )
SwH$high<- with(SwH,predMass.mean+predMass.standard.error*CI )
SwH$low<- with(SwH,predMass.mean-predMass.standard.error*CI )
SwW$high<- with(SwW,predMass.mean+predMass.standard.error*CI )
SwW$low<- with(SwW,predMass.mean-predMass.standard.error*CI )

polygon(x = c(SnH$Time, rev(SnH$Time)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SnW$Time, rev(SnW$Time)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
polygon(x = c(SwH$Time, rev(SwH$Time)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SwW$Time, rev(SwW$Time)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
mtext("Temperate",3,line=-2,at=14,cex=1.5, outer=F)
mtext(text="Total biomass (g)", outer=T, side=2, line=1, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2)



###############################################
#2by2 version
windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,6,6,6))

plotBy(predMass.mean~Time,data=NnH,legend=F,type="l",las=1,#yaxs="i",xaxs="i",
       ylim=c(0,82),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,
       ylab=expression(Total~mass~(g)),
       xlab="")
lines(predMass.mean~Time, data=NnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,82),lty=2,lwd=2)
polygon(x = c(NnH$Time, rev(NnH$Time)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NnW$Time, rev(NnW$Time)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
mtext(text="Narrow", side=3, line=0.5, cex=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n")
                      
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","a", bty="n", cex=1.2)

plotBy(predMass.mean~Time, data=NwH,col="black",legend=F,#yaxs="i",xaxs="i",
      xaxt='n', yaxt='n',ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
lines(predMass.mean~Time, data=NwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
polygon(x = c(NwH$Time, rev(NwH$Time)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NwW$Time, rev(NwW$Time)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Wide", side=3, line=0.5, cex=1.2)
legend("topright","b", bty="n", cex=1.2)

plotBy(predMass.mean~Time,data=SnH,legend=F,type="l",las=1,#yaxs="i",xaxs="i",
       ylim=c(0,82),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,xlab="")
lines(predMass.mean~Time, data=SnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,82),lty=2,lwd=2)
polygon(x = c(SnH$Time, rev(SnH$Time)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SnW$Time, rev(SnW$Time)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","c", bty="n", cex=1.2)

plotBy(predMass.mean~Time, data=SwH,col="black",legend=F, yaxt='n',#yaxs="i",xaxs="i",
      xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
lines(predMass.mean~Time, data=SwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
polygon(x = c(SwH$Time, rev(SwH$Time)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SwW$Time, rev(SwW$Time)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","d", bty="n", cex=1.2)
mtext(text="Total biomass (g)", outer=T, side=2, line=2.5, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)

text(70,y=117,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
text(70,y=22,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
########################################################################

#provenance specific version
dat<- summaryBy(predMass+Range~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error))
dat$high<-with(dat,predMass.mean+predMass.standard.error*CI)
dat$low<-with(dat,predMass.mean-predMass.standard.error*CI)
combostemp <- c("BOT","ATER","BTER","LONG","ACAM","BCAM","SMIT","CCAM")
combostrop<- c("BRA","CTER","DTER","PEL","ETER","DCAM","PLAT",
               "ECAM","FCAM")
combos<-c(combostrop,combostemp)
windows(8.27,11.69)
par(mfrow=c(6,3),mar=c(0,0,0,0),oma=c(6,6,1,1))

for (i in 1:length(combos)){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  
  with(subset(dat2,Treatment=="Home"),
       plot(predMass.mean~Time,col="black",legend=F,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0, 82),axes=F,xlab="Time",ylab="Mass"))  
  with(subset(dat2,Treatment=="Home"),
       polygon(x = c(subset(dat2,Treatment=="Home")$Time, 
                     rev(subset(dat2,Treatment=="Home")$Time)), 
               y = c(subset(dat2,Treatment=="Home")$high,
                     rev(subset(dat2,Treatment=="Home")$low)),
               col = alpha("black",0.4), border = NA))
  par(new=T)
  with(subset(dat2,Treatment=="Warmed"),
       plot(predMass.mean~Time,col="red",legend=F,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0,82),axes=F,xlab="Time",ylab="Mass"))  
  with(subset(dat2,Treatment=="Warmed"),
       polygon(x = c(subset(dat2,Treatment=="Warmed")$Time, 
                     rev(subset(dat2,Treatment=="Warmed")$Time)), 
               y = c(subset(dat2,Treatment=="Warmed")$high,
                     rev(subset(dat2,Treatment=="Warmed")$low)),
               col = alpha("red",0.4), border = NA))
  #first plot
  ifelse(dat2$Taxa %in% c("BRA","PEL","PLAT","BOT","LONG"),
         magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2),
         ifelse(dat2$Taxa %in% c("CCAM","BCAM"),
                magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2),
                ifelse(dat2$Taxa=="SMIT",
                       magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2),
                       magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2))))
  
  legend("topleft",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(Biomass~(g)),side=2,line=3,outer=T,cex=1.5)
mtext(expression(Time~(days)),side=1,line=3,outer=T,cex=1.5)
legend(x=200,y=50,legend=c("Warmed","Home"),pch=15,cex=1.5,xpd=NA,col=c(alpha("red",0.4),alpha("black",0.4)))
