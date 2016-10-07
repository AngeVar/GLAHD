#Figure 4: Relative growth rate
g.trt <- summaryBy(dydt~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

CI<-1.645 #90% CI #1.96 95% CI

windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,6,6,6))
NnH<-subset(g.trt, combotrt=="N_narrow_Home")
NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
NwH<-subset(g.trt, combotrt=="N_wide_Home")
NwW<- subset(g.trt, combotrt=="N_wide_Warmed")
SnH<-subset(g.trt, combotrt=="S_narrow_Home")
SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
SwH<-subset(g.trt, combotrt=="S_wide_Home")
SwW<- subset(g.trt, combotrt=="S_wide_Warmed")

NnH$high<- with(NnH,dydt.mean+dydt.standard.error*CI )
NnH$low<- with(NnH,dydt.mean-dydt.standard.error*CI )
NnW$high<- with(NnW,dydt.mean+dydt.standard.error*CI )
NnW$low<- with(NnW,dydt.mean-dydt.standard.error*CI )
NwH$high<- with(NwH,dydt.mean+dydt.standard.error*CI )
NwH$low<- with(NwH,dydt.mean-dydt.standard.error*CI )
NwW$high<- with(NwW,dydt.mean+dydt.standard.error*CI )
NwW$low<- with(NwW,dydt.mean-dydt.standard.error*CI )
SnH$high<- with(SnH,dydt.mean+dydt.standard.error*CI )
SnH$low<- with(SnH,dydt.mean-dydt.standard.error*CI )
SnW$high<- with(SnW,dydt.mean+dydt.standard.error*CI )
SnW$low<- with(SnW,dydt.mean-dydt.standard.error*CI )
SwH$high<- with(SwH,dydt.mean+dydt.standard.error*CI )
SwH$low<- with(SwH,dydt.mean-dydt.standard.error*CI )
SwW$high<- with(SwW,dydt.mean+dydt.standard.error*CI )
SwW$low<- with(SwW,dydt.mean-dydt.standard.error*CI )

plotBy(dydt.mean~Time,data=NnH,legend=FALSE,type="l",las=1,yaxs="i",xaxs="i",
       ylim=c(0.02,0.17),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,
       ylab=expression(Total~mass~(g)),
       xlab="")
lines(dydt.mean~Time, data=NnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=2,lwd=2)
polygon(x = c(NnH$Time, rev(NnH$Time)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NnW$Time, rev(NnW$Time)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
mtext(text="Narrow", side=3, line=0.5, cex=1.2)
legend(4,70, legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n")
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","a", bty="n", cex=1.2)
#legend(12,0.17, legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n")
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",1),alpha("black",0.6)),
       bty="n",cex=1.2)

plotBy(dydt.mean~Time, data=NwH,col="black",legend=FALSE,yaxs="i",xaxs="i",
       xaxt='n', yaxt='n',ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
lines(dydt.mean~Time, data=NwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
polygon(x = c(NwH$Time, rev(NwH$Time)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(NwW$Time, rev(NwW$Time)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Wide", side=3, line=0.5, cex=1.2)
legend("topright","b", bty="n", cex=1.2)


plotBy(dydt.mean~Time,data=SnH,legend=FALSE,type="l",las=1,yaxs="i",xaxs="i",
       ylim=c(0.02,0.17),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,xlab="")
lines(dydt.mean~Time, data=SnW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=2,lwd=2)
polygon(x = c(SnH$Time, rev(SnH$Time)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SnW$Time, rev(SnW$Time)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
#mtext(text="Temperate",side=3, line=-2,at=14,cex=1.5, outer=FALSE)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","c", bty="n", cex=1.2)

plotBy(dydt.mean~Time, data=SwH,col="black",legend=FALSE, yaxt='n',yaxs="i",xaxs="i",
       xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
lines(dydt.mean~Time, data=SwW,col="red",
      xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
polygon(x = c(SwH$Time, rev(SwH$Time)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(SwW$Time, rev(SwW$Time)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","d", bty="n", cex=1.2)
mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=3, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)

text(68,y=0.22,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
text(68,y=0.07,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)

###############################################################################


#provenance specific version
dat<- summaryBy(dydt+Range~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error))
dat$high<-with(dat,dydt.mean+dydt.standard.error*CI)
dat$low<-with(dat,dydt.mean-dydt.standard.error*CI)
combostrop<- c("BRA",NA,"CTER","DTER","PEL",NA,"ETER","DCAM","PLAT",NA,
               "ECAM","FCAM")
combostemp <- c("BOT",NA,"ATER","BTER","LONG",NA,"ACAM","BCAM","SMIT",NA,"CCAM")
combos<-c(combostrop,rep(NA,4),combostemp)
windows(8.27,11.69)
par(mfrow=c(8,8),mar=c(0,0,0,0),oma=c(6,6,3,3))
layout(matrix(c(1:28), nrow=7, ncol=4,byrow=T),
       heights=c(1,1,1,0.3,1,1,1),
       widths=c(1,0.3,1,1,1))


for (i in 1:length(combos)){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  with(subset(dat2,Treatment=="Home"),
       plot(dydt.mean~Time,col="black",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0, 0.19),axes=FALSE,xlab="Time",ylab="Mass"))  
  with(subset(dat2,Treatment=="Home"),
       polygon(x = c(subset(dat2,Treatment=="Home")$Time, 
                     rev(subset(dat2,Treatment=="Home")$Time)), 
               y = c(subset(dat2,Treatment=="Home")$high,
                     rev(subset(dat2,Treatment=="Home")$low)),
               col = alpha("black",0.4), border = NA))
  par(new=T)
  with(subset(dat2,Treatment=="Warmed"),
       plot(dydt.mean~Time,col="red",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0,0.19),axes=FALSE,xlab="Time",ylab="Mass"))  
  with(subset(dat2,Treatment=="Warmed"),
       polygon(x = c(subset(dat2,Treatment=="Warmed")$Time, 
                     rev(subset(dat2,Treatment=="Warmed")$Time)), 
               y = c(subset(dat2,Treatment=="Warmed")$high,
                     rev(subset(dat2,Treatment=="Warmed")$low)),
               col = alpha("red",0.4), border = NA))
  #first plot
  ifelse(dat2$Taxa %in% c("BRA","PEL","BOT","LONG"),
         magaxis(side=c(1,2),labels=c(0,1),frame.plot=T,las=1,cex.axis=1.2),
         ifelse(dat2$Taxa %in% c("CCAM","BCAM","ECAM","FCAM"),
                magaxis(side=c(1,2),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2),
                ifelse(dat2$Taxa %in% c("SMIT","PLAT"),
                       magaxis(side=c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2),
                       magaxis(side=c(1,2),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2))))
  
  legend("topleft",legend=dat2$Taxa[1],bty="n",cex=1.3)
}
mtext(expression(RGR~(g~g^-1~day^-1)),side=2,line=3,outer=T,cex=1.5)
mtext(expression(Time~(days)),side=1,line=3,outer=T,cex=1.5)
mtext("Narrow",side=3,line=1,outer=T,cex=1, adj=0.08)
mtext("Wide",side=3,line=1,outer=T,cex=1, adj=0.72)
text(155,y=0.9,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
text(155,y=0.3,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)

legend(x=70,y=0.15,legend=c("Warmed","Home"),cex=1.4,
       xpd=NA,fill=c(alpha("red",0.4),alpha("black",0.4)), border="black",
       bty="n")

