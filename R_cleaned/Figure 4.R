#Figure 4: Relative growth rate
#source("R_cleaned/3. Create_datasets.R")
#source("R_cleaned/gamfits_14g.R") # first

#changed all Time <40 to Time<32
g.trt <- summaryBy(dydt~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error,NROW))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

CI<-1.96 #90% CI #1.96 95% CI

windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,8,6,6))
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

with(subset(NnH,Time<33),
     plot(dydt.mean~Time,type="l",las=1,yaxs="i",xaxs="i",
       ylim=c(0.01,0.13),lty=2,lwd=2,cex.lab=2, xlim=c(0.01,34),axes=FALSE,
       ylab=expression(Total~mass~(g)),xlab=""))
with(subset(NnW,Time <33),
     lines(dydt.mean~Time, col="red",xaxt='n', ylab="", type="l",ylim=c(0.01,0.2),lty=2,lwd=2))
with(subset(NnH, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="black", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(NnW, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="red", pch=19,cex=2,
       panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(NnH,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("black",0.4), border = NA))
with(subset(NnW,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("red",0.4), border = NA))
mtext(text="Narrow", side=3, line=0.5, cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","a", bty="n", cex=1.5)

with(subset(NwH,Time<33),
     plot(dydt.mean~Time,type="l",las=1,yaxs="i",xaxs="i",
          ylim=c(0.01,0.13),lty=1,lwd=2,cex.lab=2, xlim=c(0.01,34),axes=FALSE,
          ylab=expression(Total~mass~(g)),xlab=""))
with(subset(NwW,Time <33),
     lines(dydt.mean~Time, col="red",xaxt='n', ylab="", type="l",ylim=c(0.01,0.2),lty=1,lwd=2))
with(subset(NwH, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="black", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(NwW, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="red", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(NwH,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("black",0.4), border = NA))
with(subset(NwW,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("red",0.4), border = NA))
mtext(text="Wide", side=3, line=0.5, cex=1.2)
legend(4,70, legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n")
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","b", bty="n", cex=1.5)


with(subset(SnH,Time<33),
     plot(dydt.mean~Time,type="l",las=1,yaxs="i",xaxs="i",
          ylim=c(0.01,0.13),lty=2,lwd=2,cex.lab=2, xlim=c(0.01,34),axes=FALSE,
          ylab=expression(Total~mass~(g)),xlab=""))
with(subset(SnW,Time <33),
     lines(dydt.mean~Time, col="red",xaxt='n', ylab="", type="l",ylim=c(0.01,0.2),lty=2,lwd=2))
with(subset(SnH, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="black", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(SnW, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="red", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(SnH,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("black",0.4), border = NA))
with(subset(SnW,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("red",0.4), border = NA))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","c", bty="n", cex=1.5)
legend("bottomleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",1),alpha("black",0.6)),
       bty="n",cex=1.2)


with(subset(SwH,Time<33),
     plot(dydt.mean~Time,type="l",las=1,yaxs="i",xaxs="i",
          ylim=c(0.01,0.13),lty=1,lwd=2,cex.lab=2, xlim=c(0.01,34),axes=FALSE,
          ylab=expression(Total~mass~(g)),xlab=""))
with(subset(SwW,Time <33),
     lines(dydt.mean~Time, col="red",xaxt='n', ylab="", type="l",ylim=c(0.01,0.2),lty=1,lwd=2))
with(subset(SwH, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="black", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(SwW, Time %in% unique(rate$Time)),
     points(dydt.mean~Time, col="red", pch=19,cex=2,
            panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
with(subset(SwH,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("black",0.4), border = NA))
with(subset(SwW,Time <33),
     polygon(x = c(Time, rev(Time)), y = c(high,rev(low)),col = alpha("red",0.4), border = NA))
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","d", bty="n", cex=1.5)

mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=3.5, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)

text(38,y=0.18,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
text(38,y=0.055,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)

###############################################################################
# Need to investigate:
# #where is ECAM warmed? - it is there with original gamfit, not there if only fitting GAM to <15 g.
# #Why are some completely flat?

#provenance specific version
dat<- summaryBy(dydt+Range~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error, length))
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
  dat2 <- subset(dat,Taxa==as.character(combos[i])&Time<33)
  with(subset(dat2,Treatment=="Home"),
       plot(dydt.mean~Time,col="black",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0.01,34),ylim=c(0, 0.19),axes=FALSE,xlab="Time",ylab="Mass"))
  with(subset(dat2,Treatment=="Home"),
       polygon(x = c(subset(dat2,Treatment=="Home")$Time,
                     rev(subset(dat2,Treatment=="Home")$Time)),
               y = c(subset(dat2,Treatment=="Home")$high,
                     rev(subset(dat2,Treatment=="Home")$low)),
               col = alpha("black",0.4), border = NA))
  with(subset(dat2,Treatment=="Home"& Time %in% unique(rate$Time)),
       points(dydt.mean~Time, col="black", pch=19,cex=2,
              panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
  par(new=T)
  with(subset(dat2,Treatment=="Warmed"),
       plot(dydt.mean~Time,col="red",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0.01,34),ylim=c(0,0.19),axes=FALSE,xlab="Time",ylab="Mass"))
  with(subset(dat2,Treatment=="Warmed"& Time %in% unique(rate$Time)),
       points(dydt.mean~Time, col="red", pch=19,cex=2,
              panel.first=adderrorbars(x=Time,y=dydt.mean,SE=dydt.standard.error, direction="updown",col="black",las=1)))
  
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

legend(x=34,y=0.15,legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),cex=1.4,
       xpd=NA,fill=c(alpha("red",0.4),alpha("black",0.4)), border="black",
       bty="n")
text(77,y=0.95,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
text(77,y=0.3,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)

# ###################################################################################
# ###    OVER MASS
# 
# #Figure 5: Relative growth rate
# g.trt <- summaryBy(dydt+predMass~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
# g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
# g.trt.S<- subset(g.trt, Location == "S")
# g.trt.N<- subset(g.trt, Location == "N")
# 
# CI<-1.645 #90% CI #1.96 95% CI
# 
# windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,6,6,6))
# NnH<-subset(g.trt, combotrt=="N_narrow_Home")
# NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
# NwH<-subset(g.trt, combotrt=="N_wide_Home")
# NwW<- subset(g.trt, combotrt=="N_wide_Warmed")
# SnH<-subset(g.trt, combotrt=="S_narrow_Home")
# SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
# SwH<-subset(g.trt, combotrt=="S_wide_Home")
# SwW<- subset(g.trt, combotrt=="S_wide_Warmed")
# 
# NnH$high<- with(NnH,dydt.mean+dydt.standard.error*CI )
# NnH$low<- with(NnH,dydt.mean-dydt.standard.error*CI )
# NnW$high<- with(NnW,dydt.mean+dydt.standard.error*CI )
# NnW$low<- with(NnW,dydt.mean-dydt.standard.error*CI )
# NwH$high<- with(NwH,dydt.mean+dydt.standard.error*CI )
# NwH$low<- with(NwH,dydt.mean-dydt.standard.error*CI )
# NwW$high<- with(NwW,dydt.mean+dydt.standard.error*CI )
# NwW$low<- with(NwW,dydt.mean-dydt.standard.error*CI )
# SnH$high<- with(SnH,dydt.mean+dydt.standard.error*CI )
# SnH$low<- with(SnH,dydt.mean-dydt.standard.error*CI )
# SnW$high<- with(SnW,dydt.mean+dydt.standard.error*CI )
# SnW$low<- with(SnW,dydt.mean-dydt.standard.error*CI )
# SwH$high<- with(SwH,dydt.mean+dydt.standard.error*CI )
# SwH$low<- with(SwH,dydt.mean-dydt.standard.error*CI )
# SwW$high<- with(SwW,dydt.mean+dydt.standard.error*CI )
# SwW$low<- with(SwW,dydt.mean-dydt.standard.error*CI )
# 
# xlims<-c(0.1,70)
# plotBy(dydt.mean~predMass.mean,data=NnH,legend=FALSE,type="l",las=1,yaxs="i",xaxs="i",log="x",
#        ylim=c(0.02,0.17),lty=2,lwd=2,cex.lab=2, xlim=xlims,axes=FALSE,ylab=expression(Total~mass~(g)),xlab="")
# lines(dydt.mean~predMass.mean, data=NnW,col="red",xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=2,lwd=2)
# 
# polygon(x = c(NnH$predMass.mean, rev(NnH$predMass.mean)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NnW$predMass.mean, rev(NnW$predMass.mean)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
# 
# mtext(text="Narrow", side=3, line=0.5, cex=1.2)
# legend(4,70, legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n")
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","a", bty="n", cex=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",1),alpha("black",0.6)),
#        bty="n",cex=1.2)
# rect(0.8, 0.02, 1, 0.17 ,col = alpha("black",0.2), border = NA)
# #comparison at a very small mass
# rect(19, 0.02, 25, 0.17 ,col = alpha("black",0.2), border = NA)
# 
# plotBy(dydt.mean~predMass.mean, data=NwH,col="black",legend=FALSE,yaxs="i",xaxs="i",log="x",
#        xaxt='n', yaxt='n',ylab="", type="l",ylim=c(0.02,0.17),xlim=xlims,lty=1,lwd=2)
# lines(dydt.mean~predMass.mean, data=NwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
# polygon(x = c(NwH$predMass.mean, rev(NwH$predMass.mean)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NwW$predMass.mean, rev(NwW$predMass.mean)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text="Wide", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.2)
# rect(0.8, 0.02, 1, 0.17 ,col = alpha("black",0.2), border = NA)
# rect(19, 0.02, 25, 0.17 ,col = alpha("black",0.2), border = NA)
# 
# plotBy(dydt.mean~predMass.mean,data=SnH,legend=FALSE,type="l",las=1,yaxs="i",xaxs="i",log="x",
#        ylim=c(0.02,0.17),lty=2,lwd=2,cex.lab=2, xlim=xlims,axes=FALSE,xlab="")
# lines(dydt.mean~predMass.mean, data=SnW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=2,lwd=2)
# polygon(x = c(SnH$predMass.mean, rev(SnH$predMass.mean)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SnW$predMass.mean, rev(SnW$predMass.mean)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
# #mtext(text="Temperate",side=3, line=-2,at=14,cex=1.5, outer=FALSE)
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","c", bty="n", cex=1.2)
# rect(0.8, 0.02, 1, 0.17 ,col = alpha("black",0.2), border = NA)
# rect(19, 0.02, 25, 0.17 ,col = alpha("black",0.2), border = NA)
# 
# plotBy(dydt.mean~predMass.mean, data=SwH,col="black",legend=FALSE, yaxt='n',yaxs="i",xaxs="i",log="x",
#        xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),xlim=xlims,lty=1,lwd=2)
# lines(dydt.mean~predMass.mean, data=SwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0.02,0.17),lty=1,lwd=2)
# polygon(x = c(SwH$predMass.mean, rev(SwH$predMass.mean)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SwW$predMass.mean, rev(SwW$predMass.mean)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
# 
# magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","d", bty="n", cex=1.2)
# mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=3, cex=1.2)
# mtext(text="Log Total Mass (g)", side=1, line=3, cex=1.2, adj=-0.3)
# 
# text(125,y=0.22,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
# text(130,y=0.07,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
# 
# rect(0.8, 0.02, 1, 0.17 ,col = alpha("black",0.2), border = NA)
# rect(19, 0.02, 25, 0.17 ,col = alpha("black",0.2), border = NA)
# 
#   ###############################################################################
# 
# 
# #provenance specific version
# dat<- summaryBy(dydt+Range+predMass~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error))
# dat$high<-with(dat,dydt.mean+dydt.standard.error*CI)
# dat$low<-with(dat,dydt.mean-dydt.standard.error*CI)
# combostrop<- c("BRA",NA,"CTER","DTER","PEL",NA,"ETER","DCAM","PLAT",NA,
#                "ECAM","FCAM")
# combostemp <- c("BOT",NA,"ATER","BTER","LONG",NA,"ACAM","BCAM","SMIT",NA,"CCAM")
# combos<-c(combostrop,rep(NA,4),combostemp)
# windows(8.27,11.69)
# par(mfrow=c(8,8),mar=c(0,0,0,0),oma=c(6,6,3,3))
# layout(matrix(c(1:28), nrow=7, ncol=4,byrow=T),
#        heights=c(1,1,1,0.3,1,1,1),
#        widths=c(1,0.3,1,1,1))
# 
# 
# for (i in 1:length(combos)){
#   dat2 <- subset(dat,Taxa==as.character(combos[i]))
#   with(subset(dat2,Treatment=="Home"),
#        plot(dydt.mean~predMass.mean,col="black",legend=FALSE,type="p",lty=ifelse(Range.mean == 1,2,1),log="x",
#             xlim=c(0,70),ylim=c(0, 0.19),axes=FALSE,xlab="predMass.mean",ylab="Mass"))  
#   with(subset(dat2,Treatment=="Home"),
#        polygon(x = c(subset(dat2,Treatment=="Home")$predMass.mean, 
#                      rev(subset(dat2,Treatment=="Home")$predMass.mean)), 
#                y = c(subset(dat2,Treatment=="Home")$high,
#                      rev(subset(dat2,Treatment=="Home")$low)),
#                col = alpha("black",0.4), border = NA))
#   par(new=T)
#   with(subset(dat2,Treatment=="Warmed"),
#        plot(dydt.mean~predMass.mean,col="red",legend=FALSE,type="p",lty=ifelse(Range.mean == 1,2,1),log="x",
#             xlim=c(0,70),ylim=c(0,0.19),axes=FALSE,xlab="predMass.mean",ylab="Mass"))  
#   with(subset(dat2,Treatment=="Warmed"),
#        polygon(x = c(subset(dat2,Treatment=="Warmed")$predMass.mean, 
#                      rev(subset(dat2,Treatment=="Warmed")$predMass.mean)), 
#                y = c(subset(dat2,Treatment=="Warmed")$high,
#                      rev(subset(dat2,Treatment=="Warmed")$low)),
#                col = alpha("red",0.4), border = NA))
#   #first plot
#   ifelse(dat2$Taxa %in% c("BRA","PEL","BOT","LONG"),
#          magaxis(side=c(1,2),labels=c(0,1),frame.plot=T,las=1,cex.axis=1.2),
#          ifelse(dat2$Taxa %in% c("CCAM","BCAM","ECAM","FCAM"),
#                 magaxis(side=c(1,2),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2),
#                 ifelse(dat2$Taxa %in% c("SMIT","PLAT"),
#                        magaxis(side=c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2),
#                        magaxis(side=c(1,2),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2))))
#   
#   legend("topleft",legend=dat2$Taxa[1],bty="n",cex=1.3)
#   rect(0.8, 0.02, 1, 0.17 ,col = alpha("black",0.4), border = NA)
#   rect(15, 0.02, 20, 0.17 ,col = alpha("black",0.4), border = NA)
#   
# }
# mtext(expression(RGR~(g~g^-1~day^-1)),side=2,line=3,outer=T,cex=1.5)
# mtext(expression(Total~mass~(g)),side=1,line=3,outer=T,cex=1.5)
# mtext("Narrow",side=3,line=1,outer=T,cex=1, adj=0.08)
# mtext("Wide",side=3,line=1,outer=T,cex=1, adj=0.72)
# text(160,y=0.9,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
# text(160,y=0.3,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)
# 
# legend(x=70,y=0.15,legend=c("Warmed","Home"),cex=1.4,
#        xpd=NA,fill=c(alpha("red",0.4),alpha("black",0.4)), border="black",
#        bty="n")
# 
# ###############################################################################
# 
# ###################################################################################
# ###  RGR points at measurements points OVER MASS, 40 days
# 
# #Figure 5: Relative growth rate
# g.trt <- summaryBy(dydt+predMass~Time+Treatment+Location+Range,data=subset(rate,Time<40),FUN=c(mean,standard.error))
# g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
# g.trt.S<- subset(g.trt, Location == "S")
# g.trt.N<- subset(g.trt, Location == "N")
# 
# CI<-1.645 #90% CI #1.96 95% CI
# 
# 
# NnH<-subset(g.trt, combotrt=="N_narrow_Home")
# NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
# NwH<-subset(g.trt, combotrt=="N_wide_Home")
# NwW<- subset(g.trt, combotrt=="N_wide_Warmed")
# SnH<-subset(g.trt, combotrt=="S_narrow_Home")
# SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
# SwH<-subset(g.trt, combotrt=="S_wide_Home")
# SwW<- subset(g.trt, combotrt=="S_wide_Warmed")
# 
# NnH$high<- with(NnH,dydt.mean+dydt.standard.error*CI )
# NnH$low<- with(NnH,dydt.mean-dydt.standard.error*CI )
# NnW$high<- with(NnW,dydt.mean+dydt.standard.error*CI )
# NnW$low<- with(NnW,dydt.mean-dydt.standard.error*CI )
# NwH$high<- with(NwH,dydt.mean+dydt.standard.error*CI )
# NwH$low<- with(NwH,dydt.mean-dydt.standard.error*CI )
# NwW$high<- with(NwW,dydt.mean+dydt.standard.error*CI )
# NwW$low<- with(NwW,dydt.mean-dydt.standard.error*CI )
# SnH$high<- with(SnH,dydt.mean+dydt.standard.error*CI )
# SnH$low<- with(SnH,dydt.mean-dydt.standard.error*CI )
# SnW$high<- with(SnW,dydt.mean+dydt.standard.error*CI )
# SnW$low<- with(SnW,dydt.mean-dydt.standard.error*CI )
# SwH$high<- with(SwH,dydt.mean+dydt.standard.error*CI )
# SwH$low<- with(SwH,dydt.mean-dydt.standard.error*CI )
# SwW$high<- with(SwW,dydt.mean+dydt.standard.error*CI )
# SwW$low<- with(SwW,dydt.mean-dydt.standard.error*CI )
# 
# windows(11.69,11.69);par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(6,8,6,6))
# 
# plotBy(dydt.mean~predMass.mean,data=NnH,legend=FALSE,
#        yaxs="i",xaxs="i", ylim=c(0.05,0.15),
#        cex.lab=2, axes=FALSE,xlim=c(0.1,40),
#        col="black", pch=19, log="x")
# points(dydt.mean~predMass.mean, data=NnW,col="red", pch=19)
# 
# points(dydt.mean~predMass.mean, data=NwH,col="black", pch=15)
# points(dydt.mean~predMass.mean, data=NwW,col="red", pch=15)
# 
# points(dydt.mean~predMass.mean,data=SnH, col="black", pch=1)
# points(dydt.mean~predMass.mean, data=SnW,col="red", pch=1)
# 
# points(dydt.mean~predMass.mean, data=SwH,col="black", pch=0)
# points(dydt.mean~predMass.mean, data=SwW,col="red", pch=0)
# 
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=3.5, cex=1.2)
# mtext(text="Total Mass (g)", side=1, line=3, cex=1.2)
# 
# legend("topright", pch=c(19,19,15,15,1,1,0,0), col=c("black","red"), 
#        legend=c("Tropical-Narrow-Home","Tropical-Narrow-Warmed",
#                 "Tropical-Wide-Home","Tropical-Wide-Warmed",
#                 "Temperate-Narrow-Home","Temperate-Narrow-Warmed",
#                 "Temperate-Wide-Home","Temperate-Wide-Warmed"))
# 
# 
# 
# 
# 
# ###############################################################################
# 
