#Figure 2
source("R_cleaned/3. Create_datasets.R")

#Figure 2: Mass over time
###############################################
#2by2 version on log base cut off at last dh measurement time point where masses were <15 g.

#add predmass from "rate" as means and SEs
rate.m<-summaryBy(predMass~Location+Treatment+Range+Time,data=rate, FUN=c(mean, standard.error))

g.trt <- summaryBy(predMass~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

NnH<-subset(g.trt, combotrt=="N_narrow_Home")
NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
NwH<-subset(g.trt, combotrt=="N_wide_Home")
NwW<- subset(g.trt, combotrt=="N_wide_Warmed")

SnH<-subset(g.trt, combotrt=="S_narrow_Home")
SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
SwH<-subset(g.trt, combotrt=="S_wide_Home")
SwW<- subset(g.trt, combotrt=="S_wide_Warmed")

CI<-1.96 #90% CI #1.96 95% CI

NnH$high<- with(NnH,predMass.mean+predMass.standard.error*CI )
NnH$low<- with(NnH,predMass.mean-predMass.standard.error*CI )
NnW$high<- with(NnW,predMass.mean+predMass.standard.error*CI )
NnW$low<- with(NnW,predMass.mean-predMass.standard.error*CI )
NwH$high<- with(NwH,predMass.mean+predMass.standard.error*CI )
NwH$low<- with(NwH,predMass.mean-predMass.standard.error*CI )
NwW$high<- with(NwW,predMass.mean+predMass.standard.error*CI )
NwW$low<- with(NwW,predMass.mean-predMass.standard.error*CI )


SnH$high<- with(SnH,predMass.mean+predMass.standard.error*CI )
SnH$low<- with(SnH,predMass.mean-predMass.standard.error*CI )
SnW$high<- with(SnW,predMass.mean+predMass.standard.error*CI )
SnW$low<- with(SnW,predMass.mean-predMass.standard.error*CI )
SwH$high<- with(SwH,predMass.mean+predMass.standard.error*CI )
SwH$low<- with(SwH,predMass.mean-predMass.standard.error*CI )
SwW$high<- with(SwW,predMass.mean+predMass.standard.error*CI )
SwW$low<- with(SwW,predMass.mean-predMass.standard.error*CI )


windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,8,6,6))

plot(predMass.mean~Time,data=NnH,type="l",las=1,ylim=c(0.2,130),
     lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,log="y",col="darkgrey",
     ylab=expression(Total~mass~(g)),
     xlab="")
lines(predMass.mean~Time,data=subset(NnH,predMass.mean <15),type="l",las=1,ylim=c(0.2,130),
      lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,log="y",
      ylab=expression(Total~mass~(g)),
      xlab="")

points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"), 
       col="darkgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$Time,
                                y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15), 
       col="black", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$Time,
                                y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

lines(predMass.mean~Time, data=NnW,col="lightgrey",
      xaxt='n', ylab="", type="l",lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(NnW, predMass.mean <15),col="red",
      xaxt='n', ylab="", type="l",lty=2,lwd=2)
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"), 
       col="lightgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$Time,
                                y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15), 
       col="red", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$Time,
                                y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$predMass.standard.error,
                                direction="updown",col="black",las=1))


polygon(x = c(subset(NnH, predMass.mean <15)$Time, rev(subset(NnH, predMass.mean <15)$Time)), y = c(subset(NnH, predMass.mean <15)$high,rev(subset(NnH, predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(subset(NnW, predMass.mean <15)$Time, rev(subset(NnW, predMass.mean <15)$Time)), y = c(subset(NnW, predMass.mean <15)$high,rev(subset(NnW, predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
mtext(text="Narrow", side=3, line=0.5, cex=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2, log="y", logpretty=T)
legend("topright","a", bty="n", cex=1.5)

plotBy(predMass.mean~Time, data=NwH,col="darkgrey",legend=FALSE,ylim=c(0.2,130),xlim=c(0,60),
       xaxt='n', yaxt='n',ylab="", type="l",log="y",lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(NwH,predMass.mean <15),col="black",legend=FALSE,ylim=c(0.2,130),xlim=c(0,60),
      xaxt='n', yaxt='n',ylab="", type="l",log="y",lty=1,lwd=2)

points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"), 
       col="darkgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$Time,
                                y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15), 
       col="black", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

lines(predMass.mean~Time, data=NwW,col="lightgrey",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(NwW,predMass.mean <15),col="red",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"), 
       col="lightgray", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$Time,
                                y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15), 
       col="red", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

polygon(x = c(subset(NwH,predMass.mean <15)$Time, rev(subset(NwH,predMass.mean <15)$Time)), y = c(subset(NwH,predMass.mean <15)$high,rev(subset(NwH,predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(subset(NwW,predMass.mean <15)$Time, rev(subset(NwW,predMass.mean <15)$Time)), y = c(subset(NwW,predMass.mean <15)$high,rev(subset(NwW,predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
mtext(text="Wide", side=3, line=0.5, cex=1.2)
legend("topright","b", bty="n", cex=1.5)

plotBy(predMass.mean~Time,data=SnH,legend=FALSE,type="l",las=1,ylim=c(0.2,130),col="darkgrey",
       log="y",lty=2,lwd=2,cex.lab=2, xlim=c(0,60),axes=FALSE,xlab="")
lines(predMass.mean~Time,data=subset(SnH,Time <47),legend=FALSE,type="l",las=1,ylim=c(0.2,130),
      log="y",lty=2,lwd=2,cex.lab=2, xlim=c(0,60),axes=FALSE,xlab="")
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"), 
       col="darkgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$Time,
                                y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15), 
       col="black", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

lines(predMass.mean~Time, data=SnW,col="lightgrey",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(SnW,Time<40),col="red",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=2,lwd=2)
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"), 
       col="lightgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$Time,
                                y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15), 
       col="red", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

polygon(x = c(subset(SnH,Time <47)$Time, rev(subset(SnH,Time <47)$Time)), y = c(subset(SnH,Time <47)$high,rev(subset(SnH,Time <47)$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(subset(SnW,Time <40)$Time, rev(subset(SnW,Time <40)$Time)), y = c(subset(SnW,Time <40)$high,rev(subset(SnW,Time <40)$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
legend("topright","c", bty="n", cex=1.5)

plotBy(predMass.mean~Time, data=SwH,col="darkgrey",legend=FALSE, yaxt='n',ylim=c(0.2,130),#yaxs="i",xaxs="i",
       xaxt='n', ylab="", type="l",log="y",lty=1,lwd=2,xlim=c(0,60))
lines(predMass.mean~Time, data=subset(SwH,Time <47),col="black",legend=FALSE, yaxt='n',ylim=c(0.2,130),#yaxs="i",xaxs="i",
      xaxt='n', ylab="", type="l",log="y",lty=1,lwd=2,xlim=c(0,60))
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"), 
       col="darkgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$Time,
                                y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15), 
       col="black", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

lines(predMass.mean~Time, data=SwW, col="lightgrey",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(SwW,Time <33),col="red",
      xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"), 
       col="lightgrey", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$Time,
                                y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
                                direction="updown",col="black",las=1))
points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15), 
       col="red", pch=19,cex=2,
       panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$Time,
                                y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
                                SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
                                direction="updown",col="black",las=1))

polygon(x = c(subset(SwH,Time <47)$Time, rev(subset(SwH,Time <47)$Time)), y = c(subset(SwH,Time <47)$high,rev(subset(SwH,Time <47)$low)),col = alpha("black",0.4), border = NA)
polygon(x = c(subset(SwW,Time <33)$Time, rev(subset(SwW,Time <33)$Time)), y = c(subset(SwW,Time <33)$high,rev(subset(SwW,Time <33)$low)),col = alpha("red",0.4), border = NA)
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
legend("topright","d", bty="n", cex=1.5)
mtext(text="Total biomass (g)", outer=T, side=2, line=3.5, cex=1.2)
mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)

text(70,y=2500,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
text(70,y=2,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
########################################################################
###############################################
# #2by2 version on log base
# 
# #add predmass from "rate" as means and SEs
# rate.m<-summaryBy(predMass~Location+Treatment+Range+Time,data=rate, FUN=c(mean, standard.error))
# 
# g.trt <- summaryBy(predMass~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
# g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
# g.trt.S<- subset(g.trt, Location == "S")
# g.trt.N<- subset(g.trt, Location == "N")
# 
# NnH<-subset(g.trt, combotrt=="N_narrow_Home")
# NnW<-subset(g.trt, combotrt=="N_narrow_Warmed")
# NwH<-subset(g.trt, combotrt=="N_wide_Home")
# NwW<- subset(g.trt, combotrt=="N_wide_Warmed")
# 
# SnH<-subset(g.trt, combotrt=="S_narrow_Home")
# SnW<-subset(g.trt, combotrt=="S_narrow_Warmed")
# SwH<-subset(g.trt, combotrt=="S_wide_Home")
# SwW<- subset(g.trt, combotrt=="S_wide_Warmed")
# 
# CI<-1.96 #90% CI #1.96 95% CI
# 
# NnH$high<- with(NnH,predMass.mean+predMass.standard.error*CI )
# NnH$low<- with(NnH,predMass.mean-predMass.standard.error*CI )
# NnW$high<- with(NnW,predMass.mean+predMass.standard.error*CI )
# NnW$low<- with(NnW,predMass.mean-predMass.standard.error*CI )
# NwH$high<- with(NwH,predMass.mean+predMass.standard.error*CI )
# NwH$low<- with(NwH,predMass.mean-predMass.standard.error*CI )
# NwW$high<- with(NwW,predMass.mean+predMass.standard.error*CI )
# NwW$low<- with(NwW,predMass.mean-predMass.standard.error*CI )
# 
# 
# SnH$high<- with(SnH,predMass.mean+predMass.standard.error*CI )
# SnH$low<- with(SnH,predMass.mean-predMass.standard.error*CI )
# SnW$high<- with(SnW,predMass.mean+predMass.standard.error*CI )
# SnW$low<- with(SnW,predMass.mean-predMass.standard.error*CI )
# SwH$high<- with(SwH,predMass.mean+predMass.standard.error*CI )
# SwH$low<- with(SwH,predMass.mean-predMass.standard.error*CI )
# SwW$high<- with(SwW,predMass.mean+predMass.standard.error*CI )
# SwW$low<- with(SwW,predMass.mean-predMass.standard.error*CI )
# 
# #subset(rate.m, predMass.mean<=20)
# windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,8,6,6))
# 
# plot(predMass.mean~Time,data=NnH,type="l",las=1,ylim=c(0.2,130),
#      lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,log="y",col="darkgrey",
#      ylab=expression(Total~mass~(g)),
#      xlab="")
# lines(predMass.mean~Time,data=subset(NnH,predMass.mean <15),type="l",las=1,ylim=c(0.2,130),
#      lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=F,log="y",
#      ylab=expression(Total~mass~(g)),
#      xlab="")
# 
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"), 
#        col="darkgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"&predMass.mean<15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# lines(predMass.mean~Time, data=NnW,col="lightgrey",
#       xaxt='n', ylab="", type="l",lty=2,lwd=2)
# lines(predMass.mean~Time, data=subset(NnW, predMass.mean <15),col="red",
#       xaxt='n', ylab="", type="l",lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"), 
#        col="lightgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"&predMass.mean<15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# 
# polygon(x = c(subset(NnH, predMass.mean <15)$Time, rev(subset(NnH, predMass.mean <15)$Time)), y = c(subset(NnH, predMass.mean <15)$high,rev(subset(NnH, predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(NnW, predMass.mean <15)$Time, rev(subset(NnW, predMass.mean <15)$Time)), y = c(subset(NnW, predMass.mean <15)$high,rev(subset(NnW, predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
# mtext(text="Narrow", side=3, line=0.5, cex=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n", cex=1.2)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2, log="y", logpretty=T)
# legend("topright","a", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=NwH,col="darkgrey",legend=FALSE,ylim=c(0.2,130),xlim=c(0,60),
#        xaxt='n', yaxt='n',ylab="", type="l",log="y",lty=1,lwd=2)
# lines(predMass.mean~Time, data=subset(NwH,predMass.mean <15),col="black",legend=FALSE,ylim=c(0.2,130),xlim=c(0,60),
#        xaxt='n', yaxt='n',ylab="", type="l",log="y",lty=1,lwd=2)
# 
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"), 
#        col="darkgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# lines(predMass.mean~Time, data=NwW,col="lightgrey",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# lines(predMass.mean~Time, data=subset(NwW,predMass.mean <15),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"), 
#        col="lightgray", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# polygon(x = c(subset(NwH,predMass.mean <15)$Time, rev(subset(NwH,predMass.mean <15)$Time)), y = c(subset(NwH,predMass.mean <15)$high,rev(subset(NwH,predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(NwW,predMass.mean <15)$Time, rev(subset(NwW,predMass.mean <15)$Time)), y = c(subset(NwW,predMass.mean <15)$high,rev(subset(NwW,predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# mtext(text="Wide", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time,data=SnH,legend=FALSE,type="l",las=1,ylim=c(0.2,130),col="darkgrey",
#        log="y",lty=2,lwd=2,cex.lab=2, xlim=c(0,60),axes=FALSE,xlab="")
# lines(predMass.mean~Time,data=subset(SnH,predMass.mean <15),legend=FALSE,type="l",las=1,ylim=c(0.2,130),
#        log="y",lty=2,lwd=2,cex.lab=2, xlim=c(0,60),axes=FALSE,xlab="")
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"), 
#        col="darkgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# lines(predMass.mean~Time, data=SnW,col="lightgrey",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=2,lwd=2)
# lines(predMass.mean~Time, data=subset(SnW,predMass.mean <15),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"), 
#        col="lightgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# polygon(x = c(subset(SnH,predMass.mean <15)$Time, rev(subset(SnH,predMass.mean <15)$Time)), y = c(subset(SnH,predMass.mean <15)$high,rev(subset(SnH,predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(SnW,predMass.mean <15)$Time, rev(subset(SnW,predMass.mean <15)$Time)), y = c(subset(SnW,predMass.mean <15)$high,rev(subset(SnW,predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# legend("topright","c", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=SwH,col="darkgrey",legend=FALSE, yaxt='n',ylim=c(0.2,130),#yaxs="i",xaxs="i",
#        xaxt='n', ylab="", type="l",log="y",lty=1,lwd=2,xlim=c(0,60))
# lines(predMass.mean~Time, data=subset(SwH,predMass.mean <15),col="black",legend=FALSE, yaxt='n',ylim=c(0.2,130),#yaxs="i",xaxs="i",
#        xaxt='n', ylab="", type="l",log="y",lty=1,lwd=2,xlim=c(0,60))
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"), 
#        col="darkgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# lines(predMass.mean~Time, data=SwW, col="lightgrey",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# lines(predMass.mean~Time, data=subset(SwW,predMass.mean <15),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"), 
#        col="lightgrey", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"&predMass.mean <15)$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# 
# polygon(x = c(subset(SwH,predMass.mean <15)$Time, rev(subset(SwH,predMass.mean <15)$Time)), y = c(subset(SwH,predMass.mean <15)$high,rev(subset(SwH,predMass.mean <15)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(SwW,predMass.mean <15)$Time, rev(subset(SwW,predMass.mean <15)$Time)), y = c(subset(SwW,predMass.mean <15)$high,rev(subset(SwW,predMass.mean <15)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# legend("topright","d", bty="n", cex=1.5)
# mtext(text="Total biomass (g)", outer=T, side=2, line=3.5, cex=1.2)
# mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)
# 
# text(70,y=2500,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
# text(70,y=2,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
# ########################################################################
#
#
#
#

# 
# #Wide and narrow together
# windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))
# 
# plotBy(predMass.mean~Time,data=NnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
#        ylim=c(1,70),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,
#        ylab=expression(Total~mass~(g)),
#        xlab="")
# lines(predMass.mean~Time, data=NnW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
# lines(predMass.mean~Time, data=NwH,col="black",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
# lines(predMass.mean~Time, data=NwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
# polygon(x = c(NnH$Time, rev(NnH$Time)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NnW$Time, rev(NnW$Time)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
# polygon(x = c(NwH$Time, rev(NwH$Time)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NwW$Time, rev(NwW$Time)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)
# 
# 
# legend(0,60, legend=c("Narrow Home","Narrow +3.5","Wide Home","Wide +3.5"),
#        col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
# mtext("Tropical",3,line=-2,at=12.5, cex=1.5, outer=F)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# 
# 
# plotBy(predMass.mean~Time,data=SnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
#        ylim=c(1,70),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,
#        ylab=expression(Total~mass~(g)),
#        xlab="")
# lines(predMass.mean~Time, data=SnW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
# lines(predMass.mean~Time, data=SwH,col="black",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
# lines(predMass.mean~Time, data=SwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
# polygon(x = c(SnH$Time, rev(SnH$Time)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SnW$Time, rev(SnW$Time)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
# polygon(x = c(SwH$Time, rev(SwH$Time)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SwW$Time, rev(SwW$Time)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
# 
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext("Temperate",3,line=-2,at=14,cex=1.5, outer=FALSE)
# mtext(text="Total biomass (g)", outer=T, side=2, line=1, cex=1.2)
# mtext(text="Time (Days)", side=1, line=3, cex=1.2)
#
# ###############################################
# #2by2 version without greyed-out sections
# #add predmass from "rate" as means and SEs
# rate.m<-summaryBy(predMass~Location+Treatment+Range+Time,data=rate, FUN=c(mean, standard.error))
# 
# windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,6,6,6))
# 
# plotBy(predMass.mean~Time,data=NnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
#        ylim=c(0,82),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,
#        ylab=expression(Total~mass~(g)),
#        xlab="")
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# lines(predMass.mean~Time, data=NnW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,82),lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# polygon(x = c(NnH$Time, rev(NnH$Time)), y = c(NnH$high,rev(NnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NnW$Time, rev(NnW$Time)), y = c(NnW$high,rev(NnW$low)),col = alpha("red",0.4), border = NA)
# mtext(text="Narrow", side=3, line=0.5, cex=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n", cex=1.2)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","a", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=NwH,col="black",legend=FALSE,#yaxs="i",xaxs="i",
#       xaxt='n', yaxt='n',ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# lines(predMass.mean~Time, data=NwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# polygon(x = c(NwH$Time, rev(NwH$Time)), y = c(NwH$high,rev(NwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(NwW$Time, rev(NwW$Time)), y = c(NwW$high,rev(NwW$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text="Wide", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time,data=SnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
#        ylim=c(0,82),lty=2,lwd=2,cex.lab=2, xlim=c(1,60),axes=FALSE,xlab="")
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# lines(predMass.mean~Time, data=SnW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,82),lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# polygon(x = c(SnH$Time, rev(SnH$Time)), y = c(SnH$high,rev(SnH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SnW$Time, rev(SnW$Time)), y = c(SnW$high,rev(SnW$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","c", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=SwH,col="black",legend=FALSE, yaxt='n',#yaxs="i",xaxs="i",
#       xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# lines(predMass.mean~Time, data=SwW,col="red",
#       xaxt='n', ylab="", type="l",ylim=c(0,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",add=T,las=1))
# polygon(x = c(SwH$Time, rev(SwH$Time)), y = c(SwH$high,rev(SwH$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(SwW$Time, rev(SwW$Time)), y = c(SwW$high,rev(SwW$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","d", bty="n", cex=1.5)
# mtext(text="Total biomass (g)", outer=T, side=2, line=2.5, cex=1.2)
# mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)
# 
# text(70,y=123,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
# text(70,y=33,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
# ########################################################################
# ###############################################
# #2by2 version on log base cut at 40 days
# #add predmass from "rate" as means and SEs
# rate.m<-summaryBy(predMass~Location+Treatment+Range+Time,data=rate, FUN=c(mean, standard.error))
# 
# windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,8,6,6))
# 
# plot(predMass.mean~Time,data=subset(NnH,Time <40),type="l",las=1,ylim=c(0.2,130),
#       lty=2,lwd=2,cex.lab=2, xlim=c(1,40),axes=F,log="y",
#        ylab=expression(Total~mass~(g)),
#        xlab="")
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# lines(predMass.mean~Time, data=subset(NnW, Time<40),col="red",
#       xaxt='n', ylab="", type="l",lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# polygon(x = c(subset(NnH, Time< 40)$Time, rev(subset(NnH, Time< 40)$Time)), y = c(subset(NnH, Time< 40)$high,rev(subset(NnH, Time< 40)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(NnW, Time< 40)$Time, rev(subset(NnW, Time< 40)$Time)), y = c(subset(NnW, Time< 40)$high,rev(subset(NnW, Time< 40)$low)),col = alpha("red",0.4), border = NA)
# mtext(text="Narrow", side=3, line=0.5, cex=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), lwd=2,bty="n", cex=1.2)
# magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2, log="y", logpretty=T)
# legend("topright","a", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=subset(NwH,Time<40),col="black",legend=FALSE,ylim=c(0.2,130),
#        xaxt='n', yaxt='n',ylab="", type="l",log="y",lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# lines(predMass.mean~Time, data=subset(NwW,Time<40),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="N" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# polygon(x = c(subset(NwH,Time<40)$Time, rev(subset(NwH,Time<40)$Time)), y = c(subset(NwH,Time<40)$high,rev(subset(NwH,Time<40)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(NwW,Time<40)$Time, rev(subset(NwW,Time<40)$Time)), y = c(subset(NwW,Time<40)$high,rev(subset(NwW,Time<40)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# mtext(text="Wide", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time,data=subset(SnH,Time<40),legend=FALSE,type="l",las=1,ylim=c(0.2,130),#yaxs="i",xaxs="i",
#        log="y",lty=2,lwd=2,cex.lab=2, xlim=c(1,40),axes=FALSE,xlab="")
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# lines(predMass.mean~Time, data=subset(SnW,Time<40),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=2,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="narrow" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# polygon(x = c(subset(SnH,Time<40)$Time, rev(subset(SnH,Time<40)$Time)), y = c(subset(SnH,Time<40)$high,rev(subset(SnH,Time<40)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(SnW,Time<40)$Time, rev(subset(SnW,Time<40)$Time)), y = c(subset(SnW,Time<40)$high,rev(subset(SnW,Time<40)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# legend("topright","c", bty="n", cex=1.5)
# 
# plotBy(predMass.mean~Time, data=subset(SwH,Time<40),col="black",legend=FALSE, yaxt='n',ylim=c(0.2,130),#yaxs="i",xaxs="i",
#        xaxt='n', ylab="", type="l",log="y",lty=1,lwd=2,xlim=c(1,40),)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home"), 
#        col="black", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Home")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# lines(predMass.mean~Time, data=subset(SwW,Time<40),col="red",
#       xaxt='n', ylab="", type="l",ylim=c(1,82),lty=1,lwd=2)
# points(predMass.mean~Time, data=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed"), 
#        col="red", pch=19,cex=2,
#        panel.first=adderrorbars(x=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$Time,
#                                 y=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.mean,
#                                 SE=subset(rate.m, Location =="S" & Range=="wide" & Treatment == "Warmed")$predMass.standard.error,
#                                 direction="updown",col="black",las=1))
# polygon(x = c(subset(SwH,Time<40)$Time, rev(subset(SwH,Time<40)$Time)), y = c(subset(SwH,Time<40)$high,rev(subset(SwH,Time<40)$low)),col = alpha("black",0.4), border = NA)
# polygon(x = c(subset(SwW,Time<40)$Time, rev(subset(SwW,Time<40)$Time)), y = c(subset(SwW,Time<40)$high,rev(subset(SwW,Time<40)$low)),col = alpha("red",0.4), border = NA)
# magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2,log="y", logpretty=T)
# legend("topright","d", bty="n", cex=1.5)
# mtext(text="Total biomass (g)", outer=T, side=2, line=3.5, cex=1.2)
# mtext(text="Time (Days)", side=1, line=3, cex=1.2, adj=-0.3)
# 
# text(50,y=2500,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.5)
# text(50,y=2,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.5)
# ########################################################################
# 
# #provenance specific version
# dat<- summaryBy(predMass+Range~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error))
# dat$high<-with(dat,predMass.mean+predMass.standard.error*CI)
# dat$low<-with(dat,predMass.mean-predMass.standard.error*CI)
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
#        plot(predMass.mean~Time,col="black",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
#             xlim=c(0,65),ylim=c(0, 88),axes=FALSE,xlab="Time",ylab="Mass"))
#   with(subset(dat2,Treatment=="Home"),
#        polygon(x = c(subset(dat2,Treatment=="Home")$Time,
#                      rev(subset(dat2,Treatment=="Home")$Time)),
#                y = c(subset(dat2,Treatment=="Home")$high,
#                      rev(subset(dat2,Treatment=="Home")$low)),
#                col = alpha("black",0.4), border = NA))
#   par(new=T)
#   with(subset(dat2,Treatment=="Warmed"),
#        plot(predMass.mean~Time,col="red",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
#             xlim=c(0,65),ylim=c(0,88),axes=FALSE,xlab="Time",ylab="Mass"))
#   with(subset(dat2,Treatment=="Warmed"),
#        polygon(x = c(subset(dat2,Treatment=="Warmed")$Time,
#                      rev(subset(dat2,Treatment=="Warmed")$Time)),
#                y = c(subset(dat2,Treatment=="Warmed")$high,
#                      rev(subset(dat2,Treatment=="Warmed")$low)),
#                col = alpha("red",0.4), border = NA))
#   #first plot
#     ifelse(dat2$Taxa %in% c("BRA","PEL","BOT","LONG"),
#            magaxis(side=c(1,2),labels=c(0,1),frame.plot=T,las=1,cex.axis=1.2),
#            ifelse(dat2$Taxa %in% c("CCAM","BCAM","ECAM","FCAM"),
#                   magaxis(side=c(1,2),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2),
#                   ifelse(dat2$Taxa %in% c("SMIT","PLAT"),
#                          magaxis(side=c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2),
#                          magaxis(side=c(1,2),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2))))
#
#     legend("topleft",legend=dat2$Taxa[1],bty="n",cex=1.3)
# }
# mtext(expression(Biomass~(g)),side=2,line=3,outer=T,cex=1.5)
# mtext(expression(Time~(days)),side=1,line=3,outer=T,cex=1.5)
# mtext("Narrow",side=3,line=1,outer=T,cex=1, adj=0.08)
# mtext("Wide",side=3,line=1,outer=T,cex=1, adj=0.72)
# text(155,y=415,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
# text(155,y=135,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)
#
# legend(x=70,y=60,legend=c("Warmed","Home"),cex=1.4,
#        xpd=NA,fill=c(alpha("red",0.4),alpha("black",0.4)), border="black",
#       bty="n")

#################################################################################################
#provenance specific version on a log-base
dat<- summaryBy(predMass+Range~Time+Treatment+Taxa,data=gamfits2,FUN=c(mean,standard.error))
dat$high<-with(dat,predMass.mean+predMass.standard.error*CI)
dat$low<-with(dat,predMass.mean-predMass.standard.error*CI)
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
       plot(predMass.mean~Time,col="darkgrey",type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0.2, 88),axes=FALSE,xlab="Time",ylab="Mass",log="y"))
  with(subset(dat2,Treatment=="Home" & Time %in% unique(rate$Time)),
       points(predMass.mean~Time,col="darkgrey", pch=19,cex=2,
         panel.first=adderrorbars(x=Time,y=predMass.mean,
                                  SE=predMass.standard.error,
                                  direction="updown",col="black",las=1)))
  with(subset(dat2,Treatment=="Home"&predMass.mean<15),
       lines(predMass.mean~Time,col="black",type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0.2, 88),axes=FALSE,xlab="Time",ylab="Mass",log="y"))
  with(subset(dat2,Treatment=="Home" & Time %in% unique(rate$Time)&predMass.mean<15),
       points(predMass.mean~Time,col="black", pch=19,cex=2,
              panel.first=adderrorbars(x=Time,y=predMass.mean,
                                       SE=predMass.standard.error,
                                       direction="updown",col="black",las=1)))
  with(subset(dat2,Treatment=="Home"&predMass.mean<15),
       polygon(x = c(Time,rev(Time)),y = c(high,rev(low)),col = alpha("black",0.4), border = NA))
  par(new=T)
  with(subset(dat2,Treatment=="Warmed"),
       plot(predMass.mean~Time,col="lightgrey",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0.2,88),axes=FALSE,xlab="Time",ylab="Mass",log="y"))
  with(subset(dat2,Treatment=="Warmed" & Time %in% unique(rate$Time)),
       points(predMass.mean~Time,col="lightgrey", pch=19,cex=2,
              panel.first=adderrorbars(x=Time,y=predMass.mean,
                                       SE=predMass.standard.error,
                                       direction="updown",col="black",las=1)))
  with(subset(dat2,Treatment=="Warmed"&predMass.mean<15),
       lines(predMass.mean~Time,col="red",legend=FALSE,type="l",lty=ifelse(Range.mean == 1,2,1),
            xlim=c(0,65),ylim=c(0.2,88),axes=FALSE,xlab="Time",ylab="Mass",log="y"))
  with(subset(dat2,Treatment=="Warmed" & Time %in% unique(rate$Time)& predMass.mean<15),
       points(predMass.mean~Time,col="red", pch=19,cex=2,
              panel.first=adderrorbars(x=Time,y=predMass.mean,
                                       SE=predMass.standard.error,
                                       direction="updown",col="black",las=1)))
  with(subset(dat2,Treatment=="Warmed"&predMass.mean<15),
       polygon(x = c(Time,rev(Time)),y = c(high,rev(low)),
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
mtext(expression(Biomass~(g)),side=2,line=3,outer=T,cex=1.5)
mtext(expression(Time~(days)),side=1,line=3,outer=T,cex=1.5)
mtext("Narrow",side=3,line=1,outer=T,cex=1, adj=0.1)
mtext("Wide",side=3,line=1,outer=T,cex=1, adj=0.72)
text(150,y=10^(12.3),labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
text(150,y=400,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)

legend(x=70,y=15,legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),cex=1.4,
       xpd=NA,fill=c("red","black"), border="black",
       bty="n")

