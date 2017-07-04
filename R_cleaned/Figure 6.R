#Figure 6: Photosynthesis
Asatm$Photom<-with(Asatm, (Photo*SLA/10000)*1000)#from micro to nano
acifits$Jmaxm<-with(acifits, Jmax*SLA/10000)
acifits$Vcmaxm<-with(acifits, Vcmax*SLA/10000)
acifits$JtoVm<-with(acifits,Jmaxm/Vcmaxm)

asat.tm <- summaryBy(Photo~Taxa+Treatment+Location+Range,data=Asat,FUN=mean,keep.names=T)
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)
jmax.tm <- summaryBy(Jmax~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
vcmax.tm <- summaryBy(Vcmax~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=mean,keep.names=T)


colors <- c(alpha("black",0.4),"red")
#with SMIT
windows(11.69,11.69);par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(6,9,6,9))
#Jmax
ylims=c(90,250)
boxplot(Jmax~Treatment*Range,data=subset(jmax.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, 
       pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
       bty="n",cex=1.5)

mtext(text=expression(J["max"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
legend("topright","a", bty="n", cex=1.5)

boxplot(Jmax~Treatment*Range,data=subset(jmax.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
legend("topright","b", bty="n", cex=1.5)


#Vcmax
ylims=c(60,220)
boxplot(Vcmax~Treatment*Range,data=subset(vcmax.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.5,line=3)
legend("topright","c", bty="n", cex=1.5)

boxplot(Vcmax~Treatment*Range,data=subset(vcmax.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","d", bty="n", cex=1.5)

#Rmass
ylims=c(5,19)
boxplot(Rmass~Treatment*Range,data=subset(Rmass.tm,Location=="S"),
        ylim=ylims,axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(1,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
legend("topright","e", bty="n", cex=1.5)      
mtext(text=expression(R["mass"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.12,line=3)

boxplot(Rmass~Treatment*Range,data=subset(Rmass.tm,Location=="N"),ylim=ylims,
                axes=F,las=2,col=colors)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(0,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
legend("topright","f", bty="n", cex=1.5)

######################################
windows(11.69,11.69);par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(6,9,6,9))

#Asat area
ylims=c(20,35)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
legend("topright","a", bty="n", cex=1.5)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, 
       pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
       bty="n",cex=1.5)

boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","b", bty="n", cex=1.5)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)

#Asat mass
ylims=c(350,700)
boxplot(Photom~Treatment*Range,data=subset(asatm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.5,line=3)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topright","c", bty="n", cex=1.5)

boxplot(Photom~Treatment*Range,data=subset(asatm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topright","d", bty="n", cex=1.5)

# #W/o SMIT
# jmax.tms<- subset(jmax.tm, Taxa != "SMIT")
# vcmax.tms<- subset(vcmax.tm, Taxa != "SMIT")
# asat.tms<- subset(asat.tm, Taxa != "SMIT")
# 
# windows(11.69,11.69);par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(6,9,6,9))
# #Jmax
# ylims=c(90,250)
# boxplot(Jmax~Treatment*Range,data=subset(jmax.tms,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
#        bty="n",cex=1.2)
# 
# mtext(text=expression(J["max"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
# mtext(text="Temperate", side=3, line=0.5, cex=1.2)
# legend("topright","a", bty="n", cex=1.5)
# 
# boxplot(Jmax~Treatment*Range,data=subset(jmax.tms,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text="Tropical", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.5)
# 
# 
# #Vcmax
# ylims=c(60,220)
# boxplot(Vcmax~Treatment*Range,data=subset(vcmax.tms,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.5,line=3)
# legend("topright","c", bty="n", cex=1.5)
# 
# boxplot(Vcmax~Treatment*Range,data=subset(vcmax.tms,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","d", bty="n", cex=1.5)
# 
# #Asat
# ylims=c(20,35)
# boxplot(Photo~Treatment*Range,data=subset(asat.tms,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.12,line=3)
# axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
# legend("topright","e", bty="n", cex=1.5)
# 
# boxplot(Photo~Treatment*Range,data=subset(asat.tms,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
# legend("topright","f", bty="n", cex=1.5)
# ##---------------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------------

## --- Mass based rates

#Figure 5: Photosynthesis

# asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=mean,keep.names=T)
# Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)
# jmaxm.tm <- summaryBy(Jmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
# vcmaxm.tm <- summaryBy(Vcmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
# 
# colors <- c(alpha("black",0.4),"red")
# 
# windows(11.69,11.69)
# par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(6,9,6,9))
# #Jmaxm
# ylims=c(2,4)
# boxplot(Jmaxm~Treatment*Range,data=subset(jmaxm.tm,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
#        bty="n",cex=1.2)
# 
# mtext(text=expression(J["max"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
# mtext(text="Temperate", side=3, line=0.5, cex=1.2)
# legend("topright","a", bty="n", cex=1.5)
# 
# boxplot(Jmaxm~Treatment*Range,data=subset(jmaxm.tm,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text="Tropical", side=3, line=0.5, cex=1.2)
# legend("topright","b", bty="n", cex=1.5)
# 
# 
# #Vcmaxm
# ylims=c(1,4.5)
# boxplot(Vcmaxm~Treatment*Range,data=subset(vcmaxm.tm,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.5,line=3)
# legend("topright","c", bty="n", cex=1.5)
# 
# boxplot(Vcmaxm~Treatment*Range,data=subset(vcmaxm.tm,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topright","d", bty="n", cex=1.5)
# 
# #Asat
# ylims=c(350,700)
# boxplot(Photom~Treatment*Range,data=subset(asatm.tm,Location=="S"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
# axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
# legend("topright","e", bty="n", cex=1.5)
# 
# boxplot(Photom~Treatment*Range,data=subset(asatm.tm,Location=="N"),ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
# legend("topright","f", bty="n", cex=1.5)
# 
# ##---------------------------------------------
# 
# 
# ##--- SLA test : SLA increased by warming in S but not N (treat x Loc) and in narrow not wide (treat x Range)
# 
# asatmt<- Asatm[,c("Location","Range","Taxa","Treatment","SLA")]
# 
# #from Asat
# windows(18,5.845)
# par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(4,6,2,2))
# 
# asatmt$TaxTreat<- as.factor(paste(asatmt$Taxa,asatmt$Treatment, sep="."))
# asatmt$TaxTreat<- factor(asatmt$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# 
# ylims=c(120,320)
# boxplot(SLA~TaxTreat,data=asatmt,ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression("From"~A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("SLA ("*cm^2~g^-1*")"),side=2,outer=F,cex=1,line=3)
# #axis(side=1,at=c(1:34),labels=rep(levels(asatm_sort$Taxa),2),las=1,cex.axis=1)
# axis(side=1,at=c(seq(1.5,34, by=2)),labels=c(levels(asatmt$Taxa)),las=1,cex.axis=1.5)
# 
# #taxa specific responses
# #area based
# 
# colors <- c(alpha("black",0.4),"red")
# acifits$TaxTreat<- as.factor(paste(acifits$Taxa,acifits$Treatment, sep="."))
# Asat$TaxTreat<- as.factor(paste(Asat$Taxa,Asat$Treatment, sep="."))
# Rdark$TaxTreat<- as.factor(paste(Rdark$Taxa,Rdark$Treatment, sep="."))
# 
# acifits$TaxTreat<- factor(acifits$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# Asat$TaxTreat<- factor(Asat$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# Rdark$TaxTreat<- factor(Rdark$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# 
# windows(16.53,11.69);par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(6,9,6,5))
# #Jmax
# ylims=c(90,280)
# boxplot(Jmax~TaxTreat,data=acifits,ylim=ylims,axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, 
#        pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
#        bty="n",cex=1.2)
# 
# mtext(text=expression(J["max"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
# legend("topright","a", bty="n", cex=1.5)
# 
# #Vcmax
# ylims=c(60,260)
# boxplot(Vcmax~TaxTreat,data=acifits,ylim=ylims,axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.64,line=3)
# legend("topright","b", bty="n", cex=1.5)
# 
# #Asat
# ylims=c(15,35)
# boxplot(Photo~TaxTreat,data=Asat,ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.34,line=3)
# legend("topright","c", bty="n", cex=1.5)
# 
# #Rdark
# ylims=c(0,1)
# boxplot(R~TaxTreat,data=Rdark,ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(R["dark"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
# axis(side=1,at=c(seq(1.5,34, by=2)),labels=c(levels(Asat$Taxa)),las=1,cex.axis=1.5)
# legend("topright","d", bty="n", cex=1.5)
# 
# #mass based
# Asatm$Photom<-with(Asatm, Photo*SLA/10000)
# 
# colors <- c(alpha("black",0.4),"red")
# acifits$TaxTreat<- as.factor(paste(acifits$Taxa,acifits$Treatment, sep="."))
# Asatm$TaxTreat<- as.factor(paste(Asatm$Taxa,Asatm$Treatment, sep="."))
# 
# acifits$TaxTreat<- factor(acifits$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# Asatm$TaxTreat<- factor(Asatm$TaxTreat, levels=c(
#   "ATER.Home","ATER.Warmed","BTER.Home","BTER.Warmed","ACAM.Home","ACAM.Warmed","BCAM.Home","BCAM.Warmed",
#   "CCAM.Home","CCAM.Warmed","BOT.Home","BOT.Warmed","LONG.Home","LONG.Warmed","SMIT.Home","SMIT.Warmed",
#   "CTER.Home","CTER.Warmed","DTER.Home","DTER.Warmed", "ETER.Home", "ETER.Warmed","DCAM.Home","DCAM.Warmed", 
#   "ECAM.Home","ECAM.Warmed","FCAM.Home","FCAM.Warmed","BRA.Home","BRA.Warmed","PEL.Home","PEL.Warmed",
#   "PLAT.Home","PLAT.Warmed"))
# 
# windows(16.53,11.69);par(mfrow=c(4,1),mar=c(0,0,0,0),oma=c(6,9,6,5))
# #Jmax
# ylims=c(1.5,6)
# boxplot(Jmaxm~TaxTreat,data=acifits,ylim=ylims,axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, 
#        pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
#        bty="n",cex=1.2)
# 
# mtext(text=expression(J["max"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=3)
# legend("topright","a", bty="n", cex=1.5)
# 
# #Vcmax
# ylims=c(1,5.5)
# boxplot(Vcmaxm~TaxTreat,data=acifits,ylim=ylims,axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(V["cmax"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.65,line=3)
# legend("topright","b", bty="n", cex=1.5)
# 
# #Asat
# ylims=c(0.28,0.85)
# boxplot(Photom~TaxTreat,data=Asatm,ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.35,line=3)
# legend("topright","c", bty="n", cex=1.5)
# 
# #rdark
# ylims=c(5,20)
# boxplot(Rmass~TaxTreat,data=Rdark,ylim=ylims,
#         axes=F,las=2,col=colors)
# magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext(text=expression(R["mass"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
# mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
# axis(side=1,at=c(seq(1.5,34, by=2)),labels=c(levels(Asatm$Taxa)),las=1,cex.axis=1.5)
# legend("topright","d", bty="n", cex=1.5)
