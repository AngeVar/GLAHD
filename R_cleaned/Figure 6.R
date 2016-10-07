#Figure 5: Photosynthesis
Asatm$Photom<-with(Asatm, Photo*SLA/10000)
acifits$Jmaxm<-with(acifits, Jmax*SLA/10000)
acifits$Vcmaxm<-with(acifits, Vcmax*SLA/10000)

asat.tm <- summaryBy(Photo~Taxa+Treatment+Location+Range,data=Asat,FUN=mean,keep.names=T)
asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=mean,keep.names=T)
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)
jmax.tm <- summaryBy(Jmax~Taxa+Treatment+Location+Range,data=acifit,FUN=mean,keep.names=T)
vcmax.tm <- summaryBy(Vcmax~Taxa+Treatment+Location+Range,data=acifit,FUN=mean,keep.names=T)
jmaxm.tm <- summaryBy(Jmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
vcmaxm.tm <- summaryBy(Vcmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)

colors <- c(alpha("black",0.4),"red")

windows(11.69,11.69);par(mfrow=c(3,2),mar=c(0,0,0,0),oma=c(6,9,6,9))
#Jmax
ylims=c(90,250)
boxplot(Jmax~Treatment*Range,data=subset(jmax.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
       bty="n",cex=1.2)

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

#Asat
ylims=c(20,35)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1.2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topright","e", bty="n", cex=1.5)

boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topright","f", bty="n", cex=1.5)


