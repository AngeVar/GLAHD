#Figure 7: respiration
#- Respiration
source("R_cleaned/3. Create_datasets.R") #get "Rdark" - spot measures of R and "rt" - RvT data and "predR" - modelled values of R

colors <- c(alpha("black",0.4),alpha("red",1))
windows(11.69,11.69)
par(mfrow=c(3,2),mar=c(3,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)

#Rmass
#- average across species
Rdark.m <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)

ylims=c(5,19)
boxplot(Rmass~Treatment*Range,data=subset(Rdark.m,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(1,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",1),alpha("black",0.6)),
       bty="n",cex=1.5)
legend("topright","a", bty="n", cex=1.5)

boxplot(Rmass~Treatment*Range,data=subset(Rdark.m,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(0,0),frame.plot=T,las=1)

mtext(text=expression(R["mass"]~"("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1.2,adj=0.5,line=3)

axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
legend("topright","b", bty="n", cex=1.5)





#- plot R vs. T

#- break data into bins of leaf temperatures, for averaging and plotting data points
rt$Tleaf_bin <- cut(rt$Tleaf,breaks=seq(from=10,to=65,by=2))
rt$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rt$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 
rt.m <- summaryBy(R_mass+R_area~Species+Taxa+Treatment+Tleaf_bin_mid,data=rt,FUN=c(mean,standard.error))
toplot <- subset(rt.m,Tleaf_bin_mid < 40)

CI<-1.96 #90% CI #1.96 95% CI

predR.df$Location<- as.factor(ifelse(predR.df$Taxa %in% c("BTER","BOT"), "S","N"))
predR.df$Range<- as.factor(ifelse(predR.df$Taxa %in% c("BRA","BOT"), "Narrow","Wide"))

predRmass<- summaryBy(R~T+Taxa+Treatment+Location+Range,data=predR.df,FUN=c(mean, standard.error))
SnH<- subset(predRmass,Taxa=="BOT" & Treatment =="Home")
SnH$high<- with(SnH,R.mean+R.standard.error*CI );SnH$low<- with(SnH,R.mean-R.standard.error*CI )
SnW<- subset(predRmass,Taxa=="BOT" & Treatment =="Warmed")
SnW$high<- with(SnW,R.mean+R.standard.error*CI );SnW$low<- with(SnW,R.mean-R.standard.error*CI )
SwH<- subset(predRmass,Taxa=="BTER" & Treatment =="Home")
SwH$high<- with(SwH,R.mean+R.standard.error*CI );SwH$low<- with(SwH,R.mean-R.standard.error*CI )
SwW<- subset(predRmass,Taxa=="BTER" & Treatment =="Warmed")
SwW$high<- with(SnW,R.mean+R.standard.error*CI );SwW$low<- with(SwW,R.mean-R.standard.error*CI )

NnH<- subset(predRmass,Taxa=="BRA" & Treatment =="Home")
NnH$high<- with(NnH,R.mean+R.standard.error*CI );NnH$low<- with(NnH,R.mean-R.standard.error*CI )
NnW<- subset(predRmass,Taxa=="BRA" & Treatment =="Warmed")
NnW$high<- with(NnW,R.mean+R.standard.error*CI );NnW$low<- with(NnW,R.mean-R.standard.error*CI )
NwH<- subset(predRmass,Taxa=="CTER" & Treatment =="Home")
NwH$high<- with(NwH,R.mean+R.standard.error*CI );NwH$low<- with(NwH,R.mean-R.standard.error*CI )
NwW<- subset(predRmass,Taxa=="CTER" & Treatment =="Warmed")
NwW$high<- with(SnW,R.mean+R.standard.error*CI );NwW$low<- with(NwW,R.mean-R.standard.error*CI )


ylims <- c(0,41)
xlims <- c(10,41)
cexpoints=1.8
colors <- c(alpha("black",0.6),alpha("red",0.6))
par(mar=c(0,0,0,0))

# plot south narrow
plotBy(R.mean~T,data=SnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
       lty=2,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
       ylab=expression(Total~mass~(g)),xlab="")
lines(R.mean~T,data=SnW,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
       lty=2,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
       ylab=expression(Total~mass~(g)),col="red",xlab="")

points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BOT"&Treatment=="Home"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[1],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BOT"&Treatment=="Home")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BOT"&Treatment=="Home")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BOT"&Treatment=="Home")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BOT"&Treatment=="Warmed"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[2],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BOT"&Treatment=="Warmed")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BOT"&Treatment=="Warmed")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BOT"&Treatment=="Warmed")$R_mass.standard.error,direction="updown",col="red"),
       axes=F,las=1,xlab="")

polygon(x = c(SnH$T,rev(SnH$T)),
        y = c(SnH$high,rev(SnH$low)),
        col = alpha("black",0.2), border = NA)
polygon(x = c(SnW$T,rev(SnW$T)),
        y = c(SnW$high,rev(SnW$low)),
        col = alpha("red",0.2), border = NA)

magaxis(c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1)
legend("topright","c", bty="n", cex=1.5)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),lty=c(1,1), 
       lwd=2,bty="n", cex=1.5)


# plot north narrow
plotBy(R.mean~T,data=NnH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
       lty=2,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
       ylab=expression(Total~mass~(g)),xlab="")
lines(R.mean~T,data=NnW,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
      lty=2,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
      ylab=expression(Total~mass~(g)),col="red",xlab="")

points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BRA"&Treatment=="Home"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[1],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BRA"&Treatment=="Home")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BRA"&Treatment=="Home")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BRA"&Treatment=="Home")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BRA"&Treatment=="Warmed"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[2],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BRA"&Treatment=="Warmed")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BRA"&Treatment=="Warmed")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BRA"&Treatment=="Warmed")$R_mass.standard.error,direction="updown",col="red"),
       axes=F,las=1,xlab="")

polygon(x = c(NnH$T,rev(NnH$T)),
        y = c(NnH$high,rev(NnH$low)),
        col = alpha("black",0.2), border = NA)
polygon(x = c(NnW$T,rev(NnW$T)), 
        y = c(NnW$high,rev(NnW$low)),
        col = alpha("red",0.2), border = NA)

magaxis(c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1)
legend("topright","d", bty="n", cex=1.5)


# plot south wide
plotBy(R.mean~T,data=SwH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
       lty=1,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
       ylab=expression(Total~mass~(g)),xlab="")
lines(R.mean~T,data=SwW,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
      lty=1,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
      ylab=expression(Total~mass~(g)),col="red",xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BTER"&Treatment=="Home"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[1],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BTER"&Treatment=="Home")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BTER"&Treatment=="Home")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BTER"&Treatment=="Home")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="BTER"&Treatment=="Warmed"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[2],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="BTER"&Treatment=="Warmed")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="BTER"&Treatment=="Warmed")$R_mass.mean,
                                SE=subset(toplot,Taxa=="BTER"&Treatment=="Warmed")$R_mass.standard.error,direction="updown",col="red"),
       axes=F,las=1,xlab="")
polygon(x = c(SwH$T,rev(SwH$T)),
        y = c(SwH$high,rev(SwH$low)),
        col = alpha("black",0.2), border = NA)
polygon(x = c(SwW$T,rev(SwW$T)), 
        y = c(SwW$high,rev(SwW$low)),
        col = alpha("red",0.2), border = NA)
magaxis(c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1)
legend("topright","e", bty="n", cex=1.5)

# plot north wide

plotBy(R.mean~T,data=NwH,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
       lty=1,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
       ylab=expression(Total~mass~(g)),xlab="")
lines(R.mean~T,data=NwW,legend=FALSE,type="l",las=1,#yaxs="i",xaxs="i",
      lty=1,lwd=2,cex.lab=2, axes=FALSE,xlim=xlims,ylim=ylims,
      ylab=expression(Total~mass~(g)),col="red",xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="CTER"&Treatment=="Home"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[1],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="CTER"&Treatment=="Home")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="CTER"&Treatment=="Home")$R_mass.mean,
                                SE=subset(toplot,Taxa=="CTER"&Treatment=="Home")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")
points(R_mass.mean~Tleaf_bin_mid,data=subset(toplot,Taxa=="CTER"&Treatment=="Warmed"),xlim=xlims,ylim=ylims, pch=19,cex=cexpoints,col=colors[2],legend=F,
       panel.first=adderrorbars(x=subset(toplot,Taxa=="CTER"&Treatment=="Warmed")$Tleaf_bin_mid,
                                y=subset(toplot,Taxa=="CTER"&Treatment=="Warmed")$R_mass.mean,
                                SE=subset(toplot,Taxa=="CTER"&Treatment=="Warmed")$R_mass.standard.error,direction="updown",col="red"),
       axes=F,las=1,xlab="")
polygon(x = c(NwH$T,rev(NwH$T)),
        y = c(NwH$high,rev(NwH$low)),
        col = alpha("black",0.2), border = NA)
polygon(x = c(NwW$T,rev(NwW$T)), 
        y = c(NwW$high,rev(NwW$low)),
        col = alpha("red",0.2), border = NA)
magaxis(c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1)
legend("topright","f", bty="n", cex=1.5)

mtext(text=expression(T["leaf"]~"("*degree~C*")"),side=1,cex=1.2,adj=-0.1,line=3)
text(45,y=60,labels="Narrow", xpd=NA, srt=-90, pos=2, cex=2)
text(45,y=17,labels="Wide", xpd=NA, srt=-90, pos=2, cex=2)
