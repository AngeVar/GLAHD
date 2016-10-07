#Figure 6: respiration
#- Respiration
dat <- read.csv("Data/GHS39_GLAHD_MAIN_GX-RVT_20150114-20150120_L2.csv")

dat$Taxa <- as.factor(sapply(1:nrow(dat),function(i)strsplit(as.character(dat$Code[i]),split="-")[[1]][1])) # get the "Taxa" factor out of Code. Note this could have been a loop
dat$R_area <- dat$Photo*-1

#- read in the leaf mass data, calculate R per unit leaf mass
leafmass <- read.csv("Data/GHS39_GLAHD_MAIN_BIOMASS-RVTLEAVES_20150120_L2.csv")

r <- merge(dat,leafmass,by="Code")
r$R_mass <- r$R_area/10000*r$Area/r$Leafmass*1000 #note mass was in g

rt <- subset(r,Taxa!="SMIT" & Taxa !="ACAM") #remove SMIT and ACAM
rt$Taxa <- factor(rt$Taxa,levels=c("BOT","BTER","BRA","CTER"))
rt$lnRmass <- log(rt$R_mass)
rt <- subset(rt,Code!="CTER-4" & Code!="CTER-22" & Code!="CTER-27") #CTER4 looks weird... peaks much much lower than others. CTER-22 also looks bad
rt$Code <- factor(rt$Code)
gx2 <- getRdark()
#- break data into bins of leaf temperatures, for averaging and plotting
rt$Tleaf_bin <- cut(rt$Tleaf,breaks=seq(from=10,to=65,by=2))#length=60))
rt$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rt$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 

#- average across taxa
gx2$Species <- factor(gx2$Species)
gx2.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=gx2,FUN=mean,keep.names=T)

#- make plots
colors <- c(alpha("black",0.4),alpha("red",1))
windows(11.69,11.69);par(mfrow=c(2,2),mar=c(1.5,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)

#Rmass
ylims=c(5,19)
boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(1,0),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
#legend("topright","b",bty="n",cex=1.3)
#mtext("Range size",side=1,cex=1.5,outer=T,line=-17)
legend("topleft", legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22, pt.cex=2, pt.bg=c(alpha("red",1),alpha("black",0.6)),
       bty="n",cex=1.2)
legend("topright","a", bty="n", cex=1.2)

boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
magaxis(c(2,4),minorn=5, majorn=3,labels=c(0,0),frame.plot=T,las=1)
mtext(text=expression(R["mass"]~"("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1.2,adj=0.5,line=3)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.2)
#legend("topright","a",bty="n",cex=1.3)
legend("topright","b", bty="n", cex=1.2)

#- plot R vs. T
rt.m <- summaryBy(R_mass+R_area~Species+Taxa+Treatment+Tleaf_bin_mid,data=rt,FUN=c(mean,standard.error))
toplot <- subset(rt.m,Tleaf_bin_mid < 40)
toplot.s <- subset(toplot,Taxa=="BOT" | Taxa=="BTER")
toplot.n <- subset(toplot,Taxa=="BRA" | Taxa=="CTER")

ylims <- c(0,49)
xlims <- c(10,45)
cexpoints=1.8
colors <- c(alpha("black",0.6),alpha("red",0.6))
# plot south
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.s,Taxa=="BOT"),xlim=xlims,ylim=ylims,pch=1,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.s,Taxa=="BOT")$Tleaf_bin_mid,
                                y=subset(toplot.s,Taxa=="BOT")$R_mass.mean,
                                SE=subset(toplot.s,Taxa=="BOT")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")

plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.s,Taxa=="BTER"),ylim=ylims,pch=16,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.s,Taxa=="BTER")$Tleaf_bin_mid,
                                y=subset(toplot.s,Taxa=="BTER")$R_mass.mean,
                                SE=subset(toplot.s,Taxa=="BTER")$R_mass.standard.error,direction="updown",col="black"),
       add=T,las=1,axes=F,xlab="")
legend("topright","c", bty="n", cex=1.2)

magaxis(c(1,2,3,4),labels=c(1,1,0,0),frame.plot=T,las=1)

legend("topleft",legend=c(expression(paste("Narrow","(+3.5)",degree,"C")), "Narrow Home",
                          expression(paste("Wide","(+3.5)",degree,"C")),
                          "Wide Home"),pch=c(1,1,16,16),col=c("red","black"),ncol=1,
       cex=1.2,bty="n")

# plot north
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.n,Taxa=="BRA"),xlim=xlims,ylim=ylims,pch=1,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.n,Taxa=="BRA")$Tleaf_bin_mid,
                                y=subset(toplot.n,Taxa=="BRA")$R_mass.mean,
                                SE=subset(toplot.n,Taxa=="BRA")$R_mass.standard.error,direction="updown",col="black"),
       las=1,axes=F,xlab="")
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.n,Taxa=="CTER"),ylim=ylims,pch=16,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.n,Taxa=="CTER")$Tleaf_bin_mid,
                                y=subset(toplot.n,Taxa=="CTER")$R_mass.mean,
                                SE=subset(toplot.n,Taxa=="CTER")$R_mass.standard.error,direction="updown",col="black"),
       add=T,las=1)
legend("topright","d", bty="n", cex=1.2)
magaxis(c(1,2,3,4),labels=c(1,0,0,0),frame.plot=T,las=1)
mtext(text=expression(T["leaf"]~"("*degree~C*")"),side=1,cex=1.2,adj=-0.2,line=3)

# 
# home<- subset(gx2.tm, Treatment=="Home")
# warm<- subset(gx2.tm, Treatment=="Warmed")
# names(warm)<- c("Taxa","Treatment.w","Location","Range","Rmass.w")
# gx3<- merge(home,warm, by=c("Taxa", "Location","Range"))
# gx3$acSettemp<- with(gx3, Rmass/Rmass.w)
# boxplot(gx3$acSettemp~gx3$Location*gx3$Range)
