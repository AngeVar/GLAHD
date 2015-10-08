#load some libraries
source("R/loadLibraries.R")
library(raster)
library(rgdal)
library(stringr)
library(scales)#load file
seed_atsc<- read.csv("Data/GLAHDseed.csv")

#load known distribution
dist <- read.csv("Data/spALAdata2503v2.csv")

#load data
biodat <- getData("worldclim", var="bio", res=2.5, path="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data")
#subset data
biodat10 <- subset(biodat,10)
#clip data to E Australia
YbrevRange <- extent(141.00, 154.00, -44, -10.0)
biodat.oz10 <- crop(biodat10,YbrevRange)
xy <- SpatialPoints(cbind(dist$Longitude...processed,dist$Latitude...processed))
dist$bio10 <- extract(biodat.oz10/10,xy,method="bilinear",fun=mean, buffer=15000)
#add new coordinates to brassiana
seed_atsc$lat[seed_atsc$Species == "brassiana"] <- -14.2107
seed_atsc$lon[seed_atsc$Species == "brassiana"] <- 144.1304
#Extract climate data for new coordinates
xy <- SpatialPoints(cbind(seed_atsc$lon,seed_atsc$lat))
seed_atsc$bio10 <- extract(biodat.oz10/10,xy,method="bilinear", fun=mean,buffer=15000)

dis <- subset(dist,dist$Coordinate.Uncertainty.in.Metres...parsed<5000)

#subset data per species
cam<- subset(seed_atsc, Taxa == "camaldulensis")
cam$Taxa <- factor(cam$Taxa)
camdist <- subset(dis, Species...matched == "Eucalyptus camaldulensis")
camdist$Species...matched <- factor(camdist$Species...matched)
ter<- subset(seed_atsc, Taxa == "tereticornis")
ter$Taxa <- factor(ter$Taxa)
terdist <- subset(dis, Species...matched == "Eucalyptus tereticornis")
terdist$Species...matched <- factor(terdist$Species...matched)

bra<- subset(seed_atsc, Taxa == "brassiana")
bra$Taxa <- factor(bra$Taxa)
bradist <- subset(dis, Species...matched == "Eucalyptus brassiana")
bradist$Species...matched <- factor(bradist$Species...matched)
pel<- subset(seed_atsc, Taxa == "pellita")
pel$Taxa <- factor(pel$Taxa)
peldist <- subset(dis, Species...matched == "Eucalyptus pellita")
peldist$Species...matched <- factor(peldist$Species...matched)
pla<- subset(seed_atsc, Taxa == "platyphylla")
pla$Taxa <- factor(pla$Taxa)
pladist <- subset(dis, Species...matched == "Eucalyptus platyphylla")
pladist$Species...matched <- factor(pladist$Species...matched)

bot<- subset(seed_atsc, Taxa == "botryoides")
bot$Taxa <- factor(bot$Taxa)
botdist <- subset(dis, Species...matched == "Eucalyptus botryoides")
botdist$Species...matched <- factor(botdist$Species...matched)
lon<- subset(seed_atsc, Taxa == "longifolia")
lon$Taxa <- factor(lon$Taxa)
londist <- subset(dis, Species...matched == "Eucalyptus longifolia")
londist$Species...matched <- factor(londist$Species...matched)
smi<- subset(seed_atsc, Taxa == "smithii")
smi$Taxa <- factor(smi$Taxa)
smidist <- subset(dis, Species...matched == "Eucalyptus smithii")
smidist$Species...matched <- factor(smidist$Species...matched)

wide<-rbind(cam,ter)
widedist<-rbind(camdist,terdist)
narrowS<-rbind(bot,lon,smi)
narrowN<-rbind(bra, pel, pla)
narrowSdist<-rbind(botdist,londist,smidist)
narrowNdist<-rbind(bradist, peldist, pladist)
  
#density plot
a<-density(terdist$bio10, na.rm=T)
b<-density(camdist$bio10, na.rm=T)
c<-density(pladist$bio10, na.rm=T)
d<-density(bradist$bio10, na.rm=T)
e<-density(peldist$bio10, na.rm=T)
f<-density(botdist$bio10, na.rm=T)
g<-density(londist$bio10, na.rm=T)
h<-density(smidist$bio10, na.rm=T)

i<-density(widedist$bio10, na.rm=T)
j<-density(narrowNdist$bio10, na.rm=T)
k<-density(narrowSdist$bio10, na.rm=T)
 
# plot(a, ylim=c(0,1.4))
# polygon(a, col=alpha("black",0.3), border= "black")
# polygon(b, col=alpha("white",0.3), border= "black")
# polygon(c, col=alpha("red",0.3), border= "black")
# polygon(d, col=alpha("yellow",0.3), border= "black")
# polygon(e, col=alpha("orange",0.3), border= "black")
# polygon(f, col=alpha("dodgerblue",0.3), border= "black")
# polygon(g, col=alpha("cyan",0.3), border= "black")
# polygon(h, col=alpha("purple",0.3), border= "black")

plot(i, ylim=c(0,0.6))
polygon(i, col=alpha("black",0.3), border= "black")
polygon(j, col=alpha("white",0.3), border= "black")
polygon(k, col=alpha("white",0.3), border= "black")


widehist <- hist(widedist$bio10)
widemultiplier <- widehist$counts / widehist$density
widedensity <- density(widedist$bio10, na.rm=T)
widedensity$y <- widedensity$y * widemultiplier[1]

narrowNhist <- hist(narrowNdist$bio10)
narrowNmultiplier <- narrowNhist$counts / narrowNhist$density
narrowNdensity <- density(narrowNdist$bio10, na.rm=T)
narrowNdensity$y <- narrowNdensity$y * narrowNmultiplier[1]

narrowShist <- hist(narrowSdist$bio10)
narrowSmultiplier <- narrowShist$counts / narrowShist$density
narrowSdensity <- density(narrowSdist$bio10, na.rm=T)
narrowSdensity$y <- narrowSdensity$y * narrowSmultiplier[1]

windows(width=5, height=5)
  par(mar=c(5,0,0,0))
  plot(widehist, border="white", ylim=c(0,800), xlim=c(13,33), xlab="",
       ylab="", main="",axes=F,line=1)
  axis(1, at=seq(14,32,1), labels=F, line=-0.5)
  axis(1, at=seq(15,32,5), line=-0.5)
  #axis(2, at=seq(0,1000,200), line=-1.5)
  mtext(text=expression(Mean~temperature~of~warmest~quarter~(degree*C)),side=1,outer=T,line=-3,cex=1)
  polygon(widedensity, col=alpha("black",0.4), border= "black")
  polygon(narrowNdensity, col=alpha("white",0.5), border= "black")
  polygon(narrowSdensity, col=alpha("white",0.5), border= "black")
  #axis(1, at=c(28.22286,22.04737,19.91316,20.97619,26.92,26.77), labels=F, tick=T, lwd.ticks=8, col.ticks= "black", tcl=1, line=0.5,lwd=0)
  #axis(1, at=c(28.22286,22.04737,19.91316,20.97619,26.92,26.77), labels=F, tick=T, lwd.ticks=6, col.ticks= "darkgrey", tcl=1, line=0.5,lwd=0)
#   segments(18.77274,-180,18.93623,-180, col=alpha("blue",0.5),lty=1,lwd=8,xpd=T)
#   segments(27.36636,-180,27.53607,-180, col=alpha("blue",0.5),lty=1,lwd=8,xpd=T)
#   segments(22.20515,-180,22.37852,-180, col=alpha("red",0.5),lty=1,lwd=8,xpd=T)
#   segments(30.84427,-180,31.02614,-180, col=alpha("red",0.5),lty=1,lwd=8,xpd=T)
  arrows(18.85448,30,22.29183,30,code=2, lwd=2, length=0.1, angle=20)
  arrows(27.45121,30,30.93521,30,code=2, lwd=2, length=0.1, angle=20)
legend(27.5,800,legend=c("Wide","Narrow"), fill=c(alpha("black",0.4),alpha("white",0.5)),border= "black")
  #axis(1, at=c(28.22286,22.04737,19.91316,20.97619,26.92,26.77), labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
  #axis(1, at=ter$bio10, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=pla$bio10, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=pel$bio10, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=bra$bio10, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=18.7, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=18.8, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
#   axis(1, at=smi$bio10, labels=F, tick=T, lwd.ticks=4, col.ticks= "black", tcl=1, line=0.5,lwd=0)
  
#pick random subsample and see if it makes the figure prettier
terdist<- terdist[sample(1:nrow(terdist),2000,replace=F),]
camdist<- camdist[sample(1:nrow(camdist),2000,replace=F),]
botdist<- botdist[sample(1:nrow(botdist),400,replace=F),]
londist<- londist[sample(1:nrow(londist),400,replace=F),]
smidist<- smidist[sample(1:nrow(smidist),400,replace=F),]

wide<-rbind(cam,ter)
widedist<-rbind(camdist,terdist)
# narrow<-rbind(bra, pel, pla,bot,lon,smi)
# narrowdist<-rbind(bradist, peldist, pladist,botdist,londist,smidist)

narrowS<-rbind(bot,lon,smi)
narrowN<-rbind(bra, pel, pla)
narrowSdist<-rbind(botdist,londist,smidist)
narrowNdist<-rbind(bradist, peldist, pladist)  