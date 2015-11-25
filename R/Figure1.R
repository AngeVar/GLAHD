#load some libraries
source("R/loadLibraries.R")
library(oz)
library(maps)
library(mapdata)
library(RNCEP)
library(mapplots)
library(plotrix)
library(plotBy)
library(doBy)
library(dismo)
library(rgeos)
library(raster)
library(rgdal)
library(stringr)
library(scales)#load file

biodat <- getData("worldclim", var="bio", res=2.5, path="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data")
biodat10 <- subset(biodat,10)
biodat1<-subset(biodat,1)
YbrevRange <- extent(141.00, 154.00, -40, -6.0)
biodat.oz10 <- crop(biodat10,YbrevRange)
biodat.oz1 <- crop(biodat1,YbrevRange)

seed_atsc<- read.csv("Data/GLAHDseed.csv")
dist <- read.csv("Data/ALAdata20151022Clean.csv")

cam<- subset(seed_atsc, Taxa == "camaldulensis");ter<- subset(seed_atsc, Taxa == "tereticornis")
bra<- subset(seed_atsc, Taxa == "brassiana");pel<- subset(seed_atsc, Taxa == "pellita");pla<- subset(seed_atsc, Taxa == " platyphylla")
bot<- subset(seed_atsc, Taxa == "botryoides");lon<- subset(seed_atsc, Taxa == "longifolia");smi<- subset(seed_atsc, Taxa == "smithii")

wide<-subset(seed_atsc, Taxa == "camaldulensis"|Taxa == "tereticornis")
Nn<-subset(seed_atsc, Taxa == "brassiana"|Taxa == "pellita"|Taxa == " platyphylla")
Sn<-subset(seed_atsc, Taxa == "botryoides"|Taxa == "longifolia"|Taxa == "smithii")

widedist<- subset(dist, Species...matched == "Eucalyptus tereticornis"|Species...matched == "Eucalyptus camaldulensis"& "Longitude...processed" > 141.00)
Nndist<-subset(dist, Species...matched == "Eucalyptus brassiana"|Species...matched == "Eucalyptus pellita"|Species...matched == "Eucalyptus platyphylla"& "Longitude...processed" > 141.00)
Sndist<- subset(dist, Species...matched == "Eucalyptus botryoides"|Species...matched == "Eucalyptus longifolia"|Species...matched == "Eucalyptus smithii"& "Longitude...processed" > 141.00)

coordinates(widedist) <- c("Longitude...processed","Latitude...processed")
coordinates(Nndist) <- c("Longitude...processed","Latitude...processed")
coordinates(Sndist) <- c("Longitude...processed","Latitude...processed")

projection(widedist) <- CRS('+proj=longlat') 
projection(Nndist) <- CRS('+proj=longlat') 
projection(Sndist) <- CRS('+proj=longlat') 

w <- circles(widedist, d=50000, lonlat=TRUE) 
Nn <- circles(Nndist, d=50000, lonlat=TRUE) 
Sn <- circles(Sndist, d=50000, lonlat=TRUE) 

wpol <- gUnaryUnion(w@polygons) 
Nnpol <- gUnaryUnion(Nn@polygons) 
Snpol <- gUnaryUnion(Sn@polygons) 

windows(5.845,5)
par(mfrow=c(1,3), mar=c(2,0,2,0),oma=c(2,5,2,10))
plot(biodat.oz10/10,main="Selected",cex.main=1.5 , xlim=c(141,154), ylim=c(-40,-6),legend=F, col=alpha("black",0.3))
points(x=cam$lon,y=cam$lat,col="black", bg="white",cex=2.5,pch=21,lwd=2)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=2.5,pch=21,lwd=2)
points(x=pla$lon,y=pla$lat,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=pel$lon,y=pel$lat,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=145.15,y=-16.58,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)
points(x=lon$lon,y=lon$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)
points(x=smi$lon,y=smi$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)


plot(biodat.oz10/10,main="Wide",cex.main=1.5, xlim=c(141,154), ylim=c(-40,-6),yaxt='n', legend=F, col=alpha("black",0.3))
plot(wpol,add=T,col="white")
points(x=cam$lon,y=cam$lat,col="black", bg="white",cex=2.5,pch=21,lwd=2)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=2.5,pch=21,lwd=2)
box()

plot(biodat.oz10/10,main="Narrow",cex.main=1.5 ,xlim=c(141,154), ylim=c(-40,-6),yaxt='n',legend=F, col=alpha("black",0.3))
plot(Nnpol,add=T,col="red")
plot(Snpol,add=T,col="dodgerblue")
points(x=145.15,y=-16.58,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=pla$lon,y=pla$lat,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=pel$lon,y=pel$lat,col="black", bg="red",cex=2.5,pch=21,lwd=2)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)
points(x=lon$lon,y=lon$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)
points(x=smi$lon,y=smi$lat,col="black", bg="dodgerblue",cex=2.5,pch=21,lwd=2)

mtext("Latitude", side=2, outer=T, line=2.5)
mtext("Longitude", side=1, outer=T, line=1)
text(155.5,-10, labels= "Tropical", xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-12, labels= expression(italic("E. brassiana")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-14.2, labels= expression(italic("E. pellita")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-16.2, labels= expression(italic("E. platyphylla")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-33, labels= "Temperate", xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-35, labels= expression(italic("E. botryoides")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-36.5, labels= expression(italic("E. longifolia")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-37.8, labels= expression(italic("E. smithii")), xpd=NA, cex=1.5, adj = c(0,0))

text(155.5,-22, labels= "Wide", xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-24, labels= expression(italic("E. camaldulensis")), xpd=NA, cex=1.5, adj = c(0,0))
text(155.5,-25.8, labels= expression(italic("E. tereticornis")), xpd=NA, cex=1.5, adj = c(0,0))

  #########################################################################################################################
library(raster)
library(rgdal)
library(stringr)
library(scales)

xy <- SpatialPoints(cbind(dist$Longitude...processed,dist$Latitude...processed))
dist$bio10 <- extract(biodat.oz10/10,xy,method="bilinear",fun=mean, buffer=15000)
dist$bio1 <-extract(biodat.oz1/10,xy,method="bilinear",fun=mean, buffer=15000)
#add new coordinates to brassiana
seed_atsc$lat[seed_atsc$Species == "brassiana"] <- -14.2107
seed_atsc$lon[seed_atsc$Species == "brassiana"] <- 144.1304

#subset data per species
cam<- subset(seed_atsc, Taxa == "camaldulensis")
cam$Taxa <- factor(cam$Taxa)
camdist <- subset(dist, Species...matched == "Eucalyptus camaldulensis")
camdist$Species...matched <- factor(camdist$Species...matched)
ter<- subset(seed_atsc, Taxa == "tereticornis")
ter$Taxa <- factor(ter$Taxa)
terdist <- subset(dist, Species...matched == "Eucalyptus tereticornis")
terdist$Species...matched <- factor(terdist$Species...matched)

bra<- subset(seed_atsc, Taxa == "brassiana")
bra$Taxa <- factor(bra$Taxa)
bradist <- subset(dist, Species...matched == "Eucalyptus brassiana")
bradist$Species...matched <- factor(bradist$Species...matched)
pel<- subset(seed_atsc, Taxa == "pellita")
pel$Taxa <- factor(pel$Taxa)
peldist <- subset(dist, Species...matched == "Eucalyptus pellita")
peldist$Species...matched <- factor(peldist$Species...matched)
pla<- subset(seed_atsc, Taxa == "platyphylla")
pla$Taxa <- factor(pla$Taxa)
pladist <- subset(dist, Species...matched == "Eucalyptus platyphylla")
pladist$Species...matched <- factor(pladist$Species...matched)

bot<- subset(seed_atsc, Taxa == "botryoides")
bot$Taxa <- factor(bot$Taxa)
botdist <- subset(dist, Species...matched == "Eucalyptus botryoides")
botdist$Species...matched <- factor(botdist$Species...matched)
lon<- subset(seed_atsc, Taxa == "longifolia")
lon$Taxa <- factor(lon$Taxa)
londist <- subset(dist, Species...matched == "Eucalyptus longifolia")
londist$Species...matched <- factor(londist$Species...matched)
smi<- subset(seed_atsc, Taxa == "smithii")
smi$Taxa <- factor(smi$Taxa)
smidist <- subset(dist, Species...matched == "Eucalyptus smithii")
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
       ylab="", main="",axes=F)
  axis(1, at=seq(14,32,1), labels=F, line=-0.5)
  axis(1, at=seq(15,32,5), line=-0.5)
  #axis(2, at=seq(0,1000,200), line=-1.5)
  mtext(text=expression(Mean~temperature~of~warmest~quarter~(degree*C)),side=1,outer=T,line=-3,cex=1)
  polygon(widedensity, col=alpha("black",0.4), border= "black")
  polygon(narrowNdensity, col=alpha("white",0.5), border= "black")
  polygon(narrowSdensity, col=alpha("white",0.5), border= "black")
  arrows(18.85448,30,22.29183,30,code=2, lwd=2, length=0.1, angle=20)
  arrows(27.45121,30,30.93521,30,code=2, lwd=2, length=0.1, angle=20)
legend(27.5,800,legend=c("Wide","Narrow"), fill=c(alpha("black",0.4),alpha("white",0.5)),border= "black")



#pick random subsample and see if it makes the figure prettier
narrowSdist<- narrowSdist[sample(1:nrow(narrowSdist),1000,replace=F),]
narrowNdist<- narrowNdist[sample(1:nrow(narrowNdist),1000,replace=F),]
widedist<- widedist[sample(1:nrow(widedist),4000,replace=F),]

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

windows(width=6, height=6)
par(mar=c(5,3,1,0))
plot(widehist, border="white", ylim=c(0,800), xlim=c(13,33), xlab="",
     ylab="", main="",axes=F,line=1)
axis(1, at=seq(14,32,1), labels=F, line=-0.5)
axis(1, at=seq(15,32,5), line=-0.5)
axis(2, at=seq(0,800,200),line=-1.5)
mtext(text=expression(Mean~temperature~of~warmest~quarter~(degree*C)),side=1,outer=T,line=-3,cex=1.3)
polygon(widedensity, col=alpha("white",0.4), border= "black", lwd=2)
polygon(narrowNdensity, col=alpha("red",0.5), border= "black")
polygon(narrowSdensity, col=alpha("dodgerblue",0.5), border= "black")
mtext("Frequency", side=2, line=1.3, cex=1.5)
arrows(18.85448,30,22.29183,30,code=2, lwd=2.5, length=0.15, angle=20)
arrows(27.45121,30,30.93521,30,code=2, lwd=2.5, length=0.15, angle=20)
legend("topright",legend=c("2 Wide","3 Temperate Narrow", "3 Tropical Narrow"), 
       fill=c(alpha("white",0.4),alpha("dodgerblue",0.5),alpha("red",0.5)),
       border= "black", bty="n", cex=1.2)

windows(width=5, height=5)
par(mar=c(5,3,1,0))
plot(narrowNdensity, border="white", ylim=c(0,800), xlim=c(13,33), xlab="",
     ylab="", main="",axes=F,line=1)
axis(1, at=seq(14,32,1), labels=F, line=-0.5)
axis(1, at=seq(15,32,5), line=-0.5)
axis(2, at=seq(0,800,200),line=-1.5)
mtext(text=expression(Mean~temperature~of~warmest~quarter~(degree*C)),side=1,outer=T,line=-3,cex=1)


widehist <- hist(widedist$bio1)
widemultiplier <- widehist$counts / widehist$density
widedensity <- density(widedist$bio1, na.rm=T)
widedensity$y <- widedensity$y * widemultiplier[1]

narrowNhist <- hist(narrowNdist$bio1)
narrowNmultiplier <- narrowNhist$counts / narrowNhist$density
narrowNdensity <- density(narrowNdist$bio1, na.rm=T)
narrowNdensity$y <- narrowNdensity$y * narrowNmultiplier[1]

narrowShist <- hist(narrowSdist$bio1)
narrowSmultiplier <- narrowShist$counts / narrowShist$density
narrowSdensity <- density(narrowSdist$bio1, na.rm=T)
narrowSdensity$y <- narrowSdensity$y * narrowSmultiplier[1]

windows(width=5, height=5)
par(mar=c(5,3,1,0))
plot(widehist, border="white", ylim=c(0,800), xlim=c(5,30), xlab="",
     ylab="", main="",axes=F,line=1)
axis(1, at=seq(5,30,1), labels=F, line=-0.5)
axis(1, at=seq(5,30,5), line=-0.5)
axis(2, at=seq(0,800,200),line=-1.5)
mtext(text=expression(Mean~annual~temperature~(degree*C)),side=1,outer=T,line=-3,cex=1)
polygon(widedensity, col=alpha("white",0.4), border= "black")
polygon(narrowNdensity, col=alpha("red",0.5), border= "black")
polygon(narrowSdensity, col=alpha("dodgerblue",0.5), border= "black")
mtext("Frequency of occurrance", side=2, line=1.5)
arrows(18.85448,30,22.29183,30,code=2, lwd=2.5, length=0.15, angle=20)
arrows(27.45121,30,30.93521,30,code=2, lwd=2.5, length=0.15, angle=20)
legend(24,700,legend=c("Wide","Narrow"), fill=c(alpha("black",0.4),alpha("white",0.5)),border= "black", bty="n", cex=1.5)


bothist <- hist(botdist$bio1)
botmultiplier <- bothist$counts / bothist$density
botdensity <- density(botdist$bio1, na.rm=T)
botdensity$y <- botdensity$y * botmultiplier[1]

botdist$bio1future<- botdist$bio1 +3.5

bothist2 <- hist(botdist$bio1future)
botmultiplier2 <- bothist2$counts / bothist2$density
botdensity2 <- density(botdist$bio1future, na.rm=T)
botdensity2$y <- botdensity2$y * botmultiplier2[1]

windows(width=5, height=5)
par(mar=c(5,3,2,1))
plot(bothist, border="white", ylim=c(0,1000), xlim=c(10,25), xlab="",
     ylab="", main="",axes=F,line=1)
#polygon(botdensity2, col=alpha("red",0.5), border= "black")
axis(1, at=seq(10,25,1), labels=F, line=-0.5)
axis(1, at=seq(10,25,5), line=-0.5)
axis(2, at=seq(0,1000,500),line=-1.5)
#arrows(15.5,250,19,250,code=2, lwd=2.5, length=0.15, angle=20)
mtext(text=expression(Mean~Annual~Temperature~(degree*C)),side=1,outer=T,line=-3,cex=1.2)
#mtext(text=expression(italic("E. botryoides")), side=3, line=-3,adj=0.2,cex=1)
polygon(botdensity, col=alpha("white",0.4), border= "black")
mtext("Frequency of occurrence", side=2, line=1.5, cex=1.2)

windows(3,5)
par(mfrow=c(1,1), oma=c(1,1,1,0))
plot(biodat.oz1/10,main="",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F, col= "grey")
points(x=botdist$Longitude...processed,y=botdist$Latitude...processed,col=alpha("red",0.3), bg=alpha("red",0.3),cex=1,pch=21)
mtext("Latitude", side=2, outer=F, line=2.5)
mtext("Longitude", side=1, outer=F, line=2.5)

