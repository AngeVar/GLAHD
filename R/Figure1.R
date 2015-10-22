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
YbrevRange <- extent(141.00, 154.00, -40, -6.0)
biodat.oz10 <- crop(biodat10,YbrevRange)

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

windows(5,4.5)
par(mfrow=c(1,3), mar=c(2,0,2,0),oma=c(2,4,2,0))
plot(biodat.oz10/10,main="Selected",cex.main=1.5 , xlim=c(141,154), ylim=c(-40,-6),legend=F, col=alpha("black",0.3))
points(x=cam$lon,y=cam$lat,col="black", bg="white",cex=1.5,pch=21,lwd=2)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=1.5,pch=21,lwd=2)
points(x=pla$lon,y=pla$lat,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=pel$lon,y=pel$lat,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=145.15,y=-16.58,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)
points(x=lon$lon,y=lon$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)
points(x=smi$lon,y=smi$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)


plot(biodat.oz10/10,main="Wide",cex.main=1.5, xlim=c(141,154), ylim=c(-40,-6),yaxt='n', legend=F, col=alpha("black",0.3))
plot(wpol,add=T,col="white")
points(x=cam$lon,y=cam$lat,col="black", bg="white",cex=1.5,pch=21,lwd=2)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=1.5,pch=21,lwd=2)
box()

plot(biodat.oz10/10,main="Narrow",cex.main=1.5 ,xlim=c(141,154), ylim=c(-40,-6),yaxt='n',legend=F, col=alpha("black",0.3))
plot(Nnpol,add=T,col="red")
plot(Snpol,add=T,col="dodgerblue")
points(x=145.15,y=-16.58,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=pla$lon,y=pla$lat,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=pel$lon,y=pel$lat,col="black", bg="red",cex=1.75,pch=21,lwd=2)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)
points(x=lon$lon,y=lon$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)
points(x=smi$lon,y=smi$lat,col="black", bg="dodgerblue",cex=1.75,pch=21,lwd=2)

mtext("Latitude", side=2, outer=T, line=2.5)
mtext("Longitude", side=1, outer=T, line=1)

#########################################################################################################################
library(raster)
library(rgdal)
library(stringr)
library(scales)

xy <- SpatialPoints(cbind(dist$Longitude...processed,dist$Latitude...processed))
dist$bio10 <- extract(biodat.oz10/10,xy,method="bilinear",fun=mean, buffer=15000)
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
       ylab="", main="",axes=F,line=1)
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

windows(width=5, height=5)
par(mar=c(5,3,1,0))
plot(widehist, border="white", ylim=c(0,800), xlim=c(13,33), xlab="",
     ylab="", main="",axes=F,line=1)
axis(1, at=seq(14,32,1), labels=F, line=-0.5)
axis(1, at=seq(15,32,5), line=-0.5)
axis(2, at=seq(0,800,200),line=-1.5)
mtext(text=expression(Mean~temperature~of~warmest~quarter~(degree*C)),side=1,outer=T,line=-3,cex=1)
polygon(widedensity, col=alpha("black",0.4), border= "black")
polygon(narrowNdensity, col=alpha("white",0.5), border= "black")
polygon(narrowSdensity, col=alpha("white",0.5), border= "black")
mtext("Frequency of occurrance", side=2, line=1.5)
arrows(18.85448,30,22.29183,30,code=2, lwd=2, length=0.1, angle=20)
arrows(27.45121,30,30.93521,30,code=2, lwd=2, length=0.1, angle=20)
legend(27.5,800,legend=c("Wide","Narrow"), fill=c(alpha("black",0.4),alpha("white",0.5)),border= "black")
