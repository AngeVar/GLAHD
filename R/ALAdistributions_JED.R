setwd("W:/WorkingData/GHS39/GLAHD/Share/")
#load some libraries
library(raster)
library(rgdal)
library(stringr)
library(scales)


#load file
seed_atsc<- read.csv("./Data/GLAHDseed.csv")

#load known distribution
dist <- read.csv("./Data/spALAdata.csv")

#plot final

cam<- subset(seed_atsc, Taxa == "camaldulensis")
camdist <- subset(dist, Species == "Eucalyptus camaldulensis")
ter<- subset(seed_atsc, Taxa == "tereticornis")
terdist <- subset(dist, Species == "Eucalyptus tereticornis")
                  
bra<- subset(seed_atsc, Taxa == "brassiana")
bradist <- subset(dist, Species == "Eucalyptus brassiana")
pel<- subset(seed_atsc, Taxa == "pellita")
peldist <- subset(dist, Species == "Eucalyptus pellita")
pla<- subset(seed_atsc, Taxa == " platyphylla")
pladist <- subset(dist, Species == "Eucalyptus platyphylla")

bot<- subset(seed_atsc, Taxa == "botryoides")
botdist <- subset(dist, Species == "Eucalyptus botryoides")
lon<- subset(seed_atsc, Taxa == "longifolia")
londist <- subset(dist, Species == "Eucalyptus longifolia")
smi<- subset(seed_atsc, Taxa == "smithii")
smidist <- subset(dist, Species == "Eucalyptus smithii")

#load data
biodat <- getData("worldclim", var="bio", res=2.5, path="W:/WorkingData/GHS39/GLAHD/Share/Data")

#subset data
#biodat1 <- subset(biodat,1)
#biodat4 <- subset(biodat,4)
#biodat5 <- subset(biodat,5)
#biodat6 <- subset(biodat,6)
#biodat8 <- subset(biodat,8)
#biodat9 <- subset(biodat,9)
biodat10 <- subset(biodat,10)
#biodat12 <- subset(biodat,12)
#biodat13 <- subset(biodat,13)
#biodat14 <- subset(biodat,14)
#biodat15 <- subset(biodat,15)
biodat18 <- subset(biodat,18)
#biodat19 <- subset(biodat,19)

#clip data to E Australia

YbrevRange <- extent(141.00, 154.00, -45, -30.0)
#biodat.oz1 <- crop(biodat1,YbrevRange)
#biodat.oz4 <- crop(biodat4,YbrevRange)
#biodat.oz5 <- crop(biodat5,YbrevRange)
#biodat.oz6 <- crop(biodat6,YbrevRange)
#biodat.oz8 <- crop(biodat8,YbrevRange)
#biodat.oz9 <- crop(biodat9,YbrevRange)
biodat.oz10 <- crop(biodat10,YbrevRange)
#biodat.oz12 <- crop(biodat12,YbrevRange)
#biodat.oz13 <- crop(biodat13,YbrevRange)
#biodat.oz14 <- crop(biodat14,YbrevRange)
#biodat.oz15 <- crop(biodat15,YbrevRange)
biodat.oz18 <- crop(biodat18,YbrevRange)
#biodat.oz19 <- crop(biodat19,YbrevRange)

#write rasters
#writeRaster(biodat.oz1,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz1.grd",overwrite=T)
#writeRaster(biodat.oz4,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz4.grd",overwrite=T)
#writeRaster(biodat.oz5,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz5.grd",overwrite=T)
#writeRaster(biodat.oz6,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz6.grd",overwrite=T)
#writeRaster(biodat.oz8,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz8.grd",overwrite=T)
#writeRaster(biodat.oz9,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz9.grd",overwrite=T)
#writeRaster(biodat.oz10,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz10.grd",overwrite=T)
#writeRaster(biodat.oz12,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz12.grd",overwrite=T)
#writeRaster(biodat.oz13,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz13.grd",overwrite=T)
#writeRaster(biodat.oz14,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz14.grd",overwrite=T)
#writeRaster(biodat.oz15,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz15.grd",overwrite=T)
#writeRaster(biodat.oz18,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz18.grd",overwrite=T)
#writeRaster(biodat.oz19,filename="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data/bioclim/biodat_oz19.grd",overwrite=T)

#plot provenances on map
windows(9,12)

plot(biodat.oz10/10,main="Seed Provenances",xlim=c(141,155))

points(x=bra$lon,y=bra$lat,col="black", bg="yellow",cex=1.75,pch=21)
points(x=pel$lon,y=pel$lat,col="black", bg="orange",cex=1.75,pch=21)
points(x=145.15,y=-16.58,col="black", bg="red",cex=1.75,pch=21)

points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21)
points(x=lon$lon,y=lon$lat,col="black", bg="cyan",cex=1.75,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=1.75,pch=21)

points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=1.5,pch=21)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=1.5,pch=21)

legend("topright",legend=c("E. camaldulensis","E. tereticornis", NA, "E. brassiana","E. pellita","E. platyphylla",NA,"E. botryoides","E. longifolia", "E. smithii")
       ,col=c("black","black",NA,"black","black","black",NA,"black","black","black"), 
       pt.bg=c("black","white",NA,"yellow","orange","red",NA,"dodgerblue","cyan","purple"),pch=21,cex=0.75, pt.cex=1.2, bg="white")

dev.copy2pdf(file="Seed Provenances.pdf")

#Plot individual species

windows(50,15)

par(mfrow=c(1,5), mar=c(2,2,2,1), oma=c(2,1,1,6))
plot(biodat.oz10/10,main="E. camaldulensis",cex.main=1.5,font.main=3, xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=camdist$Longitude...processed,y=camdist$Latitude...processed,col=alpha("red",0.3), bg=alpha("grey85",0.3),cex=1,pch=21)
points(x=cam$lon,y=cam$lat,col="black", bg="red",cex=2.5,pch=21)

plot(biodat.oz10/10,main="E. tereticornis",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=terdist$Longitude...processed,y=terdist$Latitude...processed,col=alpha("grey85",0.3), bg=alpha("grey85",0.3),cex=1,pch=21)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=2.5,pch=21)

plot(biodat.oz10/10,main="E. botryoides",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=botdist$Longitude...processed,y=botdist$Latitude...processed,col=alpha("dodgerblue",0.3), bg=alpha("dodgerblue",0.3),cex=1,pch=21)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=2.55,pch=21)

plot(biodat.oz10/10,main="E. longifolia",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=londist$Longitude...processed,y=londist$Latitude...processed,col=alpha("cyan",0.3), bg=alpha("cyan",0.3),cex=1,pch=21)
points(x=lon$lon,y=lon$lat,col="black", bg="cyan",cex=2.5,pch=21)

plot(biodat.oz10/10,main="E. smithii",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=smidist$Longitude...processed,y=smidist$Latitude...processed,col=alpha("purple",0.3), bg=alpha("purple",0.3),cex=1,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=2.5,pch=21)

par(mfrow=c(1, 2), mar=c(0, 0, 0, 1), oma=c(0, 0, 0, 2),new=FALSE, xpd=TRUE)
plot(biodat.oz10/10, legend.only=TRUE,legend.shrink=0.75,legend.width=1.5,
     axis.args=list(at=seq(10,30,5),labels=seq(10, 30, 5),  cex.axis=0.8),
     legend.args=list(text='Mean T, warmest quarter', side=4, font=2, line=2, cex=0.8))

dev.copy2pdf(file="Seed Provenances + distribution2.pdf")


