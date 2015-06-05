  setwd("W:/WorkingData/GHS39/GLAHD/Varhammar_A/")
  #load some libraries
  library(raster)
  library(rgdal)
  library(stringr)
  library(scales)
  
  
  #load file
  seed_atsc<- read.csv("GLAHDseed.csv")
  
  #load known distribution
  d <- read.csv("distributions_bioclimv2.csv")
  dist <- subset(d,d$Coordinate.Uncertainty.in.Metres...parsed<5000)
  #plot final
  
  cam<- subset(seed_atsc, Taxa == "camaldulensis")
  camdist <- subset(dist, Species...matched == "Eucalyptus camaldulensis")
  ter<- subset(seed_atsc, Taxa == "tereticornis")
  terdist <- subset(dist, Species...matched == "Eucalyptus tereticornis")
                    
  bra<- subset(seed_atsc, Taxa == "brassiana")
  bradist <- subset(dist, Species...matched == "Eucalyptus brassiana")
  pel<- subset(seed_atsc, Taxa == "pellita")
  peldist <- subset(dist, Species...matched == "Eucalyptus pellita")
  pla<- subset(seed_atsc, Taxa == " platyphylla")
  pladist <- subset(dist, Species...matched == "Eucalyptus platyphylla")
  
  bot<- subset(seed_atsc, Taxa == "botryoides")
  botdist <- subset(dist, Species...matched == "Eucalyptus botryoides")
  lon<- subset(seed_atsc, Taxa == "longifolia")
  londist <- subset(dist, Species...matched == "Eucalyptus longifolia")
  smi<- subset(seed_atsc, Taxa == "smithii")
  smidist <- subset(dist, Species...matched == "Eucalyptus smithii")
  
  #load data
  biodat <- getData("worldclim", var="bio", res=2.5, path="//ad.uws.edu.au/dfshare/HomesHWK$/30034792/My Documents/Projects/ARC/Seed/T data")
  
  #subset data
  #biodat1 <- subset(biodat,1)
  #biodat4 <- subset(biodat,4)
  biodat5 <- subset(biodat,5)
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
  
  YbrevRange <- extent(141.00, 154.00, -44, -10.0)
  #biodat.oz1 <- crop(biodat1,YbrevRange)
  #biodat.oz4 <- crop(biodat4,YbrevRange)
  biodat.oz5 <- crop(biodat5,YbrevRange)
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

#plot provenances on map
windows(9,12)

plot(biodat.oz10/10,main="Seed Provenances",xlim=c(144.1,144.3),ylim=c(-14.4,-14))

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
#Remove E. tereticornis ssp. mediana
terdist<- subset(terdist, Latitude...processed >-37.6333)

windows(20,33)

par(mfrow=c(3,3), mar=c(2,2,2,1), oma=c(2,1,1,6))
plot(biodat.oz10/10,xlim=c(141,155), ylim=c(-44,-11),legend=F)

points(x=bra$lon,y=bra$lat,col="black", bg="yellow",cex=1.75,pch=21)
points(x=pel$lon,y=pel$lat,col="black", bg="orange",cex=1.75,pch=21)
points(x=145.15,y=-16.58,col="black", bg="red",cex=1.75,pch=21)

points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21)
points(x=lon$lon,y=lon$lat,col="black", bg="cyan",cex=1.75,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=1.75,pch=21)

points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=1.5,pch=21)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=1.5,pch=21)

plot(biodat.oz10/10,main="E. camaldulensis",cex.main=1.5,font.main=3, xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=camdist$Longitude...processed,y=camdist$Latitude...processed,col=alpha("grey85",0.3), bg=alpha("grey85",0.3),cex=1,pch=21)
points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=1.5,pch=21)

plot(biodat.oz10/10,main="E. tereticornis",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=terdist$Longitude...processed,y=terdist$Latitude...processed,col=alpha("grey85",0.3), bg=alpha("grey85",0.3),cex=1,pch=21)
points(x=ter$lon,y=ter$lat,col="black", bg="white",cex=1.5,pch=21)

plot(biodat.oz10/10,main="E. platyphylla",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=pladist$Longitude...processed,y=pladist$Latitude...processed,col=alpha("red",0.3), bg=alpha("red",0.3),cex=1,pch=21)
points(x=145.15,y=-16.58,col="black", bg="red",cex=1.75,pch=21)

plot(biodat.oz10/10,main="E. brassiana",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=bradist$Longitude...processed,y=bradist$Latitude...processed,col=alpha("yellow",0.3), bg=alpha("yellow",0.3),cex=1,pch=21)
points(x=bra$lon,y=bra$lat,col="black", bg="yellow",cex=1.75,pch=21)

plot(biodat.oz10/10,main="E. pellita",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=bradist$Longitude...processed,y=bradist$Latitude...processed,col=alpha("orange",0.3), bg=alpha("orange",0.3),cex=1,pch=21)
points(x=pel$lon,y=pel$lat,col="black", bg="orange",cex=1.75,pch=21)

plot(biodat.oz10/10,main="E. botryoides",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=botdist$Longitude...processed,y=botdist$Latitude...processed,col=alpha("dodgerblue",0.3), bg=alpha("dodgerblue",0.3),cex=1,pch=21)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21)

plot(biodat.oz10/10,main="E. longifolia",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=londist$Longitude...processed,y=londist$Latitude...processed,col=alpha("cyan",0.3), bg=alpha("cyan",0.3),cex=1,pch=21)
points(x=lon$lon,y=lon$lat,col="black", bg="cyan",cex=1.75,pch=21)

plot(biodat.oz10/10,main="E. smithii",cex.main=1.5 ,font.main=3,xlim=c(141,155), ylim=c(-44,-11),legend=F)
points(x=smidist$Longitude...processed,y=smidist$Latitude...processed,col=alpha("purple",0.3), bg=alpha("purple",0.3),cex=1,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=1.75,pch=21)

par(mfrow=c(1, 2), mar=c(0, 0, 0, 1), oma=c(0, 0, 0, 2),new=FALSE, xpd=TRUE)
plot(biodat.oz10/10, legend.only=TRUE,legend.shrink=0.75,legend.width=1.5,
     axis.args=list(at=seq(10,30,5),labels=seq(10, 30, 5),  cex.axis=0.8),
     legend.args=list(text='Mean Temperature of Warmest Quarter', side=4, font=2, line=2, cex=0.8))

dev.copy2pdf(file="Seed Provenances + distribution2.pdf")


#plot heatwave provenances on map
windows(6,8)

plot(biodat.oz10/10,main="Seed Provenances",xlim=c(144.00, 154.00),ylim=c(-44, -32.0))

points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=1.75,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=1.75,pch=21)

points(x=146.47,y=-36.36,col="black", bg="black",cex=1.75,pch=21)
points(x=150.07,y=-35.4,col="black", bg="white",cex=1.75,pch=21)

legend("bottomright",legend=c("E. camaldulensis","E. tereticornis", NA,"E. botryoides", "E. smithii")
       ,col=c("black","black",NA,"black","black"), 
       pt.bg=c("black","white",NA,"dodgerblue","purple"),pch=21,cex=1.2, pt.cex=1.2, bg="white")

dev.copy2pdf(file="Seed Provenances Heatwave.pdf")

#Plot individual species Heatwave

windows(18,17)

par(mfrow=c(2,3), mar=c(2,2,2,1), oma=c(2,1,1,6))
plot(biodat.oz10/10,xlim=c(145,155), ylim=c(-44,-32),legend=F)
points(x=146.47,y=-36.36,col="black", bg="black",cex=2,pch=21)
points(x=150.07-0.1,y=-35.4-0.1,col="black", bg="white",cex=2,pch=21)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=2,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=2,pch=21)

plot(biodat.oz10/10,main="E. camaldulensis",cex.main=1.5,font.main=3, xlim=c(145,155), ylim=c(-44,-32),legend=F)
points(x=camdist$Longitude...processed,y=camdist$Latitude...processed,col=alpha("grey25",0.3), bg=alpha("grey",0.3),cex=1,pch=21)
points(x=146.47,y=-36.36,col="black", bg="black",cex=2,pch=21)

plot(biodat.oz10/10,main="E. tereticornis",cex.main=1.5 ,font.main=3,xlim=c(145,155), ylim=c(-44,-32),legend=F)
points(x=terdist$Longitude...processed,y=terdist$Latitude...processed,col=alpha("grey55",0.3), bg=alpha("white",0.3),cex=1,pch=21)
points(x=150.07,y=-35.4,col="black", bg="white",cex=2,pch=21)

plot(biodat.oz10/10,main="E. botryoides",cex.main=1.5 ,font.main=3,xlim=c(145,155), ylim=c(-44,-32),legend=F)
points(x=botdist$Longitude...processed,y=botdist$Latitude...processed,col=alpha("dodgerblue",0.3), bg=alpha("dodgerblue",0.3),cex=1,pch=21)
points(x=bot$lon,y=bot$lat,col="black", bg="dodgerblue",cex=2,pch=21)

plot(biodat.oz10/10,main="E. smithii",cex.main=1.5 ,font.main=3,xlim=c(145,155), ylim=c(-44,-32),legend=F)
points(x=smidist$Longitude...processed,y=smidist$Latitude...processed,col=alpha("purple",0.3), bg=alpha("purple",0.3),cex=1,pch=21)
points(x=smi$lon,y=smi$lat,col="black", bg="purple",cex=2,pch=21)

par(mfrow=c(1, 2), mar=c(0, 0, 0, 1), oma=c(0, 0, 0, 2),new=FALSE, xpd=TRUE)
plot(biodat.oz10/10, legend.only=TRUE,legend.shrink=0.75,legend.width=1.5,
     axis.args=list(at=seq(10,30,5),labels=seq(10, 30, 5),  cex.axis=0.8),
     legend.args=list(text='Mean Temperature of Warmest Quarter', side=4, font=2, line=2, cex=0.8))

dev.copy2pdf(file="Seed Provenances Heatwave + distribution.pdf")


