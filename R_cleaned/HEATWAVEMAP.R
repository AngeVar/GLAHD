#Heatwave map
source("R_cleaned/GLAHD_LoadLibraries.R")
library(raster)
library(rgdal)
library(stringr)
library(scales)
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


#load file

biodat <- getData("worldclim", var="bio", res=2.5, path="Data/Climate/T data")
biodat10 <- subset(biodat,10)
biodat1<-subset(biodat,1)
YbrevRange <- extent(141.00, 154.00, -39.2, -11.0)
biodat.oz10 <- crop(biodat10,YbrevRange)
biodat.oz1 <- crop(biodat1,YbrevRange)

#seed_atsc<- read.csv("Data/GLAHDseed.csv")
seed_atsc<-read.csv("R:/WORKING_DATA/GHS39/GLAHD/Varhammar_A/newbrascoords.csv")
dist <- read.csv("C:/Repos/GLAHD/Data/ALAdata20151022Clean.csv")

cam<- subset(seed_atsc, Code == "ACAM");ter<- subset(seed_atsc, Code == "BTER")
bot<- subset(seed_atsc, Code == "BOT");smi<- subset(seed_atsc, Code == "SMIT")

terdist<- subset(dist, Species...matched == "Eucalyptus tereticornis"& "Longitude...processed" > 141.00)
camdist<-subset(dist, Species...matched == "Eucalyptus camaldulensis"& "Longitude...processed" > 141.00)
botdist<-subset(dist, Species...matched == "Eucalyptus botryoides"& "Longitude...processed" > 141.00)
smitdist<- subset(dist, Species...matched == "Eucalyptus smithii"& "Longitude...processed" > 141.00)

coordinates(camdist) <- c("Longitude...processed","Latitude...processed")
coordinates(terdist) <- c("Longitude...processed","Latitude...processed")
coordinates(botdist) <- c("Longitude...processed","Latitude...processed")
coordinates(smitdist) <- c("Longitude...processed","Latitude...processed")

projection(camdist) <- CRS('+proj=longlat') 
projection(terdist) <- CRS('+proj=longlat') 
projection(botdist) <- CRS('+proj=longlat') 
projection(smitdist) <- CRS('+proj=longlat') 

t <- circles(terdist, d=50000, lonlat=TRUE) 
c <- circles(camdist, d=50000, lonlat=TRUE) 
b <- circles(botdist, d=50000, lonlat=TRUE) 
s <- circles(smitdist, d=50000, lonlat=TRUE) 

tpol <- gUnaryUnion(t@polygons)
cpol <- gUnaryUnion(c@polygons)
bpol <- gUnaryUnion(b@polygons)
spol <- gUnaryUnion(s@polygons)

gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

tpol2 <- gClip(tpol, YbrevRange)
cpol2 <- gClip(cpol, YbrevRange)
bpol2 <- gClip(bpol, YbrevRange)
spol2 <- gClip(spol, YbrevRange)

windows(4.5,7.2)
par(mfrow=c(2,2), mar=c(2,0,2,0), oma=c(2,4,1,6))

plot(biodat.oz10/10, xlim=c(141,154), ylim=c(-44,-11),legend=F, xaxt='n')
plot(cpol2,col=alpha("white",0.6),xlim=c(141,154), ylim=c(-44,-11), add=T)
mtext("E. camaldulensis",3,cex=0.8,font=3)
outline <- map("worldHires", regions="Australia", exact=TRUE, plot=FALSE) # returns a list of x/y coords
xrange <- range(outline$x, na.rm=TRUE) # get bounding box
yrange <- range(outline$y, na.rm=TRUE)
xbox <- xrange + c(-2, 2)
ybox <- yrange + c(-2, 2)
subset <- !is.na(outline$x)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd")
box()
points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=1.5,pch=21)
legend("topright","a", bty='n', cex=1.2,)

plot(biodat.oz10/10, xlim=c(141,154), ylim=c(-44,-11),legend=F, xaxt='n', yaxt='n')
plot(tpol2,col=alpha("yellow",0.6),xlim=c(141,154), ylim=c(-44,-11), add=T)
mtext("E. tereticornis",3, cex=0.8,font=3)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=ter$lon,y=ter$lat,col="black", bg="black",cex=1.5,pch=21)
legend("topright","b", bty='n', cex=1.2)

plot(biodat.oz10/10,xlim=c(141,154), ylim=c(-44,-11),legend=F)
plot(bpol2,col=alpha("blue",0.6),xlim=c(141,154), ylim=c(-44,-11), add=T)
mtext("E. botryoides",3,cex=0.8,font=3)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=bot$lon,y=bot$lat,col="black", bg="black",cex=1.5,pch=21)
legend("topright","c", bty='n', cex=1.2)

plot(biodat.oz10/10,xlim=c(141,154),ylim=c(-44,-11),legend=F, yaxt='n')
plot(spol2,col=alpha("red",0.6),xlim=c(141,154), ylim=c(-44,-11), add=T)
mtext("E. smithii",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=smi$lon,y=smi$lat,col="black", bg="black",cex=1.5,pch=21)
legend("topright","d", bty='n', cex=1.2,)

mtext("Latitude", side=2, outer=T, line=2)
mtext("Longitude", side=1, outer=T, line=0.5, at=0.45)

par(mfrow=c(1, 2), mar=c(0, 0, 0, 1), oma=c(0, 0, 0, 2),new=FALSE, xpd=TRUE)
plot(biodat.oz10/10, legend.only=TRUE,legend.shrink=0.75,legend.width=1.5,
     axis.args=list(at=seq(10,30,5),labels=seq(10, 30, 5),  cex.axis=0.8),
     legend.args=list(text='Mean Temperature of Warmest Quarter', side=4, font=1, line=2, cex=0.8))

dev.copy2pdf(file="heatwave_provenances.pdf")
