#load some libraries
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

biodat <- raster::getData('worldclim', var='bio', res=2.5, path="C:/Repos/GLAHD/Data/Climate/T data")
biodat1<-subset(biodat,1) #mean annual T
#biodat5 <- subset(biodat,5)#mean max T of warmest month
#biodat6<-subset(biodat,6) #mean min T of coldest month
biodat10 <- subset(biodat,10) #mean T warmest quarter
#biodat11 <- subset(biodat,11) #mean T coldest quarter

YbrevRange <- extent(100.00, 154.00, -40, -11.0)
biodat.oz1 <- crop(biodat1,YbrevRange)
# biodat.oz5 <- crop(biodat5,YbrevRange)
# biodat.oz6 <- crop(biodat6,YbrevRange)
biodat.oz10 <- crop(biodat10,YbrevRange)
# biodat.oz11 <- crop(biodat11,YbrevRange)

seed_atsc<-read.csv("W:/WORKING_DATA/GHS39/GLAHD/Varhammar_A/newbrascoords.csv")
dist <- read.csv("C:/Repos/GLAHD/Data/ALAdata20151022Clean.csv")
#dist<- subset(dist, Country...parsed == "Australia" )#only Australian records
cam<- subset(seed_atsc, Taxa == "camaldulensis");ter<- subset(seed_atsc, Taxa == "tereticornis")
bra<- subset(seed_atsc, Taxa == "brassiana");pel<- subset(seed_atsc, Taxa == "pellita");pla<- subset(seed_atsc, Taxa == "platyphylla")
bot<- subset(seed_atsc, Taxa == "botryoides");lon<- subset(seed_atsc, Taxa == "longifolia");smi<- subset(seed_atsc, Taxa == "smithii")

coordinates(dist) <- c("Longitude...processed","Latitude...processed")
projection(dist) <- CRS('+proj=longlat') 

terdist<- subset(dist, Species...matched == "Eucalyptus tereticornis"& "Longitude...processed" > 141.00)
camdist<- subset(dist, Species...matched == "Eucalyptus camaldulensis"& "Longitude...processed" > 141.00)
bradist<- subset(dist, Species...matched == "Eucalyptus brassiana"& "Longitude...processed" > 141.00)
peldist<- subset(dist, Species...matched == "Eucalyptus pellita"& "Longitude...processed" > 141.00)
pladist<- subset(dist, Species...matched == "Eucalyptus platyphylla"& "Longitude...processed" > 141.00)
botdist<- subset(dist, Species...matched == "Eucalyptus botryoides"& "Longitude...processed" > 141.00)
londist<- subset(dist, Species...matched == "Eucalyptus longifolia"& "Longitude...processed" > 141.00)
smidist<- subset(dist, Species...matched == "Eucalyptus smithii"& "Longitude...processed" > 141.00)

t<- circles(terdist, d=50000, lonlat=TRUE)
c<- circles(camdist, d=50000, lonlat=TRUE)
br<- circles(bradist, d=50000, lonlat=TRUE)
pe<- circles(peldist, d=50000, lonlat=TRUE)
pl<- circles(pladist, d=50000, lonlat=TRUE)
b<- circles(botdist, d=50000, lonlat=TRUE)
l<- circles(londist, d=50000, lonlat=TRUE)
s<- circles(smidist, d=50000, lonlat=TRUE)

tpol<- gUnaryUnion(t@polygons)
cpol<- gUnaryUnion(c@polygons)
brpol<- gUnaryUnion(br@polygons)
pepol<- gUnaryUnion(pe@polygons)
plpol<- gUnaryUnion(pl@polygons)
bpol<- gUnaryUnion(b@polygons)
lpol<- gUnaryUnion(l@polygons)
spol<- gUnaryUnion(s@polygons)
 
gClip <- function(shp, bb){
  if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
  else b_poly <- as(extent(bb), "SpatialPolygons")
  gIntersection(shp, b_poly, byid = T)
}

tpol2 <- gClip(tpol, YbrevRange)
cpol2 <- gClip(cpol, YbrevRange)
bpol2 <- gClip(bpol, YbrevRange)
spol2 <- gClip(spol, YbrevRange)
lpol2 <- gClip(lpol, YbrevRange)
brpol2 <- gClip(brpol, YbrevRange)
pepol2 <- gClip(pepol, YbrevRange)
plpol2 <- gClip(plpol, YbrevRange)

outline <- map("worldHires", regions="Australia", exact=TRUE, plot=FALSE) # returns a list of x/y coords
outline2<-subset(outline, YbrevRange)
xrange <- range(outline$x, na.rm=TRUE) # get bounding box
yrange <- range(outline$y, na.rm=TRUE)
xbox <- xrange + c(-2, 11)
ybox <- yrange + c(-2, 2)
subset <- !is.na(outline$x)

#######################################################
windows(9,5)
par(mfrow=c(2,5),mar=c(1.5,0,0,0),oma=c(6,6,6,4))
layout(matrix(c(1:10), nrow=2, ncol=5,byrow=T),
       heights=c(1,1,1,0.01,1),
       widths=c(1,1,1,0.01,1))
#windows(7.79,11.69);par(mfrow=c(2,4),mar=c(1.5,0,0,0),oma=c(6,6,6,4))
xlims<-c(130,154)
ylims<-c(-40,-11)
##############
plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F,xaxt='n',col="white")
plot(brpol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. platyphylla",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=pla$lon,y=pla$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","a", bty='n', cex=1.2)
###################
plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F,xaxt='n',yaxt='n',col="white")
plot(pepol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. pellita",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", xlim=xlims,ylim=ylims,box=F)
box()
points(x=pel$lon,y=pel$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","b", bty='n', cex=1.2)
###################

plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F,xaxt='n', yaxt='n',col="white", bty='n')
plot(brpol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. brassiana",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=144.1304,y=-14.51,col="black", bg="black",cex=2,pch=21)
legend("topright","c", bty='n', cex=1.2)
#####################
plot(NA)

##############
plot(biodat.oz10/10, xlim=xlims, ylim=ylims,legend=F, xaxt='n', col="white")
plot(cpol,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext(text="E. camaldulensis",3, cex=0.8,font=3)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","d", bty='n', cex=1.2)
###############
plot(biodat.oz10/10,xlim=xlims, ylim=ylims,legend=F, col="white")
plot(bpol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. botryoides",3,cex=0.8,font=3)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=bot$lon,y=bot$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","e", bty='n', cex=1.2)

#######################
plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F, yaxt='n',col="white")
plot(spol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. smithii",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=smi$lon,y=smi$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","f", bty='n', cex=1.2)
#######################

plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F, yaxt='n',col="white")
plot(lpol2,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext("E. longiflora",3,cex=0.8,font=3) 
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd", box=F)
box()
points(x=lon$lon,y=lon$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","g", bty='n', cex=1.2)
#######################
plot(NA)
##################
plot(biodat.oz10/10, xlim=xlims, ylim=ylims,legend=F, col="white")
plot(tpol,col=alpha("red",0.8),border=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
mtext(text="E. tereticornis",3,cex=0.8,font=3)
polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
         c(outline$y[subset], NA, rep(ybox, each=2)),
         col="white", rule="evenodd")
box()
points(x=ter$lon,y=ter$lat,col="black", bg="black",cex=2,pch=21)
legend("topright","h", bty='n', cex=1.2)

mtext("Latitude", side=2, outer=T, line=3)
mtext("Longitude", side=1, outer=T, line=2, at=0.45)
mtext("Narrow", side=3, outer=T, line=2, at=0.3)
mtext("Wide", side=3, outer=T, line=2, at=0.9)
text(100,y=0,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.7)
text(100,y=-35,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.7)
##################################


#########################################################################################################################

#extract climate for coordinates
xy <- SpatialPoints(cbind(dist$Longitude...processed,dist$Latitude...processed))
dist$bio10 <- extract(biodat.oz10/10,xy,method="bilinear",fun=mean, buffer=15000)
#dist$bio11 <- extract(biodat.oz11/10,xy,method="bilinear",fun=mean, buffer=15000)
dist$bio1 <- extract(biodat.oz1,xy,method="bilinear",fun=mean, buffer=15000) 
# dist$bio5 <- extract(biodat.oz5/10,xy,method="bilinear",fun=mean, buffer=15000)
# dist$bio6 <- extract(biodat.oz6/10,xy,method="bilinear",fun=mean, buffer=15000)
# # write out a csv
# write.csv(dist,file="dist_bioclim.csv")
# dist<-read.csv("dist_bioclim.csv")
#extracting climate for selected provenances
y <- SpatialPoints(cbind(seed_atsc$lon,seed_atsc$lat))
seed_atsc$bio10 <- extract(biodat.oz10/10,y,method="bilinear", fun=mean,buffer=15000)
# seed_atsc$bio11 <- extract(biodat.oz11/10,y,method="bilinear", fun=mean,buffer=15000)
seed_atsc$bio1 <- extract(biodat.oz1/10,y,method="bilinear", fun=mean,buffer=15000)
# seed_atsc$bio5 <- extract(biodat.oz5/10,y,method="bilinear",fun=mean, buffer=15000)
# seed_atsc$bio6 <- extract(biodat.oz6/10,y,method="bilinear",fun=mean, buffer=15000)

seed_atsc$Species...matched<- as.factor(paste(seed_atsc$Genus,seed_atsc$Taxa, sep=" "))
#removing distribution data where location uncertainty is more than 5 km
dis <- as.data.frame(subset(dist,dist$Coordinate.Uncertainty.in.Metres...parsed<5000))

#summarise climatic conditions of distributions per species

dat<-summaryBy(bio10~Species...matched,data=as.data.frame(dis), FUN=c(mean,min,max),na.rm=T)
dat$Sp<- factor(dat$Species...matched, levels = 
                  c("Eucalyptus brassiana","Eucalyptus pellita",
                    "Eucalyptus platyphylla","Eucalyptus longifolia",
                    "Eucalyptus smithii","Eucalyptus botryoides",
                    "Eucalyptus tereticornis","Eucalyptus camaldulensis"
))

seed_atsc$Sp<- factor(seed_atsc$Species...matched, levels = 
                        c("Eucalyptus brassiana","Eucalyptus pellita",
                          "Eucalyptus platyphylla","Eucalyptus longifolia",
                          "Eucalyptus smithii","Eucalyptus botryoides",
                          "Eucalyptus tereticornis","Eucalyptus camaldulensis"
                        ))
seed_atsc$bio10.mean<- seed_atsc$bio10
####################################################
#new plot:
windows(7.79,11.69)

plot1<-ggplot(data=dat,
              aes(x=Sp,y=bio10.mean,ymax=bio10.max,
                  ymin=bio10.min))+ ylim(15,33)+geom_pointrange(shape="", col='red')

#plot provenances 
# plot1.1<- plot1+geom_point(data = seed_atsc,
#              aes(x=Sp, y=bio10.mean,ymax=bio10.mean,
#                  ymin=bio10.mean),
#              color = 'black', size=3)

#flip and add lines
plot2<-plot1+coord_flip()+geom_hline(aes(yintercept=30.7), lty=3,size=1)+
  geom_hline(aes(yintercept=27.4), lty=2,size=1)+
  geom_hline(aes(yintercept=22.4), lty=3,size=1)+
  geom_hline(aes(yintercept=18.9), lty=2,size=1)+xlab("")+ylab("")

plot2+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"), 
             axis.text=element_text(size=12))+
  ylab("Temperature of Warmest Quarter")+xlab("")+
  theme(axis.title.x=element_text(margin=margin(20,0,0,0)))+
  theme(plot.margin=unit(c(1,1,1,0),"cm"))

#add density plots
cam.df<- as.data.frame(cam)
camdist.df<- as.data.frame(camdist)
ter.df<- as.data.frame(ter)
terdist.df<- as.data.frame(terdist)
bra.df<- as.data.frame(bra)
bradist.df<- as.data.frame(bradist)
pel.df<- as.data.frame(pel)
peldist.df<- as.data.frame(peldist)
pla.df<- as.data.frame(pla)
pladist.df<- as.data.frame(pladist)
bot.df<- as.data.frame(bot)
botdist.df<- as.data.frame(botdist)
lon.df<- as.data.frame(lon)
londist.df<- as.data.frame(londist)
smi.df<- as.data.frame(smi)
smidist.df<- as.data.frame(smidist)

#density plot
a<-density(terdist.df$bio10, na.rm=T)
b<-density(camdist.df$bio10, na.rm=T)
c<-density(pladist.df$bio10, na.rm=T)
d<-density(bradist.df$bio10, na.rm=T)
e<-density(peldist.df$bio10, na.rm=T)
f<-density(botdist.df$bio10, na.rm=T)
g<-density(londist.df$bio10, na.rm=T)
h<-density(smidist.df$bio10, na.rm=T)

windows(7.79,11.69)
par(mfrow=c(8,1), mar=c(1,1,1,1), oma=c(6,8,4,4))
plot(a, ylim=c(0,0.3), main="", xlab="", xlim=c(15,33), ylab="E. tereticornis", xaxt="n")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
arrows(18.85448,0.3,22.29183,0.3,code=2, lwd=2, length=0.1, angle=20, xpd=TRUE)
arrows(27.45121,0.3,30.93521,0.3,code=2, lwd=2, length=0.1, angle=20)
plot(b, ylim=c(0,0.25), main="", xlab="", xlim=c(15,33), ylab="E. camaldulensis", xaxt="n")
polygon(b, col=alpha("black",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(c, ylim=c(0,0.5), main="", xlab="", xlim=c(15,33), ylab="E. platyphylla", xaxt="n")
polygon(c, col=alpha("red",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(d, ylim=c(0,1.2), main="", xlab="", xlim=c(15,33), ylab="E. brassiana", xaxt="n")
polygon(d, col=alpha("red",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(e, ylim=c(0,0.7), main="", xlab="", xlim=c(15,33), ylab="E. pellita", xaxt="n")
polygon(e, col=alpha("red",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(f, ylim=c(0,0.4), main="", xlab="", xlim=c(15,33), ylab="E. botryoides", xaxt="n")
polygon(f, col=alpha("blue",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(g, ylim=c(0,1.2), main="", xlab="", xlim=c(15,33), ylab="E. longifolia", xaxt="n")
polygon(g, col=alpha("blue",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
plot(h, ylim=c(0,0.5), main="", xlab="", xlim=c(15,33), ylab="E. smithii")
polygon(h, col=alpha("blue",0.3), border= "black")
abline(v=c(18.9,27.4), lty=2)
abline(v=c(22.4,30.7), lty=3)
mtext(text="Temperature of Warmest Quarter", side=1, outer=T, line=3)
mtext(text="Temperate                       Tropical", side=3, outer=T, line=1)




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

# ####################################################
# #subset data per species
# cam<- subset(seed_atsc, Taxa == "camaldulensis")
# camdist <- subset(dis, Species...matched == "Eucalyptus camaldulensis")
# ter<- subset(seed_atsc, Taxa == "tereticornis")
# terdist <- subset(dis, Species...matched == "Eucalyptus tereticornis")
# 
# bra<- subset(seed_atsc, Taxa == "brassiana")
# bradist <- subset(dis, Species...matched == "Eucalyptus brassiana")
# pel<- subset(seed_atsc, Taxa == "pellita")
# peldist <- subset(dis, Species...matched == "Eucalyptus pellita")
# pla<- subset(seed_atsc, Taxa == "platyphylla")
# pladist <- subset(dis, Species...matched == "Eucalyptus platyphylla")
# 
# bot<- subset(seed_atsc, Taxa == "botryoides")
# botdist <- subset(dis, Species...matched == "Eucalyptus botryoides")
# lon<- subset(seed_atsc, Taxa == "longifolia")
# londist <- subset(dis, Species...matched == "Eucalyptus longifolia")
# smi<- subset(seed_atsc, Taxa == "smithii")
# smidist <- subset(dis, Species...matched == "Eucalyptus smithii")
# 
# #windows(7.79,11.69);par(mfrow=c(4,2),mar=c(0,0,0,0),oma=c(6,6,6,4))
# windows(7.79,11.69);par(mfrow=c(4,2),mar=c(1.5,0,0,1.5),oma=c(6,6,6,4))
# ylims=c(12,38)
# xlims=c(0,28)
# #plot species on climate space
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus tereticornis"), 
#      col=alpha("red",0.1), xlim=c(-5,25), ylim=ylims, xaxt="n",xaxs="i",
#      pch= 16, cex=1)
# points(bio10~bio11, data=ter, col="black", pch=16)
# abline(h=c(22.4 ,30.7), lty =3)
# abline(h=c(18.9 ,27.4), lty =2)
# mtext(text="E. tereticornis",3,cex=0.8,font=3)
# axis(side = 1, labels=F, tcl=0.5 )
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus camaldulensis"), 
#      pch=16, cex=1, ylim=ylims,xlim=xlims,col=alpha("red",0.1), xaxt="n", yaxt="n",xaxs="i")
# points(bio10~bio11, data=cam, col="black", pch=16)
# abline(h=c(22.4 ,30.7), lty =3)
# abline(h=c(18.9 ,27.4), lty =2)
# axis(side=4)
# mtext(text="E. camaldulensis",3,cex=0.8,font=3)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus botryoides"), ylim=ylims,
#      pch=16, cex=1,xlim=xlims,col=alpha("red",0.1), xaxt="n",xaxs="i")
# points(bio10~bio11, data=bot, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# mtext(text="E. botryoides",3,cex=0.8,font=3)
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus platyphylla"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt="n", yaxt="n",xaxs="i")
# points(bio10~bio11, data=pla, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext(text="E. platyphylla",3,cex=0.8,font=3)
# axis(side=4)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus smithii"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1),pch=16, cex=1, xaxt="n",xaxs="i")
# points(bio10~bio11, data=smi, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# mtext(text="E. smithii",3,cex=0.8,font=3)
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus pellita"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt="n", yaxt="n",xaxs="i")
# points(bio10~bio11, data=pel, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext(text="E. pellita",3,cex=0.8,font=3)
# axis(side=4)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus longifolia"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt='n',xaxs="i")
# points(bio10~bio11, data=lon, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# mtext(text="E. longifolia",3,cex=0.8,font=3)
# #magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=T, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus brassiana"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, yaxt="n", xaxt='n',xaxs="i")
# points(bio10~bio11, data=bra, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext(text="E. brassiana",3,cex=0.8,font=3)
# axis(side=4)
# axis(side = 1, labels=T, tcl=0.5 )
# #magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)
# 
# mtext(text=expression(Temperature~of~warmest~Quarter~(degree*C)),side=2, outer=T, line=3)
# mtext(text=expression(Temperature~of~coldest~Quarter~(degree*C)),side=1, outer=T, line=2, at=0.45)
# 
########################################################################################################
#with black outline on distributions
# #Combine the two
# #left, right bottom top
# split.screen(rbind(c(0.1,0.292,0.1, 0.98), c(0.312, 0.95, 0.1, 0.98)))
# screen(1)
# windows(11.69,11.69);par(mfrow=c(4,5),mar=c(1.5,0,0,0),oma=c(6,6,6,4))
# 
# 
# xlims<-c(130,154)
# ylims<-c(-40,-11)
# plot(biodat.oz10/10, xlim=xlims, ylim=ylims,legend=F, xaxt='n',
#      col="white")
# plot(tpol,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext(text="E. tereticornis",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd")
# box()
# points(x=ter$lon,y=ter$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","a", bty='n', cex=1.2)
# 
# plot(biodat.oz10/10, xlim=xlims, ylim=ylims,legend=F, yaxt='n',xaxt='n',col="white")
# plot(cpol,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext(text="E. camaldulensis",3, cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=cam$lon,y=cam$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","b", bty='n', cex=1.2)
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# 
# ylims=c(12,38)
# xlims=c(0,28)
# #plot species on climate space
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus tereticornis"),
#      col=alpha("red",0.1), xlim=c(-5,25), ylim=ylims, xaxt="n",
#      pch= 16, cex=1)
# points(bio10~bio11, data=ter, col="black", pch=16)
# abline(h=c(22.4 ,30.7), lty =3)
# abline(h=c(18.9 ,27.4), lty =2)
# axis(side = 1, labels=F, tcl=0.5 )
# mtext(text="E. tereticornis",3,cex=0.8,font=3)
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus camaldulensis"),
#      pch=16, cex=1, ylim=ylims,xlim=xlims,col=alpha("red",0.1), yaxt='n',xaxt="n")
# points(bio10~bio11, data=cam, col="black", pch=16)
# abline(h=c(22.4 ,30.7), lty =3)
# abline(h=c(18.9 ,27.4), lty =2)
# axis(side=2)
# #legend("topleft", legend="E. camaldulensis", cex=1.3, bty="n", text.font=3)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# mtext(text="E. camaldulensis",3, cex=0.8,font=3)
# 
# 
# xlims<-c(130,154)
# ylims<-c(-40,-11)
# plot(biodat.oz10/10,xlim=xlims, ylim=ylims,legend=F, xaxt='n',col="white")
# plot(bpol2,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext("E. botryoides",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=bot$lon,y=bot$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","c", bty='n', cex=1.2)
# 
# plot(biodat.oz10/10,xlim=xlims, ylim=ylims,legend=F, xaxt='n', yaxt='n',col="white")
# plot(plpol2,col="red",xlim=xlims, ylim=ylims, add=T)
# mtext("E. platyphylla",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=pla$lon,y=pla$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","d", bty='n', cex=1.2)
# 
# ylims=c(12,38)
# xlims=c(0,28)
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus botryoides"), ylim=ylims,
#      pch=16, cex=1,xlim=xlims,col=alpha("red",0.1), xaxt="n")
# points(bio10~bio11, data=bot, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# #legend("topleft", legend="E. botryoides", cex=1.3, bty="n", text.font=3)
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# mtext("E. botryoides",3,cex=0.8,font=3)
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus platyphylla"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt="n", yaxt="n")
# points(bio10~bio11, data=pla, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext("E. platyphylla",3,cex=0.8,font=3)
# axis(side=4)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# xlims<-c(130,154)
# ylims<-c(-40,-11)
# plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F, xaxt='n',col="white")
# plot(spol2,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext("E. smithii",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=smi$lon,y=smi$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","e", bty='n', cex=1.2)
# 
# plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F,  xaxt='n',yaxt='n',col="white")
# plot(pepol2,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext("E. pellita",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=pel$lon,y=pel$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","f", bty='n', cex=1.2)
# 
# ylims=c(12,38)
# xlims=c(0,28)
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus smithii"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1),pch=16, cex=1, xaxt="n")
# points(bio10~bio11, data=smi, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# mtext("E. smithii",3,cex=0.8,font=33)
# #magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus pellita"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt="n", yaxt="n")
# points(bio10~bio11, data=pel, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext("E. pellita",3,cex=0.8,font=3)
# axis(side=4)
# #magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
# axis(side = 1, labels=F, tcl=0.5 )
# 
# xlims<-c(130,154)
# ylims<-c(-40,-11)
# plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F, col="white")
# plot(lpol2,col=alpha("red",0.8)),xlim=xlims, ylim=ylims, add=T)
# mtext("E. longiflora",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=lon$lon,y=lon$lat,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","g", bty='n', cex=1.2)
# 
# plot(biodat.oz10/10,xlim=xlims,ylim=ylims,legend=F, yaxt='n',col="white", bty='n')
# plot(brpol2,col=alpha("red",0.8),xlim=xlims, ylim=ylims, add=T)
# mtext("E. brassiana",3,cex=0.8,font=3)
# polypath(c(outline$x[subset], NA, c(xbox, rev(xbox))),
#          c(outline$y[subset], NA, rep(ybox, each=2)),
#          col="white", rule="evenodd", box=F)
# box()
# points(x=144.1304,y=-14.51,col="black", bg="black",cex=1.5,pch=21)
# legend("topright","h", bty='n', cex=1.2)
# 
# ylims=c(12,38)
# xlims=c(0,28)
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus longifolia"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, xaxt='n')
# points(bio10~bio11, data=lon, col="black", pch=16)
# abline(h=22.4 , lty =3)
# abline(h=c(18.9), lty =2)
# #legend("topleft", legend="E. longifolia", cex=1.3, bty="n", text.font=3)
# #magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)
# mtext("E. longiflora",3,cex=0.8,font=3)
# 
# plot(bio10~bio11, data=subset(dis,dis$Species...matched == "Eucalyptus brassiana"), ylim=ylims,
#      xlim=xlims,col=alpha("red",0.1), pch=16, cex=1, yaxt="n", xaxt='n')
# points(bio10~bio11, data=bra, col="black", pch=16)
# abline(h=30.7, lty =3)
# abline(h=c(27.4), lty =2)
# mtext("E. brassiana",3,cex=0.8,font=3)
# axis(side=4)
# axis(side = 1, labels=T, tcl=0.5 )








cam.df<- as.data.frame(cam)
camdist.df<- as.data.frame(camdist)
ter.df<- as.data.frame(ter)
terdist.df<- as.data.frame(terdist)
bra.df<- as.data.frame(bra)
bradist.df<- as.data.frame(bradist)
pel.df<- as.data.frame(pel)
peldist.df<- as.data.frame(peldist)
pla.df<- as.data.frame(pla)
pladist.df<- as.data.frame(pladist)
bot.df<- as.data.frame(bot)
botdist.df<- as.data.frame(botdist)
lon.df<- as.data.frame(lon)
londist.df<- as.data.frame(londist)
smi.df<- as.data.frame(smi)
smidist.df<- as.data.frame(smidist)

wide<-rbind(cam.df,ter.df)
widedist<-rbind(camdist.df,terdist.df)
narrowS<-rbind(bot.df,lon.df,smi.df)
narrowN<-rbind(bra.df, pel.df, pla.df)
narrowSdist<-rbind(botdist.df,londist.df,smidist.df)
narrowNdist<-rbind(bradist.df, peldist.df, pladist.df)

#density plot
a<-density(terdist.df$bio10, na.rm=T)
b<-density(camdist.df$bio10, na.rm=T)
c<-density(pladist.df$bio10, na.rm=T)
d<-density(bradist.df$bio10, na.rm=T)
e<-density(peldist.df$bio10, na.rm=T)
f<-density(botdist.df$bio10, na.rm=T)
g<-density(londist.df$bio10, na.rm=T)
h<-density(smidist.df$bio10, na.rm=T)

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

# #_______________________________
# #Different way to plot data on a map
# library(maps)
# library(mapdata)
# 
# etdist<- subset(dist, Species...matched == "Eucalyptus tereticornis" & Longitude...processed > 141.00 & Latitude...processed < -11)
# ecdist<- subset(dist, Species...matched == "Eucalyptus camaldulensis" & Longitude...processed > 141.00& Latitude...processed < -11)
# botdist<- subset(dist, Species...matched == "Eucalyptus botryoides" & Longitude...processed > 144.00& Latitude...processed < -11)
# smdist<- subset(dist, Species...matched == "Eucalyptus smithii" & Longitude...processed > 141.00& Latitude...processed < -11)
# 
# #one remarkable outlier record of botryoides - remove
# botdist<- subset(botdist,Latitude...processed < -25)
# 
# coordinates(etdist) <- c("Longitude...processed","Latitude...processed")
# coordinates(ecdist) <- c("Longitude...processed","Latitude...processed")
# coordinates(botdist) <- c("Longitude...processed","Latitude...processed")
# coordinates(smdist) <- c("Longitude...processed","Latitude...processed")
# 
# projection(etdist) <- CRS('+proj=longlat')
# projection(ecdist) <- CRS('+proj=longlat')
# projection(botdist) <- CRS('+proj=longlat')
# projection(smdist) <- CRS('+proj=longlat')
# 
# et <- circles(etdist, d=50000, lonlat=TRUE)
# ec <- circles(ecdist, d=50000, lonlat=TRUE)
# bt <- circles(botdist, d=50000, lonlat=TRUE)
# sm<- circles(smdist, d=50000, lonlat=TRUE)
# 
# etpol <- gUnaryUnion(et@polygons)
# ecpol <- gUnaryUnion(ec@polygons)
# botpol <- gUnaryUnion(bt@polygons)
# smpol <- gUnaryUnion(sm@polygons)
# 
# windows(10,10);par(mfrow=c(1,2), mar=c(2,0,2,0),oma=c(2,5,2,2))
# map("worldHires","Australia",xlim=c(141.00, 154.00), col="gray90", fill=T)
# plot(ecpol,add=T,col=alpha("#E0B55A", 0.6), border=F)
# plot(etpol,add=T,col=alpha("#595BDE",0.6), border=F)
# 
# points(x=subset(cam, Code=="ACAM")$lon,y=subset(cam, Code=="ACAM")$lat,
#        col="black", bg=alpha("#E0B55A",0.6),cex=2,pch=21,lwd=2)
# points(x=subset(ter, Code=="BTER")$lon,y=subset(ter, Code=="BTER")$lat,
#        col="black", bg=alpha("#595BDE",0.6),cex=2,pch=21,lwd=2)
# box()
# map("worldHires","Australia",xlim=c(141.00, 154.00), col="gray90", fill=T)
# 
# plot(botpol,add=T,col=alpha("#C74444",0.6), border=F)
# plot(smpol,add=T,col=alpha("#8541B5",0.6), border=F)
# points(x=smi$lon,y=smi$lat,col="black", bg=alpha("#C74444",0.6),cex=2,pch=21,lwd=2)
# points(x=bot$lon,y=bot$lat,col="black", bg=alpha("#8541B5",0.6),cex=2,pch=21,lwd=2)
# box()
# mtext("Latitude", side=2, outer=T, line=2.5)
# mtext("Longitude", side=1, outer=T, line=1)
# legend(0,2, c("E. camaldulensis","E.tereticornis"),
#        col=c(alpha("#E0B55A",0.6),alpha("#595BDE",0.6)))
# 
# 
# 
# 
# require(spatialEco)
# require(sp)
# data(meuse)
# coordinates(meuse) = ~x+y
# sr1=Polygons(list(Polygon(cbind(c(180114, 180553, 181127, 181477, 181294, 181007, 180409,
#                                   180162, 180114), c(332349, 332057, 332342, 333250, 333558, 333676,
#                                                      332618, 332413, 332349)))),'1')
# sr2=Polygons(list(Polygon(cbind(c(180042, 180545, 180553, 180314, 179955, 179142, 179437,
#                                   179524, 179979, 180042), c(332373, 332026, 331426, 330889, 330683,
#                                                              331133, 331623, 332152, 332357, 332373)))),'2')
# sr3=Polygons(list(Polygon(cbind(c(179110, 179907, 180433, 180712, 180752, 180329, 179875,
#                                   179668, 179572, 179269, 178879, 178600, 178544, 179046, 179110),
#                                 c(331086, 330620, 330494, 330265, 330075, 330233, 330336, 330004,
#                                   329783, 329665, 329720, 329933, 330478, 331062, 331086)))),'3')
# sr4=Polygons(list(Polygon(cbind(c(180304, 180403,179632,179420,180304),
#                                 c(332791, 333204, 333635, 333058, 332791)))),'4')
# sr=SpatialPolygons(list(sr1,sr2,sr3,sr4))
# srdf=SpatialPolygonsDataFrame(sr, data.frame(row.names=c('1','2','3','4'), PIDS=1:4, y=runif(4)))
# head(srdf@data)  # polygons
# head(meuse@data) # points
# plot(srdf)
# points(meuse, pch=20)
# pts.poly <- point.in.poly(meuse, srdf)
# head(pts.poly@data)
# srdf@data$poly.ids <- 1:nrow(srdf)
# # Number of points in each polygon
# tapply(pts.poly@data$lead, pts.poly@data$PIDS, FUN=length)
# # Mean lead in each polygon
# tapply(pts.poly@data$lead, pts.poly@data$PIDS, FUN=mean)
# 
# m1<-over(srdf,meuse)
