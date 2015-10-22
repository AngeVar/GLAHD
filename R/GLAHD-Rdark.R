#------------------------------------------------------------------------------------------------------------------------------
#- Reads in and processes the leaf gas exchange data for dark respiration
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")

gx2 <- getRdark()

#---------------------------
#-- plot Rdark
windows(14,14);par(mfrow=c(1,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")
ylims=c(0.2,1.2)
boxplot(R~Treatment+Taxa,data=gx2,ylim=ylims,ylab=expression(R[dark]~(mu*mol~m^-2~s^-1)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)
dev.copy2pdf(file="Output/Rdark_area_results.pdf")


windows(14,14);par(mfrow=c(1,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")
ylims=c(4,18)
boxplot(Rmass~Treatment+Taxa,data=gx2,ylim=ylims,ylab=expression(R[dark]~(n*mol~g^-1~s^-1)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)
dev.copy2pdf(file="Output/Rdark_mass_results.pdf")



windows(14,14);par(mfrow=c(1,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")
ylims=c(100,300)
boxplot(SLA~Treatment+Taxa,data=gx2,ylim=ylims,ylab=expression(SLA~(cm^2~g^-1)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)
dev.copy2pdf(file="Output/Leaf_SLA.pdf")

#---------------------------










#---------------------------
#-- plot leaf size
windows(14,14);par(mfrow=c(1,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")
ylims=c(0,120)
boxplot(Area~Treatment+Taxa,data=gx2,ylim=ylims,ylab=expression(Leaf~size~(cm^-2)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)
dev.copy2pdf(file="Output/Leaf_size.pdf")

#---------------------------






#-----------------------------------------------------------------------------------------------------------------------------------
#- so now we know that respiration measured at a common temperature was reduced by warming across most/all provenances.
#-  but... was this homeostatic? That is, let's use the RvsT curves measured on a subset of these provenances to predict Rmass at their growth temperatures
#-  this is only possible for a subset of the taxa (ATER, BOT, CTER, BRA)

#- read in the polynomial fits output from the script "GLAHD-RvsT.R"
polyfits <- read.csv("Data/GasEx/RvsT/poly_coefs_RvsT_JED.csv")

#- merge with measured Rdark data
gx3 <- merge(polyfits,gx2,by.x=c("Taxa","Treatment"),by.y=c("Taxa","Treatment"))

#- predict Rmass at the mean nightly growth temperature. This was 15.64542 in south home, 19.07866 in south warmed,
#- 24.13100 in north home, and 27.63018 in north warmed. (see "climate_S39.R")


#- get the locaiton data in to merge with the climate data
allom <- unique(read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")[,c(3,5)])


gx4 <- merge(gx3,allom,by=c("Taxa"))
gx4$mean_airT <- ifelse(gx4$Location.x=="S" & gx4$Treatment=="Home",15.65,
                        ifelse(gx4$Location.x=="S" & gx4$Treatment=="Warmed",19.08,
                               ifelse(gx4$Location.x=="N" & gx4$Treatment=="Home",24.13,27.63)))

#- predict Rmass at the mean nightly growth temperature. Uses equation 8 from O'Sullivan (2013)
gx4$Rpred <- with(gx4,Rmass*exp(b*(mean_airT-Tleaf)+c*(mean_airT^2-Tleaf^2)))
gx4$Taxa <- factor(gx4$Taxa,levels=c("BTER","BOT","CTER","BRA"))

windows(20,30);par(mar=c(7,6,1,1),mfrow=c(2,1))
boxplot(Rmass~Treatment*Taxa,data=gx4,las=2,col=c("blue","red"),ylab="R at common T")
boxplot(Rpred~Treatment*Taxa,data=gx4,las=2,col=c("blue","red"),ylab="R at growth T")
#-----------------------------------------------------------------------------------------------------------------------------------

#Compare Rdark vs. Asat at growth temperature

#- read Asat data
sat <- read.csv(file="Data/GasEx/Asat/GHS39_GLAHD_MAIN_Asat_04122014_L1.csv")

#- get the first bit of the code (the taxa)
sat$Taxa <- unlist(strsplit(x=as.character(sat$Code),split="-"))[seq(from=1,to=2071,by=2)]
sat$Taxa <- factor(sat$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                   "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
sat$Pot <- as.numeric(unlist(strsplit(x=as.character(sat$Code),split="-"))[seq(from=2,to=2072,by=2)])
sat$Treat <- as.factor(ifelse(sat$Pot < 20, "Home","Warmed"))
sat$Location <- as.factor(ifelse(sat$Taxa == "ACAM" | sat$Taxa == "BCAM"| sat$Taxa == "CCAM"| 
                                   sat$Taxa == "ATER"| sat$Taxa == "BTER"| sat$Taxa == "LONG"| 
                                   sat$Taxa == "SMIT" | sat$Taxa == "BOT", "S","N"))
sat$Range <- as.factor(ifelse(sat$Taxa == "LONG"| sat$Taxa == "SMIT" |sat$Taxa == "BOT"| sat$Taxa == "PEL"
                             | sat$Taxa == "PLAT"|sat$Taxa =="BRA","Narrow","Wide"))
#- average over sub-replicate logs
sat2 <- summaryBy(.~Code+Taxa+Treat+Location+Range,data=sat,keep.names=T)

#merge sat and gx4
gx5 <- merge(sat,gx4,by="Taxa")
gx5$AR <- gx5$Photo.x/gx5$Rpred
gx5$Taxa <- factor(gx5$Taxa,levels=c("BTER","BOT","CTER","BRA"))
boxplot(AR~Treatment*Taxa,data=gx5,las=2,col=c("blue","red"),ylab="AR at growth T")
#Predict R at temperatures at which Asat was measured (Tleaf.x)

gx5$Rp <- with(gx5,Rmass*exp(b*(Tleaf.x-Tleaf.y)+c*(Tleaf.x^2-Tleaf.y^2)))
gx5$AsatR <- gx5$Photo.x/gx5$Rp
boxplot(AsatR~Treatment*Taxa,data=gx5,las=2,col=c("blue","red"),ylab="AR at common T")

windows(20,30);par(mar=c(7,6,1,1),mfrow=c(2,1))
boxplot(AsatR~Treatment*Taxa,data=gx5,las=2,col=c("blue","red"),ylab="Asat:R at common day growth T")
boxplot(AR~Treatment*Taxa,data=gx5,las=2,col=c("blue","red"),ylab="Day Asat:night Rd at growth T")
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/AsatRratios.pdf")


gx5.fm1<- lm(AsatR~Treatment*Range.x,data=gx5)
gx5.fm2<- lm(AR~Treatment*Range.x,data=gx5)



