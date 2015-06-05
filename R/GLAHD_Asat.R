#------------------------------------------------------------------------------------------------------------------------------
#- Reads in and processes the leaf gas exchange data
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")


#- read data
gx <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Asat/GHS39_GLAHD_MAIN_Asat_04122014_L1.csv")

#- get the first bit of the code (the taxa)
gx$Taxa <- unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=1,to=2071,by=2)]
gx$Taxa <- factor(gx$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
gx$Pot <- as.numeric(unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=2,to=2072,by=2)])
gx$Treat <- as.factor(ifelse(gx$Pot < 20, "Home","Warmed"))
gx$Location <- as.factor(ifelse(gx$Taxa == "ACAM" | gx$Taxa == "BCAM"| gx$Taxa == "CCAM"| 
                                   gx$Taxa == "ATER"| gx$Taxa == "BTER"| gx$Taxa == "LONG"| 
                                   gx$Taxa == "SMIT" | gx$Taxa == "BOT", "S","N"))
gx$Range <- as.factor(ifelse(gx$Taxa == "LONG"| gx$Taxa == "SMIT" |gx$Taxa == "BOT"| gx$Taxa == "PEL"
                              | gx$Taxa == "PLAT"|gx$Taxa =="BRA","Narrow","Wide"))
#- average over sub-replicate logs
gx2 <- summaryBy(.~Code+Taxa+Treat+Location+Range,data=gx,keep.names=T)



boxplot(VpdL~Treat+Taxa,data=gx2,ylab=expression(VPD~(kPa)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)


#---------------------------
#-- plot photosynthesis...
windows(14,14);par(mfrow=c(2,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")
#- home photo
ylims=c(18,35)
boxplot(Photo~Treat+Taxa,data=gx2,ylim=ylims,ylab=expression(A[sat]~(mu*mol~m^-2~s^-1)),axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)


#- ... and conductance
boxplot(Cond~Treat+Taxa,data=gx2,ylab=expression(g[s]~(mol~m^-2~s^-1)),axes=F,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Asat_results.pdf")

#---------------------------

L1<-lm((Photo)^2~Taxa*Treat, data=gx) # The warming effect on Asat is taxa specific
anova(L1)
summary(L1)

L2<-lm(Photo~Range*Treat*Location, data=gx) 
L2<-lm(Photo~Range+Treat, data=gx) 
anova(L2)
summary(L2)

L3<-lm((Photo)^2~Range+Treat, data=gx)
anova(L3)
summary(L3)
#No difference in Asat in North and South
#but warmed has lower Asat than home and wide has higher than narrow

#---------------------------


#####plot Asat at Tleaf ###

#average Photo and Tleaf for warmed and home for each taxa
gx3 <- summaryBy(Photo+Tleaf~Taxa+Treat, data=gx2, FUN=c(mean, standard.error))


#- split into list across all taxa for plotting
dat <- split(gx3,gx3$Taxa)

#plot mean Photo vs. mean Tleaf for home and warmed, for each taxa

windows(30,30)
par(mfrow=c(5,4), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(18,35)
xlims <- c(20,37)
palette(c("black","red"))
for (i in 1:length(dat)){
  toplot <- dat[[i]]
  plotBy(Photo.mean~Tleaf.mean|Treat,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="Asat",xlab="Tleaf",legend=F,
         panel.first=adderrorbars(x=toplot$Tleaf.mean,y=toplot$Photo.mean,
                                  SE=toplot$Photo.standard.error,direction="updown"))
  
  mtext(text=paste(toplot$Taxa),side=1,line=-1.5)
  if (i==length(dat)){
    mtext(text="Asat",side=2,outer=T,line=2,cex=2)
    mtext(text="Tleaf",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Asat_@Tleaf.pdf")

### Process and plot narrow vs. wide in N vs. S

#- average across location (N vs. S), treatment, and range (wide vs. narrow)

gx4 <- summaryBy(Photo+Tleaf~Location+Treat+Range, data=gx2, FUN=c(mean, standard.error))
dat.l <- split(gx4,as.factor(paste(gx4$Location,gx4$Range,sep="-")))

windows(30,30)
par(mfrow=c(2,2), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(18,35)
xlims <- c(20,37)
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(Photo.mean~Tleaf.mean|Treat,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="Asat",xlab="Tleaf",legend=F,
         panel.first=adderrorbars(x=toplot$Tleaf.mean,y=toplot$Photo.mean,
                                  SE=toplot$Photo.standard.error,direction="updown"))

  mtext(text=paste(toplot$Location,toplot$Range,sep="-"),side=1,line=-1.5)
  if (i==length(dat.l)){
    mtext(text="Asat",side=2,outer=T,line=2,cex=1.8)
    mtext(text="Tleaf",side=1,outer=T,line=2,cex=1.8)
    
  }
}

#get Asat values from ACI curves

aci1 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day1_compiled.csv")
aci2 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day2_compiled.csv")
aci3 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day3_compiled.csv")
aci <- rbind(aci1,aci2,aci3)
aci$Taxa <- unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=1,to=nrow(aci)*2,by=2)]
aci$Taxa <- factor(aci$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
aci$Pot <- as.numeric(unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=2,to=nrow(aci)*2,by=2)])
aci$Treat <- as.factor(ifelse(aci$Pot < 20, "ACiHome","ACiWarmed"))
aci$Location <- as.factor(ifelse(aci$Taxa == "ACAM" | aci$Taxa == "BCAM"| aci$Taxa == "CCAM"| 
                                   aci$Taxa == "ATER"| aci$Taxa == "BTER"| aci$Taxa == "LONG"| 
                                   aci$Taxa == "SMIT" | aci$Taxa == "BOT", "S","N"))
aci$Range <- as.factor(ifelse(aci$Taxa == "LONG"| aci$Taxa == "SMIT" |aci$Taxa == "BOT"| aci$Taxa == "PEL"
                              | aci$Taxa == "PLAT"|aci$Taxa =="BRA","Narrow","Wide"))

ambientco2 <- subset(aci,aci$CO2R <500 & aci$CO2R>400)

############
#try to improve subsetting Asat data
#write.csv(aci, "W:/WorkingData/GHS39/GLAHD/Varhammar_A/aci.csv")#manually mark what to fit
data<- read.csv("W:/WorkingData/GHS39/GLAHD/Varhammar_A/aci.csv")
ambientco2<- droplevels(subset(data, toFit == "1"))
#############

asci <- summaryBy(.~Code+Taxa+Treat+Location+Range,data=ambientco2,keep.names=T)

#Combine all values of photo and Tleaf per taxa, treatment and location
asc<- asci[,!(names(asci) %in% "toFit")]
asat <- rbind(asc,gx2)

asatS <- droplevels(subset(asat, asat$Location == "S"))
asS<- droplevels(subset(asatS, asatS$Treat != "ACiWarmed"))
asatN <- droplevels(subset(asat, asat$Location == "N"))
asN<- droplevels(subset(asatN, asatS$Treat != "ACiWarmed"))
#reorder factor levels
asS$Treat2 <- factor(asS$Treat, levels= c("Home","ACiHome","Warmed"))
asN$Treat2 <- factor(asN$Treat, levels= c("Home","ACiHome","Warmed"))


windows(30,15);par(mfrow=c(1,2),oma=c(6,5,2,1),mar=c(0,0,0,0), cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","white","red")
#- home photo
ylims=c(13,35)
boxplot(Photo~Treat2+Taxa,data=asS,ylim=ylims,axes=F,las=2,col=colors,
        at =c(2,3,4, 6,7,8, 10,11,12, 14,15,16, 18,19,20, 22,23,24, 26,27,28, 30,31,32))
magaxis(c(2,4),labels=c(1,0),box=T)
axis(side=1,at=seq(from=3,to=31,by=4),labels=levels(asS$Taxa),las=2)

boxplot(Photo~Treat2+Taxa,data=asN,ylim=ylims,axes=F,las=2,col=colors,
        at =c(2,3,4, 6,7,8, 10,11,12, 14,15,16, 18,19,20, 22,23,24, 26,27,28, 30,31,32, 34,35,36))
magaxis(c(2,4),labels=F,box=T)
axis(side=1,at=seq(from=3,to=35,by=4),labels=levels(asN$Taxa),las=2)
legend("topright",legend=c("Home","Warmed - no acc","Warmed - acc"),col="black", pt.bg=colors, bg="white", pch=22)
mtext(text=expression(A[sat]~(mu*mol~m^-2~s^-1)), side=2, outer=T,line=2,cex=1.5)

dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Asat+AsatfromACi_@Tleaf.pdf")


######## Ci/Ca ratios#######
windows(14,14);par(mfrow=c(2,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
boxplot(Ci.Ca~Treat+Taxa,data=gx2,ylab="Ci/Ca (Asat measurements)",axes=F,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)

boxplot(Ci.Ca~Treat+Taxa,data=ambientco2,ylab="Ci/Ca (Aci curves)",axes=F,col=colors)
magaxis(c(2,3,4),labels=c(1,0,0),box=T)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(gx2$Taxa),las=2)
abline(v=16.5)

dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/CiCaforAsat.pdf")



######Asat:Rd ratios##########

