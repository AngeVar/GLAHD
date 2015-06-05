#- to install the newest version of plantecophys
#library(devtools)
#install_bitbucket("remkoduursma/plantecophys")

#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")

#- read data from csv files. These are rawish- they are the files right off the machine, but manually processed to remove remarks
#- and code which points to exclude from the A:Ci curve fitting.
aci1 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day1_compiled.csv")
aci2 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day2_compiled.csv")
aci3 <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/Aci/Aci_day3_compiled.csv")

aci <- rbind(aci1,aci2,aci3)

#- get the first bit of the code (the taxa)
aci$Taxa <- unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=1,to=nrow(aci)*2,by=2)]
aci$Taxa <- factor(aci$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                   "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
aci$Pot <- as.numeric(unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=2,to=nrow(aci)*2,by=2)])
aci$Treat <- ifelse(aci$Pot < 20, "Home","Warmed")

#- have a look
palette(rainbow(n=length(levels(aci$Code))))
plotBy(Photo~Ci|Code,data=subset(aci,toFit==1),legend=F,type="b")

#- extract the data to fit (exclude lines that I've manually decided to remove from the curve fitting)
aci.fit <- subset(aci,toFit==1)

#- fit the aci curves for each plant
fits <- fitacis(aci.fit,group="Code",PPFD="PARi",Tcorrect=F,citransition=450)
fits.params <- coef(fits) #extract the Vcmax and Jmax parameters

#- assign other factor variables based on the "Code".
fits.params$Taxa <- unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=1,to=nrow(fits.params)*2,by=2)]
fits.params$Taxa <- factor(fits.params$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                  "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
fits.params$Pot <- as.numeric(unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=2,to=nrow(fits.params)*2,by=2)])
fits.params$Treat <- ifelse(fits.params$Pot < 20, "Home","Warmed")

#- get summary statistics
summaryBy(Vcmax+Jmax~Taxa+Treat,data=fits.params,FUN=c(mean,sd,length))


#---------------------------------------------------------------------
#make a big pdf with all of the fits for each curve.
pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Aci fits.pdf")

for (i in 1:length(fits)){
  par(xpd=F)
  plot(fits[[i]],xlim=c(0,2000),ylim=c(-2,50))
  title(paste("ch.campaign =", fits.params[i,1]),sep="")
  par(xpd=T)
  legend("topleft",legend=paste("Vcmax=",round(fits.params[i,2],1),",","Jmax=",round(fits.params[i,3],1)),
         bty="n")
}
dev.off()


#----------------------------------------------------------------------------------
#-- plot Vcmax and Jmax for each taxa
windows(16,14);par(mfrow=c(2,1),mar=c(6,5,2,3),cex.lab=1.3,cex.axis=1.3)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")

#plot Vcmax
ylims=c(10,35)
boxplot(Vcmax~Treat+Taxa,data=fits.params,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(V["c,max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(fits.params$Taxa),las=2)
abline(v=16.4)

#plot Jmax
boxplot(Jmax~Treat+Taxa,data=fits.params,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(J["max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(fits.params$Taxa),las=2)
abline(v=16.4)
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Aci_results_Vcmax_Jmax.pdf")
#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
#- fit and plot each taxa-treatment combination
head(aci)
aci$TT <- factor(paste(aci$Taxa,aci$Treat,sep="-"))
fits.TT <- fitacis(subset(aci,toFit==1),group="TT",PPFD="PARi",Tcorrect=F,citransition=450)
fits.TT.params <- coef(fits.TT) #extract the Vcmax and Jmax parameters

windows(30,20)
xlims=c(0,1900)
ylims=c(0,50)
par(mfrow=c(6,6),mar=c(0,0,0,0),oma=c(2,2,1,1))
for (i in 1:length(fits.TT)){
  plot(fits.TT[[i]],addlegend=F,axes=F,ylim=ylims,xlim=xlims)
  magaxis(1:4,labels=rep(0,4))
  mtext(text=fits.TT.params$TT[i],side=1,line=-1.5)
  
}
#----------------------------------------------------------------------------------


#Plot J:V ratios
windows(10,5)
boxplot(Jmax/Vcmax~Treat+Taxa,data=fits.params,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(J["max"]/V["c,max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(fits.params$Taxa),las=2)
abline(v=16.4)

dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/JmaxVcmaxratios.pdf")




