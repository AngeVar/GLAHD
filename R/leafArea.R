#------------------------------------------------------------------------------------------------------------------------------
#- Reads in and processes the leaf area data from the photos in the GLAHD experiment
#------------------------------------------------------------------------------------------------------------------------------

#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")


#- read data
la <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/Leaf photos/GLAHD-LeafArea.csv")
la$Date <- as.Date(la$Date,format="%d/%m/%Y")
#- get the first bit of the code (the taxa)
la$Taxa <- unlist(strsplit(x=as.character(la$Code),split="-"))[seq(from=1,to=nrow(la)*2,by=2)]
la$Taxa <- factor(la$Taxa,levels=c("BTER","ACAM","BOT","LONG","CTER","DCAM","BRA","PLAT"))
la$Pot <- as.numeric(unlist(strsplit(x=as.character(la$Code),split="-"))[seq(from=2,to=nrow(la)*2,by=2)])
la$Treat <- ifelse(la$Pot < 20, "Home","Warmed")

#- Carrie's leaves were smaller
boxplot(Area~Observer,data=la,ylab="Leaf area (cm2)")


boxplot(Area~Treat+Taxa,data=la,ylab="Leaf area (cm2)",col=c("blue","red"),axes=F,las=2)
axis(side=1,at=seq(from=1.5,to=16.5,by=2),labels=levels(la$Taxa),las=2)
magaxis(side=c(2,4),labels=c(1,1),box=T,las=1)
