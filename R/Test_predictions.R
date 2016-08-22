#Predict LAR for harvested plants and compare predicted and observed

#- load libraries from script
source("R/loadLibraries.R")
source("R/gamplotfunctions.R")
library(scales)

#- read in the data, do a few conversions
allom <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
allom$logd2h <- log10(allom$d2h)
allom$logTM <- log10(allom$Totmass)
allom$logLA <- log10(allom$Leafarea)
allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                ifelse(allom$Pot>=40,"Pre","Warmed")))
allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

allom$LAR <- with(allom, Leafarea/Totmass)
allom.2 <- subset(allom,Treatment!="Pre-treatment")
#####################

allom.2$TotpredMass <- predict_mass(Taxa=allom.2$Taxa,Treat=allom.2$Treatment,Diameter=allom.2$Diameter,Height=allom.2$Height,model_flag="complex")
allom.2$predleafArea <- predict_LA(Taxa=allom.2$Taxa,Treat=allom.2$Treatment,Diameter=allom.2$Diameter,Height=allom.2$Height)
allom.2$predLAR <- with(allom.2,predleafArea/TotpredMass)


windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(TotpredMass~Totmass,data=allom.2,pch=3,col="red",cex.lab=1.5,
     xlab="Observed mass (g)",ylab="Predicted mass (g)", xlim=c(0,125), ylim=c(0,120))
abline(0,1)
#- log scale
plot(TotpredMass~Totmass,data=allom.2,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed mass (g)",ylab="Predicted mass (g)")
abline(0,1)
summary(lm(TotpredMass~Totmass,data=allom.2)) #adj r2 - 0.9241 

windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(predleafArea~Leafarea,data=allom.2,pch=3,col="red",cex.lab=1.5,
     xlab="Observed Leafarea (cm2)",ylab="Predicted Leafarea (cm2)")
abline(0,1)
#- log scale
plot(predleafArea~Leafarea,data=allom.2,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed Leafarea (cm2)",ylab="Predicted Leafarea (cm2)")
abline(0,1)
summary(lm(predleafArea~Leafarea,data=allom.2)) #adj r2 - 0.9145 

windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(predLAR~LAR,data=allom.2,pch=3,col="red",cex.lab=1.5,
     xlab="Observed LAR",ylab="Predicted LAR")
abline(0,1)
#- log scale
plot(predLAR~LAR,data=allom.2,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed LAR",ylab="Predicted LAR")
abline(0,1)
summary(lm(predLAR~LAR,data=allom.2)) #adj r2 - 0.5761
#-------------------------------------------------------------------------------------

#Are predictions worse for some taxa? - Yes, much worse

allom.3<- subset(allom.2, Species == "BRA")

windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(predLAR~LAR,data=allom.3,pch=3,col="red",cex.lab=1.5,
     xlab="Observed LAR",ylab="Predicted LAR")
abline(0,1)
#- log scale
plot(predLAR~LAR,data=allom.3,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed LAR",ylab="Predicted LAR")
abline(0,1)
summary(lm(predLAR~LAR,data=allom.3))

#CAM 0.4541 
#TER 0.6067
#BOT 0.167
#LONG 0.1443
#SMIT -0.01554
#PEL 0.4134
#PLAT 0.6111
#BRA 0.3423