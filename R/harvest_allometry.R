#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plots the biomass data for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")
source("W:/WorkingData/GHS39/GLAHD/Share/R/gamplotfunctions.R")
library(scales)

#- read in the data, do a few conversions
dat <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Totmass <- base::rowSums(dat[,11:13]) #total mass is the sum of leaf, stem, and root mass
dat$LAR <- with(dat,Leafarea/Totmass)
dat$SLA <- with(dat,Leafarea/Leafmass)
dat$d2h <- with(dat,(Diameter/10)^2*(Height)) #calculate d2h in cm3
dat$logd2h <- log10(dat$d2h)
dat$logTM <- log10(dat$Totmass)
dat$logVI <- log10(dat$Height*dat$Diameter)
dat$Treat <- as.factor(ifelse(dat$Pot < 20, "Home",
                    ifelse(dat$Pot>=40,"Pre","Warmed")))
dat$Taxa <- factor(dat$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                   "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
dat[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

#-- allometry?
windows(12,12)

#- plot taxa
palette(rainbow(length(levels(dat$Taxa))))
colors <- (rainbow(length(levels(dat$Taxa))))

with(dat,plot(logTM~logd2h,xlim=c(-2,3.5),ylim=c(-1.5,2),col=Taxa,pch=15,
     xlab=expression(log[10]~(Diameter^2~"*"~Height)),
     ylab=expression(log[10]~Total~mass~(g))))
dat.l <- split(dat,dat$Taxa)
for (i in 1:length(dat.l)){
  fit <- lm(logTM~logd2h,data=dat.l[[i]])
  predline(fit,col=alpha(colors[i],alpha=0.5))

}
with(dat,points(logTM~logd2h,xlim=c(-2,3.5),ylim=c(-1.5,2),col=Taxa,pch=15))
legend("topleft",legend=levels(dat$Taxa),pch=15,col=colors)
title(main="All taxa")

#-----------------------------------
#-- allometry ANCOVA for all taxa
lm.taxa <- lm(logTM~logd2h*Taxa,data=dat)
anova(lm.taxa)
summary(lm.taxa)
#-----------------------------------





#--------------------------------------------------
#- evaluate the relationship between d2h and mass

pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/allom_d2h_height.pdf")

dat2 <- subset(dat,Treat!="Pre")
dat2$Treat <- factor(dat2$Treat)
dat2$TT <- as.factor(paste(dat2$Taxa,dat2$Treat,sep="-"))

dat.l <- split(dat2,dat2$TT)

#par(mfrow=c(2,1),mar=c(2,4,1,1),oma=c(4,1,2,1))
xlims <- c(0,400)
for (i in 1:length(dat.l)){
  plotBy(Height~d2h,data=dat.l[[i]],type="p",xlim=xlims,xlab="d2h",legend=F,col="black",cex=2)
  title(main=dat.l[[i]]$TT[1],outer=T,line=-3)
  
}

dev.off()

#--------------------------------------------------





#-----------------------------------
#- plot leaf vs. totalmass allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
colors <-(c("red","orange","green"))
palette(colors)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logLM <- log10(dat2$Leafmass)
  
  with(subset(dat2,Treat=="Home"),plot(logLM~logTM,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,axes=F,cex=1.5,
                xlab=expression(log[10]~(Total~mass~(g))),
                ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logLM~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logLM~logTM,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logLM~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))

  with(subset(dat2,Treat=="Pre"),points(logLM~logTM,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g)))) 
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Leaf~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0.4,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/LeafMass_TotalMass.pdf")
#-----------------------------------




#-----------------------------------
#- plot root allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logRM <- log10(dat2$Rootmass)
  
  with(subset(dat2,Treat=="Home"),plot(logRM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,axes=F,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logRM~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logRM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logRM~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  with(subset(dat2,Treat=="Pre"),points(logRM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Root~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/RootMass_TotalMass.pdf")
#-----------------------------------





#-----------------------------------
#- plot LAR vs. totalmass allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  
  with(subset(dat2,Treat=="Home"),plot(LAR~logTM,xlim=c(-1.3,2.2),ylim=c(35,200),col=Treat,pch=15,axes=F,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(LAR~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(LAR~logTM,xlim=c(-1.3,2.2),ylim=c(-2.2,1),col=Treat,pch=15,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(LAR~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Pre"),points(LAR~logTM,xlim=c(-1.3,2.2),ylim=c(-2.2,1),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(Leaf~area~ratio~(cm^2*g^-1)),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=150,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/LAR_TotalMass.pdf")
#-----------------------------------






#-----------------------------------
#- plot SLA vs. totalmass allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  
  with(subset(dat2,Treat=="Home"),plot(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(80,500),col=Treat,pch=15,axes=F,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(SLA~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(-2.2,1),col=Treat,pch=15,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(SLA~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Pre"),points(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(-2.2,1),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(Specific~leaf~area~(cm^2*g^-1)),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=380,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/SLA_TotalMass.pdf")
#-----------------------------------

#-----------------------------------
#- plot LA vs. totalleafmass allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  
  with(subset(dat2,Treat=="Home"),plot(Leafarea~Leafmass,xlim=c(0,50),ylim=c(0,9000),col=Treat,pch=15,axes=F,cex=1.5,
                                       xlab=expression((Leaf~mass~(g))),
                                       ylab=expression(Leaf~area~(cm2))))  
  fit <- lm(Leafarea~Leafmass,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(Leafarea~Leafmass,xlim=c(0,50),ylim=c(0,9000),col=Treat,pch=15,cex=1.5,
                                           xlab=expression((Leaf~mass~(g))),
                                           ylab=expression(Leaf~area~(cm2))))  
  fit <- lm(Leafarea~Leafmass,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Pre"),points(Leafarea~Leafmass,xlim=c(0,50),ylim=c(0,9000),col=Treat,pch=15,cex=1.5,
                                        xlab=expression((Leaf~mass~(g))),
                                        ylab=expression(Leaf~area~(cm2))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(Leaf~mass~(g)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(Leaf~area~(cm^2)),side=2,line=2,outer=T,cex=1.5)
legend(x=60,y=6000,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])

  dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/LATotalleafmass.pdf")





