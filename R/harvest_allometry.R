#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plots the biomass data for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")
source("R/gamplotfunctions.R")
library(scales)

#- read in the data, do a few conversions
dat <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Totmass <- base::rowSums(dat[,11:13]) #total mass is the sum of leaf, stem, and root mass
dat$LAR <- with(dat,Leafarea/Totmass)
<<<<<<< HEAD
dat$SLA <- with(dat,(Leafarea/10000)/(Leafmass))#m2/g
=======
dat$SLA <- with(dat,(Leafarea/10000)/(Leafmass/1000))#m2/g
>>>>>>> 79536e13318ec1b971c3ee32f3b73ab07bc43866
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

pdf(file="Output/allom_d2h_height.pdf")

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
dev.copy2pdf(file="Output/LeafMass_TotalMass.pdf")
#-----------------------------------
dat2 <- subset(dat,Treat!="Pre")
dat2$Treat <- factor(dat2$Treat)
dat2$TT <- as.factor(paste(dat2$Taxa,dat2$Treat,sep="-"))
dat2$logLM <- log10(dat2$Leafmass)

ancova.full <- lm(logLM~logTM*Taxa*Treat,data=dat2) # most higher order terms not significant
ancova.2 <- lm(logLM~logTM+Taxa+Treat+logTM:Taxa+logTM:Treat+Taxa:Treat,data=dat2) # drop 3-way interaction
ancova.3 <- lm(logLM~logTM+Taxa+Treat+logTM:Taxa+logTM:Treat,data=dat2)
ancova.4 <- lm(logLM~logTM+Taxa+Treat+logTM:Taxa,data=dat2)
anova(ancova.full,ancova.4)
plot(ancova.4) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.6)

windows(20,20)
colors <-(c("red","orange","green"))
palette(colors)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logLM <- log10(dat2$Leafmass)
  
  with(subset(dat2,Treat=="Home"),plot(logLM~logd2h,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,axes=F,cex=1.5,
                                       xlab=expression(log[10]~(Total~mass~(g))),
                                       ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logLM~logd2h,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logLM~logd2h,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logLM~logd2h,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  with(subset(dat2,Treat=="Pre"),points(logLM~logd2h,xlim=c(-1.3,2.5),ylim=c(-2.2,2),col=Treat,pch=15,cex=1.5,
                                        xlab=expression(log[10]~(Total~mass~(g))),
                                        ylab=expression(log[10]~Leaf~mass~(g)))) 
  fit <- lm(logLM~logd2h,data=subset(dat2,Treat=="Pre"))
  predline(fit,col=alpha("orange",alpha=0.1))
  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(logd2h~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Leaf~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0.4,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])

#-----------------------------------
dat2 <- subset(dat,Treat!="Pre")
dat2$Treat <- factor(dat2$Treat)
dat2$TT <- as.factor(paste(dat2$Taxa,dat2$Treat,sep="-"))
dat2$logLM <- log10(dat2$Leafmass)
ancova.full <- lm(logLM~logd2h*Taxa*Treat,data=dat2) # most higher order terms not significant
ancova.2 <- lm(logLM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=dat2) # drop 3-way interaction
ancova.3 <- lm(logLM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat,data=dat2)
ancova.4 <- lm(logLM~logd2h+Taxa+Treat+logd2h:Taxa,data=dat2)
ancova.5 <- lm(logLM~logd2h+Taxa+Treat,data=dat2)
ancova.6 <- lm(logLM~logd2h+Taxa,data=dat2) #Taxa specific intercepts
anova(ancova.full,ancova.6)
plot(ancova.6) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.6)


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
dev.copy2pdf(file="Output/RootMass_TotalMass.pdf")
#-----------------------------------

dat2 <- subset(dat,Treat!="Pre")
dat2$Treat <- factor(dat2$Treat)
dat2$TT <- as.factor(paste(dat2$Taxa,dat2$Treat,sep="-"))
dat2$logRM <- log10(dat2$Rootmass)
ancova.full <- lm(logRM~logTM*Taxa*Treat,data=dat2) # three-way interaction significant
anova(ancova.full, ancova.2)

ancova.2 <- lm(logRM~logTM+Taxa+Treat+logTM:Treat,data=dat2)
plot(allEffects(ancova.2)) 

combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logRM <- log10(dat2$Rootmass)
  
  with(subset(dat2,Treat=="Home"),plot(logRM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,axes=F,cex=1.5,
                                       xlab=expression(log[10]~(Total~mass~(g))),
                                       ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logRM~logd2h,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logRM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(logRM~logd2h,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  with(subset(dat2,Treat=="Pre"),points(logRM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,1.5),col=Treat,pch=15,cex=1.5,
                                        xlab=expression(log[10]~(Total~mass~(g))),
                                        ylab=expression(log[10]~Leaf~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(logd2h~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Root~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])

dat2 <- subset(dat,Treat!="Pre")
dat2$Treat <- factor(dat2$Treat)
dat2$TT <- as.factor(paste(dat2$Taxa,dat2$Treat,sep="-"))
dat2$logRM <- log10(dat2$Rootmass)
ancova.full <- lm(logRM~logd2h*Taxa*Treat,data=dat2)
ancova.2 <- lm(logRM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=dat2) # drop 3-way interaction
ancova.3 <- lm(logRM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat,data=dat2)
ancova.4 <- lm(logRM~logd2h+Taxa+Treat+logd2h:Taxa,data=dat2)
ancova.5 <- lm(logRM~logd2h*Taxa+Treat,data=dat2)
anova(ancova.full,ancova.5)
plot(ancova.5) 
anova(ancova.5)

#-----------------------------------


#-----------------------------------
#-----------------------------------
#- plot stem allometries for all taxa
combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logSM <- log10(dat2$Stemmass)
  
  with(subset(dat2,Treat=="Home"),plot(logSM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,axes=F,cex=1.5,
                                       xlab=expression(log[10]~(Total~mass~(g))),
                                       ylab=expression(log[10]~Stem~mass~(g))))  
  fit <- lm(logSM~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logSM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Stem~mass~(g))))  
  fit <- lm(logSM~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  with(subset(dat2,Treat=="Pre"),points(logSM~logTM,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,cex=1.5,
                                        xlab=expression(log[10]~(Total~mass~(g))),
                                        ylab=expression(log[10]~Stem~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Stem~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
dev.copy2pdf(file="Output/StemMass_TotalMass.pdf")
#-----------------------------------

ancova.full <- lm(logSM~logTM*Taxa*Treat,data=dat2) # most higher order terms not significant
ancova.2 <- lm(logSM~logTM+Taxa+Treat+logTM:Taxa+logTM:Treat+Taxa:Treat,data=dat2) # drop 3-way interaction
ancova.3 <- lm(logSM~logTM+Taxa+Treat+logTM:Taxa+logTM:Treat,data=dat2)
ancova.4 <- lm(logSM~logTM+Taxa+Treat+logTM:Taxa,data=dat2)
plot(ancova.full) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.full,ancova.4)

plot(allEffects(ancova.4)) 

combos <- levels(dat$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  dat2$logSM <- log10(dat2$Stemmass)
  
  with(subset(dat2,Treat=="Home"),plot(logSM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,axes=F,cex=1.5,
                                       xlab=expression(log[10]~(Total~mass~(g))),
                                       ylab=expression(log[10]~Stem~mass~(g))))  
  fit <- lm(logSM~logd2h,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(logSM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Stem~mass~(g))))  
  fit <- lm(logSM~logd2h,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  with(subset(dat2,Treat=="Pre"),points(logSM~logd2h,xlim=c(-1.3,2.2),ylim=c(-2,2),col=Treat,pch=15,cex=1.5,
                                        xlab=expression(log[10]~(Total~mass~(g))),
                                        ylab=expression(log[10]~Stem~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(logd2h~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Stem~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=0,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])

#-----------------------------------
ancova.full <- lm(logSM~logd2h*Taxa*Treat,data=dat2)
ancova.2 <- lm(logSM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=dat2) # drop 3-way interaction
ancova.3 <- lm(logSM~logd2h+Taxa+Treat+logd2h:Treat+Taxa:Treat,data=dat2)
anova(ancova.full,ancova.3)
plot(ancova.3) 


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



#plot sla
dat$Range<- as.factor(ifelse(dat$Species == "TER"|dat$Species == "CAM", "wide","narrow"))
dat2<-droplevels(subset(dat,Treat!="Pre"))
dat2$combotrt <- as.factor(paste(dat2$Location,dat2$Range,dat2$Treatment,sep="_"))

windows(20,15);par(mfrow=c(1,2),mar=c(2,0,1,0),oma=c(5,9,3,5),cex.axis=1.2)

<<<<<<< HEAD
ylims=c(0.01,0.04)
boxplot(SLA~Treatment*Range,data=subset(dat2,Location=="N"),ylim=ylims,
        axes=F,las=2,col=c("blue","red"))
=======
ylims=c(5,40)
boxplot(SLA~Treatment*Range,data=subset(dat2,Location=="N"),ylim=ylims,
        axes=F,las=2,col=c("red","blue"))
>>>>>>> 79536e13318ec1b971c3ee32f3b73ab07bc43866
legend("topleft","Tropical",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(SLA),side=2,outer=T,cex=2,line=5)
mtext(text=expression("("*m^2~g^-1*")"),side=2,outer=T,cex=1,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(dat2$Range),las=1,cex.axis=1.5)
boxplot(SLA~Treatment*Range,data=subset(dat2,Location=="S"),ylim=ylims,
<<<<<<< HEAD
        axes=F,las=2,col=c("blue","red"))
=======
        axes=F,las=2,col=c("red","blue"))
>>>>>>> 79536e13318ec1b971c3ee32f3b73ab07bc43866
legend("topleft","Temperate",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(dat2$Range),las=1,cex.axis=1.5)
mtext(text=expression(Range~size),side=1,outer=T,cex=2,line=3)


#-----------------------------------
#- plot SLA vs. totalmass allometries for all taxa
combos <- levels(dat2$Taxa)
#names(combos) <- c("Treat","Taxa")
windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(levels(dat2$Taxa))){
  dat2 <- subset(dat,Taxa==as.character(combos[i]))
  
  with(subset(dat2,Treat=="Home"),plot(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(4,40),col="green",pch=15,axes=F,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(SLA~logTM,data=subset(dat2,Treat=="Home"))
  predline(fit,col=alpha("green",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Warmed"),points(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(4,40),col="red",pch=15,cex=1.5,
                                         xlab=expression(log[10]~(Total~mass~(g))),
                                         ylab=expression(log[10]~Leaf~mass~(g))))  
  fit <- lm(SLA~logTM,data=subset(dat2,Treat=="Warmed"))
  predline(fit,col=alpha("red",alpha=0.1))
  
  
  with(subset(dat2,Treat=="Pre"),points(SLA~logTM,xlim=c(-1.3,2.2),ylim=c(4,40),col="black",pch=15,cex=1.5,
                                           xlab=expression(log[10]~(Total~mass~(g))),
                                           ylab=expression(log[10]~Leaf~mass~(g))))  
  magaxis(side=1:4,labels=c(1,1,0,0),box=T)
  legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
}
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=2,outer=T,cex=1.5)
mtext(expression(Specific~leaf~area~(m^2*g^-1)),side=2,line=2,outer=T,cex=1.5)
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


  #- plot LA vs. totalmass allometries for all taxa
  combos <- levels(dat$Taxa)
  #names(combos) <- c("Treat","Taxa")
  windows(20,20)
  par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
  for (i in 1:length(levels(dat$Taxa))){
    dat2 <- subset(dat,Taxa==as.character(combos[i]))
    
    with(subset(dat2,Treat=="Home"),plot(log(Leafarea)~logTM,xlim=c(0,2.5),ylim=c(3.5,9),col=Treat,pch=15,axes=F,cex=1.5,
                                         xlab=expression((Leaf~mass~(g))),
                                         ylab=expression(Leaf~area~(cm2))))  
    fit <- lm(log(Leafarea)~logTM,data=subset(dat2,Treat=="Home"))
    predline(fit,col=alpha("red",alpha=0.1))
    
    
    with(subset(dat2,Treat=="Warmed"),points(log(Leafarea)~logTM,xlim=c(0,2.5),ylim=c(3.5,9),col=Treat,pch=15,cex=1.5,
                                             xlab=expression((Leaf~mass~(g))),
                                             ylab=expression(Leaf~area~(cm2))))  
    fit <- lm(log(Leafarea)~logTM,data=subset(dat2,Treat=="Warmed"))
    predline(fit,col=alpha("green",alpha=0.1))
    
    
    with(subset(dat2,Treat=="Pre"),points(log(Leafarea)~logTM,xlim=c(0,2.5),ylim=c(3.5,9),col=Treat,pch=15,cex=1.5,
                                          xlab=expression((Leaf~mass~(g))),
                                          ylab=expression(Leaf~area~(cm2))))  
    magaxis(side=1:4,labels=c(1,1,0,0),box=T)
    legend("bottom",legend=dat2$Taxa[1],bty="n",cex=1.5)
  }
  mtext(expression(logTotal~mass~(g)),side=1,line=2,outer=T,cex=1.5)
  mtext(expression(logLeaf~area~(cm^2)),side=2,line=2,outer=T,cex=1.5)
  legend(x=60,y=6000,legend=c("Warmed","Pre","Home"),pch=15,cex=1.5,xpd=NA,col=colors[1:3])
  
ancova.full <- lm(log(Leafarea)~logTM*Taxa*Treat,data=dat2) 
plot(ancova.full) 
anova(ancova.full)  
plot(allEffects(ancova.full)) 


#do analysis using the mixed effects model.
dat$Range<- as.factor(ifelse(dat$Species == "TER"|dat$Species == "CAM", "wide","narrow"))
dat2<-droplevels(subset(dat,Treat!="Pre"))
dat2$combotrt <- as.factor(paste(dat2$Location,dat2$Range,dat2$Treatment,sep="_"))

dat2$Location <- factor(dat2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat2$Sp_RS_EN <- as.factor(with(dat2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat2$Prov_Sp_EN <- as.factor(with(dat2,paste(Taxa,Species)))
dat2$Sp_Loc_EN <- as.factor(with(dat2,paste(Species,Location)))


#---------------------------------------------------------------------------------------------------------------
#- allocation tests for leaf mass and area

#Leaf mass - 3way interaction with treatment and location P=0.1465
#            3way interaction with Location and range P=0.1641
fm1LM <- lme(Leafmass~Totmass*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2)#, method="ML")
plot(fm1LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LM,Leafmass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LM,Leafmass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LM$residuals[,1])
anova(fm1LM)    
summary(fm1LM) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
               #the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1LM)) 
plot(effect("Totmass:Range",fm1LM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location:Range",fm1LM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment",fm1LM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment:Location",fm1LM), multiline=TRUE) #- compares slopes (overlayed)



# #- try recentering data
# #- recenter the leaf mass data
# dat2$Totmass.c <- scale(dat2$Totmass, center = TRUE, scale = FALSE)
# fm2LM <- lme(Leafmass~Totmass.c*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2)#, method="ML")
# plot(effect("Totmass.c:Range",fm2LM), multiline=TRUE) #- compares slopes (overlayed)
# plot(effect("Totmass.c:Location",fm1LM), multiline=TRUE) #- compares slopes (overlayed)
# plot(effect("Totmass.c:Treatment",fm2LM), multiline=TRUE) #- compares slopes (overlayed)
# plot(effect("Totmass.c:Treatment:Location",fm1LM), multiline=TRUE) #- compares slopes (overlayed)


#-- REPEAT FOR LEAF AREA
fm1LA <- lme(sqrt(Leafarea)~Totmass*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2)#, method="ML")
plot(fm1LA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LA,Leafarea~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LA,Leafarea~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LA$residuals[,1])
anova(fm1LA)    
summary(fm1LA) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1LA)) 
plot(effect("Totmass:Range",fm1LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location",fm1LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment",fm1LA), multiline=TRUE) #- compares slopes (overlayed)




#-- REPEAT FOR ROOT MASS
fm1RM <- lme(Rootmass~Totmass*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2,
             weights=varFunc(~Totmass))#, method="ML")
plot(fm1RM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1RM,Rootmass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1RM,Rootmass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1RM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1RM$residuals[,1])
anova(fm1RM)    
summary(fm1RM) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1RM)) 
#plot(effect("Totmass:Range",fm1RM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location:Range",fm1RM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment",fm1RM), multiline=TRUE) #- compares slopes (overlaye  d)
plot(effect("Totmass:Treatment:Location",fm1RM), multiline=TRUE) #- compares slopes (overlayed)



#-- REPEAT FOR STEM MASS
fm1SM <- lme(Stemmass~Totmass*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2,
             weights=varFunc(~Totmass))#, method="ML")
plot(fm1SM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1SM,Stemmass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1SM,Stemmass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1SM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1SM$residuals[,1])
anova(fm1SM)    
summary(fm1SM) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1SM)) 
#plot(effect("Totmass:Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location",fm1SM),multiline=T)
plot(effect("Totmass:Treatment:Location",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location:Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)

#-- LEAF AREA OVER LEAF MASS
fm2LA <- lme((Leafarea)~Leafmass*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2)#, method="ML")
plot(fm2LA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2LA,Leafarea~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm2LA,Leafarea~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm2LA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2LA$residuals[,1])
anova(fm2LA)    
summary(fm2LA) 

plot(allEffects(fm2LA)) 
plot(effect("Leafmass:Range",fm2LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Leafmass:Location",fm2LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Leafmass:Treatment",fm2LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Leafmass:Treatment:Location",fm2LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Leafmass:Treatment:Range",fm2LA), multiline=TRUE) #- compares slopes (overlayed)
coef(summary(fm1LA))

#---------------------------------------------------------------------------------------------------------------

#Was increase in stem mass reflected in stem volume?
dat3<-subset(dat2, Code !="BCAM-2")#d2h NA

fm1d2h <- lme(log(d2h)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat3)#, method="ML")
plot(fm1d2h,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1d2h,d2h~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1d2h,d2h~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1d2h, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1d2h$residuals[,1])
anova(fm1d2h)    
summary(fm1d2h) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1d2h)) 
#plot(effect("Totmass:Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("log(Totmass):Treatment",fm1d2h), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("log(Totmass):Treatment:Location",fm1d2h),multiline=T)
plot(effect("log(Totmass):Treatment:Range",fm1d2h),multiline=T)
