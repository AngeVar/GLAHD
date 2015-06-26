#------------------------------------------------------------------------------------------------------------------------------
# This script assesses the allometry between size and mass for the GLAHD experiment
#------------------------------------------------------------------------------------------------------------------------------



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






#####################
#--- did warming change the allometric relationship between d2h and mass?

#--- No. It is justifiable to combine the home and warmed treatments in this allometry.


#- split into list across all taxa for plotting
allom.2 <- subset(allom,Treat!="Pre")
allom.2$Treat <- factor(allom.2$Treat)
allom.l <- split(allom.2,allom.2$Taxa)

#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(-0.5,2.1)
xlims <- c(-.5,3)
palette(c("black","red"))
for (i in 1:length(allom.l)){
  toplot <- allom.l[[i]]
  plotBy(logTM~logd2h|Treat,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,axes=F,
         ylab="H",xlab="",legend=F)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  
  fit.h <- lm(logTM~logd2h,data=subset(toplot,Treat=="Home"))
  predline(fit.h,col=alpha("black",alpha=0.1))
  
  fit.w <- lm(logTM~logd2h,data=subset(toplot,Treat=="Warmed"))
  predline(fit.w,col=alpha("red",alpha=0.1))
  
}
mtext(expression(log[10]~(diam^2~"*"~height)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Total~plant~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=2,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

ancova.full <- lm(logTM~logd2h*Taxa*Treat,data=allom.2) # most higher order terms not significant
ancova.2 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=allom.2) # drop 3-way interaction
ancova.3 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat,data=allom.2)
ancova.4 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Treat,data=allom.2)
ancova.5 <- lm(logTM~logd2h+Taxa+logd2h:Treat,data=allom.2) # significance of logd2h:Treat interaction goes away after dropping non-significant terms
ancova.6 <- lm(logTM~logd2h+Taxa,data=allom.2)
ancova.y <- lm(logTM~logd2h*Taxa+Taxa,data=allom.2)
anova(ancova.full,ancova.6) # full ancova has slightly lower AIC score, but is not significantly different than the full ancova.
plot(ancova.6) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.6)
#####################




#####################
#-- So, can we add the pre-treatment data then?

#-- yes, It looks like it to me

#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(-1.4,2.1)
xlims <- c(-1.6,3)
palette(c("black","red"))
allom.l <- split(allom,allom$Taxa)
for (i in 1:length(allom.l)){
  toplot <- allom.l[[i]]
  plotBy(logTM~logd2h|Treat,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,axes=F,
         ylab="H",xlab="",legend=F)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  
  fit <- lm(logTM~logd2h,data=toplot)
  predline(fit,col=alpha("black",alpha=0.1))
  
  
}
mtext(expression(log[10]~(diam^2~"*"~height)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Total~plant~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=2,legend=c("Home and Warmed","Pretreatment"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#####################
a1 <- lm(logTM~logd2h+Taxa, data=allom) #with pre gives better fit (higher R-squared, lower AIC)
a2 <- lm(logTM~logd2h+Taxa, data=allom.2) #without pre

#####################
#-- is there a good relationship between total mass and size?
#- plot each taxa on a separate panel. Results in a huge figure
#-- allometry?
windows(12,12);par(mar=c(5,6,1,1))

#- plot taxa
palette(rainbow(length(levels(allom$Taxa))))
colors <- (rainbow(length(levels(allom$Taxa))))

with(allom,plot(logTM~logd2h,xlim=c(-2,3.5),ylim=c(-2,2.5),col=Taxa,pch=15,cex.lab=1.5,
                xlab=expression(log[10]~(Diameter^2~"*"~Height)),
                ylab=expression(log[10]~Total~mass~(g))))
dat.l <- split(allom,allom$Taxa)
for (i in 1:length(dat.l)){
  fit <- lm(logTM~logd2h,data=dat.l[[i]])
  predline(fit,col=alpha(colors[i],alpha=0.5))
  
}
with(allom,points(logTM~logd2h,col=Taxa,pch=15))
legend("topleft",legend=levels(allom$Taxa),pch=15,col=colors)
title(main="All taxa")

allom.2 <- subset(allom,Treatment!="Pre-treatment")
allom.2$Treatment <- factor(allom.2$Treatment)
ancova.full <- lm(logLA~logd2h*Taxa*Treat,data=allom.2) # 3-way term significant, can't be dropped.
# this means that the warming effect on the slope is taxa-specific






#####################
#-- is there a good relationship between total leaf area and size?

#-- allometry?
windows(12,12);par(mar=c(5,6,1,1))

#- plot taxa
palette(rainbow(length(levels(allom$Taxa))))
colors <- (rainbow(length(levels(allom$Taxa))))

with(allom,plot(logLA~logd2h,xlim=c(-2,3.5),ylim=c(0.3,4),col=Taxa,pch=15,cex.lab=1.5,
              xlab=expression(log[10]~(Diameter^2~"*"~Height)),
              ylab=expression(log[10]~Total~leaf~area~(cm^2))))
dat.l <- split(allom,allom$Taxa)
for (i in 1:length(dat.l)){
  fit <- lm(logLA~logd2h,data=dat.l[[i]])
  predline(fit,col=alpha(colors[i],alpha=0.5))
  
}
with(allom,points(logLA~logd2h,col=Taxa,pch=15))
legend("topleft",legend=levels(allom$Taxa),pch=15,col=colors)
title(main="All taxa")

allom.2 <- subset(allom,Treatment!="Pre-treatment")
allom.2$Treatment <- factor(allom.2$Treatment)
ancova.full <- lm(logLA~logd2h*Taxa*Treat,data=allom.2) # 3-way term significant, can't be dropped. This means that the warming effect on the slope is taxa-specific
ancova.full1 <- update(ancova.full,.~. - Taxa:Treat)  # 2-way interaction Taxa:Treat looks non-significant, but it can't be dropped without and logLik and AIC penalty
anova(ancova.full,ancova.full1)
AIC(ancova.full,ancova.full1)
logLik(ancova.full);logLik(ancova.full1)





#- so let's visualize that effect, because I don't understand it
#- plot each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(1.5,4)
xlims <- c(-.5,3)
palette(c("black","red"))
allom.l <- split(allom.2,allom.2$Taxa)
for (i in 1:length(allom.l)){
  toplot <- allom.l[[i]]
  plotBy(logLA~logd2h|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,axes=F,col=c("black","red"),
         ylab="H",xlab="",legend=F)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  
  fit.h <- lm(logLA~logd2h,data=subset(toplot,Treatment=="Home"))
  predline(fit.h,col=alpha("black",alpha=0.1))
  
  fit.w <- lm(logLA~logd2h,data=subset(toplot,Treatment=="Warmed"))
  predline(fit.w,col=alpha("red",alpha=0.1))
  
}
mtext(expression(log[10]~(diam^2~"*"~height)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Total~leaf~area~(cm^2))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=2,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))
#####################
#####################

#Assign three pretreatment plants to Warmed and Home
pretre<- subset(allom, Treat== "Pre")
pretre$sign<- c(rep(sample(c(1,2,1,2,1,2)),16),sample(c(1,2,1,2,1)))
pretre$Treatment <- ifelse(pretre$sign == "1", "Warmed","Home")

allom.3<- droplevels(rbind(allom.2,pretre[,c(1:20)]))
###--------------------------------------------------
windows(12,12);par(mar=c(5,6,1,1))

#- plot taxa
palette(rainbow(length(levels(allom.3$Taxa))))
colors <- (rainbow(length(levels(allom.3$Taxa))))

with(allom.3,plot(logLA~logd2h,xlim=c(-2,3.5),ylim=c(0.3,4),col=Taxa,pch=15,cex.lab=1.5,
                xlab=expression(log[10]~(Diameter^2~"*"~Height)),
                ylab=expression(log[10]~Total~leaf~area~(cm^2))))
dat.l <- split(allom.3,allom.3$Taxa)
for (i in 1:length(dat.l)){
  fit <- lm(logLA~logd2h,data=dat.l[[i]])
  predline(fit,col=alpha(colors[i],alpha=0.5))
  
}
with(allom.3,points(logLA~logd2h,col=Taxa,pch=15))
legend("topleft",legend=levels(allom.3$Taxa),pch=15,col=colors)
title(main="All taxa")

ancova.full <- lm(logLA~logd2h*Taxa*Treatment,data=allom.3) # 3-way term almost significant

allom.l <- split(allom.3,allom.3$Taxa)

windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(-1.4,2.1)
xlims <- c(-1.4,3)
palette(c("black","red"))
for (i in 1:length(allom.l)){
  toplot <- allom.l[[i]]
  plotBy(logTM~logd2h|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,axes=F,
         ylab="H",xlab="",legend=F)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  
  fit.h <- lm(logTM~logd2h,data=subset(toplot,Treatment=="Home"))
  predline(fit.h,col=alpha("black",alpha=0.1))
  
  fit.w <- lm(logTM~logd2h,data=subset(toplot,Treatment=="Warmed"))
  predline(fit.w,col=alpha("red",alpha=0.1))
  
}
mtext(expression(log[10]~(diam^2~"*"~height)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Total~plant~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=2,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

#####################