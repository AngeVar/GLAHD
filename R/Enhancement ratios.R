

#- load libraries from script
source("R/loadLibraries.R")
source("R/fitGAM/derivSimulCI.R")
source("R/fitGAM/plotCIdate.R")
source("R/fitGAM/smoothplot.R")



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- read in the data, do a few conversions
dat2 <- return_size_mass_all(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1


#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

dat3$LMF<- with(dat3, leafMass/TotMass)             #Leaf Mass Fraction
dat3$SMF<- with(dat3, stemMass/TotMass)             #Stem Mass Fraction
dat3$RMF<- with(dat3, rootMass/TotMass)             #Root Mass Fraction
dat3$LMA<- with(dat3, leafMass/leafArea)            #Leaf Mass per Area
dat3$SLA<- with(dat3, leafArea/leafMass)            #Specific leaf area
dat3$LAR<- with(dat3, leafArea/TotMass)             #Leaf Area Ratio
dat3$RSR<- with(dat3, rootMass/(stemMass+leafMass)) #Root:Shoot Ratio

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------

#- plot Leaf, stema nd root mass over time

M.trt <- summaryBy(leafMass+stemMass+rootMass~Time+Treatment+Location+Range,data=dat3,FUN=c(mean,standard.error))
M.trt$combotrt <- as.factor(paste(M.trt$Location,M.trt$Range,M.trt$Treatment,sep="_"))

windows(20,10);par(mfrow=c(1,2), mar=c(0.5,3,0.5,3), oma=c(5,1,1,1))

palettS <- c("black","grey","blue","cyan"); palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(M.trt, Location == "S"));ng.trt<- droplevels(subset(M.trt, Location == "N"))

plotBy(leafMass.mean~Time|combotrt,data=sg.trt,col=palettS,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,40),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$leafMass.mean,
                                SE=sg.trt$leafMass.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, text="Leaf Mass (g)",  line=3)
mtext(text="Day", side=1, line=3)

plotBy(leafMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,40),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$leafMass.mean,
                                SE=ng.trt$leafMass.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, text="Leaf Mass (g)",  line=3)
mtext(text="Day", side=1, line=3)

windows(20,10);par(mfrow=c(1,2), mar=c(0.5,3,0.5,3), oma=c(5,1,1,1))
plotBy(stemMass.mean~Time|combotrt,data=sg.trt,col=palettS,,
       legendwhere="topleft",pch=3, ylab= "",xaxt='n', type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$stemMass.mean,
                                SE=sg.trt$stemMass.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text="Stem Mass (g)")
mtext(text="Day", side=1, line=3)

plotBy(stemMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legendwhere="topleft",pch=3, ylab= "",xaxt='n', type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$stemMass.mean,
                                SE=ng.trt$stemMass.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text="Stem Mass (g)")
mtext(text="Day", side=1, line=3)

windows(20,10);par(mfrow=c(1,2), mar=c(0.5,3,0.5,3), oma=c(5,1,1,1))
plotBy(rootMass.mean~Time|combotrt,data=sg.trt,col=palettS,
       legendwhere="topleft",pch=3,ylab= "",xaxt='n', type="l",ylim=c(0,12),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$rootMass.mean,
                                SE=sg.trt$rootMass.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text="Root Mass (g)")
mtext(text="Day", side=1, line=3)


plotBy(rootMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legendwhere="topleft",pch=3,ylab= "",xaxt='n', type="l",ylim=c(0,12),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$rootMass.mean,
                                SE=ng.trt$rootMass.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text="Root Mass (g)")
mtext(text="Day", side=1, line=3)



#Did leaf mass increase?
ber<- summaryBy(leafMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","leafMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, (leafMassWarm)/leafMass)
BER<- bio[,c(1:3,5,8)]

windows(10,10);par(mfrow=c(1,1),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(1,2), xlim=c(0,64),yaxp  = c(1, 2, 10))
#text(2,1.95,labels="A",cex=2,adj=0.5)
abline(h=1)
legend(0,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5)
mtext(text="Leaf mass enhancement",side=2,outer=F,cex=1,line=2.5)

plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main ="", ylab="",xlab="",ylim=c(0.5,1.5), xlim=c(0,64), yaxp  = c(0.5, 1.5, 10))
legend(45,1.5, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
#text(2,1.95,labels="B",cex=2,adj=0.5)
abline(h=1)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5)
mtext(text="Leaf mass enhancement",side=2,outer=F,cex=1,line=2.5)

#over mass
windows(10,10);par(mfrow=c(1,1),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
plotBy(ber~leafMass,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="",ylab="", ylim=c(0.9,2), xlim=c(0,15),xaxp  = c(0, 16, 8),yaxp  = c(0.9, 2, 11))
mtext(text="Leafmass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~leafMass,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
#text(1,1.95,labels="C",cex=2,adj=0.5)
abline(h=1)
mtext(text="Total 'Home' Mass (g)",side=1,outer=F,cex=1,line=2.5)

windows(10,10);par(mfrow=c(1,1),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
plotBy(ber~leafMass,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylab="", ylim=c(0.5,1.5),xlim=c(0,40),xaxp  = c(0, 40, 10),yaxp  = c(0.5, 1.5, 10))
points(ber~leafMass,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Leafmass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
mtext(text="Total 'Home' Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
#text(2,1.95,labels="D",cex=2,adj=0.5)

#Did stem mass increase?
ber<- summaryBy(stemMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","stemMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, stemMassWarm/stemMass)
BER<- bio[,c(1:3,5,8)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.9,2.2), xlim=c(0,64),yaxp  = c(0.9, 2.2, 13))
text(2,2.1,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.9,2.2),yaxt='n', xlim=c(0,64),yaxp  = c(0.9, 2.2, 13))
abline(h=1)
legend(45,2.2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
text(2,2.1,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~stemMass,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.9,2.2), xlim=c(0,13),xaxp  = c(0, 12, 6),yaxp  = c(0.9, 2.2, 13))
mtext(text="Stemmass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~stemMass,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,2.1,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~stemMass,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.9,2.2),yaxt='n', xlim=c(0,29),xaxp  = c(0, 28, 7),yaxp  = c(0.9, 2.2, 13))
points(ber~stemMass,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,2.1,labels="D",cex=2,adj=0.5)

#Did root mass increase?
ber<- summaryBy(rootMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","rootMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, rootMassWarm/rootMass)
BER<- bio[,c(1:3,5,8)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.8,1.9), xlim=c(0,64),yaxp  = c(0.8, 1.9, 11))
text(2,1.85,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.8,1.9),yaxt='n', xlim=c(0,64),yaxp  = c(0.8, 1.9, 11))
abline(h=1)
legend(45,1.9, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
text(2,1.85,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~rootMass,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.8,1.9), xlim=c(0,8),xaxp  = c(0, 7, 7),yaxp  = c(0.8, 1.9, 11))
mtext(text="Rootmass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~rootMass,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(0.5,1.85,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~rootMass,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.8,1.9),yaxt='n', xlim=c(0,12),xaxp  = c(0, 12, 6),yaxp  = c(0.8, 1.9, 11))
points(ber~rootMass,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(0.5,1.85,labels="D",cex=2,adj=0.5)


#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------

#Did leaf mass fraction increase?
#LMF decreased over time with warming in the south

ber<- summaryBy(LMF+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","LMFWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, LMFWarm/LMF)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.93,1.08), xlim=c(0,64),yaxp  = c(0.9, 1.05, 3))
text(2,1.07,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.93,1.08),yaxt='n', xlim=c(0,64),yaxp  = c(0.9, 1.05, 3))
abline(h=1)
legend(45,1.08, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
text(2,1.07,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.93,1.08), xlim=c(0,37),xaxp  = c(0, 35, 7),yaxp  = c(0.9, 1.05, 3))
mtext(text="LeafmassFraction Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.07,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.93,1.08),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8),yaxp  = c(0.9, 1.05, 3))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.07,labels="D",cex=2,adj=0.5)

anova(lm(ber~Time*Location*Range, data=BER))
plot(allEffects(lm(ber~Time*Location*Range, data=BER)))  

#Did stem mass fraction increase?
#No, but there is a Location Range effect.

ber<- summaryBy(SMF+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","SMFWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, SMFWarm/SMF)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.85,1.4), xlim=c(0,64))
text(2,1.35,labels="A",cex=2,adj=0.5)
abline(h=1)
legend(45,1.4, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.85,1.4),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.35,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.85,1.4), xlim=c(0,37),xaxp  = c(0, 35, 7))
mtext(text="StemmassFraction Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.35,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.85,1.4),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.35,labels="D",cex=2,adj=0.5)

summary(lm(ber~Time*Location*Range, data=BER))
plot(allEffects(lm(ber~Time*Location*Range, data=BER)))  

#Did root mass fraction increase?
ber<- summaryBy(RMF+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","RMFWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, RMFWarm/RMF)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.8,1.05), xlim=c(0,64))
text(2,1.025,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.8,1.05),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.025,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,1.05, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.8,1.05), xlim=c(0,37),xaxp  = c(0, 35, 7))
mtext(text="RootmassFraction Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.025,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.8,1.05),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.025,labels="D",cex=2,adj=0.5)

#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------


#Diameter enhancement
ber<- summaryBy(Diameter+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","DiameterWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, DiameterWarm/Diameter)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.9,1.4), xlim=c(0,64))
text(2,1.375,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.9,1.4),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.375,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,1.4, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.9,1.4), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="Diameter Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.375,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.9,1.4),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.375,labels="D",cex=2,adj=0.5)


#Height enhancement
ber<- summaryBy(Height+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","HeightWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, HeightWarm/Height)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.9,1.4), xlim=c(0,64))
text(2,1.375,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.9,1.4),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.375,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,1.4, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.9,1.4), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="Height Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.375,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.9,1.4),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.375,labels="D",cex=2,adj=0.5)



#LMA enhancement
ber<- summaryBy(LMA+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","LMAWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, LMAWarm/LMA)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,64))
text(2,1.275,labels="A",cex=2,adj=0.5)
abline(h=1)
legend(45,1.3, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.275,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)


#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="LMA Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.275,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.275,labels="D",cex=2,adj=0.5)

summary(lm(ber~Time*Location*Range, data=BER))
plot(allEffects(lm(ber~Time*Location*Range, data=BER)))  

#SLA enhancement
ber<- summaryBy(SLA+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","SLAWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, SLAWarm/SLA)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,64))
text(2,1.275,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.275,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,1.3, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="SLA Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.275,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.275,labels="D",cex=2,adj=0.5)

#LAR enhancement
ber<- summaryBy(LAR+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","LARWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, LARWarm/LAR)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,64))
text(2,1.275,labels="A",cex=2,adj=0.5)
abline(h=1)
legend(45,1.3, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.275,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)


#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.7,1.3), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="LAR Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.275,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.7,1.3),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.275,labels="D",cex=2,adj=0.5)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, min), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

summary(lm(ber~Time*Location*Range, data=BER))
plot(allEffects(lm(ber~Time*Location*Range, data=BER)))     


#Root shoot ratio enhancement
ber<- summaryBy(RSR+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","RSRWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, RSRWarm/RSR)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.8,1.05), xlim=c(0,64))
text(2,1.015,labels="A",cex=2,adj=0.5)
abline(h=1)

plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.8,1.05),yaxt='n', xlim=c(0,64))
abline(h=1)
text(2,1.015,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,1.05, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.8,1.05), xlim=c(0,33),xaxp  = c(0, 35, 7))
mtext(text="RSR Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,1.015,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.8,1.05),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.015,labels="D",cex=2,adj=0.5)

#leafArea enhancement
ber<- summaryBy(leafArea+TotMass~Time+Range+Location+Treatment, data=dat3, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","leafAreaWarm","TotMass" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, leafAreaWarm/leafArea)
BER<- bio[,c(1:3,5:6,10)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.7,2.3), xlim=c(0,64),yaxp  = c(0.7, 2.3, 4))
text(2,2.2,labels="A",cex=2,adj=0.5)
abline(h=1)

plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.7,2.3),yaxt='n', xlim=c(0,64),yaxp  = c(0.7, 2.3, 4))
abline(h=1)
text(2,2.2,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
legend(45,2.3, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3, bg="white")

#over mass
plotBy(ber~TotMass.x,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.7,2.3), xlim=c(0,33),xaxp  = c(0, 35, 7),yaxp  = c(0.7, 2.3, 4))
mtext(text="leafArea Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~TotMass.x,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(1,2.2,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~TotMass.x,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main = "", ylim=c(0.7,2.3),yaxt='n', xlim=c(0,80),xaxp  = c(0, 80, 8),yaxp  = c(0.7, 2.3, 4))
points(ber~TotMass.x,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,2.2,labels="D",cex=2,adj=0.5)

