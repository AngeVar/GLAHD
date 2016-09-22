#Figure 3
#Response ratios of biomass and RGR
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=c(mean,standard.error, length), keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment.w","predMassWarm","predMass.standard.error.w","predMass.length.w")
bio <- merge(berh,berw, by=c("Time","Range","Location"))

#Confidence intervals using Fieller's theorem (Intuitive Biostatistics H.J. Motulsky 2014 ISBN10: 0199946647)
alph <- 0.9
dfs<-bio$predMass.length-1
bio$g<- (qt(alph,df=(dfs))*(bio$predMass.standard.error/bio$predMass.mean))^2
A<-bio$predMassWarm
B<-bio$predMass.mean
bio$Q<-A/B
SEMA<-bio$predMass.standard.error.w
SEMB<-bio$predMass.standard.error

bio$SEQ<-(bio$Q/(1-bio$g))*sqrt((1-bio$g)*((SEMA^2)/(A^2))+((SEMB^2)/(B^2)))
bio$CIQhigh<-(bio$Q/(1-bio$g))+qt(alph,df=bio$predMass.length.w+bio$predMass.length-2)*bio$SEQ
bio$CIQlow<-(bio$Q/(1-bio$g))-qt(alph,df=bio$predMass.length+bio$predMass.length-2)*bio$SEQ


windows(11.69,11.69);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(6,6,6,6))

Sbion<- subset(bio, Location =="S"&Range=="narrow"); Sbion<-Sbion[order(Sbion$Time),]
Sbiow<- subset(bio, Location =="S"&Range=="wide"); Sbiow<-Sbiow[order(Sbiow$Time),]
Nbion<- subset(bio,Location =="N"&Range=="narrow"); Nbion<-Nbion[order(Nbion$Time),]
Nbiow<- subset(bio,Location =="N"&Range=="wide"); Nbiow<-Nbiow[order(Nbiow$Time),]
plotBy(Q~Time,data=Sbion,col="black",lty = 2,lwd=2,axes=F,
       legend=F,type="l", main="", ylim=c(0.75,2.2), xaxt='n',xlim=c(0,64),yaxp  = c(0.8, 2.2, 7))
lines(Q~Time,data=Sbiow,col="black",lty=1,lwd=2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
abline(h=1)
mtext("Biomass response ratio", side=2, line=3)
polygon(x = c(Sbion$Time, rev(Sbion$Time)), y = c(Sbion$CIQhigh,rev(Sbion$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Sbiow$Time, rev(Sbiow$Time)), y = c(Sbiow$CIQhigh,rev(Sbiow$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","a", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)

plotBy(Q~Time,data=Nbion,col="black",lty = 2,lwd=2,axes=F,
       legend=F,type="l", main = "", ylim=c(0.75,2.2),yaxt='n',xaxt='n', xlim=c(0,64),yaxp  = c(0.8, 2.2, 7))
lines(Q~Time,data=Nbiow,col="black",lty=1,lwd=2)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
abline(h=1)
legend("topleft", legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1.2,lwd=2, bty="n")
polygon(x = c(Nbion$Time, rev(Nbion$Time)), y = c(Nbion$CIQhigh,rev(Nbion$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Nbiow$Time, rev(Nbiow$Time)), y = c(Nbiow$CIQhigh,rev(Nbiow$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","b", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)

biod<- summaryBy(dydt~Time+Range+Location+Treatment, data=gamfits2, FUN=c(mean,standard.error, length), keep.names=T) 
biodh <- subset(biod, Treatment == "Home")
biodw <- subset (biod, Treatment == "Warmed")
names(biodw)<- c("Time","Range","Location","Treatmentw","dydtWarm", "dydt.standard.error.w","dydt.length.w" )
biod <- merge(biodh,biodw, by=c("Time","Range","Location"))

#Confidence intervals using Fieller's theorem (Intuitive Biostatistics H.J. Motulsky 2014 ISBN10: 0199946647)

dfsd<-biod$dydt.length-1
biod$g<-(qt(alph,df=(dfsd))*(biod$dydt.standard.error/biod$dydt.mean))^2
Ad<-biod$dydtWarm
Bd<-biod$dydt.mean
biod$Q<-Ad/Bd
SEMAd<-biod$dydt.standard.error.w
SEMBd<-biod$dydt.standard.error

biod$SEQ<-(biod$Q/(1-biod$g))*sqrt((1-biod$g)*((SEMAd^2)/(Ad^2))+((SEMBd^2)/(Bd^2)))
biod$CIQhigh<-(biod$Q/(1-biod$g))+qt(alph,df=biod$dydt.length.w+biod$dydt.length-2)*biod$SEQ
biod$CIQlow<-(biod$Q/(1-biod$g))-qt(alph,df=biod$dydt.length+biod$dydt.length-2)*biod$SEQ


Sbiodn<- subset(biod, Location =="S"&Range=="narrow"); Sbiodn<-Sbiodn[order(Sbiodn$Time),]
Sbiodw<- subset(biod, Location =="S"&Range=="wide"); Sbiodw<-Sbiodw[order(Sbiodw$Time),]
Nbiodn<- subset(biod,Location =="N"&Range=="narrow"); Nbiodn<-Nbiodn[order(Nbiodn$Time),]
Nbiodw<- subset(biod,Location =="N"&Range=="wide"); Nbiodw<-Nbiodw[order(Nbiodw$Time),]
plotBy(Q~Time,data=Sbiodn,col="black",lty = 2,lwd=2,axes=F,
       legend=F,type="l", main="", ylim=c(0.75,2.18), xlim=c(0,64),yaxp  = c(0.8, 2, 6))
lines(Q~Time,data=Sbiodw,col="black",lty=1,lwd=2)
abline(h=1)
mtext("RGR response ratio", side=2, line=3)
polygon(x = c(Sbiodn$Time, rev(Sbiodn$Time)), y = c(Sbiodn$CIQhigh,rev(Sbiodn$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Sbiodw$Time, rev(Sbiodw$Time)), y = c(Sbiodw$CIQhigh,rev(Sbiodw$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","c", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,las=1,cex.axis=1.2)

plotBy(Q~Time,data=Nbiodn,col="black",lty = 2,lwd=2, axes=F,
       legend=F,type="l", main = "", ylim=c(0.75,2.18),yaxt='n', xlim=c(0,64),yaxp  = c(0.8, 2, 6))
lines(Q~Time,data=Nbiodw,col="black",lty=1,lwd=2)
abline(h=1)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
polygon(x = c(Nbiodn$Time, rev(Nbiodn$Time)), y = c(Nbiodn$CIQhigh,rev(Nbiodn$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Nbiodw$Time, rev(Nbiodw$Time)), y = c(Nbiodw$CIQhigh,rev(Nbiodw$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","d", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(1,0,0),frame.plot=T,las=1,cex.axis=1.2)


###################

#provenance specific responses
ber<- summaryBy(predMass~Time+Taxa+Treatment, data=gamfits2, FUN=c(mean,standard.error, length), keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa","Treatment.w","predMassWarm","predMass.standard.error.w","predMass.length.w")
bio <- merge(berh,berw, by=c("Time","Taxa"))

#Confidence intervals using Fieller's theorem (Intuitive Biostatistics H.J. Motulsky 2014 ISBN10: 0199946647)
alph <- 0.9
dfs<-bio$predMass.length-1
bio$g<- (qt(alph,df=(dfs))*(bio$predMass.standard.error/bio$predMass.mean))^2
A<-bio$predMassWarm
B<-bio$predMass.mean
bio$Q<-A/B
SEMA<-bio$predMass.standard.error.w
SEMB<-bio$predMass.standard.error

bio$SEQ<-(bio$Q/(1-bio$g))*sqrt((1-bio$g)*((SEMA^2)/(A^2))+((SEMB^2)/(B^2)))
bio$CIQhigh<-(bio$Q/(1-bio$g))+qt(alph,df=bio$predMass.length.w+bio$predMass.length-2)*bio$SEQ
bio$CIQlow<-(bio$Q/(1-bio$g))-qt(alph,df=bio$predMass.length+bio$predMass.length-2)*bio$SEQ

combostrop<- c("BRA",NA,"CTER","DTER","PEL",NA,"ETER","DCAM","PLAT",NA,
               "ECAM","FCAM")
combostemp <- c("BOT",NA,"ATER","BTER","LONG",NA,"ACAM","BCAM","SMIT",NA,"CCAM")
combos<-c(combostrop,rep(NA,4),combostemp)
windows(8.27,11.69)
par(mfrow=c(8,8),mar=c(0,0,0,0),oma=c(6,6,3,3))
layout(matrix(c(1:28), nrow=7, ncol=4,byrow=T),
       heights=c(1,1,1,0.3,1,1,1),
       widths=c(1,0.3,1,1,1))
for (i in 1:length(combos)){
  dat2 <- subset(bio,Taxa==as.character(combos[i]))
  plotBy(Q~Time,data=dat2,col="black",lwd=2,axes=F,
       legend=F,type="l", main="",ylim=c(0.75,3.5),xaxt='n',
       xlim=c(0,64),yaxp=c(0.8,2.2,7),
       lty=ifelse(dat2$Taxa %in% c("BOT","LONG","SMIT","BRA",
                                   "PEL","PLAT"),2,1))
        abline(h=1,col=ifelse(dat2$Taxa==NA,"white","black"))
        polygon(x = c(dat2$Time, rev(dat2$Time)), y = c(dat2$CIQhigh,
                                                        rev(dat2$CIQlow)),
                col = alpha("black",0.2), border = NA)
        
        }
mtext("Biomass response ratio", side=2, line=3)
polygon(x = c(Sbion$Time, rev(Sbion$Time)), y = c(Sbion$CIQhigh,rev(Sbion$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Sbiow$Time, rev(Sbiow$Time)), y = c(Sbiow$CIQhigh,rev(Sbiow$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","a", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,las=1,cex.axis=1.2)

plotBy(Q~Time,data=Nbion,col="black",lty = 2,lwd=2,axes=F,
       legend=F,type="l", main = "", ylim=c(0.75,2.2),yaxt='n',xaxt='n', xlim=c(0,64),yaxp  = c(0.8, 2.2, 7))
lines(Q~Time,data=Nbiow,col="black",lty=1,lwd=2)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
abline(h=1)
legend("topleft", legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1.2,lwd=2, bty="n")
polygon(x = c(Nbion$Time, rev(Nbion$Time)), y = c(Nbion$CIQhigh,rev(Nbion$CIQlow)),col = alpha("black",0.2), border = NA)
polygon(x = c(Nbiow$Time, rev(Nbiow$Time)), y = c(Nbiow$CIQhigh,rev(Nbiow$CIQlow)),col = alpha("black",0.4), border = NA)
legend("topright","b", bty="n", cex=1.2)
magaxis(side=c(1,2,4),labels=c(0,0,0),frame.plot=T,las=1,cex.axis=1.2)
