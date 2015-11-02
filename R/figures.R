#Figures for GAM output
source("R/GLAHD_gamfits.R") #gives you the gamfits2 output
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
###############################################################################################################


#Figure 2: Mass over Time
g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
g.trt.S<- subset(g.trt, Location == "S")
g.trt.N<- subset(g.trt, Location == "N")

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(1,70),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.N$Time,y=g.trt.N$predMass.mean,
                                SE=g.trt.N$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)

axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)

legend(0,60, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=12.5, outer=F)

plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.S$Time,y=g.trt.S$predMass.mean,
                                SE=g.trt.S$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
mtext(text="Time (Days)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=70,by=10), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=70,by=10), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=14, outer=F)
mtext(text="Total biomass (g)", outer=T, side=2, line=1)

###############################################################################################################
#anova Mass Day 60
dat5<- subset(gamfits2,Time==60)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)  

dat5$combotrt <- as.factor(paste(dat5$Location,dat5$Range,dat5$Treatment,sep="_"))
summaryBy(predMass~combotrt, data=dat5, FUN=c(mean,standard.error))

###############################################################################################################
#anova Rgr Day 15
dat6<- subset(gamfits2,Time==15)

dat6$Location <- factor(dat6$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat6$Sp_RS_EN <- as.factor(with(dat6,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat6$Prov_Sp_EN <- as.factor(with(dat6,paste(Taxa,Species)))
dat6$Sp_Loc_EN <- as.factor(with(dat6,paste(Species,Location)))

fm1.rgr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat6)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,log(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,log(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)  

dat6$combotrt <- as.factor(paste(dat6$Location,dat6$Range,dat6$Treatment,sep="_"))
summaryBy(dydt~combotrt, data=dat6, FUN=c(mean,standard.error))

###############################################################################################################
#MaxRGR
maxdydt<-gamfits2[ gamfits2$dydt %in% tapply(gamfits2$dydt, gamfits2$Code, max), ]#max RGR per Code

maxdydt$Location <- factor(maxdydt$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
maxdydt$Sp_RS_EN <- as.factor(with(maxdydt,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
maxdydt$Prov_Sp_EN <- as.factor(with(maxdydt,paste(Taxa,Species)))
maxdydt$Sp_Loc_EN <- as.factor(with(maxdydt,paste(Species,Location)))

fm.maxdydt<- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt,log(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt,log(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt$residuals[,1])
anova(fm.maxdydt)

#Mass @maxRGR
fm.maxdydt2<- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt2,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt2,log(predMass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt2,log(predMass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt2, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt2$residuals[,1])
anova(fm.maxdydt2)  

###############################################################################################################
#- make a table and put it in Word
# extract.lme <- function(model){
#   mod.o <- anova(model)
#   ndf <- mod.o[2:nrow(mod.o),1]
#   ddf <- mod.o[2:nrow(mod.o),2]
#   df <- (paste(ndf,ddf,sep=", "))
#   Fvalue <-round(mod.o[2:nrow(mod.o),3],1)
#   Pvalue <-round(mod.o[2:nrow(mod.o),4],3)
#   Pvalue[which(Pvalue==0)] <- "<0.001"
#   return(data.frame("df"=df, "F" = Fvalue,"P" = Pvalue))
# }
# 
# table <- do.call(cbind,lapply(list(fm1.mass,fm1.rgr,fm.maxdydt,fm.maxdydt2),FUN=extract.lme))
# row.names(table) <- row.names(anova(fm1.mass))[2:8]
# 
# library(R2wd)
# wdGet()
# wdTable(gxtable,autoformat=2)
###############################################################################################################
#Figure 3:RGR over Time

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(dydt.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0.01,0.13),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.N$Time,y=g.trt.N$dydt.mean,
                                SE=g.trt.N$dydt.standard.error,direction="updown",
                                col="black",0))
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
abline(v=15, lty=2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.12,by=0.02), labels = T, tck = 0.01)
axis(side = 4, at = seq(from=0,to=0.12,by=0.02), labels = F, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)

legend(32,0.115, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=44.5, outer=F)

plotBy(dydt.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0.01,0.13),lty=2,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=g.trt.S$Time,y=g.trt.S$dydt.mean,
                                SE=g.trt.S$dydt.standard.error,direction="updown",
                                col="black",0))
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=2,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
lines(dydt.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="l",ylim=c(0,70),lty=1,lwd=2)
abline(v=15, lty=2)
mtext(text="Time (Days)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.12,by=0.02), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=0.12,by=0.02), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=45, outer=F)
mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=0.5)


################################################################################################################
#Figure 4: Bargraph of RGR at 15 days

g.trt.15<-subset(g.trt, Time=="15")
g.trt.15$combo <- as.factor(paste(g.trt.15$Location,g.trt.15$Range,sep="_"))

target <- c("S_narrow_Home", "S_narrow_Warmed","S_wide_Home","S_wide_Warmed","N_narrow_Home", "N_narrow_Warmed", "N_wide_Home", "N_wide_Warmed")
g.trt.15<-g.trt.15[match(target, g.trt.15$combotrt),]

windows(5,4)
bar<-barplot(g.trt.15$dydt.mean, space=c(0,0,1,0,1,0,1),ylim=c(0.05,0.15), xpd = F, col=c(rep(c(alpha("black",0.6),alpha("red",0.8)),4)), 
             axis.lty=0)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
error.bar(bar,g.trt.15$dydt.mean,g.trt.15$dydt.standard.error )
mtext(text=expression(RGR~(g~g^-1~day^-1)), outer=T, side=2, line=-1.5)
mtext(c("Narrow","Wide","Narrow","Wide"), side=1, line = 1, at=c(1,4,7,10))
abline(v=5.5)
box()
mtext("Temperate Eucalyptus",3,line=-1.5,at=2.5, outer=F)
mtext("Tropical Eucalyptus",3,line=-1.5,at=8.5, outer=F)
legend(-0.1,0.135, legend=c("H","W"),pch=22,pt.cex=2, pt.bg=c(alpha("black",0.6),"red"),bty="n")
################################################################################################################
#Response ratios

#Relative enhancement of biomass
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:3,5,8)]

windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(4,6,1,1),cex.axis=1)
SBERn<- subset(BER, Location =="S"&Range=="narrow"); SBERn<-SBERn[order(SBERn$Time),]
SBERw<- subset(BER, Location =="S"&Range=="wide"); SBERw<-SBERw[order(SBERw$Time),]
NBERn<- subset(BER,Location =="N"&Range=="narrow"); NBERn<-NBERn[order(NBERn$Time),]
NBERw<- subset(BER,Location =="N"&Range=="wide"); NBERw<-NBERw[order(NBERw$Time),]
plotBy(ber~Time,data=SBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main="", ylim=c(0.5,2), xaxt='n',xlim=c(0,64),yaxp  = c(0.5, 2, 3))
lines(ber~Time,data=SBERw,col="black",lty=1,lwd=2)
#text(2,1.95,labels="A",cex=2,adj=0.5)
abline(h=1)
mtext(expression(Delta~"Biomass"), side=2, line=3)
plotBy(ber~Time,data=NBERn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main = "", ylim=c(0.5,2),yaxt='n',xaxt='n', xlim=c(0,64),yaxp  = c(0.5, 2, 3))
lines(ber~Time,data=NBERw,col="black",lty=1,lwd=2)
abline(h=1)
legend(35,2, legend=c("Wide","Narrow"), lty=c(1,2), col="black",cex=1,lwd=2, bty="n")
#text(2,1.95,labels="B",cex=2,adj=0.5)

berd<- summaryBy(dydt~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berdh <- subset(berd, Treatment == "Home")
berdw <- subset (berd, Treatment == "Warmed")
names(berdw)<- c("Time","Range","Location","Treatment","dydtWarm" )
biod <- merge(berdh,berdw, by=c("Time","Range","Location"))
biod$berd<- with(biod, dydtWarm/dydt)
BERd<- biod[,c(1:3,5,8)]

SBERdn<- subset(BERd, Location =="S"&Range=="narrow"); SBERdn<-SBERdn[order(SBERdn$Time),]
SBERdw<- subset(BERd, Location =="S"&Range=="wide"); SBERdw<-SBERdw[order(SBERdw$Time),]
NBERdn<- subset(BERd,Location =="N"&Range=="narrow"); NBERdn<-NBERdn[order(NBERdn$Time),]
NBERdw<- subset(BERd,Location =="N"&Range=="wide"); NBERdw<-NBERdw[order(NBERdw$Time),]
plotBy(berd~Time,data=SBERdn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main="", ylim=c(0.5,2), xlim=c(0,64),yaxp  = c(0.5, 1.5, 2))
lines(berd~Time,data=SBERdw,col="black",lty=1,lwd=2)
#text(2,1.95,labels="C",cex=2,adj=0.5)
abline(h=1)
mtext(expression(Delta~"RGR"), side=2, line=3)
plotBy(berd~Time,data=NBERdn,col="black",lty = 2,lwd=2,
       legend=F,type="l", main = "", ylim=c(0.5,2),yaxt='n', xlim=c(0,64),yaxp  = c(0.5, 1.5, 2))
lines(berd~Time,data=NBERdw,col="black",lty=1,lwd=2)
abline(h=1)
#text(2,1.95,labels="D",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)
################################################################################################################
#Allocation
#LEAF AREA -- Cant get this to look nice...
dat <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Totmass <- base::rowSums(dat[,11:13]) #total mass is the sum of leaf, stem, and root mass
dat$Treat <- as.factor(ifelse(dat$Pot < 20, "Home",ifelse(dat$Pot>=40,"Pre","Warmed")))
dat[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
dat$Range <- as.factor(ifelse(dat$Species == "CAM"|dat$Species=="TER","wide","narrow"))
dat2<-droplevels(subset(dat,Treat!="Pre"))
dat2$combotrt <- as.factor(paste(dat2$Location,dat2$Range,dat2$Treatment,sep="_"))
dat2$Location <- factor(dat2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat2$Sp_RS_EN <- as.factor(with(dat2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat2$Prov_Sp_EN <- as.factor(with(dat2,paste(Taxa,Species)))
dat2$Sp_Loc_EN <- as.factor(with(dat2,paste(Species,Location)))

fm1LA <- lme(log(Leafarea)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat2)#, method="ML")
plot(fm1LA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LA,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LA,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LA$residuals[,1])
anova(fm1LA)    
summary(fm1LA) 
plot(effect("log(Totmass):Treatment:Location",fm1LA), multiline=TRUE) #- compares slopes (overlayed) 0.978667+0.139295

sdatw<-dat[dat$Treat=="Warmed"&dat$Location=="S",];sdath<-dat[dat$Treat=="Home"&dat$Location=="S",]
ndatw<-dat[dat$Treat=="Warmed"&dat$Location=="N",];ndath<-dat[dat$Treat=="Home"&dat$Location=="N",]
# plot((sdath$Totmass),(sdath$Leafarea),pch=16)
# points((Leafarea)~(Totmass), pch=16, data=sdatw, col="red")
# 
# abline(lm((Leafarea)~(Totmass),data=sdath), col="black",untf =T)
# abline(lm((Leafarea)~(Totmass),data=sdatw), col="red",untf =T)

windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))
plot((sdath$Totmass),(sdath$Leafarea),pch=16)
points((Leafarea)~(Totmass), pch=16, data=sdatw, col="red")
abline(lm((Leafarea)~(Totmass),data=sdath), col="black",untf =T)
abline(lm((Leafarea)~(Totmass),data=sdatw), col="red",untf =T)

plot((ndath$Totmass),(ndath$Leafarea),pch=16)
points((Leafarea)~(Totmass), pch=16, data=ndatw, col="red")
abline(lm((Leafarea)~(Totmass),data=ndath), col="black",untf =T)
abline(lm((Leafarea)~(Totmass),data=ndatw), col="red",untf =T)
#Mass fractions
#stacked bargraph Root:Stem:Leaf, Home and warmed next to eachother, for the different groups

library(reshape2)

dat2 <- return_size_mass_all(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1
#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)); names <- names(table(dat2$Code));keeps <- names[which(obs>6)]
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)

#make new total mass form stem, root and leaf mass
dat3$predTotMass<- with(dat3, leafMass+stemMass+rootMass)

#- Calculate some new variables
dat3$LMF<- with(dat3, leafMass/predTotMass)             #Leaf Mass Fraction
dat3$SMF<- with(dat3, stemMass/predTotMass)             #Stem Mass Fraction
dat3$RMF<- with(dat3, rootMass/predTotMass)             #Root Mass Fraction
dat3$LMA<- with(dat3, leafMass/leafArea)            #Leaf Mass per Area
dat3$SLA<- with(dat3, leafArea/leafMass)            #Specific leaf area
dat3$LAR<- with(dat3, leafArea/TotMass)             #Leaf Area Ratio
dat3$RSR<- with(dat3, rootMass/(stemMass+leafMass)) #Root:Shoot Ratio

rate <- merge(gamfits2,dat3,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,13:14,17:30)]
rate$ULR <- with(rate,dydt/LAR)

rate11<- subset(rate, Time ==11)
r.trt11 <- summaryBy(LMF+SMF+RMF~Treatment+Location+Range,data=rate11,FUN=mean, keep.names=T)
r.trt11$combotrt <- as.factor(paste(r.trt11$Location,r.trt11$Range,r.trt11$Treatment,sep="_"))
r.trt11$combo <- as.factor(paste(r.trt11$Location,r.trt11$Range,sep="_"))
r.trt11$Treat <- as.factor(ifelse(r.trt11$Treatment == "Home", "H","W"))

df1<- melt(r.trt11, c("combo","Treat"),c("RMF","SMF","LMF"),"MassFraction")
# target <- c(rep(c("S_narrow", "S_wide","N_narrow", "N_wide"),6))
# df1$combo<-df1[match(target, df1$combo),]

df2<- df1[order(df1$combo),]
m1 <- acast(df2, MassFraction~combo*Treat, value.var=c('value'), fill=0)

windows(5.845,4.135);par(mfrow=c(1,1),mar=c(6,4,2,2))
bar<- barplot(m1,beside=F, space=c(1,0), ylim=c(0,1), col=c(alpha("black",0.8),alpha("red",0.8),"white"), 
        xaxt='n',ylab="Biomass Fractions",xlab=" ", main=" ")
box()
abline(v=6.5)
legend("topleft", legend=c("Leaf","Stem","Root"),pch=22,pt.cex=2, pt.bg=c("white",alpha("red",0.8),alpha("black",0.8)),bg="white",ncol=3)
mtext(c("H","W","H","W","H","W","H","W"),side=1, line=0.25, at=c(1.5,2.5,4.5,5.5,7.5,8.5,10.5,11.5))
mtext(c("Wide","Narrow","Wide","Narrow"), side=1, line = 1.25, at=c(2,5,8,11))
mtext("Tropical",1,line=2.25,at=3.5, outer=F)
mtext("Temperate",1,line=2.25,at=9.5, outer=F)

r.tr <- summaryBy(LMF+SMF+RMF~Treatment,data=rate11,FUN=mean, keep.names=T)
r.tr2<- melt(r.tr, "Treatment",c("RMF","SMF","LMF"),"MassFraction")
m2 <- acast(r.tr2, MassFraction~Treatment, value.var=c('value'), fill=0)
barplot(m2, beside=F)

################################################################################################################
#Asat
Asat <- getAsat()
asat.tm <- summaryBy(Photo~Taxa+Treatment+Location+Range,data=Asat,FUN=mean,keep.names=T)
gx2 <- getRdark()
gx2$Species <- factor(gx2$Species)
gx2.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=gx2,FUN=mean,keep.names=T)

#- make plots
colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)

#Rmass
ylims=c(5,20)
boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)

mtext(text=expression(R["mass"]),side=2,outer=F,cex=1,adj=0.5,line=3.5)
mtext(text=expression("("*n*mol~g^-1~s^-1*")"),side=2,outer=F,cex=1,adj=0.5,line=2)
legend(-0.2,21.2,"Temperate", bty="n", cex=1.3)

boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,0),frame.plot=T,las=1)
legend(-0.2,21.2,"Tropical", bty="n", cex=1.3)
legend("topright",legend=c("H","W"),ncol=1, fill=colors,cex=1, bty="n")
#Asat
ylims=c(18,38)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
axis(side=1,at=c(1.5,3.5),labels=levels(gx2.tm$Range),las=1,cex.axis=1.5)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1,adj=0.5,line=3.5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=2)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(gx2.tm$Range),las=1,cex.axis=1.5)

################################################################################################################
#--------------------------------------------------------------------------------------------------------
#- Respiration
source("R/loadLibraries.R")
library(scales)
dat <- read.csv("Data/GasEx/RvsT/GLADH RvsT gasex processeed.csv",stringsAsFactors=T)
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Taxa <- as.factor(sapply(1:nrow(dat),function(i)strsplit(as.character(dat$Code[i]),split="-")[[1]][1])) # get the "Taxa" factor out of Code. Note this could have been a loop
dat$R_area <- dat$Photo*-1

#- read in the leaf mass data, calculate R per unit leaf mass
leafmass <- read.csv("C:/Repos/GLAHD/Data/GasEx/RvsT/GLAHD_leaf_mass_RvT.csv")
rt <- merge(dat,leafmass,by="Code")
rt$R_mass <- rt$R_area/10000*rt$Area/rt$mass*1000 #note mass was in g


rt <- subset(rt,Taxa!="SMIT" & Taxa !="ACAM")
rt$Taxa <- factor(rt$Taxa,levels=c("BOT","BTER","BRA","CTER"))
rt$lnRmass <- log(rt$R_mass)
rt <- subset(rt,Code!="CTER-4" & Code!="CTER-22" & Code!="CTER-27") #CTER4 looks weird... peaks much much lower than others. CTER-22 also looks bad
rt$Code <- factor(rt$Code)
gx2 <- getRdark()
#- break data into bins of leaf temperatures, for averaging and plotting
rt$Tleaf_bin <- cut(rt$Tleaf,breaks=seq(from=10,to=65,by=2))#length=60))
rt$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rt$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 

#- average across taxa
gx2$Species <- factor(gx2$Species)
gx2.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=gx2,FUN=mean,keep.names=T)

#- make plots
colors <- c(alpha("black",0.8),alpha("red",0.8))
windows(22,20);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1.2)

#Rmass
ylims=c(5,15)
boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
title(main="Tropical",line=0.2,cex.main=1.5,xpd=NA)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(R["mass"]~"("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1.5,adj=0.5,line=3)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topleft",c("Home","Warmed"),fill=colors,cex=1.3,bty="n")
legend("topright","a",bty="n",cex=1.3)
boxplot(Rmass~Treatment*Range,data=subset(gx2.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
title(main="Temperate",line=0.2,cex.main=1.5,xpd=NA)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=c("Narrow","Wide"),las=1,cex.axis=1.5)
legend("topright","b",bty="n",cex=1.3)
#mtext("Range size",side=1,cex=1.5,outer=T,line=-17)

#- plot R vs. T
rt.m <- summaryBy(R_mass+R_area~Species+Taxa+Treatment+Tleaf_bin_mid,data=rt,FUN=c(mean,standard.error))
toplot <- subset(rt.m,Tleaf_bin_mid < 40)
toplot.s <- subset(toplot,Taxa=="BOT" | Taxa=="BTER")
toplot.n <- subset(toplot,Taxa=="BRA" | Taxa=="CTER")

ylims <- c(0,40)
xlims <- c(10,45)
cexpoints=1.8
# plot north
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.n,Taxa=="BRA"),xlim=xlims,ylim=ylims,pch=1,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.n,Taxa=="BRA")$Tleaf_bin_mid,
                                y=subset(toplot.n,Taxa=="BRA")$R_mass.mean,
                                SE=subset(toplot.n,Taxa=="BRA")$R_mass.standard.error,direction="updown",col="black"),
       las=1,axes=F,xlab="")
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.n,Taxa=="CTER"),ylim=ylims,pch=16,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.n,Taxa=="CTER")$Tleaf_bin_mid,
                                y=subset(toplot.n,Taxa=="CTER")$R_mass.mean,
                                SE=subset(toplot.n,Taxa=="CTER")$R_mass.standard.error,direction="updown",col="black"),
       add=T,las=1)
magaxis(c(1,2,3,4),labels=c(1,1,0,0),frame.plot=T,las=1)
#title(main="North",line=0.2,cex.main=2,xpd=NA)
legend("topright","c",bty="n",cex=1.3)
legend("topleft",c("Narrow","Narrow Warmed","Wide","Wide Warmed"),pch=c(1,1,16,16),col=c(alpha("black",0.8),alpha("red",0.8)),ncol=1,
       cex=1.2,bty="n")

# plot south
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.s,Taxa=="BOT"),xlim=xlims,ylim=ylims,pch=1,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.s,Taxa=="BOT")$Tleaf_bin_mid,
                                y=subset(toplot.s,Taxa=="BOT")$R_mass.mean,
                                SE=subset(toplot.s,Taxa=="BOT")$R_mass.standard.error,direction="updown",col="black"),
       axes=F,las=1,xlab="")
plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=subset(toplot.s,Taxa=="BTER"),ylim=ylims,pch=16,cex=cexpoints,col=colors,legend=F,
       panel.first=adderrorbars(x=subset(toplot.s,Taxa=="BTER")$Tleaf_bin_mid,
                                y=subset(toplot.s,Taxa=="BTER")$R_mass.mean,
                                SE=subset(toplot.s,Taxa=="BTER")$R_mass.standard.error,direction="updown",col="black"),
       add=T,las=1,axes=F,xlab="")
#title(main="South",line=0.2,cex.main=2,xpd=NA)
legend("topright","d",bty="n",cex=1.3)

magaxis(c(1,2,3,4),labels=c(1,0,0,1),frame.plot=T,las=1)
mtext(text=expression(T["leaf"]~"("*degree~C*")"),side=1,outer=T,cex=1.5,adj=0.5,line=3)
#mtext(text=expression(R["mass"]~"("*n*mol~g^-1~s^-1*")"),side=2,outer=T,cex=2,adj=0.5,line=3)

################################################################################################################

#LA over Total mass
#- load libraries from script
source("R/loadLibraries.R");source("R/gamplotfunctions.R");library(scales)
#- read in the data, do a few conversions
dat <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Totmass <- base::rowSums(dat[,11:13]) #total mass is the sum of leaf, stem, and root mass
dat$Treat <- as.factor(ifelse(dat$Pot < 20, "Home",
                              ifelse(dat$Pot>=40,"Pre","Warmed")))
dat$Taxa <- factor(dat$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
dat[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
dat2 <- subset(dat,Treat!="Pre");dat2$Treat <- factor(dat2$Treat)
dat2$Range<- ifelse(dat2$Species == "CAM"|dat2$Species == "TER", "wide","narrow")


la.trt <- summaryBy(Leafarea~Totmass+Treatment+Location+Range,data=dat2,FUN=c(mean,standard.error))
la.trt$combotrt <- as.factor(paste(la.trt$Location,la.trt$Range,la.trt$Treatment,sep="_"))
la.trt.S<- subset(la.trt, Location == "S")
la.trt.N<- subset(la.trt, Location == "N")

slopes<-read.csv("Data/slopes for fig 3 LA.csv")#read in slopes from effect("Totmass:Treatment:Location",fm1LA)
slopesS<- subset(slopes, Location == "S")
slopesN<- subset(slopes, Location == "N")
sShome<-lm(Home~Totmass, data=slopesS)
sSwarmed<-lm(Warmed~Totmass, data=slopesS)
sNhome<-lm(Home~Totmass, data=slopesN)
sNwarmed<-lm(Warmed~Totmass, data=slopesN)
library(plotrix)
windows(8.27,11.69);par(mfrow=c(2,1),mar=c(0,2,0,1),oma=c(4,2,2,1))

plotBy(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0,10000),xlim=c(0,120),pch=1,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=la.trt.N$Totmass,y=la.trt.N$Leafarea.mean,
                                SE=la.trt.N$Leafarea.standard.error,direction="updown",
                                col="black",0))
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=1,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_wide_Home"),col="black",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="N_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)

axis(side = 1, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=10000,by=1000), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=10000,by=1000), labels = F, tck = 0.01)
ablineclip(sNhome, col="black", x1=0,x2=120)
ablineclip(sNwarmed, col="red", x1=0,x2=103)

legend(0,60, legend=c("Narrow","Narrow Warmed","Wide","Wide Warmed"),
       col=c("black","red","black","red"),lty=c(2,2,1,1), lwd=2,bty="n")
mtext("Tropical Eucalyptus",3,line=-1.5,at=22, outer=F)

plotBy(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0,10000),xlim=c(0,120),pch=1,lwd=2,yaxs="i",xaxs="i",
       panel.first=adderrorbars(x=la.trt.S$Totmass,y=la.trt.S$Leafarea.mean,
                                SE=la.trt.S$Leafarea.standard.error,direction="updown",
                                col="black",0))
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_narrow_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=1,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_wide_Home"),col="black",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
points(Leafarea.mean~Totmass, data=subset(la.trt, combotrt=="S_wide_Warmed"),col="red",
      xaxt='n', ylab="", type="p",ylim=c(0,10000),pch=16,lwd=2)
mtext(text="Total biomass (g)", side=1, line=3)
axis(side = 1, at = seq(from=0,to=110,by=10), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=10000,by=1000), labels = T, tck = 0.01)
axis(side = 3, at = seq(from=0,to=110,by=10), labels = F, tck = 0.01)
axis(side = 4, at = seq(from=0,to=10000,by=1000), labels = F, tck = 0.01)
mtext("Temperate Eucalyptus",3,line=-1.5,at=25, outer=F)
mtext(text="Leaf area", outer=T, side=2, line=1)
ablineclip(sShome, col="black",x1=0,x2=51)
ablineclip(sSwarmed, col="red",x1=0,x2=78)


###############################################################################################################

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

rate <- merge(gamfits2,dat2,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,18:19)]# merge with dat2 to get LA
rate$LAR <- rate$leafArea/rate$TotMass
rate$ULR <- with(rate,dydt/LAR)
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))

windows(60,40);par(mfrow=c(4,2), mar=c(0.5,6,0.5,3), oma=c(5,1,1,1))

palettS <- c("black","grey","blue","cyan"); palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|combotrt,data=sg.trt,col=palettS,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legendwhere="topleft", pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|combotrt,data=sg.trt,col=palettS,,
       legend=F,pch=3, ylab= "",xaxt='n', type="l",ylim=c(0.04,0.12),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(dydt.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3, ylab= "",xaxt='n', type="l",ylim=c(0.04,0.12),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(AGR.mean~Time|combotrt,data=sg.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",,log="y",
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))
  

plotBy(AGR.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",log="y",
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|combotrt,data=srate.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|combotrt,data=nrate.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)
 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#RGR and AGR over Mass
palettS <- c("black","grey","blue","cyan"); palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

windows(60,30);par(mfrow=c(3,2), mar=c(0.5,6,0.5,3), oma=c(5,1,1,1))

palettS <- c("black","grey","blue","cyan")
palettN <- c("black","grey","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"))
ng.trt<- droplevels(subset(g.trt, Location == "N"))

plotBy(dydt.mean~predMass.mean|combotrt,data=sg.trt,col=palettS,,
       legendwhere="topright",pch=3, ylim=c(0.04,0.12),ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=sg.trt$predMass.mean,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(dydt.mean~predMass.mean|combotrt,data=ng.trt,col=palettN,
       legendwhere="topright",pch=3, ylim=c(0.04,0.12),ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=ng.trt$predMass.mean,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))

plotBy(AGR.mean~predMass.mean|combotrt,data=sg.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=sg.trt$predMass.mean,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(AGR.mean~predMass.mean|combotrt,data=ng.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",
       panel.first=adderrorbars(x=ng.trt$predMass.mean,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~predMass.mean|combotrt,data=srate.trt,col=palettS,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=srate.trt$predMass.mean,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black","blue","grey","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Mass (g)", side=1, line=3)

plotBy(LAR.mean~predMass.mean|combotrt,data=nrate.trt,col=palettN,
       legend=F,pch=3,ylab= "",xaxt='n', type="l",ylim=c(60,150),
       panel.first=adderrorbars(x=nrate.trt$predMass.mean,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black","orange","grey","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Mass (g)", side=1, line=3)



#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#Relative enhancement of biomass
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:3,5,8)]

windows(10,10);par(mfrow=c(2,2),mar=c(2,0,2,0),oma=c(2,4,2,2),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="", ylim=c(0.9,2), xlim=c(0,64),yaxp  = c(0.9, 2, 11))
text(2,1.95,labels="A",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "", ylim=c(0.9,2),yaxt='n', xlim=c(0,64),yaxp  = c(0.9, 2, 11))
abline(h=1)
legend(45,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=1.3)
text(2,1.95,labels="B",cex=2,adj=0.5)
mtext(text="Time (Days)",side=1,outer=F,cex=1,line=2.5, at=0)

#over mass
plotBy(ber~predMass,data=subset(SBER, Range=="narrow"),col="black",pch = c(1),
       legend=F,type="p", main="", ylim=c(0.9,2), xlim=c(0,37),xaxp  = c(0, 35, 7),yaxp  = c(0.9, 2, 11))
mtext(text="Biomass Enhancement Ratio",side=2,outer=T,cex=1,adj=0.5,line=2.5)
points(ber~predMass,data=subset(SBER, Range=="wide"),col="black",pch = c(16))
text(2,1.95,labels="C",cex=2,adj=0.5)
abline(h=1)
plotBy(ber~predMass,data=subset(NBER,Range=="narrow"),col="black",pch = c(1),
legend=F,type="p", main = "", ylim=c(0.9,2),yaxt='n', xlim=c(0,75),xaxp  = c(0, 80, 8),yaxp  = c(0.9, 2, 11))
points(ber~predMass,data=subset(NBER, Range=="wide"),col="black",pch = c(16))
mtext(text="Control Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=0.5)
abline(h=1)
text(2,1.95,labels="D",cex=2,adj=0.5)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]


#Relative enhancement of biomass using rate data
ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=rate, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:3,5,8)]

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
SBER<- subset(BER, Location =="S")
NBER<- subset(BER,Location =="N")
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p", main="South", ylim=c(0.5,2), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p", main = "North", ylim=c(0.5,2),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(53,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=0.8)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

#Relative enhancement of biomass per species
ber<- summaryBy(predMass~Time+Species+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Species", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Species","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))
SBER$Species2<-factor(SBER$Species,levels(SBER$Species)[c(2,5,1,3,4)]);NBER$Species2<-factor(NBER$Species,levels(NBER$Species)[c(2,5,1,3,4)])

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Species2,data=SBER,col=c("black","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,2.5), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(54,2.5, legend=c(levels(SBER$Species2)), pch=16, col=c("black","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Species2,data=NBER,col=c("black","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,2.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(54,2.5, legend=c(levels(NBER$Species2)), pch=16, col=c("black","orange","grey","red","magenta"), cex=0.8)

SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Relative enhancement of biomass per taxa
ber<- summaryBy(predMass~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,3.5), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(54,3.7, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,3.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(54,3.7, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Relative enhancement of biomass per taxa
ber<- summaryBy(predMass~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","predMassWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, predMassWarm/predMass)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~predMass|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),log="x",
       legend=F,type="p", main="South", ylim=c(0.5,3.5), xlim=c(0,45))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(35,3.7, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~predMass|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),log="x",
       legend=F,type="p", main = "North", ylim=c(0.5,3.5),yaxt='n', xlim=c(0,90))
mtext(text="Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(70,3.7, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Relative enhancement of RGR per taxa
ber<- summaryBy(dydt~Time+Taxa+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")
names(berw)<- c("Time","Taxa", "Range","Location","Treatment","dydtWarm" )
bio <- merge(berh,berw, by=c("Time","Taxa","Range","Location"))
bio$ber<- with(bio, dydtWarm/dydt)
BER<- bio[,c(1:4,6,9)]
SBER<- droplevels(subset(BER, Location =="S")); NBER<- droplevels(subset(BER,Location =="N"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(ber~Time|Taxa,data=SBER,col=c("black","black","blue","blue","blue","grey","cyan","purple"),
       legend=F,type="p", main="South", ylim=c(0.5,1.8), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
legend(50,1.8, legend=c(levels(SBER$Taxa)), pch=16, col=c("black","black","blue","blue","blue","grey","cyan","purple"), cex=0.8)

plotBy(ber~Time|Taxa,data=NBER,col=c("black","black","black","orange","orange","orange","grey","red","magenta"),
       legend=F,type="p", main = "North", ylim=c(0.5,1.8), xlim=c(0,67),yaxt='n')
mtext(text="Mass (g)",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)
legend(50,1.8, legend=c(levels(NBER$Taxa)), pch=16, col=c("black","black","black","orange","orange","orange","grey","red","magenta"), cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(40,40);
split.screen( figs = c( 2, 1 ) )
split.screen( figs = c( 1, 2 ), screen = 2 )
palett <- c("red","black","blue","green","orange","cyan","grey","yellow")

screen(1)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(predMass.mean~Time|combotrt,data=g.trt,col=palett,
       legendwhere="topleft", cex=0.5, pch=3, xaxt='n', ylab="", type="o",xlab="",
       panel.first=adderrorbars(x=g.trt$Time,y=g.trt$predMass.mean,
                                SE=g.trt$predMass.standard.error,direction="updown",
                                col=c("red","blue","orange","grey","black","green","cyan","yellow")))
mtext(side=2, text="Total biomass (g)",  line=3)
mtext(side=1, text="Time (Days)",  line=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)

screen(3)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(ber~Time,data=SBER,col="black",pch = c(1,16)[as.numeric(SBER$Range)],
       legend=F,type="p",  ylim=c(0.8,2), xlim=c(0,67), ylab="",)
mtext(text=expression(M[W]:M[H]),side=2,line=3)
mtext(side=3, text="South")
abline(h=1)

screen(4)
par(mar=c(2,2,2,0),oma=c(2,2,2,2))
plotBy(ber~Time,data=NBER,col="black",pch = c(1,16)[as.numeric(NBER$Range)],
       legend=F,type="p",  ylim=c(0.8,2), xlim=c(0,67), yaxt='n')
mtext(side=3, text="North")
abline(h=1)
legend(45,2, legend=c("Wide","Narrow"), pch=c(16,1), col="black", cex=0.8)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#Mass, RGR and AGR over Time per species
g.trt <- summaryBy(dydt+predMass+AGR~Species+Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Species+Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))

palettS <- c("black","blue","green","cyan","purple"); palettN <- c("black","orange","tan4","red","magenta")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|Species,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|Species,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|Species,data=sg.trt,col=palettS,,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(dydt.mean~Time|Species,data=ng.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))
plotBy(AGR.mean~Time|Species,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))


plotBy(AGR.mean~Time|Species,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|Species,data=srate.trt,col=palettS,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|Species,data=nrate.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------


#Mass, RGR and AGR over Time per taxa
g.trt <- summaryBy(dydt+predMass+AGR~Taxa+Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
rate.trt <- summaryBy(LAR+ULR+predMass+dydt+AGR~Taxa+Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
rate.trt$combotrt <- as.factor(paste(rate.trt$Location,rate.trt$Range,rate.trt$Treatment,sep="_"))


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))

palettS <- c("black","black","blue","blue","blue","green","cyan","purple"); palettN <- c("black","black","black","orange","orange","orange","tan4","red","magenta")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
srate.trt<- droplevels(subset(rate.trt, Location == "S"));nrate.trt<- droplevels(subset(rate.trt, Location == "N"))

plotBy(predMass.mean~Time|Taxa,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(predMass.mean~Time|Taxa,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,90),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

plotBy(dydt.mean~Time|Taxa,data=sg.trt,col=palettS,,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(dydt.mean~Time|Taxa,data=ng.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(0.04,0.15),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(RGR~(g~g^-1~day^-1)))
mtext(text="Day", side=1, line=3)


windows(60,40);par(mfrow=c(2,2), mar=c(0.5,6,0.5,0.5), oma=c(5,1,1,1))
plotBy(AGR.mean~Time|Taxa,data=sg.trt,col=palettS,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$AGR.mean,
                                SE=sg.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))


plotBy(AGR.mean~Time|Taxa,data=ng.trt,col=palettN,
       legendwhere="topleft", xaxt='n', ylab="", type="p",ylim=c(0,5),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$AGR.mean,
                                SE=ng.trt$AGR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = F, tck = 0.01)
mtext(side=2, line=3, text=expression(AGR~(g~~day^-1)))

plotBy(LAR.mean~Time|Taxa,data=srate.trt,col=palettS,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(sg.trt$Treatment)))],
       panel.first=adderrorbars(x=srate.trt$Time,y=srate.trt$LAR.mean,
                                SE=srate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)

plotBy(LAR.mean~Time|Taxa,data=nrate.trt,col=palettN,
       legend=F, xaxt='n', ylab="", type="p",ylim=c(50,170),pch = c(1,16)[as.numeric(as.factor(as.character(ng.trt$Treatment)))],
       panel.first=adderrorbars(x=nrate.trt$Time,y=nrate.trt$LAR.mean,
                                SE=nrate.trt$LAR.standard.error,direction="updown",
                                col=c("black"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, line=3, text=expression(LAR~(g~~day^-1)))
mtext(text="Day", side=1, line=3)


# ltys = c("22", "44", "13", "1343", "73", "2262",
#          "12223242", "F282", "F4448444", "224282F2", "F1")

windows(8,8);par(mfrow=c(1,1), mar=c(0.5,3,0.5,0.5), oma=c(3,1,1,1))
plotBy(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
       legend=F, xaxt='n', ylab="", type="l",ylim=c(0,70),lty="22",lwd=3,
       panel.first=adderrorbars(x=g.trt$Time,y=g.trt$predMass.mean,
                                SE=g.trt$predMass.standard.error,direction="updown",
                                col="black",0))
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="22",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="44",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="44",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="13",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="13",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="solid",lwd=3)
lines(predMass.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),
       xaxt='n', ylab="", type="l",ylim=c(0,70),lty="solid",lwd=3)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)
mtext(text="Time (Days)", side=1, line=2)
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
legend("topleft", legend=c("Tropical Narrow","Tropical Narrow Warmed","Tropical Wide","Tropical Wide Warmed","Temperate Narrow","Temperate Narrow Warmed",
                           "Temperate Wide","Temperate Wide Warmed"), 
       col=c("black","red","black","red",alpha("black",0.6),alpha("red",0.6),alpha("black",0.6),alpha("red",0.6)),
       lty=c("22","22","44","44","13","13","solid","solid"), lwd=3, cex=0.8)

#Size over time

g.trt <- summaryBy(d2h+TotMass~Time+Treatment+Location+Range,data=dat3,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(18,12);par(mfrow=c(1,2), mar=c(2,0,2,0),oma=c(2,4,2,2))
  plotBy(d2h.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",
         legend=F, xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3,
         panel.first=adderrorbars(x=g.trt$Time,y=g.trt$d2h.mean,
                                  SE=g.trt$d2h.standard.error,direction="updown",
                                  col="black",0))
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  lines(d2h.mean~Time, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = FALSE, tck = 0.01)
  mtext(side=2, text=expression(paste('d'^'2','h',sep=' ')), line=2.5)
  mtext(text="Time (Days)", side=1, line=2)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = T, tck = 0.01)
  legend("topleft", ncol=2, legend=c("Tropical Narrow","Tropical Narrow Warmed","Tropical Wide","Tropical Wide Warmed","Temperate Narrow","Temperate Narrow Warmed",
                                     "Temperate Wide","Temperate Wide Warmed"), 
         col=c("black","red","black","red",alpha("black",0.6),alpha("red",0.6),alpha("black",0.6),alpha("red",0.6)),
         lty=c("22","22","44","44","13","13","solid","solid"), lwd=3, cex=0.8)
  
  #Size over mass
  plotBy(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_narrow_Home"),col="black",yaxt='n',
         legend=F, xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3,
         panel.first=adderrorbars(x=g.trt$TotMass.mean,y=g.trt$d2h.mean,
                                  SE=g.trt$d2h.standard.error,direction="updown",
                                  col="black",0))
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_narrow_Warmed"),col="red",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="22",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_wide_Home"),col="black",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="N_wide_Warmed"),col="red",yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="44",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_narrow_Home"),col=alpha("black",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_narrow_Warmed"),col=alpha("red",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="13",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_wide_Home"),col=alpha("black",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  lines(d2h.mean~TotMass.mean, data=subset(g.trt, combotrt=="S_wide_Warmed"),col=alpha("red",0.6),yaxt='n',
        xaxt='n', ylab="", type="l",ylim=c(0,250),lty="solid",lwd=3)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = FALSE, tck = 0.01)
  mtext(text="TotMass (g)", side=1, line=2)
  axis(side = 1, at = seq(from=0,to=90,by=5), labels = T, tck = 0.01)
#par(xpd=T)

