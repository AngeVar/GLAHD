#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment.
#  - NOTE- you can avoid the long model fitting process. Just skip to the bottom and read in the csv file.
#
#  - This script fits the full 3-parameter model using an evolutionary algorithm to fully explore the parameter space and 
#     (hopefully) avoid local minima to find the true global minimum (find the "true" parameter values).
#     
#  - The power function has a discontinuity at values of 
#    beta just above 1.0 (say, 1.001), where the predicted mass values are undefined.
#    This means that the evolutionary algorithm must be restricted out of that space, 
#    so I set an upper limit to the allowed beta values at 0.999 . This prevents our plants
#    from growing at rates that exceed the true exponential, which is fine.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")
library(DEoptim)
library(nlme)
library(effects)
#read in CSV file instead of re-running curve fit

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-7)) #finds first date and labels it as Time 7 i.e. 07112014 is Day 7

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)


set.seed(1234) # set seed for repeatablility
DEoptim.control <- list(VTR = -Inf, strategy = 2,itermax=500,trace=TRUE,CR=0.9)

#----------------------------------------------------------------------------------------------------------------
#- Loop over each plant and estimate the parameters by an evolutionary algorithm. This is slow (5-10 min?).
#- For testing purposes, reduce NP in DE.optim.control to 30. This gives approximately the same answers but is
#-   much much faster. Judging from my system monitor, DEoptim must use parallel processing under the hood.
#----------------------------------------------------------------------------------------------------------------

lower <- log(c(0.0001,0.4,0.001))
upper <- log(c(0.6,0.999,0.5))
dat.l <- split(dat3,dat3$Code)
output <- list()
outpar <- data.frame("Code"=rep(NA,length(dat.l)),"M0"=rep(NA,length(dat.l)),"beta"=rep(NA,length(dat.l)),"r"=rep(NA,length(dat.l)))

#- set up progress bar
pb <- txtProgressBar(min = 0, max = length(dat.l), style = 3)

for(i in 1:length(dat.l)){
  #print(i)
  tofit <- dat.l[[i]]
  output[[i]] <- DEoptim(fit.power3,lower,upper,Time=dat.l[[i]]$Time,TotMass=dat.l[[i]]$TotMass,
                    control=list(VTR = -Inf, strategy = 2,NP=300,itermax=500,trace=F,CR=0.9))
  outpar$Code[i] <- as.character(tofit$Code[1])
  outpar[i,2:4] <- exp(unname(output[[i]]$optim$bestmem))
  setTxtProgressBar(pb, i)
  
}
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------
outpar$Code <-as.factor(outpar$Code)

#- merge it all the treatment code information from dat3
DEfits.all <- merge(outpar,dat3,by="Code")

#DEfits <- merge(outpar,subset(dat3,Date==as.Date("2014-11-17")),by="Code")
#- write out a csv for future work, as the evolutionary algorithm takes a long time.
# DEfits.all$predmass <- with(DEfits.all,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
# DEfits.all$AGR <- with(DEfits.all,r*(M0^(1-beta) + r*Time*(1-beta))^(beta/(1-beta)))
# DEfits.all$RGR <- DEfits.all$AGR/DEfits.all$predmass
#write.csv(DEfits,row.names=F,file="C:/Repos/GLAHD/R/DE_power_fits.csv")
#write.csv(DEfits.all,row.names=F,file="C:/Repos/GLAHD/Output/DE_power_fits_all_complex.csv")

# plot(TotMass~predmass,data=DEfits.all,log="xy")
# abline(0,1)
# summary(lm(TotMass~predmass,data=DEfits.all))

#######-------------------------------------------------------------------------------------------------------
DEfits.all <- read.csv("C:/Repos/GLAHD/Output/DE_power_fits_all_simple.csv")
DEfits.all$Date<- as.Date(DEfits.all$Date)

#Predict daily mass, AGR and RGR
taxa<- as.data.frame(rep(unique(DEfits.all$Code),each=70));colnames(taxa) <- c("Code"); taxa$Time<- rep(seq(1:70))
code<- merge(taxa, unique(DEfits.all[,c("Code","r","beta","M0", "Treatment","Location","Range","Taxa", "Species")], by="Code"))

code$predmass <- with(code,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
code$AGR <- with(code,r*(M0^(1-beta) + r*Time*(1-beta))^(beta/(1-beta)))
code$RGR <- code$AGR/code$predmass

##-------------------------------------------------------------------------------------------------------
#Biomass enhancement ratio (mean of all data)
ber<- summaryBy(predmass~Time+Range+Location+Treatment, data=code, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]/bio[,5]
BER<- bio[,c(6:11)]
bersum<- summaryBy(ber~Range+Location+Time, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~Time|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,2.5), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(ber~Time|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0,2.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Biomass Enhancement Ratio",side=3,outer=T,cex=1,adj=0.5,line=1)

ber<- summaryBy(predmass+RGR~Time+Range+Location+Treatment, data=code, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,12]/bio[,6]
BER<- bio[,c(1:6,13)]
bersum<- summaryBy(ber~Range+Location+RGR, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~RGR|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,4), xlim=c(0,0.18))
mtext(text="Change in RGR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(ber~RGR|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0,2),yaxt='n', xlim=c(0,0.5))
mtext(text="Home RGR",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Biomass Enhancement Ratio",side=3,outer=T,cex=1,adj=0.5,line=1)


##-------------------------------------------------------------------------------------------------------
#Plot RGR, AGR and total mass over time
S<- summaryBy(predmass+RGR+AGR~Time+Treatment+Range,data=subset(code, Location == "S"),FUN=mean,keep.names=T)
N<- summaryBy(predmass+RGR+AGR~Time+Treatment+Range,data=subset(code, Location == "N"),FUN=mean,keep.names=T)
S$combotrt <- as.factor(paste(S$Range,S$Treatment,sep="_"))
N$combotrt <- as.factor(paste(N$Range,N$Treatment,sep="_"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(RGR~Time|combotrt,data=S,col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,0.5), xlim=c(0,67))
mtext(text="RGR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(RGR~Time|combotrt,data=N,col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0,0.5),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)

plotBy(AGR~Time|combotrt,data=S,col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,4), xlim=c(0,67))
mtext(text="AGR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(AGR~Time|combotrt,data=N,col=c("black","red","blue","orange"),
       legend=F, type="o", main = "North", ylim=c(0,4),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)

plotBy(predmass~Time|combotrt,data=S,col=c("black","red","blue","orange"),
       legendwhere="topleft",type="o", main="South", ylim=c(0,90), xlim=c(0,67))
mtext(text="Total Mass",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(predmass~Time|combotrt,data=N,col=c("black","red","blue","orange"),
       legend=F, type="o", main = "North", ylim=c(0,90),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)

##-------------------------------------------------------------------------------------------------------

#Plot total mass for each taxa for every day of the experinment

codesum <- summaryBy (predmass~Time+Taxa+Treatment+Location, data=code, FUN=mean, keep.names=T)
code.l <- split(codesum,codesum$Taxa)

#- plot each taxa on a separate panel
windows(40,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,80)
xlims <- c(0,70)
palette(c("black","red"))
for (i in 1:length(code.l)){
  toplot <- code.l[[i]]
  plotBy(predmass~Time|Treatment,data=toplot,type="p",pch=15,xlim=xlims,ylim=ylims,axes=F,
         ylab="H",xlab="",legend=F)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-7)

}
mtext("Time",side=1,line=2,outer=T,cex=1.5)
mtext(expression(Total~plant~mass~(g)),side=2,line=2,outer=T,cex=1.5)
legend(x=100,y=60,legend=c("Home","Warmed"),pch=15,cex=1.5,xpd=NA,col=c("black","red"))

##-------------------------------------------------------------------------------------------------------

#get total mass on specific day and test for statistical effects

fmass<- code[code$Time==35,]
fmass$Location <- factor(log(fmass)$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
fmass$Sp_RS_EN <- as.factor(with(fmass,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
fmass$Prov_Sp_EN <- as.factor(with(fmass,paste(Taxa,Species)))
fmass$Sp_Loc_EN <- as.factor(with(fmass,paste(Species,Location)))

fm1FM <- lme(log(predmass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=fmass)
plot(fm1FM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1FM,log(predmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1FM,log(predmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1FM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1FM$residuals[,1])
anova(fm1FM)                 #3-way interaction almost there

plot(allEffects(fm1FM))      #- try to make sense of Treatment:Location:Range interaction

#####----------------------------------------------------------------------------------------------------
# Predict AGR and RGR for every mass

taxa2<- as.data.frame(rep(unique(DEfits.all$Code),each=80));colnames(taxa2) <- c("Code"); taxa2$Mass<- rep(seq(1:80))
code2<- merge(taxa2, unique(DEfits.all[,c("Code","r","beta","M0", "Treatment","Location","Range","Taxa", "Species")], by="Code"))
code2$RGR <- with(code2,r*Mass^(beta-1))
code2$AGR <- code2$RGR*code2$Mass

S2<- summaryBy(RGR+AGR~Mass+Treatment+Range,data=subset(code2, Location == "S"),FUN=mean,keep.names=T)
N2<- summaryBy(RGR+AGR~Mass+Treatment+Range,data=subset(code2, Location == "N"),FUN=mean,keep.names=T)
S2$combotrt <- as.factor(paste(S2$Range,S2$Treatment,sep="_"))
N2$combotrt <- as.factor(paste(N2$Range,N2$Treatment,sep="_"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(RGR~Mass|combotrt,data=S2,col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,0.16), xlim=c(0,82))
mtext(text="RGR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(RGR~Mass|combotrt,data=N2,col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0,0.16),yaxt='n', xlim=c(0,82))
mtext(text="Mass",side=1,outer=T,cex=1,adj=0.5,line=1)

plotBy(AGR~Mass|combotrt,data=S2,col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,1.6), xlim=c(0,22))
mtext(text="AGR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(AGR~Mass|combotrt,data=N2,col=c("black","red","blue","orange"),
       legend=F,type="o", main = "North", ylim=c(0,1.6),yaxt='n', xlim=c(0,22))
mtext(text="Mass",side=1,outer=T,cex=1,adj=0.5,line=1)

#get ARG and RGR at specific mass
spmass<- code2[code2$Mass == 10,]
spmass$Location <- factor(spmass$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
spmass$Sp_RS_EN <- as.factor(with(spmass,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
spmass$Prov_Sp_EN <- as.factor(with(spmass,paste(Taxa,Species)))
spmass$Sp_Loc_EN <- as.factor(with(spmass,paste(Species,Location)))

fm1RGRg <- lme(RGR~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=spmass, method="REML")
plot(fm1RGRg,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1RGRg,RGR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1RGRg,RGR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1RGRg, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1RGRg$residuals[,1])
anova(fm1RGRg)                 #3-way interaction

plot(allEffects(fm1RGRg))      #- try to make sense of Treatment:Location:Range interaction

fm2RGR2g <- lme(RGR~Treatment+Location+Range+Treatment:Location+Treatment:Range,
                   ,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=spmass, method="ML")
anova(fm2RGR2g)
anova(fm1RGRg,fm2RGR10g)
plot(allEffects(fm2RGR10g))

#####----------------------------------------------------------------------------------------------------
#Data on measurement days to estimate LAR and ULR from allometries

DEfits.all$predmass <- with(DEfits.all,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
DEfits.all$AGR <- with(DEfits.all,r*(M0^(1-beta) + r*Time*(1-beta))^(beta/(1-beta)))
DEfits.all$RGR <- DEfits.all$AGR/DEfits.all$predmass
DEfits.all$leafarea <- predict_LA(Taxa=DEfits.all$Taxa,Treat=DEfits.all$Treatment,Diameter=DEfits.all$Diameter,Height=DEfits.all$Height)
DEfits.all$LAR <- with(DEfits.all,leafarea/TotMass)
DEfits.all$ULR <- with(DEfits.all,RGR/LAR)


dat4.m <- summaryBy(TotMass+RGR+AGR+LAR+ULR~Treatment+Location+Range+Time,FUN=mean,data=DEfits.all, keep.names=T)

windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(DEfits.all,Location=="S")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.44))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.44))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0042))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0042))
mtext(side=1,"Plant mass",line=3)

#And North
windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(DEfits.all,Location=="N")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.44))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.44))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="LAR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="LAR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="ULR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0042))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="ULR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0042))
mtext(side=1,"Plant mass",line=3)

#-----------------------------------------------------------------------------------------------


#Extract data from day 17 (the second harvest) i.e. Day 14 of the experiment

dat4<- subset(DEfits.all,Date==as.Date("2014-11-17"))

dat4$Location <- factor(dat4$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat4$Sp_RS_EN <- as.factor(with(dat4,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat4$Prov_Sp_EN <- as.factor(with(dat4,paste(Taxa,Species)))
dat4$Sp_Loc_EN <- as.factor(with(dat4,paste(Species,Location)))

fm1r.l <- lme(log(r)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1r.l,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r.l,log(r)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r.l,log(r)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r.l, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r.l$residuals[,1])
anova(fm1r.l)                 #3-way interaction almost there

plot(allEffects(fm1r.l))      #- try to make sense of Treatment:Location:Range interaction

fm1b <- lme(beta~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1b,resid(.,type="p")~fitted(.) | Treatment,abline=0)       #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1b,beta~fitted(.)|Species,abline=c(0,1))                   #predicted vs. fitted for each species
plot(fm1b,beta~fitted(.),abline=c(0,1))                           #overall predicted vs. fitted
qqnorm(fm1b, ~ resid(., type = "p"), abline = c(0, 1))            #qqplot to assess normality of residuals
hist(fm1b$residuals[,1])
anova(fm1b)                   # three-way interaction significant
plot(allEffects(fm1b))

fm1agr <- lme(log(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1agr,log(AGR)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1agr,log(AGR)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1agr, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1agr$residuals[,1])
anova(fm1agr)
plot(allEffects(fm1agr))

fm1rgr <- lme(RGR~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1rgr,RGR~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1rgr,RGR~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1rgr, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1rgr$residuals[,1])
anova(fm1rgr)
plot(allEffects(fm1rgr))    


fm1Tmass <- lme(TotMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1Tmass,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1Tmass,TotMass~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1Tmass,TotMass~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1Tmass, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1Tmass$residuals[,1])
anova(fm1Tmass)
plot(allEffects(fm1Tmass)) 

fm1LAR <- lme(log(LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4,,method="REML")
plot(fm1LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LAR,log(LAR)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1LAR,log(LAR)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1LAR, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1LAR$residuals[,1])
anova(fm1LAR)
plot(allEffects(fm1LAR)) 

fm2LAR <- lme(log(LAR)~Treatment+Location+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4,,method="ML")
anova(fm2LAR)
anova(fm1LAR,fm2LAR)
plot(allEffects(fm2LAR)) 

fm1M0 <- lme(sqrt(M0)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4,method="REML")
plot(fm1M0,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1M0,sqrt(M0)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1M0,sqrt(M0)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1M0, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1M0$residuals[,1])
anova(fm1M0)

fm2M0 <- lme(sqrt(M0)~Treatment+Location+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4,method="ML")
anova(fm2M0)
plot(allEffects(fm2M0)) 
anova(fm1M0,fm2M0) #warmed trees had lower initial mass

#Make an anova table
extract.lme <- function(model){
  mod.o <- anova(model)
  
  ndf <- mod.o[2:nrow(mod.o),1]
  ddf <- mod.o[2:nrow(mod.o),2]
  df <- (paste(ndf,ddf,sep=", "))
  Fvalue <-round(mod.o[2:nrow(mod.o),3],1)
  Pvalue <-round(mod.o[2:nrow(mod.o),4],3)
  Pvalue[which(Pvalue==0)] <- "<0.001"
  return(data.frame("df"=df, "F" = Fvalue,"P" = Pvalue))
}


#- make a big table
growthtable <- do.call(cbind,lapply(list(fm1r.l,fm1b,fm1agr,fm1rgr,fm1Tmass),FUN=extract.lme))
row.names(growthtable) <- row.names(anova(fm1r.l))[2:8]
#headings<-c("Terms", "log(r)", "b","AGR Day14", "RGR Day14", "Mass Day14")

library(R2wd)
wdGet()
wdTable(growthtable,autoformat=2)

###--------------------------------------------------------------------------------------------------------
#Does differences in m0 matter?
Time<- rep(seq(0,50),6)
m0 <- rep(seq(0,0.5, length.out=6),each=51)
r <- rep(0.08, 306)
beta <- rep(0.88, 306)
data<- data.frame(Time,m0,r,beta)
data$mass<- with(data, (m0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
data$m0<- as.factor(data$m0)
col<- c("red","black","blue","green","orange","cyan","grey")
plotBy(mass~Time|m0, data=data, col= col,legendwhere="topleft")
mtext(text="r=0.08   beta= 0.88", side=3, line = 1)

###--------------------------------------------------------------------------------------------------------

#LAR
S3<- summaryBy(LAR~Time+Treatment+Range,data=subset(DEfits.all, Location == "S"),FUN=mean,keep.names=T)
N3<- summaryBy(LAR~Time+Treatment+Range,data=subset(DEfits.all, Location == "N"),FUN=mean,keep.names=T)
S3$combotrt <- as.factor(paste(S3$Range,S3$Treatment,sep="_"))
N3$combotrt <- as.factor(paste(N3$Range,N3$Treatment,sep="_"))

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(LAR~Time|combotrt,data=S3,col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0,350), xlim=c(0,70))
mtext(text="LAR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(LAR~Time|combotrt,data=N3,col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0,350),yaxt='n', xlim=c(0,70))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)

S4<- summaryBy(LAR~TotMass+Treatment+Range,data=subset(DEfits.all, Location == "S"),FUN=mean,keep.names=T)
N4<- summaryBy(LAR~TotMass+Treatment+Range,data=subset(DEfits.all, Location == "N"),FUN=mean,keep.names=T)
S4$combotrt <- as.factor(paste(S4$Range,S4$Treatment,sep="_"))
N4$combotrt <- as.factor(paste(N4$Range,N4$Treatment,sep="_"))
windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)
plotBy(LAR~TotMass|combotrt,data=S4,col=c("black","red","blue","orange"),
       legend=F,type="p", main="South", ylim=c(0,350), xlim=c(0,80))
mtext(text="LAR",side=2,outer=T,cex=1,adj=0.5,line=3)
plotBy(LAR~TotMass|combotrt,data=N4,col=c("black","red","blue","orange"),
       legendwhere="topright",type="p", main = "North", ylim=c(0,350),yaxt='n', xlim=c(0,80))
mtext(text="Mass",side=1,outer=T,cex=1,adj=0.5,line=1)


# #LAR=SLA (LA/LM) * LMF (LM/TotM)
# DEfits.all
# allom <- read.csv("C:/Repos/GLAHD/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
# allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
# allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
# allom$logd2h <- log10(allom$d2h)
# allom$logLA <- log10(allom$Leafarea)
# allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
#                                 ifelse(allom$Pot>=40,"Pre","Warmed")))
# allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
#                                          "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
# 
# linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
# allom <- merge(allom,linkdf,by="Species")
# allom$Date <- as.Date(allom$Date, format = "%d/%mm/%yy")
# allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
# allom$LMF <- allom$Leafmass/allom$Totmass
# allom$SLA <- (10^allom$logLA)/allom$Leafmass
# 
# dat5<- subset(allom,Date==as.Date("2014-11-17"))
# dat6<- subset(allom,Date==as.Date("2014-11-07"))

#Some other plots.
windows(20,12);par(mar=c(8,7,1,1))

DEfits.all$Taxa <- factor(DEfits.all$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=DEfits.all,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=DEfits.all,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)
boxplot(M0~Treatment*Taxa,data=DEfits.all,col=c("blue","red"),las=2,ylab="initial mass",cex.lab=2)

plotBy(r~M0|Treatment,data=dat4,legendwhere="topright")
plotBy(r~beta|Treatment,data=dat4,legendwhere="topright")

#- dat4 compare predicted with measured biomass

DEfits.all$predmass <- with(DEfits.all,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))

windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(predmass~TotMass,data=DEfits.all,pch=3,col="red",cex.lab=1.5,
     xlab="Observed mass (g)",ylab="Predicted mass (g)")
abline(0,1)

#- log scale
plot(predmass~TotMass,data=DEfits.all,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed mass (g)",ylab="Predicted mass (g)")
abline(0,1)

summary(lm(predmass~TotMass,data=DEfits.all))


rates.trt <- summaryBy(TotMass+RGR+AGR~Time+Treatment+Location+Range,data=DEfits.all,FUN=mean,keep.names=T)
rates.trt$combotrt <- as.factor(paste(rates.trt$Location,rates.trt$Range,rates.trt$Treatment,sep="_"))
plotBy(RGR~Time|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright", type="o")
plotBy(AGR~Time|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft", type="o")
plotBy(TotMass~Time|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft", type="o")

