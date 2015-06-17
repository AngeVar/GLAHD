#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment.
#  - NOTE- you can avoid the long model fitting process. Just skip to the bottom and read in the csv file.
#
#  - This script fits the full 3-parameter model using an evolutionary algorithm to fully explore the parameter space and 
#     (hopefully) avoid local minima to find the true global minimum (find the "true"
#     parameter values).
#  - The power function has a discontinuity at values of 
#    beta just above 1.0 (say, 1.001), where the predicted mass values are undefined.
#    This means that the evolutionary algorithm must be restricted out of that space, 
#    so I set an upper limit to the allowed beta values at 0.999 . This prevents our plants
#    from growing at rates that exceed the true exponential, which is fine.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")
library(DEoptim)

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
#DEfits <- merge(outpar,subset(dat3,Date==as.Date("2014-11-17")),by="Code")

#- write out a csv for future work, as the evolutionary algorithm takes a long time.
#write.csv(DEfits,row.names=F,file="C:/Repos/GLAHD/R/DE_power_fits.csv")


DEfits.all <- merge(outpar,dat3,by="Code")
DEfits.all$predmass <- with(DEfits.all,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
DEfits.all$AGR <- with(DEfits.all,r*(M0^(1-beta) + r*Time*(1-beta))^(beta/(1-beta)))
DEfits.all$RGR <- DEfits.all$AGR/DEfits.all$predmass
DEfits.all$leafarea <- predict_LA(Taxa=DEfits.all$Taxa,Treat=DEfits.all$Treatment,Diameter=DEfits.all$Diameter,Height=DEfits.all$Height)
DEfits.all$LAR <- with(DEfits.all,leafarea/TotMass)
DEfits.all$ULR <- with(DEfits.all,RGR/LAR)

#write.csv(DEfits.all,row.names=F,file="C:/Repos/GLAHD/Output/DE_power_fits_all.csv")
# plot(TotMass~predmass,data=DEfits,log="xy")
# abline(0,1)
# summary(lm(TotMass~predmass,data=DEfits.all))

#---- start from here-------------------------------------------------------------------------------------------

#DEfits.all <- read.csv("C:/Repos/GLAHD/Output/DE_power_fits_all.csv")
DEfits.all$Date<- as.Date(DEfits.all$Date)
#Extract data from day 17 (the second harvest) i.e. Day 14 of the experiment

dat4<- subset(DEfits.all,Date==as.Date("2014-11-17"))

library(nlme)
library(effects)
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

fm1agr <- lme(log(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4,method="REML")
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


#extract data from mass 10 g?


#some plots over mass
windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(DEfits.all,Location=="S")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5), cols = c("black", "red"))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
mtext("South",side=3, line=2, outer=T, cex=1.4)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5), cols = c("black", "red"))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.4), cols = c("black", "red"))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.4), cols = c("black", "red"))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375), cols = c("black", "red"))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375), cols = c("black", "red"))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0032), cols = c("black", "red"))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0032), cols = c("black", "red"))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_south.pdf")


windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(DEfits.all,Location=="N")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5), cols = c("black", "red"))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
mtext("North",side=3, line=2, outer=T, cex=1.4)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5), cols = c("black", "red"))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.4), cols = c("black", "red"))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.4), cols = c("black", "red"))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375), cols = c("black", "red"))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375), cols = c("black", "red"))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0032), cols = c("black", "red"))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0032), cols = c("black", "red"))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_north.pdf")

#Biomass enhancement ratios

warm<-subset(summaryBy(TotMass~Treatment+Date+Range+Location, data=DEfits.all, FUN=mean),Treatment=="Warmed")
home<-subset(summaryBy(TotMass~Treatment+Date+Range+Location, data=DEfits.all, FUN=mean),Treatment=="Home")
BER<- cbind(warm,home)
BER$ber <- BER[,5]/BER[,10]
bioen<-BER[,7:11]
windows(30,30);par(mfrow=c(2,1), oma=c(1,1,1,1), mai=c(0,1,0.5,1))
plotBy(ber~Date|Range,data=subset(bioen, Location == "N"), type="l", lwd=2, col= c("Black","red"), ylab="Biomass enhancement factor", legendwhere= "topleft",ylim=c(0.75,1.3), main="North")
plotBy(ber~Date|Range,data=subset(bioen, Location == "S"), type="l", lwd=2, col= c("Black","red"), ylab="Biomass enhancement factor", legend=F,ylim=c(1,2), main = "South")

#initial mass
allom <- read.csv("C:/Repos/GLAHD/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
pre <- subset(allom, Treatment == "Pre-treatment")
seed_mass<- summaryBy(Totmass~Taxa, FUN=min, data=pre, var.names = "seedmass", keep.names=T) # tried the minimum instead of the mean?
seed_mass<- rename(seed_mass, c("Taxa"="Taxa", "seedmass"="M0"))
seed_mass$Treatment <- as.factor("seed")
predinit<- dat4[,c("M0","Treatment","Taxa")]
init <- rbind(predinit, seed_mass)

windows(20,12);par(mar=c(8,7,1,1))
boxplot(M0~Treatment*Taxa,data=init,col=c("blue","red"),las=2,ylab="initial mass",cex.lab=2)
abline(v=24.5)
legend("topleft", legend=c("Home","Warmed","Min pre-harvest Mass"),colors(levels(init$Treatment)))

#Some other plots.
windows(20,12);par(mar=c(8,7,1,1))

output$Taxa <- factor(output$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=dat4,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=dat4,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)
boxplot(M0~Treatment*Taxa,data=dat4,col=c("blue","red"),las=2,ylab="initial mass",cex.lab=2)

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

