#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment
#------------------------------------------------------------------------------------------------------------------------------
# Currently has an issue with fitting the power law function

#- load libraries from script
source("R/loadLibraries.R")

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)

# from Paine et al. 2012
###################################
# 2-parameter Power-law fit, initial mass fixed
###################################

seed.mass <- 0.4 #(45 mg, Hautier et al 2010, J. Ecology) # change this to suit!
fmla.pow2 <- as.formula(paste("~(", seed.mass, "^(1-beta) + r*x*(1-beta))^(1/(1-beta))", sep = ""))
SS.pow2   <- selfStart(fmla.pow2, initial = Init.pow2, parameters = c("r", "beta"))
xvals <- seq(1,70,by=0.5) # establish xvalues of time to use in the model fitting below

source("R/Paine.R")
##############################################################################################
##############################################################################################
#- fit growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
xlims <- c(0,70)
ylims <- c(0,100)

pdf(file="C:/Repos/GLAHD/Output/Growth_power_allplants.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)

output <- data.frame(Species="",Treatment="",Location="",Taxa="",Code="",Range="",r="",beta="",AGR10g="",AGR10d="",RGR10d="",RGR10g="",stringsAsFactors=F)
full.output <- list()
for (i in 1:length(growth.l)){
  print(i)
  tofit <- growth.l[[i]]
  tofit <- subset(tofit,TotMass>seed.mass)
  #fit
  tmp.pow2 <- getInitial(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  fit.pow2 <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  out.pow2 <- output.pow2.gnls(fit.pow2,xvals, CI = T, M0 = seed.mass)
  
  #plot
  plot(TotMass~Time,data=tofit,pch=20,cex=3,ylab="Total plant mass (g)",xlab="Time (days)",xlim=xlims,ylim=ylims)
  lines(M~times,data=out.pow2$rates)
  mtext(text=tofit$Code[1],side=3,line=-1)
  legend("topleft",paste("r = ",round(coef(fit.pow2)[1],digits=4)),bty="n")
  legend("left",paste("beta = ",round(coef(fit.pow2)[2],digits=4)),bty="n")
  
  #compile parameters
  output[i,1:6] <- c(as.character(tofit$Species[1]),as.character(tofit$Treatment[1]),as.character(tofit$Location[1]),
                     as.character(tofit$Taxa[1]),as.character(tofit$Code[1]),as.character(tofit$Range[1]))
  output[i,"r"] <- coef(fit.pow2)[1]
  output[i,"beta"] <- coef(fit.pow2)[2]
  index_time <- findInterval(10,out.pow2$rates$times) # this was altered from "M" to "times" to get AGR and RGR at 10 days, rather than 10 grams
  index_mass <- findInterval(10,out.pow2$rates$M) # this was altered from "M" to "times" to get AGR and RGR at 10 days, rather than 10 grams
  
  output[i,"AGR10g"] <- mean(out.pow2$rates$AGR[index_mass],out.pow2$rates$AGR[index_mass+1]) # get AGR at 10 grams
  output[i,"AGR10d"] <- out.pow2$rates$AGR[index_time]           # get AGR at 10 days
  output[i,"RGR10d"] <- out.pow2$rates$RGRt[index_time]         # get RGT at 10 days
  output[i,"RGR10g"] <- mean(out.pow2$rates$RGRm[index_mass],out.pow2$rates$RGRm[index_mass+1])# get RGT at 10 grams
  
  #- grabs all output of the "rates' dataframe
  full.output[[i]] <- out.pow2$rates
  full.output[[i]]$Species <- tofit$Species[1]
  full.output[[i]]$Treatment <- tofit$Treatment[1]
  full.output[[i]]$Location <- tofit$Location[1]
  full.output[[i]]$Taxa <- tofit$Taxa[1]
  full.output[[i]]$Pot <- tofit$Pot[1]
  full.output[[i]]$Code <- tofit$Code[1]
  full.output[[i]]$Range <- tofit$Range[1]
  
  
}
dev.off()
output[,7:12] = apply(output[,7:12], 2, function(x) as.numeric(as.character(x)))
output$Species <- as.factor(output$Species)
output$Treatment <- as.factor(output$Treatment)
output$Location <- as.factor(output$Location)
output$Taxa <- as.factor(output$Taxa)
output$Code <- as.factor(output$Code)
output$Range <- as.factor(output$Range)

#- merge full.output into a dataframe
#- use this to plot response ratio of biomass to warming over time? As Miko did?
full.output.df <- do.call(rbind,full.output)

#- get the estimate of final mass to add to this output dataframe
dat.f <- subset(dat2,Date==as.Date("2015-01-05"))[,c(1,3,4,5,7,13)]
names(dat.f)[6] <- "TotMass.f"
 
output2 <- merge(output,dat.f,by=c("Species","Treatment","Location","Taxa","Code"))

#- write out a csv
#write.csv(output2,file="Data/GLAHD_output_growth_parameters_15052015.csv",row.names=F)
output <- read.csv("Data/GLAHD_output_growth_parameters_15052015.csv")

#Some plots
windows(20,12);par(mar=c(8,7,1,1))

output$Taxa <- factor(output$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)
boxplot(AGR10d~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="AGR",cex.lab=2)
abline(v=16.5)
boxplot(RGR10g~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="RGR",cex.lab=2)
abline(v=16.5)

##############################################################################################
##############################################################################################
##############################################################################################

#Statistical testing of growth parameters

library(nlme)
library(effects)
output2$Location <- factor(output2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
output2$Sp_RS_EN <- as.factor(with(output2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
output2$Prov_Sp_EN <- as.factor(with(output2,paste(Taxa,Species)))
output2$Sp_Loc_EN <- as.factor(with(output2,paste(Species,Location)))
#-----------------------------------------------
##model from meeting: fm1 <- lme(TotMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat.f) 

fm1r.l <- lme(log(r)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
plot(fm1r.l,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r.l,log(r)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1r.l,log(r)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1r.l, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1r.l$residuals[,1])
anova(fm1r.l)    

plot(allEffects(fm1r.l))      #- try to make sense of Treatment:Location:Range interaction
#------------------------------------------------

#------------------------------------------------
fm1b <- lme(beta~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1b,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1b,beta~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1b,beta~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1b, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1b$residuals[,1])
anova(fm1b)# three-way interaction significant

#------------------------------------------------

#------------------------------------------------
fm1agr10 <- lme(log(AGR10g)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2,method="REML")
#look at model diagnostics
plot(fm1agr10,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1agr10,log(AGR10g)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1agr10,log(AGR10g)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1agr10, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1agr10$residuals[,1])

anova(fm1agr10)# 3way interaction not significant? - drop?

fm2agr10 <- lme(log(AGR10g)~Treatment+Location+Range+Treatment:Location+Treatment:Range+Location:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
                method="ML",data=output2)
anova(fm1agr10,fm2agr10)
plot(allEffects(fm7agr10))
#------------------------------------------------

#------------------------------------------------
fm1rgr10 <- lme(RGR10d~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1rgr10,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1rgr10,RGR10d~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1rgr10,RGR10d~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1rgr10, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1rgr10$residuals[,1])
anova(fm1rgr10)
plot(allEffects(fm1rgr10))     
summary(fm1rgr10)
#------------------------------------------------

#------------------------------------------------
fm1Tmass <- lme(TotMass.f~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1Tmass,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1Tmass,TotMass.f~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1Tmass,TotMass.f~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1Tmass, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1Tmass$residuals[,1])
anova(fm1Tmass)
plot(allEffects(fm1Tmass))     
summary(fm1Tmass)

#------------------------------------------------

#------------------------------------------------

fm1agr10d <- lme(log(AGR10d)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2,method="REML")
#look at model diagnostics
plot(fm1agr10d,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1agr10d,log(AGR10d)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1agr10d,log(AGR10d)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1agr10d, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1agr10$residuals[,1])

anova(fm1agr10d)# 3way interaction significant?

#------------------------------------------------

#------------------------------------------------
fm1rgr10g <- lme(log(RGR10g)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1rgr10g,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1rgr10g,log(RGR10g)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1rgr10g,log(RGR10g)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1rgr10g, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1rgr10g$residuals[,1])
anova(fm1rgr10g)
plot(allEffects(fm1rgr10g))     
summary(fm1rgr10g)



##############################################################################################
##############################################################################################
# RGR decomposition
##############################################################################################
##############################################################################################
#-- merge r and beta estimates for each plant with the size dataframe, in attempt to 
#-- decompose RGR into its components
dat4 <- merge(dat3,output2,by=c("Species","Treatment","Location","Taxa","Code","Range"))
dat4$RGR <- with(dat4,r*TotMass^(beta-1)) # calculate RGR on mass basis. See Paine et al. 2012.
dat4$AGR <- with(dat4,r*(0.5^(1-beta)+r*Time*(1-beta))^(beta/(1-beta))  )
# use function to predict total plant leaf area from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
dat4$leafarea <- predict_LA(Taxa=dat4$Taxa,Treat=dat4$Treatment,Diameter=dat4$Diameter,Height=dat4$Height)
dat4$LAR <- with(dat4,leafarea/TotMass)
#dat4$ULR <- with(dat4,AGR/leafarea)
dat4$ULR <- with(dat4,RGR/LAR)

dat4.m <- summaryBy(TotMass+RGR+AGR+LAR+ULR~Treatment+Location+Range+Time,FUN=c(mean,standard.error),data=dat4)
##############################################################################################


##############################################################################################
#-- plot AGR, RGR, LAR, and ULR over mass
windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(dat4,Location=="S")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.15))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.15))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0022))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0022))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_south2.pdf")
##############################################################################################

windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(dat4,Location=="N")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.3))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.3))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,425))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,425))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0035))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0035))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_north.pdf")

##############################################################################################
##############################################################################################
##############################################################################################
dat5 <- merge(dat3,full.output.df,by=c("Species","Treatment","Location","Taxa","Code","Range"))


dat4 <- merge(dat3,output2,by=c("Species","Treatment","Location","Taxa","Code","Range"))
dat4$RGR <- with(dat4,r*TotMass^(beta-1)) # calculate RGR on mass basis. See Paine et al. 2012.
dat4$AGR <- with(dat4,r*(0.5^(1-beta)+r*Time*(1-beta))^(beta/(1-beta))  )
# use function to predict total plant leaf area from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
dat4$leafarea <- predict_LA(Taxa=dat4$Taxa,Treat=dat4$Treatment,Diameter=dat4$Diameter,Height=dat4$Height)
dat4$LAR <- with(dat4,leafarea/TotMass)
#dat4$ULR <- with(dat4,AGR/leafarea)
dat4$ULR <- with(dat4,RGR/LAR)

dat4.m <- summaryBy(TotMass+RGR+AGR+LAR+ULR~Treatment+Location+Range+Time,FUN=c(mean,standard.error),data=dat4)

full.output.df$leafarea <- predict_LA(Taxa=full.output.df$Taxa,Treat=full.output.df$Treatment,Diameter=full.output.df$Diameter,Height=full.output.df$Height)

