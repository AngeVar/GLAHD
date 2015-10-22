#make big table with MASS AGR RGR at different timepoints

source("R/GLAHD_gamfits.R")
library(effects)

gamfits2$Location <- factor(gamfits2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2$Sp_RS_EN <- as.factor(with(gamfits2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2$Prov_Sp_EN <- as.factor(with(gamfits2,paste(Taxa,Species)))
gamfits2$Sp_Loc_EN <- as.factor(with(gamfits2,paste(Species,Location)))

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
dat3$LMA<- with(dat3, leafArea/leafMass)            #Leaf Mass per Area
dat3$SLA<- with(dat3, leafMass/leafArea)            #Specific leaf area
dat3$LAR<- with(dat3, leafArea/TotMass)             #Leaf Area Ratio
dat3$RSR<- with(dat3, rootMass/(stemMass+leafMass)) #Root:Shoot Ratio

dat3$Location <- factor(dat3$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat3$Sp_RS_EN <- as.factor(with(dat3,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat3$Prov_Sp_EN <- as.factor(with(dat3,paste(Taxa,Species)))
dat3$Sp_Loc_EN <- as.factor(with(dat3,paste(Species,Location)))
#-------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==7)

fm1.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass <- lme(sqrt(predMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr)) 

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)                 
fm1.rgr <- lme((dydt)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr)) 

ratelar<- subset(rate,Time==1)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
fm1.LAR <- lme((LAR)~Treatment+Range+Location+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(allEffects(fm1.LAR))


dat3<- subset(dat3,Time==1)
fm1.Height <- lme((Height)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat3)
plot(fm1.Height,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.Height,Height~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.Height,Height~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.Height, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.Height$residuals[,1])
anova(fm1.Height)                 
fm1.rgr <- lme((Height)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.Height)) 


# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==10)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass <- lme(sqrt(predMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment+Location+Range+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr)) 

fm1.rgr <- lme(sqrt(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,sqrt(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,sqrt(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme(sqrt(dydt)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

ratelar<- subset(rate,Time==11)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
fm1.LAR <- lme((LAR)~Treatment+Location+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==15)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass2 <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme(sqrt(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)
fm1.agr2 <- lme((AGR)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr2)) 

fm1.rgr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme((dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==20)

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
fm1.mass2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr2 <- lme((AGR)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr2)) 

fm1.rgr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme((dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

ratelar<- subset(rate,Time==20)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==25)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass) 
fm1.mass2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2))

fm1.agr <- lme(log(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)
plot(allEffects(fm1.agr))

fm1.rgr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme((dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

ratelar<- subset(rate,Time==25)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
fm1.LAR <- lme((LAR)~Treatment+Location+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(allEffects(fm1.LAR))
# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==30)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)  
plot(allEffects(fm1.mass)) 

fm1.agr <- lme(log(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
plot(allEffects(fm1.agr)) 

fm1.rgr <- lme(1/(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme(log(dydt)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

ratelar<- subset(rate,Time==32)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==35)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
plot(allEffects(fm1.mass)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr2 <- lme((AGR)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr2)) 

fm1.rgr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme(log(dydt)~Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==40)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
plot(allEffects(fm1.mass)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr2 <- lme((AGR)~Treatment+Location+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr2)) 

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
fm1.rgr2 <- lme(log(dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.rgr2)) 

ratelar<- subset(rate,Time==39)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme((LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==45)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme((AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr)) 

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
plot(allEffects(fm1.rgr)) 

ratelar<- subset(rate,Time==46)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme(log(LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==50)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme(log(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr)) 

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
plot(allEffects(fm1.rgr)) 

# -------------------------------------------------------------------------------------

dat5<- subset(gamfits2,Time==55)

dat5$Location <- factor(dat5$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat5$Sp_RS_EN <- as.factor(with(dat5,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat5$Prov_Sp_EN <- as.factor(with(dat5,paste(Taxa,Species)))
dat5$Sp_Loc_EN <- as.factor(with(dat5,paste(Species,Location)))

fm1.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.mass2 <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme(sqrt(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr))

fm1.rgr <- lme(sqrt(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
plot(allEffects(fm1.rgr))

ratelar<- subset(rate,Time==53)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme(log(LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,LAR~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,LAR~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))

# -------------------------------------------------------------------------------------

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
fm1.mass2 <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.mass2)) 

fm1.agr <- lme(sqrt(AGR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.agr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.agr,(AGR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.agr,(AGR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.agr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.agr$residuals[,1])
anova(fm1.agr)                 
fm1.agr <- lme((AGR)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.agr))

fm1.rgr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,dydt~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,dydt~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)
plot(allEffects(fm1.rgr))

ratelar<- subset(rate,Time==60)
ratelar$Location <- factor(ratelar$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
ratelar$Sp_RS_EN <- as.factor(with(ratelar,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
ratelar$Prov_Sp_EN <- as.factor(with(ratelar,paste(Taxa,Species)))
ratelar$Sp_Loc_EN <- as.factor(with(ratelar,paste(Species,Location)))

fm1.LAR <- lme(log(LAR)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ratelar)
plot(fm1.LAR,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.LAR,log(LAR)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.LAR,log(LAR)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.LAR, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.LAR$residuals[,1])
anova(fm1.LAR) 
plot(allEffects(fm1.LAR))
