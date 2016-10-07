
###########################################################

#Did warming affect biomass or RGR?
source("R_cleaned/Table 2.R")
source("R_cleaned/Figure 2.R")
source("R_cleaned/Figure 3.R")
source("R_cleaned/Figure 4.R")

############################################################
#Maximum biomass enhancement

rat<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, 
                FUN=c(mean,standard.error,length), keep.names=T)
rath <- subset(rat, Treatment == "Home")
ratw <- subset (rat, Treatment == "Warmed")
names(ratw)<- c("Time","Range","Location","Treatment","predMassWarm","predMass.se.warm","predMass.length")
resp <- merge(rath,ratw, by=c("Time","Range","Location"))
resp$ber<- with(resp, predMassWarm/predMass.mean)
SBER<- subset(resp,Location =="S")
NBER<- subset(resp,Location =="N")
SBER[ SBER$ber %in% tapply(SBER$ber, SBER$Range, max), ]
NBER[ NBER$ber %in% tapply(NBER$ber, NBER$Range, max), ]

#RGRmax (Maximum RGR was increased by warming, moreso in temperate taxa)
maxdydt<-gamfits2[ gamfits2$dydt %in% tapply(gamfits2$dydt, gamfits2$Code, max), ]#max RGR per Code

maxdydt$Location <- factor(maxdydt$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
maxdydt$Sp_RS_EN <- as.factor(with(maxdydt,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
maxdydt$Prov_Sp_EN <- as.factor(with(maxdydt,paste(Taxa,Species)))
maxdydt$Sp_Loc_EN <- as.factor(with(maxdydt,paste(Species,Location)))

fm.maxdydt<- lme(log(dydt)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt,log(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt,log(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt$residuals[,1])
anova(fm.maxdydt)

plot(effect("Treatment",fm.maxdydt))                      #<.0001 increased with warming
plot(effect("Location",fm.maxdydt))                       #0.0188 N higher than S
plot(effect("Treatment:Location",fm.maxdydt),multiline=T) # <.0001 increase more in S than N
(exp(-2.176796)-exp(-2.190548))/(exp(-2.190548))  # 1.3 %  N
(exp(-2.223111)-exp(-2.403854))/(exp(-2.403854))  # 19.8 %  S

#MassRGRmax (the mass at which the maxRGR occurred decreased among tropical taxa, i.e. the peak occurred earlier)
fm.maxdydt2<- lme(log(predMass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt2,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt2,log(predMass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt2,log(predMass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt2, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt2$residuals[,1])
anova(fm.maxdydt2)  
plot(effect("Treatment",fm.maxdydt2))                 #0.0001 decreased with warming
plot(effect("Treatment:Location",fm.maxdydt2))        #0.0024 decreased more in N than S
(exp(0.0927346)-exp(0.8952484))/(exp(0.8952484))  # -55 %  N
(exp(0.3357711)-exp(0.4649718))/(exp(0.4649718))  # -12 %  S

#Did warming affect RGR at common size?
summary(gamfits2$predMass) #mean is 14.81g
gamfits15<-subset(gamfits2,predMass>14.67462&predMass<14.94554)#mean plus and minus SE
rgr15<- summaryBy(dydt~Range+Location+Treatment, data=gamfits15, 
                FUN=c(mean,standard.error,length), keep.names=T)

fm.maxdydt2<- lme(log10(dydt)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=gamfits15)
plot(fm.maxdydt2,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt2,(dydt)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt2,(dydt)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt2, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt2$residuals[,1])
anova(fm.maxdydt2)  
plot(effect("Location",fm.maxdydt2))                 #0.0025 Higher in N than S




