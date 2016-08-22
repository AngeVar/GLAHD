###------------------------------------------------------------------------------------------------------------------
#output from GAM model:"Code" "Time" "dydt" "predMass" "Species" "Treatment" "Location" "Taxa" "Range" "AGR" 
source("R/GLAHD_gamfits.R") 
###------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------
#output from allometries: "Species" "Date" "Treatment" "Location" "Taxa" "Pot" "Code" Height" "Diameter" "Comment" "ToFit" "d2h" "Range" "TotMass" "leafArea" "leafMass" "stemMass" "rootMass" "Time" 
dat2 <- return_size_mass_all(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)); names <- names(table(dat2$Code));keeps <- names[which(obs>6)]
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)

#- Calculate some new variables
dat3$LMF<- with(dat3, leafMass/TotMass)             #Leaf Mass Fraction
dat3$SMF<- with(dat3, stemMass/TotMass)             #Stem Mass Fraction
dat3$RMF<- with(dat3, rootMass/TotMass)             #Root Mass Fraction
dat3$LMA<- with(dat3, leafMass/leafArea)            #Leaf Mass per Area
dat3$SLA<- with(dat3, leafArea/leafMass)            #Specific leaf area
dat3$LAR<- with(dat3, leafArea/TotMass)             #Leaf Area Ratio
dat3$RSR<- with(dat3, rootMass/(stemMass+leafMass)) #Root:Shoot Ratio

rate <- merge(gamfits2,dat3,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))[,c(1:10,13:14,17:30)]
rate$ULR <- with(rate,dydt/LAR)
###------------------------------------------------------------------------------------------------------------------

###------------------------------------------------------------------------------------------------------------------
#Harvested data:"Species Date Treatment Location Taxa Pot Code Height Diameter Comment ToFit d2h Range TotMass leafArea leafMass stemMass rootMass Time

data <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
data$Date <- as.Date(data$Date,format="%d/%m/%Y")
data$Totmass <- base::rowSums(data[,11:13]) #total mass is the sum of leaf, stem, and root mass
data$d2h <- with(data,(Diameter/10)^2*(Height)) #calculate d2h in cm3
data$Treat <- as.factor(ifelse(data$Pot < 20, "Home",
                               ifelse(data$Pot>=40,"Pre","Warmed")))
data$Taxa <- factor(data$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                       "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
data[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
data$Range<- as.factor(ifelse(data$Species == "TER"|data$Species == "CAM", "wide","narrow"))
data2<-droplevels(subset(data,Treat!="Pre"))

###------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------
#Gasexchange data

Rdark<-getRdark()                                  
Asat<- getAsat()                                   
source("R/GLAHD_ACIfit.R") #this takes a while     

SLARd<- Rdark[,c("Code","SLA")]                    #extract the codes we have leaf SLA for
acifits<-merge(acifit,SLARd, by="Code")            #give SLA to the codes that we have Jmax and Vcmax data for
Asatm<-merge(Asat,SLARd, by="Code")                #give SLA to the codes that we have Asat data for
###------------------------------------------------------------------------------------------------------------------
###------------------------------------------------------------------------------------------------------------------

#add nesting for mixed model
gamfits2$Location <- factor(gamfits2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2$Sp_RS_EN <- as.factor(with(gamfits2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2$Prov_Sp_EN <- as.factor(with(gamfits2,paste(Taxa,Species)))
gamfits2$Sp_Loc_EN <- as.factor(with(gamfits2,paste(Species,Location)))

rate$Location <- factor(rate$Location,levels=c("S","N")) 
rate$Sp_RS_EN <- as.factor(with(rate,paste(Species,Range)))   
rate$Prov_Sp_EN <- as.factor(with(rate,paste(Taxa,Species)))
rate$Sp_Loc_EN <- as.factor(with(rate,paste(Species,Location)))

data2$Location <- factor(data2$Location,levels=c("S","N")) 
data2$Sp_RS_EN <- as.factor(with(data2,paste(Species,Range)))   
data2$Prov_Sp_EN <- as.factor(with(data2,paste(Taxa,Species)))
data2$Sp_Loc_EN <- as.factor(with(data2,paste(Species,Location)))

###########################################################

#Did warming affect biomass?
T60<-subset(gamfits2, Time==60) 
#see daily stats.R for analysis of other time points

fm1m60 <- lme((predMass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
plot(fm1m60,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m60,(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m60,(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m60, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m60$residuals[,1])
anova(fm1m60)    
plot(effect("Treatment",fm1m60))                    #<.0001 biomass increased with warming
plot(effect("Location",fm1m60))                     #<.0001 biomass higher in N
plot(effect("Treatment:Location",fm1m60), multiline = T) #<.0001 warming increased biomass more in S than N
(66.11438-68.46033)/(68.46033)  #-3.4 % N
(47.92213-32.21786)/(32.21786)  # 48.7% S

plot(effect("Treatment:Range",fm1m60), multiline = T)
(59.29974-50.69512)/(50.69512)  #17% wide
(51.41457-46.68100)/(46.68100)  # 10.1% narrow

plot(effect("Treatment:Range:Location",fm1m60), multiline = T) #<.0001 warming increased biomass more in S than N but range size didn't matter (p=0.28)

(68.63780-69.81766)/(69.81766)  #-1.7% N wide
(61.59325-66.02844)/(66.02844)  # -6.7% Nnarrow
(51.01386-33.72724)/(33.72724)  # 51.3% S wide
(42.38279-29.51355)/(29.51355)  # 43.6% S narrow

############################################################
#biomass response ratios - find a way to estimate CIs

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

# #Fieller
# ratio.CI<-function(SEA,SEB,A,B,Q){
#   seq1<-Q*sqrt(((SEA^2)/(A^2))+((SEB)^2)/(B^2))
#   CI95low <- Q-2.06*seq1
#   CI95high <- Q+2.06*seq1
#   return(c(CI95low,CI95high))
# }
# 
# ratio.CI(0.04100491,0.04786403,0.3405541,0.3644084,1.0127511)
# 
# hist(log10(subset(gamfits2,Treatment=="Home"&Location=="S"&Range == "narrow"&Time==4)$predMass))

##############################################################

#Did warming affect RGR?
T15<-subset(gamfits2, Time==30) #change day to test other timepoints
fm1.rgr <- lme(log10(dydt)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T15)
plot(fm1.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.rgr,sqrt(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1.rgr,sqrt(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1.rgr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.rgr$residuals[,1])
anova(fm1.rgr)  
plot(effect("Treatment",fm1.rgr))                        #<.0001 increased
plot(effect("Location",fm1.rgr),multiline = T)           #0.0053 N higher than S
plot(effect("Treatment:Location",fm1.rgr),multiline = T) #<.0001 up in S but not in N
((0.3281428^2)-(0.3289495^2))/(0.3289495^2)              #-0.4 % in N
((0.3207495^2)-(0.2887872^2))/(0.2887872^2)              #23.4 % in S

plot(effect("Treatment:Range",fm1.rgr))                  #0.0382 up more in wide than narrow
((0.3277034^2)-(0.3070671^2))/(0.3070671^2)              #13.9 % in wide
((0.3179942^2)-(0.3087487^2))/(0.3087487^2)              #6.1 % in narrow

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




