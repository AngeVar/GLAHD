#######################
# GLAHD Growth Analysis
#######################

#######################
# Load all the data
#######################


###------------------------------------------------------------------------------------------------------------------
#output from GAM model:"Code" "Time" "dydt" "predMass" "Species" "Treatment" "Location" "Taxa" "Range" "AGR" 
source("R/GLAHD_gamfits.R") 
library(effects)
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

#gamfits2 has data interpolated for every day

#rate has allometrically estimated data only for the days on which we measured Height and Diameter: 
#Day 1 11 20 25 32 39 46 53 60 plus dydt (RGR), AGR and predmass from the GAM 

#data2 has harvested, actually measured data
#Rdark has Rmass and SLA for gax exchange leaves
#Asat has Asat at growth T
#acifit has Vcmax, Jmax and JtoV
#acifits has Vcmax and Jmax on a mass basis

###########################################################

#Did warming affect Diameter and Height?
#Generally Tropical taxa responded very modestly to warming with slight enhancement of height and diameter early and slight reductions towards the end.

T60<-subset(rate, Time==60) #change time to check differences at other time points

fm1size <- lme((d2h)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
plot(fm1size,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1size,(d2h)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1size,(d2h)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1size, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1size$residuals[,1])
anova(fm1size)    
plot(effect("Treatment",fm1size))                            #<.0001 increases by warming
plot(effect("Range",fm1size))                                #0.0572 Wide more than narrow
plot(effect("Location",fm1size))                             #0.0029 N more than S
plot(effect("Treatment:Range",fm1size),multiline=T)          #0.0037 increased more by warming in wide
plot(effect("Treatment:Location",fm1size),multiline=T)       #<.0001 increased more by warming in the S
plot(effect("Treatment:Range:Location",fm1size),multiline=T) #0.0689 Size increased more by warming in the S

effect("Treatment:Range:Location",fm1size)
#(Warmed-Home)/(Home)        
(229.4104-236.7054)/(236.7054)  #-3.08 %  Nw
(164.5839-186.6049)/(186.6049)  #-11.8 %  Nn
(205.2237-120.5232)/(120.5232)  # 70.3 %  Sw
(78.91821-51.43743)/(51.43743)  # 53.4 %  Sn

###########################################################

#Did warming affect biomass?
T60<-subset(gamfits2, Time==60)

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

plot(effect("Treatment:Range:Location",fm1m60), multiline = T) #<.0001 warming increased biomass more in S than N

(68.63780-69.81766)/(69.81766)  #-1.7% N wide
(61.59325-66.02844)/(66.02844)  # -6.7% Nnarrow
(51.01386-33.72724)/(33.72724)  # 51.3% S wide
(42.38279-29.51355)/(29.51355)  # 43.6% S narrow


#mass at start
T1<-subset(gamfits2, Time==1)
fm1m1 <- lme(sqrt(predMass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T1)#, method="ML")
plot(fm1m1,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m1,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m1,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m1, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m1$residuals[,1])
anova(fm1m1)
plot(effect("Location",fm1m1))                      #0.0586 slightly higher biomass in N
#plot(effect("Location", fm1m1, transformation=list(link=sqrt, inverse=function(x) x^2)))#changes scale on plot
((0.7742488^2)-(0.6423764^2))/(0.6423764^2)  #45.2 % higher in N

#mass at 15 days
T15<-subset(gamfits2, Time==15)
fm1m15 <- lme(sqrt(predMass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T15)#, method="ML")
plot(fm1m15,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m15,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m15,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m15, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m15$residuals[,1])
anova(fm1m15)
plot(effect("Treatment", fm1m15))                             #<.0001 higher in warmed
plot(effect("Location", fm1m15))                              # 0.0010 higher in N
plot(effect("Treatment:Location",fm1m15),multiline = T)       #0.0358 increased more in S than N
((1.639580^2)-(1.544822^2))/(1.544822^2)     #12.6 % N
((1.252935^2)-(1.062812^2))/(1.062812^2)     #39 % S

###########################################################

#Did warming affect RGR?
fm1.rgr <- lme(sqrt(dydt)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T15)
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

#RGRmax
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

#MassRGRmax
fm.maxdydt2<- lme(log(predMass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=maxdydt)
plot(fm.maxdydt2,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.maxdydt2,log(predMass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.maxdydt2,log(predMass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.maxdydt2, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.maxdydt2$residuals[,1])
anova(fm.maxdydt2)  
plot(effect("Treatment",fm.maxdydt2))                 #0.0001 decreased with warming
plot(effect("Treatment:Location",fm.maxdydt2))        #0.0024 decreased more in N then s
(exp(0.0927346)-exp(0.8952484))/(exp(0.8952484))  # -55 %  N
(exp(0.3357711)-exp(0.4649718))/(exp(0.4649718))  # -12 %  S

###########################################################

#Response ratio of Biomass (% increase (W-H)/H)

ber<- summaryBy(predMass~Time+Range+Location+Treatment, data=gamfits2, FUN=mean, keep.names=T)
berw<-ber[ber$Treatment=="Warmed",]; berh<-ber[ber$Treatment=="Home",]
berw$RR<-(berw$predMass-berh$predMass)/berh$predMass
SBER<-berw[berw$Location=="S",];NBER<-berw[berw$Location=="N",]

SBER[ SBER$RR %in% tapply(SBER$RR, SBER$Range, max), ] #Max ratio 96.8 on Day 37 Wide, 47.9 on Day 45 Narrow
NBER[ NBER$RR %in% tapply(NBER$RR, NBER$Range, max), ] #Max ratio 13.3 on Day 16 Wide, 10.6 on Day 12 Narrow

###########################################################

#Did warming change allocation change? (Harvested data over total mass, slope is mass fraction, 
#i.e. Totmass:Treatment, Totmass:treatment:Location, Totmass:treatment:range, Totmass:treatment:Location:range)

#Leaf mass 

fm1LM <- lme(log(Leafmass)~log(Totmass)*Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
plot(fm1LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LM,log(Leafmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LM,log(Leafmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LM$residuals[,1])
anova(fm1LM)    
plot(effect("Treatment",fm1LM))                                     #P=0.0001 Warmed trees had less Leaf mass from the start
plot(effect("log(Totmass):Location",fm1LM),multiline=T)             #P=<.0001 N taxa allocated more to leaf mass early
plot(effect("log(Totmass):Treatment:Location",fm1LM,transformation=
              list(link=log, inverse=exp)),multiline=T)             #P=0.0411 S warmed allocated less to leaf mass 

#Stem mass
fm1SM <- lme(log(Stemmass)~log(Totmass)*Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#,
             #weights=varFunc(~Totmass))#, method="ML")
plot(fm1SM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1SM,log(Stemmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1SM,log(Stemmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1SM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1SM$residuals[,1])
anova(fm1SM)    
plot(effect("Treatment",fm1SM))                                  #P=<.0001 Warmed trees had more stem mass from the start
plot(effect("log(Totmass):Location",fm1SM),multiline=T)          #P=0.0004 N taxa allocated more to stem mass 

#Root mass
fm1RM <- lme(log(Rootmass)~log(Totmass)*Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#,
             #weights=varFunc(~Totmass))#, method="ML")
plot(fm1RM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1RM,log(Rootmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1RM,log(Rootmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1RM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1RM$residuals[,1])
anova(fm1RM)    

plot(effect("Treatment",fm1RM))                                        #P=0.0008 Warmed trees had less root mass from the start
plot(effect("Location",fm1RM))                                         #P=0.0005 N had less root mass from the start
plot(effect("log(Totmass):Treatment",fm1RM),multiline=T)               #P=0.0381 Warmed trees had lower allocation to root mass
plot(effect("log(Totmass):Treatment:Location",fm1RM),multiline=T)      #P=0.0002 N taxa allocated less to root mass whereas S allocated more 



#LEAF AREA (i.e. LAR) #no treatment effects
fm1LA <- lme((Leafarea)~(Totmass)*Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
plot(fm1LA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LA,(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LA,(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LA$residuals[,1])
anova(fm1LA)

plot(effect("Totmass:Range",fm1LA), multiline=TRUE) #P=0.0657 wide had lower leaf area than narrow
plot(effect("Totmass:Location ",fm1LA), multiline=TRUE) #P=0.0882 N had slightly lower leaf area than S towards the end

###########################################################

#Did SLA change?
#Canopy SLA
fm1SLA <- lme(log(Leafarea)~log(Leafmass)*Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
plot(fm1SLA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1SLA,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1SLA,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1SLA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1SLA$residuals[,1])
anova(fm1SLA)

plot(effect("Treatment",fm1SLA), multiline=TRUE)                         #P=0.0929 Warmed had slightly higher SLA
plot(effect("log(Leafmass):Location",fm1SLA), multiline=TRUE)            #P=0.0082 N had slightly lower SLA than S towards the end
plot(effect("Treatment:Location",fm1SLA), multiline=TRUE)                #P=0.0039 WS had lower sla than HS, HN had higher sla than WN
plot(effect("log(Leafmass):Treatment:Location",fm1SLA), multiline=TRUE)  #P=0.0018 Warming increased SLA in S but not N early
plot(effect("log(Leafmass):Range:Location",fm1SLA), multiline=TRUE)      #P=0.0781 SW had lower SLA than SN

#Leaf level gasexchange SLA
fm1SLAleaf <- lme(SLA~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdark)#, method="ML")
plot(fm1SLAleaf,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1SLAleaf,SLA~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1SLAleaf,SLA~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1SLAleaf, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1SLAleaf$residuals[,1])
anova(fm1SLAleaf)

plot(effect("Treatment",fm1SLAleaf), multiline=TRUE)                    #P=0.0055 Warmed had higher SLA
plot(effect("Treatment:Range",fm1SLAleaf), multiline=TRUE)              #P=0.0025 Warming increased SLA in Narrow but not Wide 
                                                                        #if SMIT is removed this decreases to P=0.132
(184.9480-182.1307)/182.1307         #SLA increased by 1.5% in wide
(206.6620-181.4769)/181.4769         #SLA increased by 13.9% in narrow

plot(effect("Treatment:Location",fm1SLAleaf), multiline=TRUE)           #P=<.0001 WS had higher sla than HS, WN had lower sla than NW
(189.1062-195.6220)/195.6220         #SLA decreased by -3.3 % in the N with warming
(196.5071-167.2035)/167.2035         #SLA increased by 17.5 % in the South with warming

plot(effect("Treatment:Range:Location",fm1SLAleaf), multiline=TRUE) #0.3250
(185.5669-197.2897)/197.2897        #-5.9%  NW
(195.5060-192.6063)/192.6063        #1.5%   NN
(184.2853-165.8999)/165.8999        #11.1%  SW
(218.6067-169.5606)/169.5606        #28.9%  SN


# windows(14,14);par(mfrow=c(1,1),mar=c(6,5,2,1),cex.lab=1.3,cex.axis=1.3)
# #colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
# colors <- c("blue","red")
# ylims=c(0,500)#SLA
# #ylims=c(0,150)#Area
# #ylims=c(0,800)#leafDW
# boxplot(SLA~Treatment+Taxa,data=Rdark,ylim=ylims,ylab=expression(SLA~(cm^2~g^-1)),axes=F,las=2,col=colors)
# magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T)
# axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(Rdark$Taxa),las=2)
# abline(v=16.5)


###########################################################

#Photosynthsis

#Get massbased rates
acifits$Jmaxm<-with(acifits,Jmax*SLA/10000)
acifits$Vcmaxm<-with(acifits,Vcmax*SLA/10000)
Asatm$Photom<-with(Asatm,Photo*SLA/10000)

#- Vcmax
fm.vcmax <- lme(log(Vcmax)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifit)
plot(fm.vcmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmax,log(Vcmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmax,log(Vcmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmax$residuals[,1])
anova(fm.vcmax)    
plot(effect("Treatment", fm.vcmax))              #<.0001 warming reduced Vcmax
plot(effect("Location", fm.vcmax))               #<.0001 N had higher Vcmax than S
plot(effect("Range", fm.vcmax))                  #0.0128 Wide had higher Vcmax than narrow
(exp(4.913126)-exp(4.793873))/(exp(4.793873))    #- Vcmax is 12.7% lower in narrow ranged species

plot(effect("Treatment:Location",fm.vcmax),multiline =T) #-0.0002 warming reduced Vcmax more in S than N
(exp(4.613478)-exp(4.441802))/(exp(4.441802))    # 15.9% reduction in Vcmax in S
(exp(5.187537)-exp(5.166953))/(exp(5.166953))    # 2.1 % reduction in Vcmax in N

plot(effect("Treatment:Range:Location",fm.vcmax),multiline =T) 
(exp(5.188676)-exp(5.218933))/(exp(5.218933))    # 2.9  % reduction in Vcmax in Nw
(exp(5.126655)-exp(5.129296))/(exp(5.129296))    # 0.2  % reduction in Vcmax in Nn
(exp(4.482517)-exp(4.685262))/(exp(4.685262))    # 18   %reduction in Vcmax in Sw
(exp(4.366272)-exp(4.480313))/(exp(4.480313))    # 10.8 % reduction in Vcmax in Sn

acifits2<-droplevels(subset(acifits, Species!="SMIT"))#only an effect of treatment (Warmed<Home) and location (N>S) for V and J

#- Vcmax per unit mass
fm.vcmaxm <- lme(log(Vcmaxm)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.vcmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmaxm,log(Vcmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmaxm,log(Vcmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmaxm$residuals[,1])
anova(fm.vcmaxm)                                  

plot(effect("Treatment", fm.vcmaxm))              #0.1500 warming did not reduce Vcmax on a mass basis
plot(effect("Location", fm.vcmaxm))               #<.0001 N still had higher Vcmax per g than S
plot(effect("Range", fm.vcmaxm))                  #0.0799 Wide had slightly higher Vcmax per g than narrow

plot(effect("Treatment:Location",fm.vcmaxm),multiline =T) #0.5164 warming did not reduce Vcmax per g more in S than N
plot(effect("Treatment:Range",fm.vcmaxm),multiline =T)    #0.0105 but did reduce Vcmax more in wide than narrow (Vcmax in Narrows was increased)

(exp(0.8487089)-exp(0.9395010))/(exp(0.9395010))    # 8.7% reduction in Vcmax in wides
(exp(0.8706312)-exp(0.8126672))/(exp(0.8126672))    # 6 % increase in Vcmax in narrows

plot(effect("Treatment:Range:Location",fm.vcmaxm),multiline =T) 
(exp(1.189035)-exp(1.291463))/(exp(1.291463))    # 9.7  % reduction in Vcmax in Nw
(exp(1.181733)-exp(1.192434))/(exp(1.192434))    # 1  % reduction in Vcmax in Nn
(exp(0.4795326)-exp(0.5702245))/(exp(0.5702245))    # 8.7 % reduction in Vcmax in Sw
(exp(0.4918060)-exp(0.3843307))/(exp(0.3843307))    # 11.3 % increase in Vcmax in Sn

#- Jmax
fm.jmax <- lme(log(Jmax)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifit)
plot(fm.jmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmax,log(Jmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmax,log(Jmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmax$residuals[,1])
anova(fm.jmax)    

plot(effect("Treatment", fm.jmax))                     #<.0001 warming reduced Jmax
plot(effect("Treatment:Location",fm.jmax),multiline=T) #0.0018 warming reduced Jmax more in S than N

(exp(5.010598)-exp(5.104039))/(exp(5.104039))    #  8.9% reduction in Jmax in N
(exp(4.922801)-exp(5.153175))/(exp(5.153175))    # 20.6% reduction in Jmax in S

plot(effect("Treatment:Range:Location",fm.jmax),multiline =T) 
(exp(5.045225)-exp(5.144893))/(exp(5.144893))    # 9.5  % reduction in Jmax in Nw
(exp(4.946363)-exp(5.028251))/(exp(5.028251))    # 7.9  % reduction in Jmax in Nn
(exp(4.936188)-exp(5.190911))/(exp(5.190911))    # 22.4 % reduction in Jmax in Sw
(exp(4.897968)-exp(5.083172))/(exp(5.083172))    # 16.9 % reduction in Jmax in Sn

#- Jmax per unit mass
fm.jmaxm <- lme(sqrt(Jmaxm)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits2)
plot(fm.jmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmaxm$residuals[,1])
anova(fm.jmaxm)                                  

plot(effect("Treatment", fm.jmaxm))                        #0.0002 warming reduced Jmax per g
plot(effect("Location", fm.jmaxm))                         #0.0774 higher Jmax per g in N
plot(effect("Treatment:Location",fm.jmaxm),multiline=T)    #0.2217 warming did not reduce Jmax per g more in S than N
plot(effect("Treatment:Range",fm.jmaxm),multiline=T)       #0.0273 warming did however reduce Jmax per g more in wide than narrow
(exp(1.785975)-exp(1.652841))/(exp(1.785975))    # 12.4% reduction in Jmax per g in wide
(exp(1.689670)-exp(1.667141))/(exp(1.689670))    #  2.2% reduction in Jmax per g in narrow

plot(effect("Treatment:Range:Location",fm.jmaxm),multiline =T) 
((1.696783^2)-(1.846167^2))/((1.846167^2))    # 15.5  % reduction in Jmax in Nw
((1.655308^2)-(1.730862^2))/((1.730862^2))    # 8.5  % reduction in Jmax in Nn
((1.605078^2)-(1.720548^2))/((1.720548^2))    # 13 % reduction in Jmax in Sw
((1.680002^2)-(1.644897^2))/((1.644897^2))    # 4.3 % increase in Jmax in Sn


#- Jmax/Vcmax
fm.jtov <- lme(log(JtoV)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifit)
plot(fm.jtov,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtov,log(JtoV)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtov,log(JtoV)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtov, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtov$residuals[,1])
anova(fm.jtov)    

plot(effect("Treatment",fm.jtov))                     #- <.0001 warming reduced Jmax/Vcmax overall
(exp(0.1478458)-exp(0.2138038))/(exp(0.2138038))      # 6.4% reduction in Jmax/Vcmax with warming

plot(effect("Location",fm.jtov))                      #- <.0001 Jmax/Vcmax much higher in S than north. No range interactions
(exp(0.5101642)-exp(-0.1200575))/(exp(-0.1200575))     # 87.8% lower Jmax/Vcmax in N compared to S

#Asat
fm.Asat <- lme((Photo^3)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asat)
plot(fm.Asat,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asat,(Photo^3)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asat,(Photo^3)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asat, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asat$residuals[,1])
anova(fm.Asat)    

plot(effect("Range",fm.Asat))                    #-0.0475 Asat higher in wide
((24067.19^(1/3))-(17730.40^(1/3)))/(17730.40^(1/3) )     # 10.7% higher Asat in wide
plot(effect("Treatment",fm.Asat))                #-<.0001 Asat reduced by warming
((19919.04^(1/3))-(23818.58^(1/3)))/(23818.58^(1/3) )     # 5.8% reduction in Asat with warming
plot(effect("Treatment*Range*Location",fm.Asat),multiline=T) 
((21593.35^(1/3))-(26792.92^(1/3)))/(26792.92^(1/3) ) #-6.9%
((15898.25^(1/3))-(20228.68^(1/3)))/(20228.68^(1/3) ) #-7.7%
((22778.92^(1/3))-(25158.33^(1/3)))/(25158.33^(1/3) ) #-3.2%
((15615.14^(1/3))-(19193.06^(1/3)))/(19193.06^(1/3) ) #-6.6%

Asatm2<- subset(Asatm, Species !="SMIT")#treatment:Location effect still there
#Asat on mass basis
fm.Asatm <- lme((Photom)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asatm)
plot(fm.Asatm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asatm,(Photom)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asatm,(Photom)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asatm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asatm$residuals[,1])
anova(fm.Asatm)    

plot(effect("Treatment",fm.Asatm))                #-0.6956 Asat per g was not reduced by warming
plot(effect("Location",fm.Asatm))                #-0.0274 Asat per g was however higher in N than S
plot(effect("Range",fm.Asatm))                    #-0.5461 Asat per g was not higher in wide
plot(effect("Treatment:Location",fm.Asatm),multiline=T)#-0.0001 Asat per g was increased by warming in S, decreased in N
((0.4988395)-(0.5521695))/(0.5521695 )     # 9.7 % reduction in Asat with warming in N
((0.5177445)-(0.4699496))/(0.4699496 )     # 10.2 % increase in Asat with warming in S
plot(effect("Treatment:Range",fm.Asatm),multiline=T)         # 0.0821 Asat decreased with warming in wide, increased in narrow
((0.5096694)-(0.5291048))/(0.5291048 )      # 3.7 % reduction in Asat with warming in wide
((0.5044685)-(0.4827984))/(0.4827984 )      # 4.5 % increase in Asat with warming in narrow 
plot(effect("Treatment*Range*Location",fm.Asatm),multiline=T)         
((0.5071643)-(0.5749812))/(0.5749812 )      # -11.8 % NW
((0.4832917)-(0.5095652))/(0.5095652 )      # -5.1% NN
((0.5124169)-(0.4787887))/(0.4787887 )      # +7  % SW
((0.5276947)-(0.4534412))/(0.4534412 )      # +16 % SN


###########################################################
Rdark2<- subset(Rdark, Species != "SMIT")#Only a location effect (N>S) and treatment:range effect (wide decreased > narrow)
#Respiration
fm.Rmass <- lme(log(Rmass)~Treatment*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdark)
plot(fm.Rmass,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rmass,log(Rmass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rmass,log(Rmass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rmass, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rmass$residuals[,1])
anova(fm.Rmass)
plot(effect("Treatment",fm.Rmass))                                 #-<.0001 Rmass was reduced by warming
plot(effect("Location",fm.Rmass))                                  #-0.0070 Rmass was higher in N than S
plot(effect("Treatment:Range",fm.Rmass), multiline = T)            #-0.0039 Rmass was reduced more in wide than narrow
(exp(2.078187)-exp(2.374139))/(exp(2.374139))      # 25% reduction in Rmass with warming of wide species
(exp(2.196815)-exp(2.339230))/(exp(2.339230))      # 13.2% reduction in Rmass with warming of narrow species

plot(effect("Treatment:Range:Location",fm.Rmass), multiline = T)   #-0.0236 Narrow from south don't acclimate, the rest do

(exp(2.426073)-exp(2.174630))/(exp(2.426073))      # 22% reduction in Rmass with warming of wide species in the north
(exp(2.361253)-exp(2.144013))/(exp(2.361253))      # 19.5% reduction in Rmass with warming of narrow species in the north
(exp(2.318532)-exp(1.974925))/(exp(2.318532))      # 29% reduction in Rmass with warming of wide species in the south
(exp(2.315649)-exp(2.253351))/(exp(2.315649))      # 6% reduction in Rmass with warming of narrow species in the south


#Respiration Temperature sensitivity curves: See GLAHD-RvsT.R