#######################
# GLAHD Growth Analysis
######################


#output from GAM model:"Code" "Time" "dydt" "predMass" "Species" "Treatment" "Location" "Taxa" "Range" "AGR" 
source("R/GLAHD_gamfits.R") 
library(effects)

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

data <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
data$Date <- as.Date(data$Date,format="%d/%m/%Y")
data$Totmass <- base::rowSums(data[,11:13]) #total mass is the sum of leaf, stem, and root mass
data$LAR <- with(data,Leafarea/Totmass)
data$SLA <- with(data,Leafarea/Leafmass)
data$d2h <- with(data,(Diameter/10)^2*(Height)) #calculate d2h in cm3
data$Treat <- as.factor(ifelse(data$Pot < 20, "Home",
                              ifelse(data$Pot>=40,"Pre","Warmed")))
data$Taxa <- factor(data$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
data[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
data$Range<- as.factor(ifelse(data$Species == "TER"|data$Species == "CAM", "wide","narrow"))
data2<-droplevels(subset(data,Treat!="Pre"))

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
#rate has data only for the days on which we measured Height and Diameter: Day 1 11 20 25 32 39 46 53 60 based on allometries
#data2 has harvested, actually measured data

###########################################################

#Did warming affect Diameter and Height?
#Generally Tropical taxa responded very modestly to warming with slight enhancement of height and diameter early and slight reductions towards the end.

T60<-subset(rate, Time==53) #change time to check differences at other time points

fm1size <- lme(sqrt(d2h)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
plot(fm1size,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1size,sqrt(d2h)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1size,sqrt(d2h)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1size, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1size$residuals[,1])
anova(fm1size)    
summary(fm1size)

plot(effect("Treatment:Location:Range",fm1size)) #Size increased more by warming in the S. P= <0.0001

Dw<-subset(summaryBy(Diameter~Treatment+Location+Range, data=T60, keep.names = T, FUN = mean),Treatment=="Warmed")
Dh<-subset(summaryBy(Diameter~Treatment+Location+Range, data=T60, keep.names = T, FUN = mean),Treatment=="Home")
Dw$warming<- Dw$Diameter/Dh$Diameter

fm1.mass <- lme(Height~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)
anova(fm1.mass)  
plot(effect("Treatment:Range",fm1.mass)) #Height is increased more by warming in the S than in the N. P= <0.0001

Hw<-subset(summaryBy(Height~Treatment+Location+Range, data=T60, keep.names = T, FUN = mean),Treatment=="Warmed")
Hh<-subset(summaryBy(Height~Treatment+Location+Range, data=T60, keep.names = T, FUN = mean),Treatment=="Home")
Hw$warming<- Hw$Height/Hh$Height


###########################################################

#Did warming affect allocation? (Totmass:Treatment,Totmass:Treatment:Range,Totmass:Treatment:Location, Totmass:Treatment:Range:Location)

fm1LM <- lme((Leafmass)~(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
plot(fm1LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LM,(Leafmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LM,(Leafmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LM$residuals[,1])
anova(fm1LM)    
summary(fm1LM) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1LM)) #P=0.1106
plot(effect("Treatment:Location:Range",fm1LM), multiline=TRUE) #P=0.1314
plot(effect("(Totmass):Treatment:Location",fm1LM), multiline=TRUE) #P=0.0412
plot(effect("(Totmass):Location ",fm1LM), multiline=TRUE) #P=0.0001

plot(effect("Treatment:Location",fm1LM), multiline=TRUE) #P=0.337





#-- REPEAT FOR STEM MASS
fm1SM <- lme((Stemmass)~(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2,
             weights=varFunc(~Totmass))#, method="ML")
plot(fm1SM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1SM,(Stemmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1SM,(Stemmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1SM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1SM$residuals[,1])
anova(fm1SM)    
summary(fm1SM) #The slope of LM~TM varies with treatment (P=0.0566) where warmed taxa have lower allocation to leaves
#the lm~TM:Treatment effect varies among the two locations (P=0.1576)

plot(allEffects(fm1SM)) 
plot(effect("Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location",fm1SM),multiline=T)
plot(effect("Treatment:Location:Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)


plot(effect("(Totmass):Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("(Totmass):Location",fm1SM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("(Totmass):Treatment:Range",fm1SM), multiline=TRUE) #- compares slopes (overlayed)


#-- REPEAT FOR ROOT MASS
fm1RM <- lme(sqrt(Rootmass)~(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2,
             weights=varFunc(~Totmass))#, method="ML")
plot(fm1RM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1RM,sqrt(Rootmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1RM,sqrt(Rootmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1RM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1RM$residuals[,1])
anova(fm1RM)    
summary(fm1RM) 

plot(allEffects(fm1RM)) 
plot(effect("Totmass:Treatment",fm1RM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location",fm1RM), multiline=TRUE)
plot(effect("Treatment:Location",fm1RM), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Range",fm1RM), multiline=TRUE) #- compares slopes (overlaye  d)
plot(effect("Totmass:Treatment:Location",fm1RM), multiline=TRUE) #- compares slopes (overlayed)




#-- REPEAT FOR LEAF AREA
fm1LA <- lme((Leafarea)^(2/3)~(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
plot(fm1LA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LA,(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LA,(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LA$residuals[,1])
anova(fm1LA)    
summary(fm1LA) 

plot(allEffects(fm1LA)) 
plot(effect("Totmass",fm1LA), multiline=TRUE)
plot(effect("Totmass:Range",fm1LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Location",fm1LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Treatment:Location",fm1LA), multiline=TRUE) #- compares slopes (overlayed)
plot(effect("Totmass:Treatment:Location",fm1LA), multiline=TRUE) #- compares slopes (overlayed)

#figure of RGR over time and mass
g.trt <- summaryBy(dydt+predMass~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))

windows(10,10);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,5,2,2),cex.axis=1)

palettS <- c("black","purple","blue","cyan"); palettN <- c("black","brown","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
sr.trt<- droplevels(subset(r.trt, Location == "S"));nr.trt<- droplevels(subset(r.trt, Location == "N"))

plotBy(dydt.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0.025,0.15),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=10), labels = F, tck = 0.01)

plotBy(dydt.mean~predMass.mean|combotrt,data=ng.trt,col=palettN,
       legendwhere="topright", pch=3, xaxt='n', ylab="", yaxt='n',type="l",ylim=c(0.025,0.15),xlim=c(0.3,80),log="x",
       panel.first=adderrorbars(x=ng.trt$predMass.mean,y=ng.trt$dydt.mean,
                                SE=ng.trt$dydt.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=80,by=1), labels = F, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.2,by=0.02), labels = F, tck = 0.01)

plotBy(dydt.mean~Time|combotrt,data=sg.trt,col=palettS,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0.025,0.15),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
axis(side = 1, at = seq(from=0,to=80,by=10), labels = T, tck = 0.01)
mtext(side=1, text="Time (Days)", line=3)

plotBy(dydt.mean~predMass.mean|combotrt,data=sg.trt,col=palettS,
       legendwhere="topright", pch=3, xaxt='n', ylab="", yaxt='n',type="l",ylim=c(0.025,0.15),xlim=c(0.3,80),log="x",
       panel.first=adderrorbars(x=sg.trt$predMass.mean,y=sg.trt$dydt.mean,
                                SE=sg.trt$dydt.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
mtext(side=2, text=expression(RGR~(g~g^-1~day^-1)), line=3, outer=T)
axis(side = 1, at = seq(from=0,to=80,by=1), labels = T, tck = 0.01)
axis(side = 2, at = seq(from=0,to=0.2,by=0.02), labels = F, tck = 0.01)
mtext(side=1, text="logMass (g)", line=3)
#Did  warming increase biomass?


#Figure of biomass over time
g.trt <- summaryBy(dydt+predMass+AGR~Time+Treatment+Location+Range,data=gamfits2,FUN=c(mean,standard.error))
g.trt$combotrt <- as.factor(paste(g.trt$Location,g.trt$Range,g.trt$Treatment,sep="_"))
r.trt <- summaryBy(leafMass+stemMass+rootMass~Time+Treatment+Location+Range,data=rate,FUN=c(mean,standard.error))
r.trt$combotrt <- as.factor(paste(r.trt$Location,r.trt$Range,r.trt$Treatment,sep="_"))

windows(50,40);par(mfrow=c(4,2), mar=c(0.5,6,0.5,1), oma=c(5,1,1,1))

palettS <- c("black","purple","blue","cyan"); palettN <- c("black","brown","orange","red")
sg.trt<- droplevels(subset(g.trt, Location == "S"));ng.trt<- droplevels(subset(g.trt, Location == "N"))
sr.trt<- droplevels(subset(r.trt, Location == "S"));nr.trt<- droplevels(subset(r.trt, Location == "N"))

plotBy(predMass.mean~Time|combotrt,data=sg.trt,col=palettS,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=sg.trt$Time,y=sg.trt$predMass.mean,
                                SE=sg.trt$predMass.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)

legend("topleft",)

plotBy(predMass.mean~Time|combotrt,data=ng.trt,col=palettN,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,70),
       panel.first=adderrorbars(x=ng.trt$Time,y=ng.trt$predMass.mean,
                                SE=ng.trt$predMass.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Total biomass (g)",  line=3)
#-------------------------------------------------------------------------------------
plotBy(leafMass.mean~Time|combotrt,data=sr.trt,col=palettS,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=sr.trt$Time,y=sr.trt$leafMass.mean,
                                SE=sr.trt$leafMass.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Leaf  (g)",  line=3)

plotBy(leafMass.mean~Time|combotrt,data=nr.trt,col=palettN,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=nr.trt$Time,y=nr.trt$leafMass.mean,
                                SE=nr.trt$leafMass.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Leaf  (g)",  line=3)
#-------------------------------------------------------------------------------------
plotBy(stemMass.mean~Time|combotrt,data=sr.trt,col=palettS,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=sr.trt$Time,y=sr.trt$stemMass.mean,
                                SE=sr.trt$stemMass.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Stem  (g)",  line=3)

plotBy(stemMass.mean~Time|combotrt,data=nr.trt,col=palettN,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,35),
       panel.first=adderrorbars(x=nr.trt$Time,y=nr.trt$stemMass.mean,
                                SE=nr.trt$stemMass.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = FALSE, tck = 0.01)
mtext(side=2, text="Stem  (g)",  line=3)
#-------------------------------------------------------------------------------------
plotBy(rootMass.mean~Time|combotrt,data=sr.trt,col=palettS,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,12),
       panel.first=adderrorbars(x=sr.trt$Time,y=sr.trt$rootMass.mean,
                                SE=sr.trt$rootMass.standard.error,direction="updown",
                                col=c("black","blue","purple","cyan"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, text="Root  (g)",  line=3)
mtext(text="Time (days)", side=1, line=3)

plotBy(rootMass.mean~Time|combotrt,data=nr.trt,col=palettN,
       legend=F, pch=3, xaxt='n', ylab="", type="l",ylim=c(0,12),
       panel.first=adderrorbars(x=nr.trt$Time,y=nr.trt$rootMass.mean,
                                SE=nr.trt$rootMass.standard.error,direction="updown",
                                col=c("black","orange","brown","red"),0))
axis(side = 1, at = seq(from=0,to=60,by=5), labels = T, tck = 0.01)
mtext(side=2, text="Root  (g)",  line=3)
mtext(text="Time (days)", side=1, line=3)

#add nested components
gamfits2$Location <- factor(gamfits2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2$Sp_RS_EN <- as.factor(with(gamfits2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2$Prov_Sp_EN <- as.factor(with(gamfits2,paste(Taxa,Species)))
gamfits2$Sp_Loc_EN <- as.factor(with(gamfits2,paste(Species,Location)))

rate$Location <- factor(rate$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
rate$Sp_RS_EN <- as.factor(with(rate,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rate$Prov_Sp_EN <- as.factor(with(rate,paste(Taxa,Species)))
rate$Sp_Loc_EN <- as.factor(with(rate,paste(Species,Location)))


#Extract data from specific time points and look for statistically significant differences in biomass
dat5<- subset(gamfits2,Time==5)

fm1.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm1.mass,sqrt(predMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.1.mass <- lme(sqrt(predMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==10)

fm2.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm2.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2.mass,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm2.mass,sqrt(predMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm2.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2.mass$residuals[,1])
anova(fm2.mass)    
fm2.1.mass <- lme(sqrt(predMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm2.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==15)

fm3.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm3.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3.mass$residuals[,1])
anova(fm3.mass)    
fm3.1.mass <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm3.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==20)

fm4.mass <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm4.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm4.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm4.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm4.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm4.mass$residuals[,1])
anova(fm4.mass)    
fm4.1.mass <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm4.1.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==25)

fm5.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm5.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm5.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm5.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm5.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm5.mass$residuals[,1])
anova(fm5.mass) 
fm5.1.mass <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm5.1.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==30)

fm6.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm6.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm6.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm6.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm6.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm6.mass$residuals[,1])
anova(fm6.mass)  
plot(allEffects(fm6.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==35)

fm7.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm7.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm7.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm7.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm7.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm7.mass$residuals[,1])
anova(fm7.mass)    
plot(allEffects(fm7.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==40)

fm8.mass <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm8.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm8.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm8.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm8.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm8.mass$residuals[,1])
anova(fm8.mass)    
plot(allEffects(fm8.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==45)

fm9.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm9.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm9.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm9.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm9.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm9.mass$residuals[,1])
anova(fm9.mass)    
fm9.1.mass <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm9.1.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==50)

fm10.mass <- lme(predMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm10.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm10.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm10.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm10.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm10.mass$residuals[,1])
anova(fm10.mass)    
fm10.1.mass <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm10.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==55)

fm11.mass <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm11.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm11.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm11.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm11.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm11.mass$residuals[,1])
anova(fm11.mass)    
fm11.1mass <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm11.1mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(gamfits2,Time==60)

fm12.mass <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm12.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm12.mass,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm12.mass,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm12.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm12.mass$residuals[,1])
anova(fm12.mass)  
fm12.1.mass <- lme((predMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm12.1.mass))


#Same for Leaf mass
dat5<- subset(rate,Time==1)

fm1.mass <- lme(sqrt(leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(leafMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm1.mass,sqrt(leafMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.1.mass <- lme(sqrt(leafMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==11)

fm2.mass <- lme(log(leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm2.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2.mass,(leafMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm2.mass,(leafMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm2.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2.mass$residuals[,1])
anova(fm2.mass)    
fm2.1.mass <- lme((leafMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm2.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==20)

fm3.mass <- lme((leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm3.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3.mass$residuals[,1])
anova(fm3.mass)    
fm3.1.mass <- lme((leafMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm3.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==25)

fm4.mass <- lme(log(leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm4.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm4.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm4.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm4.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm4.mass$residuals[,1])
anova(fm4.mass)    
fm4.1.mass <- lme((leafMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm4.1.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==32)

fm5.mass <- lme(sqrt(leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm5.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm5.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm5.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm5.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm5.mass$residuals[,1])
anova(fm5.mass) 
plot(allEffects(fm5.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==39)

fm6.mass <- lme((leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm6.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm6.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm6.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm6.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm6.mass$residuals[,1])
anova(fm6.mass)  
fm6.1.mass <- lme((leafMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location+Location:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm6.1.mass))
 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==46)

fm7.mass <- lme((leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm7.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm7.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm7.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm7.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm7.mass$residuals[,1])
anova(fm7.mass)    
fm7.1.mass <- lme((leafMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location+Location:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm7.1.mass))


#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==53)

fm8.mass <- lme((leafMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm8.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm8.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm8.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm8.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm8.mass$residuals[,1])
anova(fm8.mass)  

fm8.1.mass <- lme((leafMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm8.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==60)

fm9.mass <- lme(leafMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm9.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm9.mass,leafMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm9.mass,leafMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm9.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm9.mass$residuals[,1])
anova(fm9.mass)    
fm9.1.mass <- lme((leafMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm9.1.mass))


#Same for stem mass
dat5<- subset(rate,Time==1)

fm1.mass <- lme(sqrt(stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(stemMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm1.mass,sqrt(stemMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.1.mass <- lme(sqrt(stemMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==11)

fm2.mass <- lme(sqrt(stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm2.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2.mass,(stemMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm2.mass,(stemMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm2.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2.mass$residuals[,1])
anova(fm2.mass)    
plot(allEffects(fm2.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==20)

fm3.mass <- lme((stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm3.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3.mass$residuals[,1])
anova(fm3.mass)    
fm3.1.mass <- lme((stemMass)~Treatment+Range+Location+Treatment:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm3.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==25)

fm4.mass <- lme(log(stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm4.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm4.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm4.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm4.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm4.mass$residuals[,1])
anova(fm4.mass)    
plot(allEffects(fm4.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==32)

fm5.mass <- lme((stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm5.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm5.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm5.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm5.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm5.mass$residuals[,1])
anova(fm5.mass) 
plot(allEffects(fm5.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==39)

fm6.mass <- lme((stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm6.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm6.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm6.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm6.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm6.mass$residuals[,1])
anova(fm6.mass)  
plot(allEffects(fm6.mass))


#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==46)

fm7.mass <- lme((stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm7.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm7.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm7.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm7.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm7.mass$residuals[,1])
anova(fm7.mass)    
plot(allEffects(fm7.mass))


#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==53)

fm8.mass <- lme((stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm8.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm8.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm8.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm8.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm8.mass$residuals[,1])
anova(fm8.mass)  

plot(allEffects(fm8.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==60)

fm9.mass <- lme(sqrt(stemMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm9.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm9.mass,stemMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm9.mass,stemMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm9.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm9.mass$residuals[,1])
anova(fm9.mass)    
fm9.1.mass <- lme(sqrt(stemMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm9.1.mass))

#Same for root mass
dat5<- subset(rate,Time==1)

fm1.mass <- lme(sqrt(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm1.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1.mass,sqrt(rootMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm1.mass,sqrt(rootMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm1.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1.mass$residuals[,1])
anova(fm1.mass)    
fm1.1.mass <- lme(sqrt(rootMass)~Treatment+Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm1.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==11)

fm2.mass <- lme(sqrt(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm2.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2.mass,(rootMass)~fitted(.)|Species,abline=c(0,1))       #predicted vs. fitted for each species
plot(fm2.mass,(rootMass)~fitted(.),abline=c(0,1))               #overall predicted vs. fitted
qqnorm(fm2.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2.mass$residuals[,1])
anova(fm2.mass)    
fm2.1.mass <- lme(sqrt(rootMass)~Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm2.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==20)

fm3.mass <- lme(log(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm3.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3.mass$residuals[,1])
anova(fm3.mass)    
fm3.1.mass <- lme((rootMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm3.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==25)

fm4.mass <- lme(log(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm4.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm4.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm4.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm4.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm4.mass$residuals[,1])
anova(fm4.mass)    
fm4.1.mass <- lme((rootMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm4.1.mass)) 


#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==32)

fm5.mass <- lme(sqrt(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm5.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm5.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm5.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm5.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm5.mass$residuals[,1])
anova(fm5.mass) 
plot(allEffects(fm5.mass))

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==39)

fm6.mass <- lme(log(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm6.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm6.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm6.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm6.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm6.mass$residuals[,1])
anova(fm6.mass)  
fm6.1.mass <- lme((rootMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm6.1.mass)) 



#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==46)

fm7.mass <- lme((rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm7.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm7.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm7.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm7.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm7.mass$residuals[,1])
anova(fm7.mass)    
plot(allEffects(fm7.mass))


#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==53)

fm8.mass <- lme((rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm8.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm8.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm8.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm8.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm8.mass$residuals[,1])
anova(fm8.mass)  
fm8.1.mass <- lme((rootMass)~Treatment+Range+Location+Treatment:Range+Treatment:Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm8.1.mass)) 

#-------------------------------------------------------------------------------------
dat5<- subset(rate,Time==60)

fm9.mass <- lme(sqrt(rootMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(fm9.mass,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm9.mass,rootMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm9.mass,rootMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm9.mass, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm9.mass$residuals[,1])
anova(fm9.mass)    
fm9.1.mass <- lme(sqrt(rootMass)~Treatment*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat5)
plot(allEffects(fm9.1.mass))


