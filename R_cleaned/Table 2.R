#Daily stats

source("R_cleaned/1. GLAHD_LoadLibraries.R")
source("R_cleaned/2. downloadData.R")
source("R_cleaned/3. Create_datasets.R")

# BIOMASS

T1<-subset(gamfits2, Time==1) 
fm1m1 <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T1, method="ML")
plot(fm1m1,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m1,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m1,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m1, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m1$residuals[,1])
anova(fm1m1)    
plot(effect("Location",fm1m1))                     #0.07 biomass higher in N

fm1m1.2 <- lme(sqrt(predMass)~Location,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T1, method="ML")
plot(fm1m1.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m1.2,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m1.2,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m1.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m1.2$residuals[,1])
anova(fm1m1.2)    
anova(fm1m1.2,fm1m1)  #0.08 biomass higher in N


T10<-subset(gamfits2, Time==10) 
fm1m10 <- lme(sqrt(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T10, method="ML")
plot(fm1m10,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m10,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m10,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m10, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m10$residuals[,1])
anova(fm1m10)    
plot(effect("Treatment",fm1m10))                    #<.0001 biomass increased with warming
plot(effect("Location",fm1m10))                     #0.004 biomass higher in N

fm1m10.2 <- lme(sqrt(predMass)~Treatment+Location,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T10, method="ML")
plot(fm1m10.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m10.2,sqrt(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m10.2,sqrt(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m10.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m10.2$residuals[,1])
anova(fm1m10.2,fm1m10)    



T20<-subset(gamfits2, Time==20) 
fm1m20 <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T20, method="ML")
plot(fm1m20,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m20,(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m20,(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m20, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m20$residuals[,1])
anova(fm1m20)    
plot(effect("Treatment",fm1m20))                    #<.0001 biomass increased with warming
plot(effect("Location",fm1m20))                     #0.0003 biomass higher in N
plot(effect("Treatment:Location",fm1m20), multiline = T) #0.015 warming increased biomass more in S than N
((effect("Treatment:Location",fm1m20)[[5]][4,1]-effect("Treatment:Location",fm1m20)[[5]][3,1])/
  effect("Treatment:Location",fm1m20)[[5]][3,1])*100 #+10 % N
((effect("Treatment:Location",fm1m20)[[5]][2,1]-effect("Treatment:Location",fm1m20)[[5]][1,1])/
  effect("Treatment:Location",fm1m20)[[5]][1,1])*100 #+51.2% S
plot(effect("Treatment:Range",fm1m20), multiline = T) #0.033 warming increased biomass more in wide than narrow
((effect("Treatment:Range",fm1m20)[[5]][4,1]-effect("Treatment:Range",fm1m20)[[5]][3,1])/
  effect("Treatment:Range",fm1m20)[[5]][3,1])*100 #+28.6 % wide
((effect("Treatment:Range",fm1m20)[[5]][2,1]-effect("Treatment:Range",fm1m20)[[5]][1,1])/
  effect("Treatment:Range",fm1m20)[[5]][1,1])*100 #+14.2.2% narrow

fm1m20.2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T20, method="ML")
plot(fm1m20.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m20.2,(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m20.2,(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m20.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m20.2$residuals[,1])
anova(fm1m20.2,fm1m20)    



T30<-subset(gamfits2, Time==30) 
fm1m30 <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T30, method="ML")
plot(fm1m30,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m30,log(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m30,log(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m30, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m30$residuals[,1])
anova(fm1m30)    
plot(effect("Treatment",fm1m30))                    #<.0001 biomass increased with warming
plot(effect("Location",fm1m30))                     #<.0001 biomass higher in N
plot(effect("Treatment:Location",fm1m30), multiline = T) #<.0001 warming increased biomass more in S than N
((exp(effect("Treatment:Location",fm1m30)[[5]][4,1])-exp(effect("Treatment:Location",fm1m30)[[5]][3,1]))/
  exp(effect("Treatment:Location",fm1m30)[[5]][3,1]))*100 #+4.06 % N
((exp(effect("Treatment:Location",fm1m30)[[5]][2,1])-exp(effect("Treatment:Location",fm1m30)[[5]][1,1]))/
  exp(effect("Treatment:Location",fm1m30)[[5]][1,1]))*100 #+73.3% S
plot(effect("Treatment:Location:Range",fm1m30), multiline = T) #0.08 warming increased biomass more in S wide
((exp(effect("Treatment:Location:Range",fm1m30)[[5]][8,1])-exp(effect("Treatment:Location:Range",fm1m30)[[5]][7,1]))/
  exp(effect("Treatment:Location:Range",fm1m30)[[5]][7,1]))*100 #+6.8 % N wide
((exp(effect("Treatment:Location:Range",fm1m30)[[5]][4,1])-exp(effect("Treatment:Location:Range",fm1m30)[[5]][3,1]))/
  exp(effect("Treatment:Location:Range",fm1m30)[[5]][3,1]))*100 #-0.65% N narrow
((exp(effect("Treatment:Location:Range",fm1m30)[[5]][6,1])-exp(effect("Treatment:Location:Range",fm1m30)[[5]][5,1]))/
  exp(effect("Treatment:Location:Range",fm1m30)[[5]][5,1]))*100 #+95.8 % S wide
((exp(effect("Treatment:Location:Range",fm1m30)[[5]][2,1])-exp(effect("Treatment:Location:Range",fm1m30)[[5]][1,1]))/
  exp(effect("Treatment:Location:Range",fm1m30)[[5]][1,1]))*100 #+39.2% S narrow

fm1m30.2 <- lme(log(predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T30, method="ML")
plot(fm1m30.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1m30.2,log(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1m30.2,log(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1m30.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1m30.2$residuals[,1])
anova(fm1m30.2)    
((exp(effect("Treatment:Location",fm1m30.2)[[5]][4,1])-exp(effect("Treatment:Location",fm1m30)[[5]][3,1]))/
    exp(effect("Treatment:Location",fm1m30.2)[[5]][3,1]))*100 #+0 % N
((exp(effect("Treatment:Location",fm1m30.2)[[5]][2,1])-exp(effect("Treatment:Location",fm1m30)[[5]][1,1]))/
    exp(effect("Treatment:Location",fm1m30.2)[[5]][1,1]))*100 #+76% S
anova(fm1m30.2,fm1m30)
AIC(fm1m30.2,fm1m30)#Simpler model is a bit better P=0.053

# T40<-subset(gamfits2, Time==40) 
# fm1m40 <- lme(log(predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
# plot(fm1m40,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1m40,log(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1m40,log(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1m40, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1m40$residuals[,1])
# anova(fm1m40)    
# plot(effect("Treatment",fm1m40))                    #<.0001 biomass increased with warming
# plot(effect("Location",fm1m40))                     #<.0001 biomass higher in N
# plot(effect("Treatment:Location",fm1m40), multiline = T) #<.0001 warming increased biomass more in S than N
# ((exp(effect("Treatment:Location",fm1m40)[[5]][4,1])-exp(effect("Treatment:Location",fm1m40)[[5]][3,1]))/
#   exp(effect("Treatment:Location",fm1m40)[[5]][3,1]))*100 #+1.9 % N
# ((exp(effect("Treatment:Location",fm1m40)[[5]][2,1])-exp(effect("Treatment:Location",fm1m40)[[5]][1,1]))/
#   exp(effect("Treatment:Location",fm1m40)[[5]][1,1]))*100 #+75.6% S
# plot(effect("Treatment:Location:Range",fm1m40), multiline = T) #0.08 warming increased biomass more in S wide
# ((exp(effect("Treatment:Location:Range",fm1m40)[[5]][8,1])-exp(effect("Treatment:Location:Range",fm1m40)[[5]][7,1]))/
#   exp(effect("Treatment:Location:Range",fm1m40)[[5]][7,1]))*100 #+3.6 % N wide
# ((exp(effect("Treatment:Location:Range",fm1m40)[[5]][4,1])-exp(effect("Treatment:Location:Range",fm1m40)[[5]][3,1]))/
#   exp(effect("Treatment:Location:Range",fm1m40)[[5]][3,1]))*100 #-1.13% N narrow
# ((exp(effect("Treatment:Location:Range",fm1m40)[[5]][6,1])-exp(effect("Treatment:Location:Range",fm1m40)[[5]][5,1]))/
#   exp(effect("Treatment:Location:Range",fm1m40)[[5]][5,1]))*100 #+96 % S wide
# ((exp(effect("Treatment:Location:Range",fm1m40)[[5]][2,1])-exp(effect("Treatment:Location:Range",fm1m40)[[5]][1,1]))/
#   exp(effect("Treatment:Location:Range",fm1m40)[[5]][1,1]))*100 #+44.17% S narrow
# 
# fm1m40.2 <- lme(log(predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range+Location:Range,
#                 random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
# plot(fm1m40.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1m40.2,log(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1m40.2,log(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1m40.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1m40.2$residuals[,1])
# anova(fm1m40.2)    
# anova(fm1m40.2,fm1m40) 
# plot(effect("Treatment:Range",fm1m40.2), multiline = T) #0.0071 warming increased biomass more in S than N
# ((exp(effect("Treatment:Range",fm1m40.2)[[5]][4,1])-exp(effect("Treatment:Range",fm1m40.2)[[5]][3,1]))/
#     exp(effect("Treatment:Range",fm1m40.2)[[5]][3,1]))*100 #+45 % Wide
# ((exp(effect("Treatment:Range",fm1m40.2)[[5]][2,1])-exp(effect("Treatment:Range",fm1m40.2)[[5]][1,1]))/
#     exp(effect("Treatment:Range",fm1m40.2)[[5]][1,1]))*100 #+20.1% Narrow
# 
# plot(effect("Treatment:Location",fm1m40.2), multiline = T) #<.0001 warming increased biomass more in S than N
# ((exp(effect("Treatment:Location",fm1m40.2)[[5]][4,1])-exp(effect("Treatment:Location",fm1m40.2)[[5]][3,1]))/
#     exp(effect("Treatment:Location",fm1m40.2)[[5]][3,1]))*100 #+1.5 % N
# ((exp(effect("Treatment:Location",fm1m40.2)[[5]][2,1])-exp(effect("Treatment:Location",fm1m40.2)[[5]][1,1]))/
#     exp(effect("Treatment:Location",fm1m40.2)[[5]][1,1]))*100 #+75.1% S
# 
# fm1m40.3 <- lme(log(predMass)~Treatment+Location+Range+Treatment:Location,
#                 random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
# anova(fm1m40.2,fm1m40.3) #significantly worse 0.0064
# 
# fm1m40.4 <- lme(log(predMass)~Treatment+Location+Range+Treatment:Range,
#                 random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
# anova(fm1m40.2,fm1m40.4) #significantly worse <.0001
# 
# fm1m40.5 <- lme(log(predMass)~Treatment+Location+Range,
#                 random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
# anova(fm1m40.5,fm1m40) #significantly worse <.0001
# anova(fm1m40.5,fm1m40.2)#significantly worse <.0001
# anova(fm1m40.5,fm1m40.3)#significantly worse <.0001
# anova(fm1m40.5,fm1m40.4)#significantly worse 0.04

# T50<-subset(gamfits2, Time==50) 
# fm1m50 <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T50, method="ML")
# plot(fm1m50,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1m50,(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1m50,(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1m50, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1m50$residuals[,1])
# anova(fm1m50)    
# plot(effect("Treatment",fm1m40))                    #<.0001 biomass increased with warming
# plot(effect("Location",fm1m40))                     #<.0001 biomass higher in N
# plot(effect("Treatment:Location",fm1m40), multiline = T) #<.0001 warming increased biomass more in S than N
# ((effect("Treatment:Location",fm1m40)[[5]][4,1]-effect("Treatment:Location",fm1m40)[[5]][3,1])/
#   effect("Treatment:Location",fm1m40)[[5]][3,1])*100 #+0.6 % N
# ((effect("Treatment:Location",fm1m40)[[5]][2,1]-effect("Treatment:Location",fm1m40)[[5]][1,1])/
#   effect("Treatment:Location",fm1m40)[[5]][1,1])*100 #+26.46% S
# plot(effect("Treatment:Range",fm1m40), multiline = T) #0.02 warming increased biomass more in wide than narrow
# ((effect("Treatment:Range",fm1m40)[[5]][4,1]-effect("Treatment:Range",fm1m40)[[5]][3,1])/
#   effect("Treatment:Range",fm1m40)[[5]][3,1])*100 #+14.1 % wide
# ((effect("Treatment:Range",fm1m40)[[5]][2,1]-effect("Treatment:Range",fm1m40)[[5]][1,1])/
#   effect("Treatment:Range",fm1m40)[[5]][1,1])*100 #+7.4% narow
# 
# fm1m50.2 <- lme((predMass)~Treatment+Location+Range+Treatment:Location+Treatment:Range,
#                 random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T50, method="ML")
# plot(fm1m50.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1m50.2,predMass~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1m50.2,predMass~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1m50.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1m50.2$residuals[,1])
# anova(fm1m50.2)    
# anova(fm1m50.2,fm1m50) 
# 
# 
# T60<-subset(gamfits2, Time==60) 
# 
# fm1m60 <- lme((predMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
# plot(fm1m60,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1m60,(predMass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1m60,(predMass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1m60, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1m60$residuals[,1])
# anova(fm1m60)    
# plot(effect("Treatment",fm1m60))                    #<.0001 biomass increased with warming
# plot(effect("Location",fm1m60))                     #<.0001 biomass higher in N
# plot(effect("Treatment:Location",fm1m60), multiline = T) #<.0001 warming increased biomass in S not N
# 
# ((effect("Treatment:Location",fm1m60)[[5]][4,1]-effect("Treatment:Location",fm1m60)[[5]][3,1])/
#   effect("Treatment:Location",fm1m60)[[5]][3,1])*100 #-3.4 % N
# ((effect("Treatment:Location",fm1m60)[[5]][2,1]-effect("Treatment:Location",fm1m60)[[5]][1,1])/
#   effect("Treatment:Location",fm1m60)[[5]][1,1])*100 #+48.7% S
# 


###-----------------------------
#RGR over time
T1<-subset(gamfits2, Time==1) 
fm1r1 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T1, method="ML")
plot(fm1r1,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r1,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r1,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r1, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r1$residuals[,1])
anova(fm1r1)    


T10<-subset(gamfits2, Time==10) 
fm1r10 <- lme(sqrt(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T10, method="ML")
plot(fm1r10,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r10,sqrt(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r10,sqrt(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r10, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r10$residuals[,1])
anova(fm1r10)    
plot(effect("Treatment",fm1r10))                    #<.0001 RGR increased with warming
plot(effect("Location",fm1r10))                     #0.0002 RGR higher in N
plot(effect("Treatment:Location",fm1r10), multiline = T) #0.0001 warming increased RGR more in S than N
(((effect("Treatment:Location",fm1r10)[[5]][4,1])^2-(effect("Treatment:Location",fm1r10)[[5]][3,1])^2)/
  (effect("Treatment:Location",fm1r10)[[5]][3,1])^2)*100 #+7 % N
(((effect("Treatment:Location",fm1r10)[[5]][2,1])^2-(effect("Treatment:Location",fm1r10)[[5]][1,1])^2)/
  (effect("Treatment:Location",fm1r10)[[5]][1,1])^2)*100 #+25.2% S

fm1r10.2 <- lme(sqrt(dydt)~Treatment+Location+Treatment:Location,
              random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T10)#, method="ML")
plot(fm1r10.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r10.2,sqrt(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r10.2,sqrt(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r10.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r10.2$residuals[,1])
anova(fm1r10.2)    
anova(fm1r10,fm1r10.2)

plot(effect("Treatment:Location",fm1r10.2), multiline = T) #0.0001 warming increased RGR more in S than N
(((effect("Treatment:Location",fm1r10.2)[[5]][4,1])^2-(effect("Treatment:Location",fm1r10.2)[[5]][3,1])^2)/
    (effect("Treatment:Location",fm1r10.2)[[5]][3,1])^2)*100 #+7 % N
(((effect("Treatment:Location",fm1r10.2)[[5]][2,1])^2-(effect("Treatment:Location",fm1r10.2)[[5]][1,1])^2)/
    (effect("Treatment:Location",fm1r10.2)[[5]][1,1])^2)*100 #+25.4% S

T20<-subset(gamfits2, Time==20) 
fm1r20 <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T20, method="ML")
plot(fm1r20,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r20,log(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r20,log(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r20, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r20$residuals[,1])
anova(fm1r20)    
plot(effect("Treatment",fm1r20))                    #<.0001 RGR increased with warming
plot(effect("Location",fm1r20))                     #0.03 RGR higher in N
plot(effect("Treatment:Location",fm1r20), multiline = T) #<.0001 warming increased RGR in S not N
((exp(effect("Treatment:Location",fm1r20)[[5]][4,1])-exp(effect("Treatment:Location",fm1r20)[[5]][3,1]))/
  exp(effect("Treatment:Location",fm1r20)[[5]][3,1]))*100 #-5.5 % N
((exp(effect("Treatment:Location",fm1r20)[[5]][2,1])-exp(effect("Treatment:Location",fm1r20)[[5]][1,1]))/
  exp(effect("Treatment:Location",fm1r20)[[5]][1,1]))*100 #+19.7% S
plot(effect("Treatment:Range",fm1r20), multiline = T) #0.01 warming increased RGR more in wide than narrow
((exp(effect("Treatment:Range",fm1r20)[[5]][4,1])-exp(effect("Treatment:Range",fm1r20)[[5]][3,1]))/
  exp(effect("Treatment:Range",fm1r20)[[5]][3,1]))*100 #+10.3 % wide
((exp(effect("Treatment:Range",fm1r20)[[5]][2,1])-exp(effect("Treatment:Range",fm1r20)[[5]][1,1]))/
  exp(effect("Treatment:Range",fm1r20)[[5]][1,1]))*100 #+1.6% narrow

fm1r20.2 <- lme(log(dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T20)#, method="ML")
plot(fm1r20.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r20.2,log(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r20.2,log(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r20.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r20.2$residuals[,1])
anova(fm1r20.2)    
anova(fm1r20.2,fm1r20)

plot(effect("Treatment:Location",fm1r20.2), multiline = T) #<.0001 warming increased RGR in S not N
((exp(effect("Treatment:Location",fm1r20.2)[[5]][4,1])-exp(effect("Treatment:Location",fm1r20.2)[[5]][3,1]))/
    exp(effect("Treatment:Location",fm1r20.2)[[5]][3,1]))*100 #-5.7 % N
((exp(effect("Treatment:Location",fm1r20.2)[[5]][2,1])-exp(effect("Treatment:Location",fm1r20.2)[[5]][1,1]))/
    exp(effect("Treatment:Location",fm1r20.2)[[5]][1,1]))*100 #+19.6% S
plot(effect("Treatment:Range",fm1r20.2), multiline = T) #0.015 warming increased RGR more in wide than narrow
((exp(effect("Treatment:Range",fm1r20.2)[[5]][4,1])-exp(effect("Treatment:Range",fm1r20.2)[[5]][3,1]))/
    exp(effect("Treatment:Range",fm1r20.2)[[5]][3,1]))*100 #+10.3 % wide
((exp(effect("Treatment:Range",fm1r20.2)[[5]][2,1])-exp(effect("Treatment:Range",fm1r20.2)[[5]][1,1]))/
    exp(effect("Treatment:Range",fm1r20.2)[[5]][1,1]))*100 #+1.6% narrow

T30<-subset(gamfits2, Time==30) 
fm1r30 <- lme(1/(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T30)#, method="ML")
plot(fm1r30,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r30,1/(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r30,1/(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r30, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r30$residuals[,1])
anova(fm1r30)    
plot(effect("Treatment:Location",fm1r30), multiline = T) #<.0001 warming increased RGR more in S than N
(((1/effect("Treatment:Location",fm1r30)[[5]][4,1])-(1/effect("Treatment:Location",fm1r30)[[5]][3,1]))/
  (1/effect("Treatment:Location",fm1r30)[[5]][3,1]))*100 #-6 % N
(((1/effect("Treatment:Location",fm1r30)[[5]][2,1])-(1/effect("Treatment:Location",fm1r30)[[5]][1,1]))/
  (1/effect("Treatment:Location",fm1r30)[[5]][1,1]))*100 #+8 % S

fm1r30.2 <- lme(1/(dydt)~Treatment+Location+Treatment:Location,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T30, method="ML")
plot(fm1r30.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r30.2,1/(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r30.2,1/(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r30.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r30.2$residuals[,1])
anova(fm1r30.2)  
anova(fm1r30.2,fm1r30)  


T40<-subset(gamfits2, Time==40) 
fm1r40 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
plot(fm1r40,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r40,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r40,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r40, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r40$residuals[,1])
anova(fm1r40)    
plot(effect("Location",fm1r40))                     #0.007 RGR lower in N
plot(effect("Treatment:Location",fm1r40), multiline = T) #0.08 warming increased RGR more in N than S
((effect("Treatment:Location",fm1r40)[[5]][4,1]-effect("Treatment:Location",fm1r40)[[5]][3,1])/
  effect("Treatment:Location",fm1r40)[[5]][3,1])*100 #+1.7 % N
((effect("Treatment:Location",fm1r40)[[5]][2,1]-effect("Treatment:Location",fm1r40)[[5]][1,1])/
  effect("Treatment:Location",fm1r40)[[5]][1,1])*100 #-5.5% S
plot(effect("Treatment:Range",fm1r40), multiline = T) #0.005 warming increased RGR more in wide than narrow
((exp(effect("Treatment:Range",fm1r40)[[5]][4,1])-exp(effect("Treatment:Range",fm1r40)[[5]][3,1]))/
  exp(effect("Treatment:Range",fm1r40)[[5]][3,1]))*100 #-0.4 % wide
((exp(effect("Treatment:Range",fm1r40)[[5]][2,1])-exp(effect("Treatment:Range",fm1r40)[[5]][1,1]))/
  exp(effect("Treatment:Range",fm1r40)[[5]][1,1]))*100 #+0.33% narrow

fm1r40.2 <- lme((dydt)~Treatment+Location+Range+Treatment:Location+Treatment:Range,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40)#, method="ML")
plot(fm1r40.2,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r40.2,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r40.2,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r40.2, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r40.2$residuals[,1])
anova(fm1r40.2)
anova(fm1r40.2,fm1r40)

fm1r40.3 <- lme((dydt)~Treatment+Location+Range+Treatment:Range,
                random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40, method="ML")
anova(fm1r40.2,fm1r40.3)

plot(effect("Treatment:Range",fm1r40.2), multiline = T) #0.005 warming increased RGR more in wide than narrow
((exp(effect("Treatment:Range",fm1r40.2)[[5]][4,1])-exp(effect("Treatment:Range",fm1r40.2)[[5]][3,1]))/
    exp(effect("Treatment:Range",fm1r40.2)[[5]][3,1]))*100 #-0.4 % wide
((exp(effect("Treatment:Range",fm1r40.2)[[5]][2,1])-exp(effect("Treatment:Range",fm1r40.2)[[5]][1,1]))/
    exp(effect("Treatment:Range",fm1r40.2)[[5]][1,1]))*100 #+0.33% narrow

# 
# T50<-subset(gamfits2, Time==50) 
# fm1r50 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T50)#, method="ML")
# plot(fm1r50,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r50,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r50,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r50, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r50$residuals[,1])
# anova(fm1r50)    
# plot(effect("Treatment",fm1r50))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fm1r50))                     #0.001 RGR lower in N
# plot(effect("Treatment:Location",fm1r50), multiline = T) #0.006 warming decreased RGR more in S than N
# ((effect("Treatment:Location",fm1r50)[[5]][4,1]-effect("Treatment:Location",fm1r50)[[5]][3,1])/
#   effect("Treatment:Location",fm1r50)[[5]][3,1])*100 #-2.9 % N
# ((effect("Treatment:Location",fm1r50)[[5]][2,1]-effect("Treatment:Location",fm1r50)[[5]][1,1])/
#   effect("Treatment:Location",fm1r50)[[5]][1,1])*100 #-14% S
# plot(effect("Treatment:Range",fm1r50), multiline = T) #0.0008 warming decreased RGR in wide not Narrow
# ((effect("Treatment:Range",fm1r50)[[5]][4,1]-effect("Treatment:Range",fm1r50)[[5]][3,1])/
#   effect("Treatment:Range",fm1r50)[[5]][3,1])*100 #-14.7 % wide
# ((effect("Treatment:Range",fm1r50)[[5]][2,1]-effect("Treatment:Range",fm1r50)[[5]][1,1])/
#   effect("Treatment:Range",fm1r50)[[5]][1,1])*100 #+0.4% narrow
# plot(effect("Treatment:Location:Range",fm1r50), multiline = T) #0.005 warming increased RGR more in S wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][8,1]-effect("Treatment:Location:Range",fm1r50)[[5]][7,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][7,1])*100 #-3.4 % N wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][4,1]-effect("Treatment:Location:Range",fm1r50)[[5]][3,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][3,1])*100 #-2% N narrow
# ((effect("Treatment:Location:Range",fm1r50)[[5]][6,1]-effect("Treatment:Location:Range",fm1r50)[[5]][5,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][5,1])*100 #-22.7 % S wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][2,1]-effect("Treatment:Location:Range",fm1r50)[[5]][1,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][1,1])*100 #+2.1% S narrow
# 
# T60<-subset(gamfits2, Time==60) 
# 
# fm1r60 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
# plot(fm1r60,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r60,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r60,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r60, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r60$residuals[,1])
# anova(fm1r60)    
# plot(effect("Treatment",fm1r60))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fm1r60))                     #<.0001 RGR lower in N
# plot(effect("Treatment:Range",fm1r60), multiline = T) #0.053 warming decreased RGR more in wide than narrow
# ((effect("Treatment:Range",fm1r60)[[5]][4,1]-effect("Treatment:Range",fm1r60)[[5]][3,1])/
#   effect("Treatment:Range",fm1r60)[[5]][3,1])*100 #-16.9 % Wide
# ((effect("Treatment:Range",fm1r60)[[5]][2,1]-effect("Treatment:Range",fm1r60)[[5]][1,1])/
#   effect("Treatment:Range",fm1r60)[[5]][1,1])*100 #+5.8% Narrow
# 
# plot(effect("Treatment:Location:Range",fm1r60), multiline = T) #0.0168 warming increased RGR more in S wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][8,1]-effect("Treatment:Location:Range",fm1r60)[[5]][7,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][7,1])*100 #-9.4 % N wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][4,1]-effect("Treatment:Location:Range",fm1r60)[[5]][3,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][3,1])*100 #-13.9 N narrow
# ((effect("Treatment:Location:Range",fm1r60)[[5]][6,1]-effect("Treatment:Location:Range",fm1r60)[[5]][5,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][5,1])*100 #-22.5 % S wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][2,1]-effect("Treatment:Location:Range",fm1r60)[[5]][1,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][1,1])*100 #+0.51% S narrow
# 

###-----------------------------
# #RGR at common mass - 
# # fit GAM lnTotMass ~ s(TotMass, k=kgam) # is this really the same RGR?
# 
# T10<-subset(gamfits2mass, Mass==10) 
# fm1r10 <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T10)#, method="ML")
# plot(fm1r10,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r10,log(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r10,log(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r10, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r10$residuals[,1])
# anova(fm1r10)    
# plot(effect("Treatment",fm1r10))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fm1r10))                     #0.0001 RGR lower in N
# plot(effect("Treatment:Location",fm1r10), multiline = T) #<0.0001 warming increased RGR more in S than N
# ((exp(effect("Treatment:Location",fm1r10)[[5]][4,1])-exp(effect("Treatment:Location",fm1r10)[[5]][3,1]))/
#   exp(effect("Treatment:Location",fm1r10)[[5]][3,1]))*100 # +5% N
# ((exp(effect("Treatment:Location",fm1r10)[[5]][2,1])-exp(effect("Treatment:Location",fm1r10)[[5]][1,1]))/
#   exp(effect("Treatment:Location",fm1r10)[[5]][1,1]))*100 #-28.6% S
# 
# T20<-subset(gamfits2mass, Mass==20) 
# fm1r20 <- lme(log(dydt+1)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T20)#, method="ML")
# plot(fm1r20,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r20,log(dydt+1)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r20,log(dydt+1)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r20, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r20$residuals[,1])
# anova(fm1r20)    
# plot(effect("Treatment",fm1r20))                    #0.0003 RGR decreased with warming
# plot(effect("Location",fm1r20))                     #0.0002 RGR lower in N
# plot(effect("Treatment:Location",fm1r20), multiline = T) #<.0001 warming increased RGR in S not N
# ((exp(effect("Treatment:Location",fm1r20)[[5]][4,1]-1)-exp(effect("Treatment:Location",fm1r20)[[5]][3,1]-1))/
#   exp(effect("Treatment:Location",fm1r20)[[5]][3,1]-1))*100 #-0.14 % N
# ((exp(effect("Treatment:Location",fm1r20)[[5]][2,1]-1)-exp(effect("Treatment:Location",fm1r20)[[5]][1,1]-1))/
#   exp(effect("Treatment:Location",fm1r20)[[5]][1,1]-1))*100 #-4.4% S
# 
# T30<-subset(gamfits2mass, Mass==30) 
# fm1r30 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T30)#, method="ML")
# plot(fm1r30,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r30,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r30,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r30, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r30$residuals[,1])
# anova(fm1r30)    
# plot(effect("Treatment:Location",fm1r30), multiline = T) #<.0001 warming increased RGR more in S than N
# (((1/effect("Treatment:Location",fm1r30)[[5]][4,1])-(1/effect("Treatment:Location",fm1r30)[[5]][3,1]))/
#   (1/effect("Treatment:Location",fm1r30)[[5]][3,1]))*100 #-6 % N
# (((1/effect("Treatment:Location",fm1r30)[[5]][2,1])-(1/effect("Treatment:Location",fm1r30)[[5]][1,1]))/
#   (1/effect("Treatment:Location",fm1r30)[[5]][1,1]))*100 #+8 % S
# 
# T40<-subset(gamfits2mass, Mass==40) 
# fm1r40 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T40)#, method="ML")
# plot(fm1r40,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r40,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r40,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r40, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r40$residuals[,1])
# anova(fm1r40)    
# plot(effect("Location",fm1r40))                     #0.007 RGR lower in N
# plot(effect("Treatment:Location",fm1r40), multiline = T) #0.08 warming increased RGR more in N than S
# ((effect("Treatment:Location",fm1r40)[[5]][4,1]-effect("Treatment:Location",fm1r40)[[5]][3,1])/
#   effect("Treatment:Location",fm1r40)[[5]][3,1])*100 #+1.7 % N
# ((effect("Treatment:Location",fm1r40)[[5]][2,1]-effect("Treatment:Location",fm1r40)[[5]][1,1])/
#   effect("Treatment:Location",fm1r40)[[5]][1,1])*100 #-5.5% S
# plot(effect("Treatment:Range",fm1r40), multiline = T) #0.005 warming increased RGR more in wide than narrow
# ((exp(effect("Treatment:Range",fm1r40)[[5]][4,1])-exp(effect("Treatment:Range",fm1r40)[[5]][3,1]))/
#   exp(effect("Treatment:Range",fm1r40)[[5]][3,1]))*100 #-0.4 % wide
# ((exp(effect("Treatment:Range",fm1r40)[[5]][2,1])-exp(effect("Treatment:Range",fm1r40)[[5]][1,1]))/
#   exp(effect("Treatment:Range",fm1r40)[[5]][1,1]))*100 #+0.33% narrow
# 
# T50<-subset(gamfits2mass, Mass==50) 
# fm1r50 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T50)#, method="ML")
# plot(fm1r50,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r50,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r50,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r50, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r50$residuals[,1])
# anova(fm1r50)    
# plot(effect("Treatment",fm1r50))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fm1r50))                     #0.001 RGR lower in N
# plot(effect("Treatment:Location",fm1r50), multiline = T) #0.006 warming decreased RGR more in S than N
# ((effect("Treatment:Location",fm1r50)[[5]][4,1]-effect("Treatment:Location",fm1r50)[[5]][3,1])/
#   effect("Treatment:Location",fm1r50)[[5]][3,1])*100 #-2.9 % N
# ((effect("Treatment:Location",fm1r50)[[5]][2,1]-effect("Treatment:Location",fm1r50)[[5]][1,1])/
#   effect("Treatment:Location",fm1r50)[[5]][1,1])*100 #-14% S
# plot(effect("Treatment:Range",fm1r50), multiline = T) #0.0008 warming decreased RGR in wide not Narrow
# ((effect("Treatment:Range",fm1r50)[[5]][4,1]-effect("Treatment:Range",fm1r50)[[5]][3,1])/
#   effect("Treatment:Range",fm1r50)[[5]][3,1])*100 #-14.7 % wide
# ((effect("Treatment:Range",fm1r50)[[5]][2,1]-effect("Treatment:Range",fm1r50)[[5]][1,1])/
#   effect("Treatment:Range",fm1r50)[[5]][1,1])*100 #+0.4% narrow
# plot(effect("Treatment:Location:Range",fm1r50), multiline = T) #0.005 warming increased RGR more in S wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][8,1]-effect("Treatment:Location:Range",fm1r50)[[5]][7,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][7,1])*100 #-3.4 % N wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][4,1]-effect("Treatment:Location:Range",fm1r50)[[5]][3,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][3,1])*100 #-2% N narrow
# ((effect("Treatment:Location:Range",fm1r50)[[5]][6,1]-effect("Treatment:Location:Range",fm1r50)[[5]][5,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][5,1])*100 #-22.7 % S wide
# ((effect("Treatment:Location:Range",fm1r50)[[5]][2,1]-effect("Treatment:Location:Range",fm1r50)[[5]][1,1])/
#   effect("Treatment:Location:Range",fm1r50)[[5]][1,1])*100 #+2.1% S narrow
# 
# T60<-subset(gamfits2mass, Mass==60) 
# 
# fm1r60 <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=T60)#, method="ML")
# plot(fm1r60,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm1r60,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm1r60,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm1r60, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm1r60$residuals[,1])
# anova(fm1r60)    
# plot(effect("Treatment",fm1r60))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fm1r60))                     #<.0001 RGR lower in N
# plot(effect("Treatment:Range",fm1r60), multiline = T) #0.053 warming decreased RGR more in wide than narrow
# ((effect("Treatment:Range",fm1r60)[[5]][4,1]-effect("Treatment:Range",fm1r60)[[5]][3,1])/
#   effect("Treatment:Range",fm1r60)[[5]][3,1])*100 #-16.9 % Wide
# ((effect("Treatment:Range",fm1r60)[[5]][2,1]-effect("Treatment:Range",fm1r60)[[5]][1,1])/
#   effect("Treatment:Range",fm1r60)[[5]][1,1])*100 #+5.8% Narrow
# 
# plot(effect("Treatment:Location:Range",fm1r60), multiline = T) #0.0168 warming increased RGR more in S wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][8,1]-effect("Treatment:Location:Range",fm1r60)[[5]][7,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][7,1])*100 #-9.4 % N wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][4,1]-effect("Treatment:Location:Range",fm1r60)[[5]][3,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][3,1])*100 #-13.9 N narrow
# ((effect("Treatment:Location:Range",fm1r60)[[5]][6,1]-effect("Treatment:Location:Range",fm1r60)[[5]][5,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][5,1])*100 #-22.5 % S wide
# ((effect("Treatment:Location:Range",fm1r60)[[5]][2,1]-effect("Treatment:Location:Range",fm1r60)[[5]][1,1])/
#   effect("Treatment:Location:Range",fm1r60)[[5]][1,1])*100 #+0.51% S narrow
# 
# 
####-----------------------------------------------------------------------------------
## Average RGR at common mass
M1<- subset(gampreds2, Mass==8)#2,4,8,12
fmr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=M1)#, method="ML")
plot(fmr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fmr,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fmr,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fmr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fmr$residuals[,1])
anova(fmr)
plot(effect("Treatment",fmr))                    #<.0001 RGR increased with warming
plot(effect("Treatment:Location",fmr))                     #<.0001 RGR increased in S not in N
plot(effect("Treatment:Range",fmr))

((effect("Treatment:Location",fmr)[[5]][4,1]-effect("Treatment:Location",fmr)[[5]][3,1])/
  effect("Treatment:Location",fmr)[[5]][3,1])*100 #-1.1 % N
((effect("Treatment:Location",fmr)[[5]][2,1]-effect("Treatment:Location",fmr)[[5]][1,1])/
  effect("Treatment:Location",fmr)[[5]][1,1])*100 #+17.1% S


# 
# #get the highest minimum mass of each treatment group
# max(summaryBy(predMass~Taxa+Treatment, data=gamfits2, FUN=c(min))$predMass.min) #DCAM HOME 0.6431564
# 
# summaryBy(predMass~Taxa+Treatment, data=subset(gamfits2, predMass > 0.8 & predMass< 0.8+0.2), 
#           FUN=length) #find range that encompasses at least 10 reps of each Treatment combo
# 
# RGRm<-droplevels(subset(gamfits2, predMass >0.8 & predMass<0.8+0.2))
# 
# fmr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=RGRm)#, method="ML")
# plot(fmr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fmr,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fmr,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fmr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fmr$residuals[,1])
# anova(fmr)    
# plot(effect("Treatment",fmr))                    #<.0001 RGR decreased with warming
# plot(effect("Location",fmr))                     #0.0044 RGR lower in N
# plot(effect("Treatment:Range",fmr), multiline = T) #0.0393 warming increased RGR more in wide than narrow
# 
# 
# #get the smallest maximum mass of each treatment group
# min(summaryBy(predMass~Taxa+Treatment, data=gamfits2, FUN=c(max))$predMass.max) #SMIT HOME 33.0923
# 
# summaryBy(predMass~Taxa+Treatment, data=subset(gamfits2, predMass >12 & predMass<14), 
#           FUN=length) #find range that encompasses at least 10 reps of each Treatment combo
# 
# RGRm<-droplevels(subset(gamfits2, predMass >12 & predMass<14))
# 
# fmr <- lme((dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=RGRm)#, method="ML")
# plot(fmr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fmr,(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fmr,(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fmr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fmr$residuals[,1])
# anova(fmr)    
# plot(effect("Treatment",fmr))                    #<.0001 RGR increased with warming
# plot(effect("Treatment:Location",fmr))                     #<.0001 RGR increased in S not in N
# plot(effect("Treatment:Range",fmr)) 
# 
# ((effect("Treatment:Location",fmr)[[5]][4,1]-effect("Treatment:Location",fmr)[[5]][3,1])/
#   effect("Treatment:Location",fmr)[[5]][3,1])*100 #-1.1 % N
# ((effect("Treatment:Location",fmr)[[5]][2,1]-effect("Treatment:Location",fmr)[[5]][1,1])/
#   effect("Treatment:Location",fmr)[[5]][1,1])*100 #+17.1% S
# 
# #get a middle value as well
# 
# summaryBy(predMass~Taxa+Treatment, data=subset(gamfits2, predMass >=6.5 & predMass<=7.5), 
#           FUN=length) #find range that encompasses at least 10 reps of each Treatment combo
# 
# RGRm<-droplevels(subset(gamfits2, predMass >=6.5 & predMass<=7.5))
# 
# 
# fmr <- lme(log(dydt)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=RGRm)#, method="ML")
# plot(fmr,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fmr,log(dydt)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fmr,log(dydt)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fmr, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fmr$residuals[,1])
# anova(fmr)    
# plot(effect("Treatment",fmr))                    #<.0001 RGR increased with warming
# plot(effect("Treatment:Location",fmr))                     #<.0001 RGR increased in S not in N
# plot(effect("Treatment:Range",fmr))              #0.0063 warming increased RGR more in wide than narrow
# 
# ((exp(effect("Treatment:Location",fmr)[[5]][4,1])-exp(effect("Treatment:Location",fmr)[[5]][3,1]))/
#     exp(effect("Treatment:Location",fmr)[[5]][3,1]))*100 #-5.8 % N
# ((exp(effect("Treatment:Location",fmr)[[5]][2,1])-exp(effect("Treatment:Location",fmr)[[5]][1,1]))/
#     exp(effect("Treatment:Location",fmr)[[5]][1,1]))*100 #+18.3% S
