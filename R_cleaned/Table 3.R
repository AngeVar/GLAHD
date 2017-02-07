# table 3
source("R_cleaned/1. GLAHD_LoadLibraries.R")
source("R_cleaned/2. downloadData.R")
source("R_cleaned/3. Create_datasets.R")

Asatm$Photom<-with(Asatm, (Photo*SLA/10000)*1000)#from micro to nano
acifits$Jmaxm<-with(acifits, Jmax*SLA/10000)
acifits$Vcmaxm<-with(acifits, Vcmax*SLA/10000)
acifits$JtoVm<-with(acifits,Jmaxm/Vcmaxm)


###
# SMIT really changes things.
# 
# For area based rates: the only difference is that the range effect on Vcmax becomes insignificant and rates are 
# no longer higher in wide than narrow (goes from being 14% higher to only 11 % higher than narrow).
# 
# For mass based rates: If SMIT is removed the Treatment * Range effects disappear for both Vcmax and Jmax as the 
# warming effect becomes more negative. For Vcmax narrow species go from +2.4% to -3.6 % with warming. For Jmax 
# narrow species change from being only -3.5 % to being -10 % reduced by warming.
# 
# I.e. SMIT did not acclimate as well as the other narrow species and it was not just because of SLA differences 
# as it still makes a big difference whether or not including SMIT in mass based rates.
###

## area based results with SMIT
fm.Asat <- lme(Photo~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asat)
plot(fm.Asat,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asat,Photo~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asat,Photo~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asat, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asat$residuals[,1])
anova(fm.Asat)    
effect("Treatment",fm.Asat)
(26.72886-28.35335)/28.35335 #-5.7% with warming

fm.vcmax <- lme(log(Vcmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.vcmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmax,log(Vcmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmax,log(Vcmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmax$residuals[,1])
anova(fm.vcmax)    
plot(effect("Range",fm.vcmax))      #- Vcmax higher in wide than narrow
(exp(4.911595)-exp(4.783875))/(exp(4.783875))    # 13.6% higher in wide

plot(effect("Treatment:Location",fm.vcmax))      #- warming reduced Vcmax in the south but not in the north. No range interactions
(exp(4.433354)-exp(4.612889))/(exp(4.612889))    # -16.4 % in Vcmax in S
(exp(5.164627)-exp(5.195303))/(exp(5.195303))    #  -3%  in Vcmax in N 


fm.jmax <- lme(log(Jmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmax,log(Jmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmax,log(Jmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmax$residuals[,1])
anova(fm.jmax)    
plot(effect("Treatment:Location",fm.jmax))
(exp(4.914717)-exp(5.153855))/(exp(5.153855))    # -21.3% in Jmax in S
(exp(5.006654)-exp(5.108184))/(exp(5.108184))    #  -9.7%  in Jmax in N 


fm.jtov <- lme(log(JtoV)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jtov,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtov,log(JtoV)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtov,log(JtoV)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtov, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtov$residuals[,1])
anova(fm.jtov)
effect("Treatment",fm.jtov)
(exp(0.1468196)-exp(0.2121764))/(exp(0.2121764))      # -6.3% in Jmax/Vcmax with warming
effect("Location",fm.jtov)
(exp(-0.1228216)-exp(0.5110372))/(exp(0.5110372))      # 46.9% lower Jmax/Vcmax in N compared to S


# --- Area based results w/o SMIT  ------------------------------------------------------------

#SMIT may be a big influencer of results - remove SMIT
Asats<- subset(Asat, Taxa != "SMIT")
acifitss<- subset(acifits, Taxa != "SMIT")


fm.Asat <- lme(Photo~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asats)
plot(fm.Asat,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asat,Photo~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asat,Photo~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asat, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asat$residuals[,1])
anova(fm.Asat)   
effect("Treatment",fm.Asat)
(26.82011-28.43675)/28.43675  #-5.7% with warming

fm.vcmax <- lme(log(Vcmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.vcmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmax,log(Vcmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmax,log(Vcmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmax$residuals[,1])
anova(fm.vcmax)
plot(effect("Range",fm.vcmax))      #- Vcmax higher in wide than narrow
(exp(4.931121)-exp(4.827046))/(exp(4.827046))    # 10.9% higher in wide

plot(effect("Treatment:Location",fm.vcmax))      #- warming reduced Vcmax in the south but not in the north. No range interactions
(exp(4.440326)-exp(4.649033))/(exp(4.649033))    # -18.8 % in S
(exp(5.167493)-exp(5.199305))/(exp(5.199305))    #  -3%  in Vcmax in N 


fm.jmax <- lme(log(Jmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.jmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmax,log(Jmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmax,log(Jmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmax$residuals[,1])
anova(fm.jmax)    
plot(effect("Treatment:Location",fm.jmax))
(exp(4.912706)-exp(5.176771))/(exp(5.176771))    # -23.2%  in Jmax in S 
(exp(5.011273)-exp(5.113089))/(exp(5.113089))    #  -9.7%  in Jmax in N 


fm.jtov <- lme(log(JtoV)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.jtov,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtov,log(JtoV)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtov,log(JtoV)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtov, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtov$residuals[,1])
anova(fm.jtov)
effect("Treatment",fm.jtov)
(exp(0.1235163)-exp(0.1869347))/(exp(0.1869347))      # -6.1%  in Jmax/Vcmax with warming
effect("Location",fm.jtov)
(exp(-0.1216975)-exp(0.4997415))/(exp(0.4997415))      # -46.3% in N

##----------------------------------
##----------------------------------
##----------------------------------
#    mass based rates

fm.Asatm <- lme(sqrt(Photom)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asatm)
plot(fm.Asatm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asatm,sqrt(Photom)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asatm,sqrt(Photom)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asatm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asatm$residuals[,1])
anova(fm.Asatm)    #no effects, still no effects if SMIT is removed
plot(effect("Treatment:Location",fm.Asatm))
((22.64264)^2-(21.61941)^2)/((21.61941)^2)    # +9.7% increase in Asat in S # 
((22.27189)^2-(23.38046)^2)/((23.38046)^2) #  -9.3% decrease in Asat in N

fm.vcmaxm <- lme(sqrt(Vcmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.vcmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmaxm,sqrt(Vcmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmaxm,sqrt(Vcmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmaxm$residuals[,1])
anova(fm.vcmaxm)    
plot(effect("Treatment:Range",fm.vcmaxm))      #- warming reduced Vcmax in the wide but not in the narrow. No location interactions
((1.554706^2)-(1.536480^2))/((1.536480^2))    # +2.4%  in Vcmax in narrow # 
((1.560650^2)-(1.639508^2))/((1.639508^2)) #  -9.4%  in Vcmax in wide

#- fit and interpret Jmax
fm.jmaxm <- lme(sqrt(Jmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmaxm$residuals[,1])
anova(fm.jmaxm)    
plot(effect("Treatment:Range",fm.jmaxm))      #- warming reduced Jmax in the wide but not in the narrow No location interactions
((1.660084)^2-(1.689630)^2)/((1.689630)^2)    # -3.5%  in Jmax in narrow #
((1.653057)^2-(1.786284)^2)/((1.786284)^2)    #  -14.4%  in Jmax in wide


#- fit and interpret Jmax/Vcmax
fm.jtovm <- lme(log(JtoVm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jtovm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtovm,log(JtoVm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtovm,log(JtoVm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtovm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtovm$residuals[,1])
anova(fm.jtovm)
effect("Treatment",fm.jtovm)
(exp(0.1468196)-exp(0.2121764))/(exp(0.2121764))      # -6.3%  in Jmax/Vcmax with warming
effect("Location",fm.jtovm)
(exp(-0.1228216)-exp(0.5110372))/(exp(0.5110372))      # -46.9% in N


#-------------------------------
#   mass based rates w/o SMIT

#SMIT may be a big influencer of results - remove SMIT
Asatms<- subset(Asatm, Taxa != "SMIT")
acifitss<- subset(acifits, Taxa != "SMIT")

fm.Asatm <- lme(sqrt(Photom)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asatms)
plot(fm.Asatm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asatm,sqrt(Photom)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asatm,sqrt(Photom)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asatm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asatm$residuals[,1])
anova(fm.Asatm)    #no effects, still no effects if SMIT is removed
plot(effect("Treatment:Location",fm.Asatm))
((22.35356)^2-(21.55693)^2)/((21.55693)^2)    # +7.5% increase in Asat in S # 
((22.28974)^2-(23.43371)^2)/((23.43371)^2) #  -9.5% decrease in Asat in N

fm.vcmaxm <- lme(sqrt(Vcmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.vcmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmaxm,sqrt(Vcmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmaxm,sqrt(Vcmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmaxm$residuals[,1])
anova(fm.vcmaxm)    
plot(effect("Treatment:Range",fm.vcmaxm))      #- warming reduced Vcmax in the wide but not in the narrow. No location interactions
((1.548087)^2-(1.577019)^2)/((1.577019)^2)    # -3.6% in Vcmax in narrow # 
((1.577773)^2-(1.658058)^2)/((1.658058)^2) #  -9.4% in Vcmax in N

#- fit and interpret Jmax
fm.jmaxm <- lme(sqrt(Jmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.jmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmaxm,sqrt(Jmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmaxm$residuals[,1])
anova(fm.jmaxm)    
plot(effect("Treatment:Range",fm.jmaxm))      #- warming reduced Jmax in the wide but not in the narrow No location interactions
((1.621423)^2-(1.710628)^2)/((1.710628)^2)    # -10.2% reduction in Jmax in narrow #
((1.655986)^2-(1.790320)^2)/((1.790320)^2) #  -14.4% decrease in Jmax in wide


#- fit and interpret Jmax/Vcmax
fm.jtovm <- lme(log(JtoVm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifitss)
plot(fm.jtovm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtovm,log(JtoVm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtovm,log(JtoVm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtovm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtovm$residuals[,1])
anova(fm.jtovm)
effect("Treatment",fm.jtovm)
(exp(0.1235163)-exp(0.1869347))/(exp(0.1869347))      # -6.1%  in Jmax/Vcmax with warming
effect("Location",fm.jtovm)
(exp(-0.1216975)-exp(0.4997415))/(exp(0.4997415))      # -46.2% in N
