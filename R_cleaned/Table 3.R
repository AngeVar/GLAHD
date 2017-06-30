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


fm.Rdark <- lme(log(R)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdark)
plot(fm.Rdark,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rdark,log(R)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rdark,log(R)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rdark, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rdark$residuals[,1])
anova(fm.Rdark)    
effect("Treatment:Location",fm.Rdark)
(exp(-0.8860410)-exp(-0.4899468))/(exp(-0.4899468)) # -32.7% in S
(exp(-0.7644012)-exp(-0.5523989))/(exp(-0.5523989)) # -19.1% in N

# --- Area based results w/o SMIT  ------------------------------------------------------------

#SMIT may be a big influencer of results - remove SMIT
Asats<- subset(Asat, Taxa != "SMIT")
acifitss<- subset(acifits, Taxa != "SMIT")
Rdarks<- subset(Rdark, Taxa != "SMIT")

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

fm.Rdark <- lme(log(R)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdarks)
plot(fm.Rdark,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rdark,R~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rdark,R~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rdark, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rdark$residuals[,1])
anova(fm.Rdark)    
effect("Treatment:Location",fm.Rdark)
(exp(-0.9234425)-exp(-0.5232258))/(exp(-0.5232258)) # -33.0% in S
(exp(-0.7610750)-exp(-0.5512533))/(exp(-0.5512533)) # -18.9% in N

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
plot(effect("Location",fm.vcmaxm)) 
((1.287201)^2-(1.848172)^2)/((1.848172)^2) #51.5 % lower in S than N
plot(effect("Treatment:Location",fm.vcmaxm)) 
((1.281249)^2-(1.293215)^2)/((1.293215)^2) #-1.8% in S
((1.811152)^2-(1.885578)^2)/((1.885578)^2) #-7.7% in N
((1.558556)^2-(1.603208)^2)/((1.603208)^2) #-5.5% in warmed than Home

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

fm.Rdark <- lme(log(Rmass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdark)
plot(fm.Rdark,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rdark,log(Rmass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rdark,log(Rmass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rdark, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rdark$residuals[,1])
anova(fm.Rdark)    
plot(effect("Treatment:Location:Range",fm.Rdark), multiline=T)
(exp(2.253351)-exp(2.315649))/(exp(2.315649)) # -6% in S narrow
(exp(2.144013)-exp(2.361253))/(exp(2.361253)) # -19.5% N narrow

(exp(1.974925)-exp(2.318532))/(exp(2.318532)) # -29.1% in S wide
(exp(2.174630)-exp(2.426073))/(exp(2.426073)) # -22.2% in N wide
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

fm.Rdark <- lme((Rmass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdarks)
plot(fm.Rdark,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rdark,(Rmass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rdark,(Rmass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rdark, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rdark$residuals[,1])
anova(fm.Rdark)    
plot(effect("Treatment:Range",fm.Rdark), multiline=T)
((8.305792)-(10.098900))/((10.098900)) # -17.8% in narrow
((8.199173)-(11.074920))/((11.074920)) # -26.0% in wide

plot(effect("Treatment:Location:Range",fm.Rdark), multiline=T)
((7.965108)-(9.286351))/((9.286351)) # -14.2% in S narrow
((8.58541)-(10.76580))/((10.76580)) # -20.3% N narrow

((7.344985)-(10.445108))/((10.445108)) # -29.7% in S wide
((8.900251)-(11.591842))/((11.591842)) # -23.2% in N wide

plot(effect("Treatment",fm.Rdark), multiline=T)
((8.232871)-(10.766437))/((10.766437)) # -23.5% in warmed

plot(effect("Location",fm.Rdark), multiline=T)
((10.032980)-(8.777051))/((8.777051)) # -14.3 % higher in N

###########################################################################################################

#Test R at 25 from temperature sensitivity curves

fm2<- aov((R18.5)~Treatment*Taxa, data=params)
hist(fm2$residuals)
anova(fm2)#0.2,0.001,
TukeyHSD(fm2)
plot(allEffects(fm2), multiline=T) #Bot had higher rates than all others 

fm1<- aov((R25)~Treatment*Taxa, data=params)
hist(fm1$residuals)
anova(fm1)#0.06,0.02 #T effect is marginal, Taxa differs at 25 C
TukeyHSD(fm1) #BOT had higher rates than BRA
plot(allEffects(fm1), multiline=T) 

fm3<- aov((R28.5)~Treatment*Taxa, data=params)
hist(fm3$residuals)
anova(fm3)#0.03,0.11 #Warming reduced Rates, taxa did not differ.
TukeyHSD(fm3)
plot(allEffects(fm3), multiline=T) 
effect("Treatment",fm3)
(13.41084-15.10342)/15.10342 #-11.2%

fm2<- aov((Q10_18.5)~Treatment*Taxa, data=params)
hist(fm2$residuals)
anova(fm2)#0.08,0.03, #taxa differed in Q10
TukeyHSD(fm2)#CTER and BOT
plot(allEffects(fm2), multiline=T)  #Bot had lower Q10 than CTER

fm1<- aov(log(Q10_25)~Treatment*Taxa, data=params)
hist(fm1$residuals)
anova(fm1)#0.1,0.01 #taxa differed in Q10
TukeyHSD(fm1) #CTER and BOT
plot(allEffects(fm1), multiline=T) #Bot had lower Q10 than CTER

fm3<- aov(log(Q10_28.5)~Treatment*Taxa, data=params)
hist(fm3$residuals)
anova(fm3)#0.1,0.008
TukeyHSD(fm3)#CTER and BOT
plot(allEffects(fm3), multiline=T) #Bot had lower Q10 than CTER

#I have cleaned up this analysis from the old script so that it now uses the data that is on HIEv and uses
#the modified Q10 model (the model with the lowest overall AIC) to predict rates at other temperatures 
#instead of the polynomial model as in the old script.
#I cannot reproduce the results John got way back whan he did it the first time.
#But I find that:
#The taxa differ in their Rmass at 25 (P=0.024), 18.5 (P=0.001) but not 28.5 (P=0.11)
#Tukeys HSD suggests that it is due to BRA differing from BOT (narrow Tropical vs. Temperate)
#There is little effect of warming at 25 (0.064), 18.5 (0.24) but at 28.5 (0.032) rates were significantly reduced -11.2%
#For Q10, Taxa differed in their Q10 at 25 (0.013), 18.5 (0.03) and 28.5 (0.008)
#Tukeys HSD suggests it is due to BOT and BTER differing (Temperate narrow vs. Wide)
#There was little effect of temperature at 25 (0.1), 18.5 (0.08) and 28.5 (0.1)

###########################################################################################################

#calculate response ratios
asat.tm <- summaryBy(Photo~Taxa+Treatment+Location+Range,data=Asat,FUN=mean,keep.names=T)
R.tm <- summaryBy(R~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)
jmax.tm <- summaryBy(Jmax~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
vcmax.tm <- summaryBy(Vcmax~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)

j.taxa.l<- split(jmax.tm,as.factor(jmax.tm$Taxa))
v.taxa.l<- split(vcmax.tm,as.factor(vcmax.tm$Taxa))
a.taxa.l<- split(asat.tm,as.factor(asat.tm$Taxa))
r.taxa.l<- split(Rmass.tm,as.factor(R.tm$Taxa))

output<- list()
for(i in 1:length(j.taxa.l)){
  rat<-(j.taxa.l[[i]]$Jmax[2]-j.taxa.l[[i]]$Jmax[1])/j.taxa.l[[i]]$Jmax[1]
  
  output[[i]] <- data.frame(Taxa=j.taxa.l[[i]]$Taxa[1],Location=j.taxa.l[[i]]$Location[1],Range=j.taxa.l[[i]]$Range[1],JRatio=rat)
}
j.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(v.taxa.l)){
  rat<-(v.taxa.l[[i]]$Vcmax[2]-v.taxa.l[[i]]$Vcmax[1])/v.taxa.l[[i]]$Vcmax[1]
  
  output[[i]] <- data.frame(Taxa=v.taxa.l[[i]]$Taxa[1],Location=v.taxa.l[[i]]$Location[1],Range=v.taxa.l[[i]]$Range[1],vRatio=rat)
}
v.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(a.taxa.l)){
  rat<-(a.taxa.l[[i]]$Photo[2]-a.taxa.l[[i]]$Photo[1])/a.taxa.l[[i]]$Photo[1]
  
  output[[i]] <- data.frame(Taxa=a.taxa.l[[i]]$Taxa[1],Location=a.taxa.l[[i]]$Location[1],Range=a.taxa.l[[i]]$Range[1],aRatio=rat)
}
a.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(r.taxa.l)){
  rat<-(r.taxa.l[[i]]$Rmass[2]-r.taxa.l[[i]]$Rmass[1])/r.taxa.l[[i]]$Rmass[1]
  
  output[[i]] <- data.frame(Taxa=r.taxa.l[[i]]$Taxa[1],Location=r.taxa.l[[i]]$Location[1],Range=r.taxa.l[[i]]$Range[1],rRatio=rat)
}
r.df <- do.call(rbind,output)


asatm.tm <- summaryBy(Photom~Taxa+Treatment+Location+Range,data=Asatm,FUN=mean,keep.names=T)
Rmass.tm <- summaryBy(Rmass~Taxa+Treatment+Location+Range,data=Rdark,FUN=mean,keep.names=T)
jmaxm.tm <- summaryBy(Jmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
vcmaxm.tm <- summaryBy(Vcmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)

jm.taxa.l<- split(jmaxm.tm,as.factor(jmaxm.tm$Taxa))
vm.taxa.l<- split(vcmaxm.tm,as.factor(vcmaxm.tm$Taxa))
am.taxa.l<- split(asatm.tm,as.factor(asatm.tm$Taxa))
rm.taxa.l<- split(Rmass.tm,as.factor(Rmass.tm$Taxa))

output<- list()
for(i in 1:length(j.taxa.l)){
  rat<-(jm.taxa.l[[i]]$Jmaxm[2]-jm.taxa.l[[i]]$Jmaxm[1])/jm.taxa.l[[i]]$Jmaxm[1]
  
  output[[i]] <- data.frame(Taxa=jm.taxa.l[[i]]$Taxa[1],Location=jm.taxa.l[[i]]$Location[1],Range=jm.taxa.l[[i]]$Range[1],JmRatio=rat)
}
jm.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(v.taxa.l)){
  rat<-(vm.taxa.l[[i]]$Vcmaxm[2]-vm.taxa.l[[i]]$Vcmaxm[1])/vm.taxa.l[[i]]$Vcmaxm[1]
  
  output[[i]] <- data.frame(Taxa=v.taxa.l[[i]]$Taxa[1],Location=v.taxa.l[[i]]$Location[1],Range=v.taxa.l[[i]]$Range[1],vmRatio=rat)
}
vm.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(am.taxa.l)){
  rat<-(am.taxa.l[[i]]$Photom[2]-am.taxa.l[[i]]$Photom[1])/am.taxa.l[[i]]$Photom[1]
  
  output[[i]] <- data.frame(Taxa=am.taxa.l[[i]]$Taxa[1],Location=am.taxa.l[[i]]$Location[1],Range=am.taxa.l[[i]]$Range[1],amRatio=rat)
}
am.df <- do.call(rbind,output)

output<- list()
for(i in 1:length(rm.taxa.l)){
  rat<-(rm.taxa.l[[i]]$Rmass[2]-rm.taxa.l[[i]]$Rmass[1])/rm.taxa.l[[i]]$Rmass[1]
  
  output[[i]] <- data.frame(Taxa=rm.taxa.l[[i]]$Taxa[1],Location=rm.taxa.l[[i]]$Location[1],Range=rm.taxa.l[[i]]$Range[1],rmRatio=rat)
}
rm.df <- do.call(rbind,output)

RR<-Reduce(function(x, y) merge(x, y, all=TRUE), list(j.df,v.df,a.df,r.df,jm.df,vm.df,am.df,rm.df))



#--- compare SLA

f.AsatSLA <- lme((SLA)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asatm)
plot(f.AsatSLA,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(f.AsatSLA,(SLA)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(f.AsatSLA,(SLA)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(f.AsatSLA, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(f.AsatSLA$residuals[,1])
anova(f.AsatSLA)    #no effects, still no effects if SMIT is removed
plot(effect("Treatment:Location",f.AsatSLA))
((194.474)-(167.919))/((167.919))    # +15.8% increase in SLA in S # 
((189.5643)-(193.3252))/((193.3252)) #  -1.9% decrease in SLA in N
plot(effect("Treatment:Range",f.AsatSLA))
((206.2998)-(182.5428))/((182.5428))    # +13.0% increase in SLA in Narrow # 
((184.1988)-(180.4939))/((180.4939)) #  +2.1% increase in SLA in Wide


#-------------------------------------

#test after fit with TPU using Dushan's script
dg3<-merge(dg2,SLARd, by="Code")            #give SLA to the codes that we have Jmax and Vcmax data for
dg3$Jmaxm<-with(dg3, Jmax*SLA/10000)
dg3$Vcmaxm<-with(dg3, Vcmax*SLA/10000)
dg3$JtoVm<-with(dg3,Jmaxm/Vcmaxm)

#not much difference in vcmax by fitting TPU
fm.vcmax <- lme(log(Vcmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dg3)
plot(fm.vcmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmax,(Vcmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmax,(Vcmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmax$residuals[,1])
anova(fm.vcmax)    
plot(effect("Range",fm.vcmax))      #- Vcmax higher in wide than narrow
(exp(4.956607)-exp(4.810580))/(exp(4.810580))    # 15.7% higher in wide

plot(effect("Treatment:Location",fm.vcmax))      #- warming reduced Vcmax in the south but not in the north. No range interactions
(exp(4.484568)-exp(4.632290))/(exp(4.632290))    # -13.7 % in Vcmax in S
(exp(5.213920)-exp(5.228832))/(exp(5.228832))    #  -0.15%  in Vcmax in N 

#some missing Jmax values - 1 BOT home, 1 BOT warmed, 1 BRA Home, 1 BTER Home, 1 PLAT home
#adding TPU limitation tot he model does not change results much. 
#The main effect of location appears where N>S, but this is trumped by the interaction with treatment anyway.
fm.jmax <- lme(log(Jmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dg3[complete.cases(dg3[, c("Jmax")]), ])
plot(fm.jmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmax,(Jmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmax,(Jmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmax$residuals[,1])
anova(fm.jmax)    
plot(effect("Treatment:Location",fm.jmax))
(exp(4.978871)-exp(5.210569))/(exp(5.210569))    # -20.6% in Jmax in S
(exp(5.178931)-exp(5.255176))/(exp(5.255176))    #  -7.3%  in Jmax in N 

#the JV ratio declined with warming in all groups except narrow N

fm.jtov <- lme(log(JtoV)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dg3[complete.cases(dg3[, c("Jmax")]), ])
plot(fm.jtov,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtov,(JtoV)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtov,(JtoV)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtov, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtov$residuals[,1])
anova(fm.jtov)
plot(effect("Treatment",fm.jtov))
(exp(0.2109301)-exp(0.2827220))/(exp(0.2827220))      # -6.9% in Jmax/Vcmax with warming
plot(effect("Location",fm.jtov))
(exp(-0.008336661)-exp(0.529045483))/(exp(0.529045483))      # 41.6% lower Jmax/Vcmax in N compared to S
plot(effect("Treatment:Location:Range",fm.jtov), multiline=T)
(exp(0.5211215)-exp(0.6297080))/(exp(0.6297080))          #-10.3 % in narrow S
(exp(0.4679126)-exp(0.5451491))/(exp(0.5451491))          #-7.4 % in wide S
(exp(-0.010326059)-exp(-0.003742756))/(exp(-0.003742756)) #-0.7 % in narrow N
(exp(-0.04982495)-exp(0.03365585))/(exp(0.03365585))      #-8 % in wide N

#only 14 replicates had no TPU limitation; ATER(1), BTER(1), CCAM(3), ECAM(1), LONG(3)
fm.tpu <- lme(log(TPU)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dg3[complete.cases(dg3[, c("TPU")]), ])
plot(fm.tpu,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.tpu,(TPU)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.tpu,(TPU)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.tpu, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.tpu$residuals[,1])
anova(fm.tpu)
plot(effect("Treatment:Location",fm.tpu))
(exp(2.250428)-exp(2.455829))/(exp(2.455829))    # -18.6% in TPU in S
(exp(2.254111)-exp(2.351137))/(exp(2.351137))    #  -9.2%  in TPU in N 
