#------------------------------------------------------------------------------------------------------------------------------
# This script provides an example statistical analysis of the final size at the end of the experiment
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")

#- read in the data, do a few conversions
dat2 <- return_size_mass()

#- extract just the final estimate of size, on the last date
dat.f <- subset(dat2,Date==as.Date("2015-01-05"))
dat.f$Species <- factor(dat.f$Species)                     # this factor had 19 levels, which was wrong. It should have 8.
dat.f$Location <- factor(dat.f$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat.f$Sp_RS_EN <- as.factor(with(dat.f,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat.f$Prov_Sp_EN <- as.factor(with(dat.f,paste(Taxa,Species)))
dat.f$Sp_Loc_EN <- as.factor(with(dat.f,paste(Species,Location)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize


#- check structure
str(dat.f) # notice that Treatment, Location, Range, Sp_RS_EN, and Prov_Sp_EN are all factors with the right number of levels


#---- let's think about the experimental design. It's not quite a simple 3-way factorial between Treatment, Location, and Range size,
#           as the observations are not independent of each other, as observations are grouped by species.
#           that is, species and range size or species and location are not independent; there is a restriction of randomization.
#           so the location and range effects will be tested with a different error term than the Treatment effect.

# 5 species in north, 5 species in south. Species{location} error term is 2(5-1)-1 = 7
table(as.factor(with(dat.f,paste(Species,Location))))

# 2 species in wide, 5 species in narrow. Species{range} error term is harder to calculate, but will be smaller than the Species{location} error term
table(as.factor(with(dat.f,paste(Species,Range))))
#----




library(nlme)
library(effects)

#- this model fits the 3-way interaction straight up. All effects are (improperly) tested against the residual error with 261 degrees of freedom
fm0 <- lm(TotMass~Treatment*Location*Range,data=dat.f)
anova(fm0)
boxplot(TotMass~Range,data=dat.f) #- main effect of Range was significant, but just looking at the data this is hard to believe

#- these models account for the fact that trees are not independent, but are nested by species and provenance
fm1 <- lme(TotMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat.f)
fm2 <- lme(TotMass~Treatment*Location*Range,random=list(~1|Species,~1|Taxa),data=dat.f)   # R can determine the random term structure itself (most of the time!)

#look at model diagnostics
plot(fm2,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2,TotMass~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm2,TotMass~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm2, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm2$residuals[,1])

anova(fm2,type="marginal") #- type III, effects assessed by "leave one out"
anova(fm2)                 #- type I, effects added to the model in the order specified
anova(fm1)                 #- exactly the same as the second fitted model
plot(allEffects(fm2))      #- try to make sense of Treatment:Location:Range interaction
plot(effect("Treatment:Location",fm2)) #- warming increases mass in the south, but not in the north


#- refit log-transformed TotMass in an attempt to stabilize the variance
fm2.l <- lme(log(TotMass)~Treatment*Location*Range,random=list(~1|Species,~1|Taxa),data=dat.f)   
plot(fm2.l,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. This looks much better
plot(fm2.l,log(TotMass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm2.l,log(TotMass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm2.l, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals


anova(fm2.l,type="marginal") #- type III, effects assessed by "leave one out"
anova(fm2.l)                 #- type I, effects added to the model in the order specified



#--------------------------------------------------------------------
#- compare models with and without provenance random term. can it be removed?

fm2.full <- lme(log(TotMass)~Treatment*Location*Range,random=list(~1|Species,~1|Taxa),data=dat.f,method="ML")   
fm2.notaxa <- lme(log(TotMass)~Treatment*Location*Range,random=list(~1|Species),data=dat.f,method="ML")   
anova(fm2.full,fm2.notaxa)

#- it looks like provenance could be removed, but this messes up the error term tests in the ANOVA.
fm3 <-  lme(log(TotMass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN),data=dat.f,method="REML")
anova(fm3)

#- 
#--------------------------------------------------------------------




#--------------------------------------------------------------------
#--- try to fit an alterative model with species as a specific factor. This doens't work because species cannot be crossed by location

fm3 <- lme(TotMass~Treatment*Location*Species,random=list(~1|Sp_RS_EN),data=dat.f) #- "singularity at backsolve" indicates that the model is not estimable.
                                                                                   #- that is, Species cannot be crossed with location because most levels of
                                                                                   #- species are not present in both locations
fm4 <- lme(TotMass~Treatment*Location+Species+Treatment:Species,random=list(~1|Sp_RS_EN),data=dat.f) #- if you drop the Species:Location interaction, the model will fit
anova(fm4) #- note how even the main effect of species is not estimable, because R can't calculate the right error term.
#--------------------------------------------------------------------