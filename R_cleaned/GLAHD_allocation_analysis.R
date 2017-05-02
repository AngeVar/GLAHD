# allometry analysis and plots

#use data2: harvest dataset
#as.Date("2014-11-07") + 40 - use only data from the first 40 days before plants got pot bound


#Leaf mass over Total mass = LMF

fm1LM <- lme(log(Leafmass)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Date < as.Date("2014-12-17")))#, method="ML")
plot(fm1LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LM,log(Leafmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LM,log(Leafmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LM$residuals[,1])
anova(fm1LM)    

#The slope of log(LM)~log(TM) was lower in N and decreased more with warming in the S, -3.4% vs. +0.99% in N
plot(effect("log(Totmass):Treatment:Location",fm1LM), multiline=T)       # 0.0412

#When only looking at 40 day data, there are no significant effects on the slope.

#######################################
#Leaf area - over Leaf mass = SLA
#            
fm2LM <- lme(log(Leafarea)~log(Leafmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Date < as.Date("2014-12-17")))#, method="ML")
plot(fm2LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2LM,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm2LM,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm2LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2LM$residuals[,1])
anova(fm2LM)    
# The slope increased in N (+8.85%) and was reduced (-2.7%) in the south

plot(effect("log(Leafmass):Treatment:Location",fm2LM), multiline=T)       # 0.0024 no difference in S, up in N

#When only looking at 40 day data, there are no significant effects on the slope.

#######################################
#Leaf area - over Total mass = LAR
#            
fm3LM <- lme(log(Leafarea)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Date < as.Date("2014-12-17")))#, method="ML")
plot(fm3LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3LM,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3LM,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3LM$residuals[,1])
anova(fm3LM)    
#+10 in N -6.6% in S
plot(effect("log(Totmass):Treatment:Location",fm3LM), multiline=T)       # 0.0001 up in N down in S

#When only looking at 40 day data, there are no significant effects on the slope.

windows(10,7)
par(mfrow=c(1,2),mar=c(3,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)
plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "S"&Date < as.Date("2014-12-17")), 
       axes=F, ylim= c(3,8.5), xlim=c(-1,4))
lm1<- lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Date < as.Date("2014-12-17")))
abline(lm1, col="red")
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Home"&Date < as.Date("2014-12-17"))), col="black",
       xlim=c(-1,3))

magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)

plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "N"&Date < as.Date("2014-12-17")), 
       ylim= c(3,8.5), xlim=c(-1,4), legend =F,axes=F)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Date < as.Date("2014-12-17"))), col="red")
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Home"&Date < as.Date("2014-12-17"))), col="black")
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
mtext(text="ln(Leaf area)",side=2,outer=T,cex=1.2,adj=0.5,line=3)
mtext(text="ln(Total biomass)",side=1,outer=T,cex=1.2,adj=0.5,line=1)



# #RGR in the same way does not show any significant effects with time and treatment.
# data2$Time<- as.numeric(data2$Date)
# #Total mass - over Time = RGR
# #
# fm3LM <- lme(log(Totmass)~Time*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=data2)#, method="ML")
# plot(fm3LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
# plot(fm3LM,log(Totmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
# plot(fm3LM,log(Totmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
# qqnorm(fm3LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
# hist(fm3LM$residuals[,1])
# anova(fm3LM)
# #
# plot(effect("Time:Treatment:Location",fm3LM), multiline=T)       #





#Leaf mass vs total mass allocation among taxa
ancova.full <- lm(log(Leafmass)~log(Totmass)*Treat*Taxa,data=data2) # most higher order terms not significant
ancova.2 <- lm(log(Leafmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa+log(Totmass):Treat+Taxa:Treat,data=data2) # drop 3-way interaction
ancova.3 <- lm(log(Leafmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa+log(Totmass):Treat,data=data2)
ancova.4 <- lm(log(Leafmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa,data=data2)
anova(ancova.full,ancova.4)
plot(ancova.4) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.full, ancova.4) #full model not significanlty better
#the slope of Leaf mass over total mass differ among taxa (P=0.036)
plot(effect("Treat",ancova.4)) #leaf mass was 4.3% lower in warmed seedlings but warming did not change the relative leaf allocation

#stem mass
ancova.full <- lm(log(Stemmass)~log(Totmass)*Treat*Taxa,data=data2) # most higher order terms not significant
ancova.2 <- lm(log(Stemmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa+log(Totmass):Treat+Taxa:Treat,data=data2) # drop 3-way interaction
ancova.3 <- lm(log(Stemmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa+log(Totmass):Treat,data=data2)
ancova.4 <- lm(log(Stemmass)~log(Totmass)+Taxa+Treat+log(Totmass):Taxa,data=data2)
plot(ancova.full) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.full,ancova.4)
plot(effect("log(Totmass):Taxa",ancova.4)) #Taxa and treatments differed in stem mass but no effect on relative stem allocation

#rootmass
ancova.full <- lm(log(Rootmass)~log(Totmass)*Treat*Taxa,data=data2) # three-way interaction significant
anova(ancova.full, ancova.2)

ancova.2 <- lm(log(Rootmass)~log(Totmass)+Treat+Taxa+log(Totmass):Treat,data=data2)
plot(allEffects(ancova.full))
plot(effect("log(Totmass):Treat",ancova.full), multiline=T) 
#relative allocation to roots (slope of roots over total mass) decreased with warming in a manner that differed among taxa.




