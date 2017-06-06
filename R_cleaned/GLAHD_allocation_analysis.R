# allometry analysis and plots

#use data2: harvest dataset
#as.Date("2014-11-07") + 30 - use only data from the first 30 days before plants got pot bound
summaryBy(Pot~Treatment+Taxa, data=subset(data2, Totmass <=14), FUN =length)

#Leaf mass over Total mass = LMF

fm1LM <- lme(log(Leafmass)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Totmass <=14), method="ML")
plot(fm1LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1LM,log(Leafmass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1LM,log(Leafmass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1LM$residuals[,1])
anova(fm1LM)    

fm1.1LM <- lme(log(Leafmass)~log(Totmass)+Treatment+Location+Treatment:Location,
             random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Totmass <=14))#, method="ML")
anova(fm1.1LM,fm1LM)
plot(effect("Treatment:Location",fm1.1LM), multiline=T)
(exp(1.071134)-exp(1.056171))/exp(1.056171) #+1.5% in S
(exp(1.145197)-exp(1.206949))/exp(1.206949) #-6% in N
# The slopes (i.e. allocation) did not change.
# But there is a significant treatment by location interaction (0.0258) where warming reduced LMF in N (-6%) but 
# not S (+1.5%)


windows(10,7)
par(mfrow=c(1,2),mar=c(3,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)
plotBy(log(Leafmass)~log(Totmass)|Treatment, data=subset(data2,Location == "S"&Totmass <=14), 
       axes=F, ylim= c(-1.2,2.5), xlim=c(-1,3))
lm1<- lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass <=14))
abline(lm1, col="red")
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass <=14)), col="black")

magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)

plotBy(log(Leafmass)~log(Totmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(-1.2,2.5), xlim=c(-1,3), legend =F,axes=F)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red")
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black")
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
mtext(text="ln(Leaf mass)",side=2,outer=T,cex=1.2,adj=0.5,line=3)
mtext(text="ln(Total mass)",side=1,outer=T,cex=1.2,adj=0.5,line=1)


#######################################
#Leaf area - over Leaf mass = SLA
#            
fm2LM <- lme(log(Leafarea)~log(Leafmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Totmass <=14), method="ML")
plot(fm2LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm2LM,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm2LM,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm2LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm2LM$residuals[,1])
anova(fm2LM)    

fm2.1LM <- lme(log(Leafarea)~log(Leafmass)+Treatment+Location+
               Treatment:Location,
             random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Totmass <=14))#, method="ML")
anova(fm2.1LM)
anova(fm2LM,fm2.1LM)
plot(effect("Treatment:Location",fm2.1LM), multiline=T)
(exp(6.472118)-exp(6.332865))/exp(6.332865) #+14.9% in S
(exp(6.370302)-exp(6.522173))/exp(6.522173) #-14.1% in N

# The slopes (i.e. allocation) did not change.
# But there is a significant treatment by location interaction (0.0014) where warming reduced SLA in N (-11.3%) 
# and increased SLA in S (+15%)

windows(10,7)
par(mfrow=c(1,2),mar=c(3,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)
plotBy(log(Leafarea)~log(Leafmass)|Treatment, data=subset(data2,Location == "S"&Totmass<=14), 
       axes=F, ylim= c(3.5,8), xlim=c(-1.5,2.5))
lm1<- lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass<=14))
abline(lm1, col="red")
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass<=14)), col="black")

magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)

plotBy(log(Leafarea)~log(Leafmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(3.5,8), xlim=c(-1.5,2.5), legend =F,axes=F)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red")
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black")
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
mtext(text="ln(Leaf mass)",side=2,outer=T,cex=1.2,adj=0.5,line=3)
mtext(text="ln(Total Leaf mass)",side=1,outer=T,cex=1.2,adj=0.5,line=1)


#######################################
#Leaf area - over Total mass = LAR
#            
fm3LM <- lme(log(Leafarea)~log(Totmass)*Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
             data=subset(data2,Totmass<=14), method="ML")
plot(fm3LM,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm3LM,log(Leafarea)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm3LM,log(Leafarea)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm3LM, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm3LM$residuals[,1])
anova(fm3LM)    

fm3.1LM <- lme(log(Leafarea)~log(Totmass)+Treatment+Location+
                     Treatment:Location,
                   random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
               data=subset(data2,Totmass<=14))#, method="ML")
anova(fm3.1LM)
anova(fm3LM,fm3.1LM)
plot(effect("Treatment:Location",fm3.1LM), multiline=T)
(exp(6.438823)-exp(6.282291))/exp(6.282291) #-+16.9% in S
(exp(6.406582)-exp(6.615842))/exp(6.615842) #-18.8% in N

# The slopes (i.e. allocation) did not change.
# But there is a significant treatment by location interaction (0.0001) where warming reduced SLA in N (-17.2%) 
# and increased SLA in S (+15.4%)

windows(10,7)
par(mfrow=c(1,2),mar=c(3,0,1.5,0),oma=c(6,7,6,7),cex.axis=1.2)
plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "S"&Totmass<=14), 
       axes=F, ylim= c(3.5,8), xlim=c(-0.5,3))
lm1<- lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass<=14))
abline(lm1, col="red")
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass<=14)), col="black",
       xlim=c(-1,3))

magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)

plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(3.5,8), xlim=c(-0.5,3), legend =F,axes=F)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red")
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black")
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
mtext(text="ln(Leaf area)",side=2,outer=T,cex=1.2,adj=0.5,line=3)
mtext(text="ln(Total mass)",side=1,outer=T,cex=1.2,adj=0.5,line=1)

###
#plot them all together
windows(10,10)

par(mfrow=c(3,2),mar=c(3,0,1.5,0),oma=c(6,10,6,7),cex.axis=1.2)

#LAR
plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "S"&Totmass<=14), 
       axes=F, ylim= c(3.5,8.5), xlim=c(-0.5,3),pch=19,cex=2, legend =F)
legend(-0.5,8, legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),col=c("red","black"),pch=19,bty="n", cex=1.2, pt.cex =2)
lm1<- lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass<=14))
abline(lm1, col="red",lwd=2)
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass<=14)), col="black",
       xlim=c(-1,3),lwd=2)
magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
mtext(text="Temperate", side=3, line=0.5, cex=1.2)
legend("topright","a", bty="n", cex=1.5)

plotBy(log(Leafarea)~log(Totmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(3.5,8.5), xlim=c(-0.5,3), legend =F,axes=F,pch=19,cex=2)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red",lwd=2)
abline(lm(log(Leafarea)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black",lwd=2)
mtext(text="Tropical", side=3, line=0.5, cex=1.2)
mtext(text="ln(Leaf area (cm2))",side=2,outer=T,cex=1,adj=0.92,line=2.5)
mtext(text="ln(Total mass (g))",side=1,outer=T,cex=1,line=-41.5)
mtext(text="LAR",side=2,outer=T,cex=1.5,adj=0.88,line=5)
legend("topright","b", bty="n", cex=1.5)

#SLA

plotBy(log(Leafarea)~log(Leafmass)|Treatment, data=subset(data2,Location == "S"&Totmass<=14), 
       axes=F, ylim= c(3.5,8.5), xlim=c(-1.5,2.5), legend=F,pch=19,cex=2)
lm1<- lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass<=14))
abline(lm1, col="red",lwd=2)
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass<=14)), col="black",lwd=2)
magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","c", bty="n", cex=1.5)

plotBy(log(Leafarea)~log(Leafmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(3.5,8.5), xlim=c(-1.5,2.5), legend =F,axes=F,pch=19,cex=2)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red",lwd=2)
abline(lm(log(Leafarea)~log(Leafmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black",lwd=2)
mtext(text="ln(Leaf area (cm2))",side=2,outer=T,cex=1,adj=0.5,line=2.5)
mtext(text="ln(Total Leaf mass (g))",side=1,outer=T,cex=1,adj=0.5,line=-21)
mtext(text="SLA",side=2,outer=T,cex=1.5,adj=0.5,line=5)
legend("topright","d", bty="n", cex=1.5)

#LMF
plotBy(log(Leafmass)~log(Totmass)|Treatment, data=subset(data2,Location == "S"&Totmass<=14), 
       axes=F, ylim= c(-1.2,3), xlim=c(-1,3), legend=F,pch=19,cex=2)
lm1<- lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Warmed"&Totmass<=14))
abline(lm1, col="red",lwd=2)
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "S"& Treatment == "Home"&Totmass<=14)), col="black",lwd=2)
magaxis(c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
legend("topright","e", bty="n", cex=1.5)

plotBy(log(Leafmass)~log(Totmass)|Treatment, data=subset(data2,Location == "N"&Totmass<=14), 
       ylim= c(-1.2,3), xlim=c(-1,3), legend =F,axes=F,pch=19,cex=2)
magaxis(c(1,4),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2)
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Warmed"&Totmass<=14)), col="red",lwd=2)
abline(lm(log(Leafmass)~log(Totmass), data=subset(data2,Location == "N"& Treatment == "Home"&Totmass<=14)), col="black",lwd=2)
mtext(text="ln(Leaf mass (g))",side=2,outer=T,cex=1,adj=0.1,line=2.5)
mtext(text="ln(Total mass (g))",side=1,outer=T,cex=1,adj=0.5,line=-1)
mtext(text="LMF",side=2,outer=T,cex=1.5,adj=0.15,line=5)
legend("topright","f", bty="n", cex=1.5)

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




