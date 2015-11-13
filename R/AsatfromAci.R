#Asat extracted from ACi curves

data<- read.csv("R:/WORKING_DATA/GHS39/GLAHD/Varhammar_A/aci.csv")
asc<- droplevels(subset(data, toFit == "1"))#extract datapoints from ACi curves that we can use

#Combine all values of photo and Tleaf per taxa, treatment and location
asc$Species <- ifelse(asc$Taxa =="ACAM"|asc$Taxa =="BCAM"|asc$Taxa =="CCAM"|asc$Taxa =="DCAM"|asc$Taxa =="ECAM"|asc$Taxa =="FCAM","CAM",
                      ifelse(asc$Taxa == "ATER"|asc$Taxa == "BTER"|asc$Taxa == "CTER"|asc$Taxa == "DTER"|asc$Taxa == "ETER", "TER", asc$Taxa ))
#asc<-asc[-128,]#very low Asat
asat.tm <- summaryBy(Photo~Taxa+Treat+Location+Range,data=asc,FUN=mean,keep.names=T)
colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)


#Asat from ACI curves
ylims=c(20,35)
boxplot(Photo~Treat*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.9,line=2.5)
legend(-0,36,"Temperate", bty="n", cex=1.3)
boxplot(Photo~Treat*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
legend(-0,36,"Tropical", bty="n", cex=1.3)

asc$Location <- factor(asc$Location,levels=c("S","N")) 
asc$Sp_RS_EN <- as.factor(with(asc,paste(Species,Range)))   
asc$Prov_Sp_EN <- as.factor(with(asc,paste(Taxa,Species)))
asc$Sp_Loc_EN <- as.factor(with(asc,paste(Species,Location)))

asc2<-subset(asc, Photo >5) #get rid of really low outlier

fm.Asataci <- lme((Photo)~Treat*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=asc2)
plot(fm.Asataci,resid(.,type="p")~fitted(.) | Treat,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asataci,(Photo)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asataci,(Photo)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asataci, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asataci$residuals[,1])
anova(fm.Asataci)

plot(effect("Range",fm.Asataci))                    #-0.0224 Asat higher in wide
((24.8168)-(21.8118))/(21.8118 )                 # 13.8% higher Asat in wide #same as in spot measures

plot(effect("Treat",fm.Asataci))                 #-<.0001 Asat reduced by warming
((22.85776)-(24.68673))/(24.68673 )              # 7.4% reduction in Asat with warming #same as spot measures

plot(effect("Treat*Location",fm.Asataci),multiline=T) #decreased more in S than N with warming #not same as spotmeasures


plot(effect("Treat*Range*Location",fm.Asataci),multiline=T) #asat spotmeasure% asat Aci curves
((23.65699)-(24.02423))/(24.02423 ) #-6.9% -1.5
((21.06331)-(22.27713))/(22.27713 ) #-7.7% -5.4
((24.17614)-(27.65193))/(27.65193 ) #-3.2% -12.6
((20.83565)-(23.18467))/(23.18467 ) #-6.6% -10.1


#massbasis
ascm<-merge(asc,SLARd, by="Code")#getSLA
ascm$Photom<-with(ascm,Photo*SLA/10000)

ascm2<- subset(ascm, Photom >0.1)
ascm3<-droplevels(subset(ascm2, Species!=17)) #when SMIT is removed there is only a treatment effect.

fm.Asataci <- lme((Photom)~Treat*Range*Location,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=ascm3)
plot(fm.Asataci,resid(.,type="p")~fitted(.) | Treat,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asataci,(Photom)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Asataci,(Photom)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Asataci, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Asataci$residuals[,1])
anova(fm.Asataci)

plot(effect("Treat",fm.Asataci))                #-Asat per g was not reduced by warming#same as spotmeasures
plot(effect("Location",fm.Asataci))                #- Asat per g was not higher in N than S#different to spot measure
plot(effect("Range",fm.Asataci))                    #- Asat per g was not higher in wide#same as spotmeasure
plot(effect("Treat:Location",fm.Asataci),multiline=T)#-0.0258 Asat per g was increased by warming in S, decreased in N#same as spot
((0.4235843)-(0.4619698))/(0.4619698 )     # -9.7 % in N -8.3%
((0.4443223)-(0.4297173))/(0.4297173 )     # 10.2 % in S +3.4%
plot(effect("Treat:Range:Location",fm.Asataci),multiline=T)
((0.4235843)-(0.4619698))/(0.4619698 )     # -9.7 % in N -8.3%
((0.4443223)-(0.4297173))/(0.4297173 )     # 10.2 % in S +3.4%
((0.4235843)-(0.4619698))/(0.4619698 )     # -9.7 % in N -8.3%
((0.4443223)-(0.4297173))/(0.4297173 )     # 10.2 % in S +3.4%

ascm.tm <- summaryBy(Photom~Taxa+Treat+Location+Range,data=ascm,FUN=mean,keep.names=T)

colors <- c(alpha("black",0.6),"red")
windows(5.845,4.135);par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(5,9,3,5),cex.axis=1.2)


#Asat
ylims=c(0.25,0.65)
boxplot(Photom~Treat*Range,data=subset(ascm.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(1,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=F,cex=1,adj=0.5,line=4)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.2,line=2.5)
#legend(-0,36,"Temperate", bty="n", cex=1.3)
axis(side=1,at=c(1.5,3.5),labels=levels(ascm.tm$Range),las=1,cex.axis=1.5)
boxplot(Photom~Treat*Range,data=subset(ascm.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,4),labels=c(0,1),frame.plot=T,las=1)
#legend(-0,36,"Tropical", bty="n", cex=1.3)
axis(side=1,at=c(1.5,3.5),labels=levels(ascm.tm$Range),las=1,cex.axis=1.5)
