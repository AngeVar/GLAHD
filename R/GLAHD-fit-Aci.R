#- to install the newest version of plantecophys
# library(devtools)
# install_bitbucket("remkoduursma/plantecophys")



#- load libraries from script
source("C:/Repos/GLAHD/R/loadLibraries.R")
library(nlme)
library(effects)

#- read data from csv files. These are rawish- they are the files right off the machine, but manually processed to remove remarks
#- and code which points to exclude from the A:Ci curve fitting.
aci1 <- read.csv(file="C:/Repos/GLAHD/Data/GasEx/Aci/Aci_day1_compiled.csv")
aci2 <- read.csv(file="C:/Repos/GLAHD/Data/GasEx/Aci/Aci_day2_compiled.csv")
aci3 <- read.csv(file="C:/Repos/GLAHD/Data/GasEx/Aci/Aci_day3_compiled.csv")

aci <- rbind(aci1,aci2,aci3)

#- get the first bit of the code (the taxa)
aci$Taxa <- unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=1,to=nrow(aci)*2,by=2)]
aci$Taxa <- factor(aci$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                   "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
aci$Pot <- as.numeric(unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=2,to=nrow(aci)*2,by=2)])
aci$Treat <- ifelse(aci$Pot < 20, "Home","Warmed")

#- sort by code
aci <- aci[with(aci,order(Code)),]

#- have a look
palette(rainbow(n=length(levels(aci$Code))))
plotBy(Photo~Ci|Code,data=subset(aci,toFit==1),legend=F,type="b")

#- extract the data to fit (exclude lines that I've manually decided to remove from the curve fitting)
aci.fit <- subset(aci,toFit==1)




#--- attempt to use parallel processing to speed up the A:Ci curve fitting process.
p.flag <- F # use parallel processing?
if(p.flag==T){
  library(doParallel)
  cl <- makeCluster(2)
  registerDoParallel(cl)
}
#example bootstrapping from doParallel example pdf. Using 2 cores sped up the bootstrapping from 30.93 to 21.57 seconds.
#- that's a speedup of about 23%.
# x <- iris[which(iris[,5] != "setosa"), c(1,5)]
# trials <- 10000
# ptime <- system.time({
#    r <- foreach(icount(trials), .combine=cbind) %do% {
#      ind <- sample(100, 100, replace=TRUE)
#      result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
#      coefficients(result1)
#      }
#   })[3]
# ptime
if(p.flag==T){
  starttime <- Sys.time()
  fits.p <- list()
  aci.fit.list <- split(aci.fit,aci.fit$Code)
  crap <- foreach(i=1:length(aci.fit.list)) %dopar% 
    plantecophys::fitaci(aci.fit.list[[i]],PPFD="PARi",Tcorrect=F,citransition=450)
  
  Sys.time()-starttime # took 55 seconds to fit all curves, rather than 1.6 mintues on 1 core.
  fits.params <- as.data.frame(do.call(rbind,lapply(crap,FUN=coef)))
  fits.params$Code <- levels(aci.fit$Code)
}#-----


#- fit the aci curves for each plant
if(p.flag==F){
  starttime <- Sys.time()
  fits <- fitacis(aci.fit,group="Code",PPFD="PARi",Tcorrect=F,citransition=450)
  fits.params <- coef(fits) #extract the Vcmax and Jmax parameters
  Sys.time()-starttime
  
}

#- assign other factor variables based on the "Code".
fits.params$Taxa <- unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=1,to=nrow(fits.params)*2,by=2)]
fits.params$Taxa <- factor(fits.params$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                  "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
fits.params$Pot <- as.numeric(unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=2,to=nrow(fits.params)*2,by=2)])
fits.params$Treat <- ifelse(fits.params$Pot < 20, "Home","Warmed")



# 
# #---------------------------------------------------------------------
# #make a big pdf with all of the fits for each curve.
# pdf(file="C:/Repos/GLAHD/Output/Aci fits.pdf")
# 
# for (i in 1:length(fits)){
#   par(xpd=F)
#   plot(fits[[i]],xlim=c(0,2000),ylim=c(-2,50))
#   title(paste("ch.campaign =", fits.params[i,1]),sep="")
#   par(xpd=T)
#   legend("topleft",legend=paste("Vcmax=",round(fits.params[i,2],1),",","Jmax=",round(fits.params[i,3],1)),
#          bty="n")
# }
# dev.off()
# #---------------------------------------------------------------------
# 







#----------------------------------------------------------------------------------
#- process Aci fits for statistical analysis

#- get the growth data, mostly just for the treatment codes
growth <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
growth2 <- summaryBy(d2h+TotMass+leafArea~Species+Treatment+Location+Taxa+Code+Range,keep.names=T,data=subset(growth,Date >= as.Date("2014-12-8") & Date <=as.Date("2014-12-20")))

#- merge size totalmass and leafarea data into dataframe with aci values
acifit <- merge(fits.params,growth2,by=c("Code","Taxa"))
acifit$JtoV <- with(acifit,Jmax/Vcmax)
acifit$Location <- factor(acifit$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
acifit$Sp_RS_EN <- as.factor(with(acifit,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
acifit$Prov_Sp_EN <- as.factor(with(acifit,paste(Taxa,Species)))


#- does Vcmax or Jmax change with plant size? note that there will be a huge N vs. S effect in Vcmax
pairs(acifit[,c("Vcmax","Jmax","TotMass","d2h","leafArea")])


#Get SLA
Rdark<-getRdark()
SLARd<- Rdark[,c("Code","SLA")]
acifits<-merge(acifit,SLARd, by="Code")

acifits$Jmaxm<-with(acifits,Jmax*SLA/10000)
acifits$Vcmaxm<-with(acifits,Vcmax*SLA/10000)
acifits$JtoVm<-with(acifits,Jmaxm/Vcmaxm)






#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#- make a nice plot of Vcmax and Jmax in the N vs. S for a manuscript


Asat <- getAsat()

#- average across taxa
acifits$Species <- factor(acifits$Species)
acifits.tm <- summaryBy(Vcmax+Jmax~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
asat.tm <- summaryBy(Photo~Taxa+Treatment+Location+Range,data=Asat,FUN=mean,keep.names=T)

#- make plots
colors <- c("blue","red")
windows(20,20);par(mfrow=c(3,2),mar=c(2,0,1,0),oma=c(5,9,3,5),cex.axis=1.2)

#Vcmax
ylims=c(60,220)
boxplot(Vcmax~Treatment*Range,data=subset(acifits.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","a",bty="n",cex=1.5,inset=-0.05)
title(main="North",line=0.2,cex.main=2,xpd=NA)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(V["c,max"]),side=2,outer=T,cex=2,adj=0.9,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.92,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
legend("bottomleft",c("Home","Warmed"),fill=colors,cex=1.3)
boxplot(Vcmax~Treatment*Range,data=subset(acifits.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","b",bty="n",cex=1.5,inset=-0.05)
title(main="South",line=0.2,cex.main=2,xpd=NA)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

#Jmax
ylims=c(100,210)
boxplot(Jmax~Treatment*Range,data=subset(acifits.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","c",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(J["max"]),side=2,outer=T,cex=2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.52,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
boxplot(Jmax~Treatment*Range,data=subset(acifits.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","d",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

#Asat
ylims=c(20,35)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","e",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=T,cex=2,adj=0.15,line=5)
mtext(text=expression("("*mu*mol~m^-2~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
boxplot(Photo~Treatment*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","f",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

mtext(text=expression(Range~size),side=1,outer=T,cex=2,line=3)
dev.copy2pdf(file="C:/Repos/GLAHD/Output/Photo_figure_Vcmax_Jmax_Asat.pdf")
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------



#- fit and interpret Vcmax
fm.vcmax <- lme(log(Vcmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.vcmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmax,log(Vcmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmax,log(Vcmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmax$residuals[,1])
anova(fm.vcmax)    

plot(effect("Treatment:Location",fm.vcmax))      #- warming reduced Vcmax in the south but not in the north. No range interactions
effect("Treatment:Location",fm.vcmax)
(exp(4.61284)-exp(4.4414))/(exp(4.61284))    # 15.7% reduction in Vcmax in S
(exp(5.187)-exp(5.1623))/(exp(5.187)) #  2% reduction in Vcmax in N
plot(effect("Range",fm.vcmax))                   #- narrow species have, on average, lower Vcmax
effect("Range",fm.vcmax)
(exp(4.9089)-exp(4.78961))/(exp(4.9089))         #- Vcmax is 11% lower in narrow ranged species



#- fit and interpret Jmax
fm.jmax <- lme(log(Jmax)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jmax,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmax,log(Jmax)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmax,log(Jmax)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmax, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmax$residuals[,1])
anova(fm.jmax)    

plot(effect("Treatment:Location",fm.jmax))      #- warming reduced Jmax to a larger degree in the south than the north. No range interactions
effect("Treatment:Location",fm.jmax)
(exp(5.153)-exp(4.923))/(exp(5.153))    # 20.5% reduction in Jmax in S
(exp(5.103)-exp(5.004))/(exp(5.103))    #  9% reduction in Jmax in N




#- fit and interpret Jmax/Vcmax
fm.jtov <- lme(log(JtoV)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jtov,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtov,log(JtoV)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtov,log(JtoV)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtov, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtov$residuals[,1])
anova(fm.jtov)    

plot(effect("Location",fm.jtov))                #- Jmax/Vcmax much higher in S than north. No range interactions
plot(effect("Treatment",fm.jtov))               #- warming reduced Jmax/Vcmax overall
effect("Treatment",fm.jtov)
(exp(0.2124)-exp(0.1459))/(exp(0.2124))      # 6% reduction in Jmax/Vcmax with warming
effect("Location",fm.jtov)
(exp(0.5103)-exp(-0.1208))/(exp(0.5103))      # 47% lower Jmax/Vcmax in N compared to S

#----------------------------------------------------------------------------------
#On a mass basis
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#- make a nice plot of Vcmax and Jmax in the N vs. S for a manuscript
Asat<- getAsat()
Rdark<- getRdark()
asatshort<-Asat[,c("Code","Species","Taxa","Treatment", "Location", "Range","Photo")]
rdarkshort<-Rdark[,c("Code","Species","Taxa","Treatment", "Location", "Range","SLA")]
gasex<-merge(asatshort,rdarkshort,by=c("Code", "Species","Taxa","Treatment","Location","Range"))


gasex$photomass<-with(gasex, Photo*SLA/10000)

#- average across taxa
acifits$Species <- factor(acifits$Species)
acifits.tm <- summaryBy(Vcmaxm+Jmaxm~Taxa+Treatment+Location+Range,data=acifits,FUN=mean,keep.names=T)
asat.tm <- summaryBy(photomass~Taxa+Treatment+Location+Range,data=gasex,FUN=mean,keep.names=T)

#- make plots
colors <- c("blue","red")
windows(20,20);par(mfrow=c(3,2),mar=c(2,0,1,0),oma=c(5,9,3,5),cex.axis=1.2)

#Vcmax
ylims=c(0,5)
boxplot(Vcmaxm~Treatment*Range,data=subset(acifits.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","a",bty="n",cex=1.5,inset=-0.05)
title(main="North",line=0.2,cex.main=2,xpd=NA)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(V["c,max"]),side=2,outer=T,cex=2,adj=0.9,line=5)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.92,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
boxplot(Vcmaxm~Treatment*Range,data=subset(acifits.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","b",bty="n",cex=1.5,inset=-0.05)
legend("topright",c("Home","Warmed"),fill=colors,cex=1.3)
title(main="South",line=0.2,cex.main=2,xpd=NA)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

#Jmax
ylims=c(0,5)
boxplot(Jmaxm~Treatment*Range,data=subset(acifits.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","c",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(J["max"]),side=2,outer=T,cex=2,adj=0.5,line=5)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.52,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
boxplot(Jmaxm~Treatment*Range,data=subset(acifits.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","d",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

#Asat
ylims=c(0,1)
boxplot(photomass~Treatment*Range,data=subset(asat.tm,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","e",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=T,cex=2,adj=0.15,line=5)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,adj=0.1,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)
boxplot(photomass~Treatment*Range,data=subset(asat.tm,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","f",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(acifits.tm$Range),las=1,cex.axis=1.5)

mtext(text=expression(Range~size),side=1,outer=T,cex=2,line=3)
dev.copy2pdf(file="C:/Repos/GLAHD/Output/Photo_figure_Vcmax_Jmax_Asat_mass.pdf")


#- fit and interpret Vcmax
fm.vcmaxm <- lme((Vcmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.vcmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.vcmaxm,log(Vcmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.vcmaxm,log(Vcmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.vcmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.vcmaxm$residuals[,1])
anova(fm.vcmaxm)    

plot(effect("Treatment",fm.vcmaxm))               #Warming increased Vcmax
plot(effect("Location",fm.vcmaxm))                #Tropical taxa had higher Vcmax than temperate
plot(effect("Treatment:Location",fm.vcmaxm))      #Warming increased Vcmax in the south but not in the north. No range interactions
plot(effect("Treatment:Location:Range",fm.vcmaxm))

(exp(4.399115)-exp(4.483795))/(exp(4.483795))    #  8.1 % reduction in Vcmax in North Wide
(exp(4.268271)-exp(4.114689))/(exp(4.114689))    #  16.6 % increase in Vcmax in North Narrow
(exp(3.688584)-exp(3.274452))/(exp(3.274452))    #  51.3 % Increase in Vcmax in South Wide
(exp(3.719271)-exp(3.632618))/(exp(3.632618))    #  9 % increase in Vcmax in South Narrow

#- fit and interpret Jmax
fm.jmaxm <- lme(log(Jmaxm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jmaxm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jmaxm,log(Jmaxm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jmaxm,log(Jmaxm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jmaxm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jmaxm$residuals[,1])
anova(fm.jmaxm)    

plot(effect("Treatment",fm.jmaxm))               #Warming did not increase Jmax
plot(effect("Location",fm.jmaxm))                #Tropical taxa had higher Jmax than temperate
plot(effect("Treatment:Location",fm.jmaxm))      #Warming increased Jmax in the south, decreased in the north
plot(effect("Location:Range",fm.jmaxm))          #N wide higher than S wide. N narrow lower than S narrow.
plot(effect("Treatment:Location:Range",fm.jmaxm))

(exp(4.257304)-exp(4.404200))/(exp(4.404200))    #  13.7 % reduction in Jmax in North Wide
(exp(4.087280)-exp(4.014411))/(exp(4.014411))    #  7.6 % increase in Jmax in North Narrow
(exp(4.139475)-exp(3.779167))/(exp(3.779167))    #  43.4 % Increase in Jmax in South Wide
(exp(4.251939)-exp(4.235751))/(exp(4.235751))    #  1.6 % increase in Jmax in South Narrow




#- fit and interpret Jmax/Vcmax
fm.jtovm <- lme(log(JtoVm)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=acifits)
plot(fm.jtovm,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.jtovm,log(JtoVm)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.jtovm,log(JtoVm)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.jtovm, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.jtovm$residuals[,1])
anova(fm.jtovm)    

plot(effect("Location",fm.jtovm))                #- Jmax/Vcmax much higher in S than north. No range interactions
plot(effect("Treatment",fm.jtovm))               #- warming reduced Jmax/Vcmax overall
effect("Treatment",fm.jtovm)
(exp(0.1495140)-exp(0.2136865))/(exp(0.2136865))      # 6.2% reduction in Jmax/Vcmax with warming
effect("Location",fm.jtovm)
(exp(0.5109005)-exp(-0.1213559))/(exp(0.5109005))      # 46.9% lower Jmax/Vcmax in N compared to S






#----------------------------------------------------------------------------------
#-- plot Vcmax and Jmax for each taxa
windows(16,14);par(mfrow=c(3,1),mar=c(6,5,2,3),cex.lab=1.3,cex.axis=1.3,las=1)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")

#plot Vcmax
ylims=c(10,35)
boxplot(Vcmax~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),frame.plot=T,las=1)
title(ylab=expression(V["c,max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(acifits$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)

#plot Jmax
boxplot(Jmax~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),frame.plot=T,las=1)
title(ylab=expression(J["max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(acifits$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)

#plot Jmax:Vcmax
boxplot(JtoV~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),frame.plot=T,las=1)
title(ylab=expression(J["max"]~"/"~V["cmax"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(acifits$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)
dev.copy2pdf(file="C:/Repos/GLAHD/Output/Aci_results_Vcmax_Jmax_atMeasuredLeafR.pdf")
#----------------------------------------------------------------------------------






#----------------------------------------------------------------------------------
#- fit the models for Asat and Rdark
Asat <- getAsat()
fm.Asat <- lme(Photo~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Asat)
anova(fm.Asat)
plot(effect("Treatment",fm.Asat))                #- Asat reduced by warming
plot(effect("Range",fm.Asat))               #- Asat higher in wide range species
effect("Treatment",fm.Asat)
(28.35-26.73)/28.35      # 5% reduction in Asat with warming


Rdark <- getRdark()
fm.Rmass <- lme(log(Rmass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=Rdark)
plot(fm.Rmass,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Rmass,log(Rmass)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.Rmass,log(Rmass)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.Rmass, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.Rmass$residuals[,1])
anova(fm.Rmass)
plot(effect("Treatment:Location:Range",fm.Rmass))                #- 3-way interaction for Rmass. Narrow from south don't acclimate, the rest do
effect("Treatment:Location:Range",fm.Rmass)
(exp(2.426073)-exp(2.174630))/(exp(2.426073))      # 22% reduction in Rmass with warming of wide species in the north
(exp(2.361253)-exp(2.144013))/(exp(2.361253))      # 19.5% reduction in Rmass with warming of narrow species in the north
(exp(2.318532)-exp(1.974925))/(exp(2.318532))      # 29% reduction in Rmass with warming of wide species in the south
(exp(2.315649)-exp(2.253351))/(exp(2.315649))      # 6% reduction in Rmass with warming of narrow species in the south


plot(effect("Treatment:Range",fm.Rmass))                         #- wide range species exhibit more substantial acclimation than narrows


#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------
# #- create ANOVA table of results, embed that table into a word document. 
# #- the table needs cleaning in word to embed the degrees of freedom and label the columns, etc.
# vcmax.o <- anova(fm.vcmax)
# 
# extract.lme <- function(model){
#   mod.o <- anova(model)
#   
#   ndf <- mod.o[2:nrow(mod.o),1]
#   ddf <- mod.o[2:nrow(mod.o),2]
#   df <- (paste(ndf,ddf,sep=", "))
#   Fvalue <-round(mod.o[2:nrow(mod.o),3],1)
#   Pvalue <-round(mod.o[2:nrow(mod.o),4],3)
#   Pvalue[which(Pvalue==0)] <- "<0.001"
#   return(data.frame("df"=df, "F" = Fvalue,"P" = Pvalue))
# }
# 
# 
# #- make a big table
# gxtable <- do.call(cbind,lapply(list(fm.vcmax,fm.jmax,fm.jtov,fm.Asat,fm.Rmass),FUN=extract.lme))
# row.names(gxtable) <- row.names(anova(fm.vcmax))[2:8]
# 
# library(R2wd)
# wdGet()
# wdTable(gxtable,autoformat=2)
#----------------------------------------------------------------------------------

################################################################################################################

#Is SLA different? warming increased SLA but only in the South
Asat<- getAsat()
Rdark<- getRdark()
asatshort<-Asat[,c("Code","Species","Taxa","Treatment", "Location", "Range","Photo")]
rdarkshort<-Rdark[,c("Code","Species","Taxa","Treatment", "Location", "Range","SLA")]
gasex<-merge(asatshort,rdarkshort,by=c("Code", "Species","Taxa","Treatment","Location","Range"))
#gasex$SLA<- with(gasex,(leafArea/10000)/(leafDW/1000))#m2 per g. This was wrong, the gx-leaf SLA values were already in Rdark (JED)
 
gasex$Location <- factor(gasex$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gasex$Sp_RS_EN <- as.factor(with(gasex,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gasex$Prov_Sp_EN <- as.factor(with(gasex,paste(Taxa,Species)))
gasex$Sp_Loc_EN <- as.factor(with(gasex,paste(Species,Location)))

fm.SLA <- lme((SLA)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=gasex)
plot(fm.SLA,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.SLA,(SLA)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm.SLA,(SLA)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm.SLA, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm.SLA$residuals[,1])

anova(fm.SLA) #SLA is significantly different 
plot(effect("Treatment",fm.SLA))                     #- SLA increased by warming P<.0001
plot(effect("Treatment:Location",fm.SLA))            #- SLA only increased with warming in the South P<.0001
plot(effect("Treatment:Location:Range",fm.SLA))      #particularly in south wide P=0.0049
effect("Treatment:Location:Range",fm.SLA)

(185.0789-193.7101)/193.7101     # 4.4 % reduction in SLA with warming of wide species in the north
(197.9415-192.6063)/192.6063     # 2.7 % increase in SLA with warming of narrow species in the north
(183.2336-165.9987)/165.9987     # 10.3% increase in SLA with warming of wide species in the south
(215.4669-171.5054)/171.5054     # 25.6% increase in SLA with warming of narrow species in the south

# (exp(-0.7940645)-exp(-0.7644813))/(exp(-0.7644813))      # 2.9% reduction in SLA with warming of wide species in the north
# (exp(-0.8559911)-exp(-1.0351236))/(exp(-1.0351236))      # 19.6% increase in SLA with warming of narrow species in the north
# (exp(-0.8233543)-exp(-1.4017938))/(exp(-1.4017938))      # 78.3% increase in SLA with warming of wide species in the south
# (exp(-0.6261996)-exp(-0.8563303))/(exp(-0.8563303))      # 25.9% increase in SLA with warming of narrow species in the south


#-- the following probably need changing (JED)
(exp(-0.7940645)-exp(-0.7644813))/(exp(-0.7644813))      # 2.9% reduction in SLA with warming of wide species in the north
(exp(-0.8559911)-exp(-1.0351236))/(exp(-1.0351236))      # 19.6% increase in SLA with warming of narrow species in the north
(exp(-0.8233543)-exp(-1.4017938))/(exp(-1.4017938))      # 78.3% increase in SLA with warming of wide species in the south
(exp(-0.6261996)-exp(-0.8563303))/(exp(-0.8563303))      # 25.9% increase in SLA with warming of narrow species in the south


windows(20,15);par(mfrow=c(1,2),mar=c(2,0,1,0),oma=c(5,9,3,5),cex.axis=1.2)

ylims=c(0,250)
boxplot(SLA~Treatment*Range,data=subset(gasex,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","Tropical",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(SLA),side=2,outer=T,cex=2,line=5)
mtext(text=expression("("*m^2~g^-1*")"),side=2,outer=T,cex=1,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(gasex$Range),las=1,cex.axis=1.5)
boxplot(SLA~Treatment*Range,data=subset(gasex,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","Temperate",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(gasex$Range),las=1,cex.axis=1.5)
mtext(text=expression(Range~size),side=1,outer=T,cex=2,line=3)

################################################################################################################
#Is Asat on mass basis different? Yes - three way interaction 
gasex$photomass<-with(gasex, Photo*(leafArea/leafDW)/10000*1000) # convert to umol CO2 g-1 s-1
asat.mass <- summaryBy(photomass~Taxa+Treatment+Location+Range,data=gasex,FUN=mean,keep.names=T)

windows(20,15);par(mfrow=c(1,2),mar=c(2,0,1,0),oma=c(5,9,3,5),cex.axis=1.2)

ylims=c(0,25)
boxplot(photomass~Treatment*Range,data=subset(gasex,Location=="N"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","Tropical",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(1,0,0),frame.plot=T,las=1)
mtext(text=expression(A["sat"]),side=2,outer=T,cex=2,line=5)
mtext(text=expression("("*mu*mol~g^-1~s^-1*")"),side=2,outer=T,cex=1,line=3)
axis(side=1,at=c(1.5,3.5),labels=levels(gasex$Range),las=1,cex.axis=1.5)
boxplot(photomass~Treatment*Range,data=subset(gasex,Location=="S"),ylim=ylims,
        axes=F,las=2,col=colors)
legend("topleft","Temperate",bty="n",cex=1.5,inset=-0.05)
magaxis(c(2,3,4),labels=c(0,0,1),frame.plot=T,las=1)
axis(side=1,at=c(1.5,3.5),labels=levels(gasex$Range),las=1,cex.axis=1.5)
mtext(text=expression(Range~size),side=1,outer=T,cex=2,line=3)



fm.Asatm <- lme(log(photomass)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=gasex)

plot(fm.Asatm,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.Asatm,log(photomass)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm.Asatm,log(photomass)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm.Asatm, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm.Asatm$residuals[,1])

anova(fm.Asatm)
plot(effect("Treatment",fm.Asatm))                #- Asatmass increased by warming
plot(effect("Location",fm.Asatm))                 #- Asatmass is higher in North than South
plot(effect("Treatment:Location",fm.Asatm))       #- Asatmass only increased with warming in the South
plot(effect("Treatment:Location:Range",fm.Asatm)) #particularly in south wide
effect("Treatment:Location:Range",fm.Asatm)
(exp(2.518584)-exp(2.625855))/(exp(2.625855))      # 10.2% reduction in Asatmass with warming of wide species in the north
(exp(2.340503)-exp(2.234243))/(exp(2.234243))      # 11.2% increase in Asatmass with warming of narrow species in the north
(exp(2.501859)-exp(1.962196))/(exp(1.962196))      # 71.5% increase in Asatmass with warming of wide species in the south
(exp(2.567042)-exp(2.412422))/(exp(2.412422))      # 16.7% increase in Asatmass with warming of narrow species in the south


