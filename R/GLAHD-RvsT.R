#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plots foliar respiraiton vs. Temperature curves for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("C:/Repos/GLAHD/R/loadLibraries.R")



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- read in the data, do a few conversions
dat <- read.csv("C:/Repos/GLAHD/Data/GasEx/RvsT/GLADH RvsT gasex processeed.csv",stringsAsFactors=T)
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Taxa <- as.factor(sapply(1:nrow(dat),function(i)strsplit(as.character(dat$Code[i]),split="-")[[1]][1])) # get the "Taxa" factor out of Code. Note this could have been a loop
dat$R_area <- dat$Photo*-1

#- read in the leaf mass data, calculate R per unit leaf mass
leafmass <- read.csv("C:/Repos/GLAHD/Data/GasEx/RvsT/GLAHD_leaf_mass_RvT.csv")
rt <- merge(dat,leafmass,by="Code")
rt$R_mass <- rt$R_area/10000*rt$Area/rt$mass*1000 #note mass was in g


rt <- subset(rt,Taxa!="SMIT" & Taxa !="ACAM")
rt$Taxa <- factor(rt$Taxa,levels=c("BTER","BOT","CTER","BRA"))
rt$lnRmass <- log(rt$R_mass)
rt <- subset(rt,Code!="CTER-4" & Code!="CTER-22") #CTER4 looks weird... peaks much much lower than others. CTER-22 also looks bad
rt$Code <- factor(rt$Code)

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
## quick plot- have a look

#- break data into bins of leaf temperatures, for averaging and plotting
rt$Tleaf_bin <- cut(rt$Tleaf,breaks=seq(from=10,to=65,length=60))
rt$Tleaf_bin_mid <- sapply(strsplit(gsub("^\\W|\\W$", "", rt$Tleaf_bin), ","), function(x)sum(as.numeric(x))/2) 


#- average across taxa, treatment, Tbin levels
rt.m <- summaryBy(R_mass+R_area~Species+Taxa+Treatment+Tleaf_bin_mid,data=rt,FUN=c(mean,standard.error))
rt.m.l <- split(rt.m,rt.m$Taxa)

#- plot respiration per unit leaf mass
windows(20,20);par(mfrow=c(2,2),mar=c(3,2,1,1),oma=c(4,4,1,1))
for (i in 1:length(rt.m.l)){
  toplot <- subset(rt.m.l[[i]],Tleaf_bin_mid < 40)
  plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=toplot,ylim=c(0,40),pch=15,
         panel.first=adderrorbars(x=toplot$Tleaf_bin_mid,y=toplot$R_mass.mean,SE=toplot$R_mass.standard.error,direction="updown",col="black"))
  title(main=toplot$Taxa[i],line=0.2)
  
}
title(xlab=expression(Leaf~temperature~(degree*C)),ylab=expression(R[mass]~(nmol~g^-1~s^-1)),outer=T,cex.lab=1.4,line=1)
dev.copy2pdf(file="C:/Repos/GLAHD/Output/Rmass_vs_T_JED.pdf")

#- plot respiration per unit leaf area
windows(20,20);par(mfrow=c(2,2),mar=c(3,2,1,1),oma=c(4,4,1,1))
for (i in 1:length(rt.m.l)){
  toplot <- subset(rt.m.l[[i]],Tleaf_bin_mid < 40)
  plotBy(R_area.mean~Tleaf_bin_mid|Treatment,data=toplot,ylim=c(0,4),pch=15,
         panel.first=adderrorbars(x=toplot$Tleaf_bin_mid,y=toplot$R_area.mean,SE=toplot$R_area.standard.error,direction="updown",col="black"))
  title(main=toplot$Taxa[i],line=0.2)
  
}
title(xlab=expression(Leaf~temperature~(degree*C)),ylab=expression(R[area]~(mu*mol~m^-2~s^-1)),outer=T,cex.lab=1.4,line=1)

#-----------------------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

#loop over each curve, fit each model, and plot the data and model predictions

q10mods <- list() #simple q10
arrmods <- list() #simple arrhenius equation
mtmods <- list()  #modified Q10 to allow Q10 to decline with measurement temperature. Mark's preferred model?
polymods <- list() #polynomial equation popularized by Owen Atkin. Basically an arrhenius where the activation energy can decline with temperature


pdf(file="C:/Repos/GLAHD/Output/RvsT_curves.pdf")
#par(mfrow=c(6,2),mar=c(0,0,0,0),oma=c(6,6,2,4))
par(cex=1,cex.axis=1,cex.lab=1)
xvals <- seq(13,45,length=101)
rt.l <- split(rt,rt$Code)
for (i in 1:length(rt.l)){

  dat <- subset(rt.l[[i]],Tleaf<40)
  print(i)
  #fit the models
  q10mods[[i]] <- nls(R_mass~ Rref*Q10^((Tleaf-15)/10),data=dat,start=list(Rref=10,Q10=2))
  arrmods[[i]] <- nls(R_mass~Rref*exp((Ea/(0.008314*(273.15+15))*(1-((273.15+15))/(273.15+Tleaf)))),data=dat,start=list(Rref=10,Ea=45))
  mtmods[[i]] <- nls(R_mass~Rref*((x-y*(Tleaf-15)/2))^((Tleaf-15)/10),data=dat,start=list(x=3,y=0.03,Rref=10))
  polymods[[i]] <- lm(lnRmass~Tleaf+I(Tleaf^2),data=dat)
  
  #plot data and predictions
  plot(R_mass~Tleaf,data=dat,type="p",pch=16,col="black",ylim=c(0,50),xlim=c(12,47),axes=F)
  lines(xvals,predict(q10mods[[i]],data.frame(Tleaf= xvals)),lwd=3,col="black")
  lines(xvals,predict(arrmods[[i]],data.frame(Tleaf= xvals)),lwd=3,col="grey")
  lines(xvals,predict(mtmods[[i]],data.frame(Tleaf= xvals)),lwd=3,col="red")
  pred.poly <- unname(predict(polymods[[i]],data.frame(Tleaf= xvals)))
  lines(xvals,exp(pred.poly),lwd=3,col="forestgreen")
  magaxis(side=1:4,labels=c(1,1,0,1))
  legend("topleft",legend=c("Q10","Arr","MT-Q10","Poly"),lty=1,lwd=3,cex=1.2,col=c("black","grey","red","forestgreen"))
  legend("top",legend=paste(dat$Code[1]),bty="n",cex=1.5)
}
title(xlab="leaf temperature (deg C)",ylab=expression(R~(n~mol~g^-1~s^-1)),outer=T,cex.lab=1.4)
dev.off()

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------






#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
## evaluate model outputs

#get coefs
q10.params <- as.data.frame(do.call(rbind,lapply(q10mods,coef)))
arr.params <- as.data.frame(do.call(rbind,lapply(arrmods,coef)))
mt.params <- as.data.frame(do.call(rbind,lapply(mtmods,coef)))
poly.params <- as.data.frame(do.call(rbind,lapply(polymods,coef)))

#rename parameter columns
names(q10.params) <- paste("Q10.",names(q10.params),sep="")
names(arr.params) <- paste("Arr.",names(arr.params),sep="")
names(mt.params) <- paste("mt.",names(mt.params),sep="")
names(poly.params) <- c("poly.a","poly.b","poly.c")

#combine fitted parameters, add chambers and treatments
params <- cbind(q10.params,arr.params,mt.params,poly.params)
params$Code  <- as.factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Code[1]))) # get the Code. Note this could have been a loop
params$Treatment  <- as.factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Treatment[1]))) # get the Code. Note this could have been a loop
params$Taxa  <- factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Taxa[1])),levels=c("BTER","BOT","CTER","BRA")) # get the Code. Note this could have been a loop

#get AIC values
q10.AIC <- abs(as.data.frame(do.call(rbind,lapply(q10mods,AIC))))
arr.AIC <- abs(as.data.frame(do.call(rbind,lapply(arrmods,AIC))))
mt.AIC <- abs(as.data.frame(do.call(rbind,lapply(mtmods,AIC))))
poly.AIC <- abs(as.data.frame(do.call(rbind,lapply(polymods,AIC))))

#rename AIC columns
names(q10.AIC) <- "Q10.AIC"
names(arr.AIC) <- "arr.AIC"
names(mt.AIC) <- "mt.AIC"
names(poly.AIC) <- "poly.AIC"

#put AIC output together with parameter values.
# On average, the MT function (variable Q10) has the lowest |AIC|, but this 
#   changes from chamber to chamber
fitted <- cbind(params,q10.AIC,arr.AIC,mt.AIC,poly.AIC)
fitted <- fitted[,11:16] # Note that Mark's modified Q10 tends to have the lowest AIC for each curve, but not always

# mark's modified Q10 has the lowest overall AIC.
summaryBy(Q10.AIC+arr.AIC+mt.AIC+poly.AIC~1,FUN=sum,data=fitted)

# parameter estimates
params.m <- summaryBy(mt.x + mt.y + mt.Rref~Taxa+Treatment,FUN=c(mean,standard.error),data=params)

#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------








#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#-- visualize polynomial model results
params$R25 <- with(params,exp(poly.a+poly.b*25+poly.c*25^2))
params$Q10 <- with(params,exp(10*(poly.b+2*poly.c*25)))

windows(20,30);par(mar=c(7,6,1,1),mfrow=c(2,1))
boxplot(R25~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="R25")
boxplot(Q10~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="Q10 at 25 C")

lmR <- lm(R25~Treatment*Taxa,data=params)
lmQ10 <- lm(Q10~Treatment*Taxa,data=params)

coefs <- summaryBy(poly.a+poly.b+poly.c~Taxa*Treatment,data=params)
names(coefs)[3:5] <- c("a","b","c")
write.csv(coefs,row.names=F,file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/RvsT/poly_coefs_RvsT_JED.csv")
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- statistically analyze R25 and Q10 estimates
library(effects)

#- R25
lm.R25 <- lm(R25~Taxa*Treatment,data=params)

#look at model diagnostics
plot(lm.R25)
hist(lm.R25$residuals)
anova(lm.R25)
plot(allEffects(lm.R25))     



#- Q10
lm.Q10 <- lm(Q10~Taxa*Treatment,data=params)

#look at model diagnostics
plot(lm.Q10) # point 37 may need to be removed
hist(lm.Q10$residuals)
anova(lm.Q10)
plot(allEffects(lm.Q10))     
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
