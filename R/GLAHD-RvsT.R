#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plots foliar respiraiton vs. Temperature curves for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- read in the data, do a few conversions
dat <- read.csv("Data/GasEx/RvsT/GLADH RvsT gasex processeed.csv",stringsAsFactors=T)
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$Taxa <- as.factor(sapply(1:nrow(dat),function(i)strsplit(as.character(dat$Code[i]),split="-")[[1]][1])) # get the "Taxa" factor out of Code. Note this could have been a loop
dat$R_area <- dat$Photo*-1

#- read in the leaf mass data, calculate R per unit leaf mass
leafmass <- read.csv("C:/Repos/GLAHD/Data/GasEx/RvsT/GLAHD_leaf_mass_RvT.csv")
rt <- merge(dat,leafmass,by="Code")
rt$R_mass <- rt$R_area/10000*rt$Area/rt$mass*1000 #note mass was in g


rt <- subset(rt,Taxa!="SMIT" & Taxa !="ACAM")
rt$Taxa <- factor(rt$Taxa,levels=c("BOT","BTER","BRA","CTER"))
rt$lnRmass <- log(rt$R_mass)
rt <- subset(rt,Code!="CTER-4" & Code!="CTER-22" & Code!="CTER-27") #CTER4 looks weird... peaks much much lower than others. CTER-22 also looks bad
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
  plotBy(R_mass.mean~Tleaf_bin_mid|Treatment,data=toplot,ylim=c(0,40),pch=15,col=c("blue","red"),legend=ifelse(i==1,T,F),
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
params$R18.5 <- with(params,exp(poly.a+poly.b*18.5+poly.c*18.5^2))
params$R28.5 <- with(params,exp(poly.a+poly.b*28.5+poly.c*28.5^2))
params$Q10_25 <- with(params,exp(10*(poly.b+2*poly.c*25)))
params$Q10_18.5 <- with(params,exp(10*(poly.b+2*poly.c*18.5)))
params$Q10_28.5 <- with(params,exp(10*(poly.b+2*poly.c*28.55)))
# 
windows(20,30);par(mar=c(7,6,1,1),mfrow=c(2,1))
boxplot(R25~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="R25")
boxplot(Q10_25~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="Q10 at 25 C")

boxplot(R18.5~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="R18.5")
boxplot(Q10_18.5~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="Q10 at 18.5 C")

boxplot(R28.5~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="R28.5")
boxplot(Q10_28.5~Treatment*Taxa,data=params,las=2,col=c("blue","red"),ylab="Q10 at 28.5 C")
# 
# lmR <- lm(R25~Treatment*Taxa,data=params)
# lmQ10 <- lm(Q10~Treatment*Taxa,data=params)
# 
# coefs <- summaryBy(poly.a+poly.b+poly.c~Taxa*Treatment,data=params)
# names(coefs)[3:5] <- c("a","b","c")
# write.csv(coefs,row.names=F,file="W:/WorkingData/GHS39/GLAHD/Share/Data/GasEx/RvsT/poly_coefs_RvsT_JED.csv")
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------------------
#- predict R and Q10 across all temperatures for all plants
tvals <- seq(from=15,to=40,by=1)
params.list <- split(params,params$Code)
predR <- list()
for(i in 1:length(params.list)){
  predR[[i]] <- data.frame("Code"=params.list[[i]]$Code[1],"Treatment"=params.list[[i]]$Treatment[1],
                            "Taxa"=params.list[[i]]$Taxa[1],"T"=tvals,
                           "R"=with(params.list[[i]],exp(poly.a+poly.b*tvals+poly.c*tvals^2)),
                           "Q10"=with(params.list[[i]],exp(10*(poly.b+2*poly.c*tvals))))
}
predR.df <- do.call(rbind,predR)
predR.df$logR <- log(predR.df$R)
plotBy(R~T|Treatment,data=predR.df)
plotBy(Q10~T|Treatment,data=predR.df)


#-------------------------------------------------------------------------------------------------------------------------
#- plot R with confidence intervals for all taxa overlapping
predR.df$Taxa_Treat <- as.factor(paste(predR.df$Taxa,predR.df$Treatment,sep="_"))
dat.l <- split(predR.df,predR.df$Taxa_Treat)

windows(18,25);par(mfrow=c(2,1),mar=c(4,5,1,1))
fit.sp <- list()
newdat <- list()
plot(Q10~T,data=dat.l[[1]],type="p",col="white",axes=F,xlab="",ylab="",xpd=NA,ylim=c(1,3.5))
colors <- c("grey","grey","green","green","yellow","yellow","cyan","cyan")
fills <- rep(c(0.3,0.9),4)
for (i in 1:length(dat.l)){
  dat.temp <- dat.l[[i]]
  fit.sp[[i]] <-  lm(Q10~T,data=dat.l[[i]])

  newdat[[i]] <- expand.grid(T=seq(from=15,to=40,length.out=99),lower=NA,upper=NA)
  newdat[[i]]$wpred <- predict(fit.sp[[i]],newdat[[i]],se.fit=T)$fit
  newdat[[i]]$se <- predict(fit.sp[[i]],newdat[[i]],se.fit=T)$se.fit
  newdat[[i]]$upper <- with(newdat[[i]],wpred+se)
  newdat[[i]]$lower <- with(newdat[[i]],wpred-se)
  
#   b <- bootCase(fit.sp[[i]])
#   for(j in 1:nrow(newdat[[i]])){
#     b02 <- SSasymp(newdat[[i]]$TDR[j],Asym=b[,"Asym"],lrc=b[,"lrc"],R0=b[,"R0"])
#     newdat[[i]]$lower[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[1]
#     newdat[[i]]$upper[j] <- unname(quantile(b02,probs=c(0.025,0.975)))[2]
#     
#   }
  lines(wpred~T,data=newdat[[i]])
  polygon(x = c(newdat[[i]]$T, rev(newdat[[i]]$T)), y = c(newdat[[i]]$lower, rev(newdat[[i]]$upper)), col = alpha(colors[i],fills[i]), border = NULL,xpd=F)
}
magaxis(side=1:4,labels=c(1,1,0,0),las=1)
title(ylab=expression(Q[10]),xlab=expression(Leaf~Temperature~(degree*C)),cex.lab=1.5,las=1)

#-----------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
#- statistically analyze R25 and Q10 estimates
library(effects)

#- R25
lm.R18.5 <- lm(R18.5~Taxa*Treatment,data=params)
lm.R25 <- lm(R25~Taxa*Treatment,data=params)
lm.R28.5 <- lm(R28.5~Taxa*Treatment,data=params)

#look at model diagnostics
plot(lm.R25)
hist(lm.R25$residuals)
anova(lm.R18.5)
anova(lm.R25)
anova(lm.R28.5)
plot(allEffects(lm.R25))     

effect("Treatment",lm.R25)
(11-9.636)/11

#- Q10
lm.Q10_18.5 <- lm(Q10_18.5~Taxa*Treatment,data=params)
lm.Q10_25 <- lm(Q10_25~Taxa*Treatment,data=params)
lm.Q10_28.5 <- lm(Q10_28.5~Taxa*Treatment,data=params)

#look at model diagnostics
plot(lm.Q10) # point 37 may need to be removed
hist(lm.Q10$residuals)
anova(lm.Q10_18.5)
anova(lm.Q10_25)
anova(lm.Q10_28.5)
plot(effect("Taxa:Treatment",lm.Q10_25))
plot(allEffects(lm.Q10)) 
#-----------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
