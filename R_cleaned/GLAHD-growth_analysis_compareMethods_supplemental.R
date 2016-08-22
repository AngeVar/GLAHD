#- load libraries from script
source("R/loadLibraries.R")
source("R/fitGAM/derivSimulCI.R")
source("R/fitGAM/plotCIdate.R")
source("R/fitGAM/smoothplot.R")

# for GAM
library(mgcv)



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- DATA

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1


#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)


#- Calculate interval-based RGR for each plant
dat3 <- dat3[with(dat3,order(Code,Date)),] #- make sure the observations for each plant are in order
growth.l <- split(dat3,dat3$Code)          #- split into a list for each plant
for (i in 1:length(growth.l)){
  #- use diff() to calculate RGR. Note that I've added a trailing "NA", or else the vector would be too short to fit
  #-   the data frame.
  growth.l[[i]]$RGR <- c(diff(growth.l[[i]]$lnTotMass),NA)/c(unname(diff(growth.l[[i]]$Time)),NA) 
  growth.l[[i]]$dDiameter <- c(diff(growth.l[[i]]$Diameter),NA) 
  
}
dat4 <- do.call(rbind,growth.l)
dat4$RGR_time <- dat4$Time+3

INTERVAL <- dat4
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
















#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- GAM

#- fit gam model to each plant one at a time, predict the derivatives
growth.l <- split(dat3,dat3$Code)
log3fits <- output <- data.out <- list()
kgam=5
gamfits <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  
  #- fit the gam
  g <- gam(lnTotMass ~ s(Time, k=kgam), data=tofit)
  
  #- plot fit
  #smoothplot(Time, lnTotMass, data=tofit, kgam=kgam)
  #title(main=tofit$Code[1])
  
  #- create a vector of "dates" on which to estimate the derivative 
  dates <- seq(min(tofit$Time), max(tofit$Time), by=1)
  
  #- extract the derivative 
  fd <- derivSimulCI(g, samples = 10000, n=length(dates))
  dydt <- fd[[1]]$deriv[,1]
  
  #- put derivatives into a list of dataframes
  gamfits[[i]] <- data.frame(Code=tofit$Code[1],Time=dates,dydt=dydt)
  
  #- get the predicted mass
  newDF <- data.frame(Time=dates) ## needs to be a data frame for predict
  X0 <- exp(predict(g, newDF))    ## exp() needed to convert from ln(mass) to mass
  
  #- put mass into the dataframe
  gamfits[[i]]$predMass <- X0
}

#- merge dataframes, combine with treatment key
gamfits.df <- do.call(rbind,gamfits)
key <- unique(subset(dat3,select=c("Code","Species","Treatment","Location","Taxa","Range")))
gamfits2 <- merge(gamfits.df,key,by=c("Code"),all.x=T)
gamfits2$AGR <- gamfits2$dydt*gamfits2$predMass

GAM <- gamfits2
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------










#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- POWER

#- fit POWER model. Or, skip fitting it (it's slow) and read in the CSV below
library(DEoptim)

set.seed(1234) # set seed for repeatablility
DEoptim.control <- list(VTR = -Inf, strategy = 2,itermax=500,trace=TRUE,CR=0.9)

#----------------------------------------------------------------------------------------------------------------
#- Loop over each plant and estimate the parameters by an evolutionary algorithm. This is slow (5-10 min?).
#- For testing purposes, reduce NP in DE.optim.control to 30. This gives approximately the same answers but is
#-   much much faster. Judging from my system monitor, DEoptim must use parallel processing under the hood.
#----------------------------------------------------------------------------------------------------------------

lower <- log(c(0.0001,0.4,0.001))
upper <- log(c(0.6,0.999,0.5))
dat.l <- split(dat3,dat3$Code)
output <- list()
outpar <- data.frame("Code"=rep(NA,length(dat.l)),"M0"=rep(NA,length(dat.l)),"beta"=rep(NA,length(dat.l)),"r"=rep(NA,length(dat.l)))

#- set up progress bar
pb <- txtProgressBar(min = 0, max = length(dat.l), style = 3)

for(i in 1:length(dat.l)){
  #print(i)
  tofit <- dat.l[[i]]
  output[[i]] <- DEoptim(fit.power3,lower,upper,Time=dat.l[[i]]$Time,TotMass=dat.l[[i]]$TotMass,
                         control=list(VTR = -Inf, strategy = 2,NP=300,itermax=500,trace=F,CR=0.9))
  outpar$Code[i] <- as.character(tofit$Code[1])
  outpar[i,2:4] <- exp(unname(output[[i]]$optim$bestmem))
  setTxtProgressBar(pb, i)
  
}
outpar$Code <-as.factor(outpar$Code)

#- merge it all the treatment code information from dat3
DEfits.all <- merge(outpar,dat3,by="Code")

#- write out a csv, read it in
write.csv(DEfits.all,row.names=F,file="C:/Repos/GLAHD/Output/DE_power_fits_all.csv")
DEfits.all <- read.csv("C:/Repos/GLAHD/Output/DE_power_fits_all.csv")
DEfits.all$Date<- as.Date(DEfits.all$Date)

#Predict daily mass, AGR and RGR
taxa<- as.data.frame(rep(unique(DEfits.all$Code),each=70));colnames(taxa) <- c("Code"); taxa$Time<- rep(seq(1:70))
code<- merge(taxa, unique(DEfits.all[,c("Code","r","beta","M0", "Treatment","Location","Range","Taxa", "Species")], by="Code"))

code$predmass <- with(code,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))
code$AGR <- with(code,r*(M0^(1-beta) + r*Time*(1-beta))^(beta/(1-beta)))
code$RGR <- code$AGR/code$predmass

POWER <- code
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------








#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- POLY

#- fit log-3polynomial growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
log3fits <- log4fits <- output <- data.out <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  log3fits[[i]] <- lm(log(TotMass)~Time+I(Time^2)+I(Time^3),data=tofit) #- fit log-polynomial
  log4fits[[i]] <- lm(log(TotMass)~Time+I(Time^2)+I(Time^3)+I(Time^4),data=tofit) #- fit log-polynomial
  
  #- third order
  #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
  output[[i]] <- output.log_3order(X=tofit$Time,Y=tofit$TotMass,
                                   params=unname(coef(log3fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$rates
  data.out[[i]] <- output.log_3order(X=tofit$Time,Y=tofit$TotMass,
                                     params=unname(coef(log3fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$data
  data.out[[i]]$Code <- tofit$Code[1]
  
#   #- fourth-order
#   #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
#   output[[i]] <- output.log_4order(X=tofit$Time,Y=tofit$TotMass,
#                                    params=unname(coef(log4fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$rates
#   data.out[[i]] <- output.log_4order(X=tofit$Time,Y=tofit$TotMass,
#                                      params=unname(coef(log4fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$data
#   data.out[[i]]$Code <- tofit$Code[1]
}
#- put rates dataframe together. This has the timecourse of mass, AGR, and RGR for each plant
rates.df <-do.call(rbind,output)
rates.df2 <- merge(rates.df,subset(dat2,Date==as.Date("2014-11-17")),by="Code")[,c(1:6,8:10,17)]#,16)] # merge with dat2 to get Treatment, Location, etc

POLY <- rates.df2[with(rates.df2,order(Code,times)),]
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
















#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#- PLOT

#- compare the gam, power, and log-polynomial fits
#- plot RGR via the segment method for each treatment combination. Overlay model outputs.

INTERVAL$combotrt <- as.factor(paste(INTERVAL$Location,INTERVAL$Range,INTERVAL$Treatment,sep=" "))
INTERVAL.mean <- summaryBy(RGR~RGR_time+Location+Range+Treatment+combotrt,data=INTERVAL,FUN=c(mean,standard.error),na.rm=T)
names(INTERVAL.mean)[1] <- "Time"
INTERVAL.mean$titles<- ifelse(INTERVAL.mean$combotrt == "N narrow Home"|INTERVAL.mean$combotrt == "N narrow Warmed", "Tropical-Narrow",
                          ifelse(INTERVAL.mean$combotrt == "N wide Home"|INTERVAL.mean$combotrt == "N wide Warmed", "Tropical-Wide",
                             ifelse(INTERVAL.mean$combotrt == "S narrow Home"|INTERVAL.mean$combotrt == "S narrow Warmed", "Temperate-Narrow",
                             "Temperate-Wide")))
INTERVAL.l <- split(INTERVAL.mean,INTERVAL.mean$combotrt)


GAM$combotrt <- as.factor(paste(GAM$Location,GAM$Range,GAM$Treatment,sep="_"))
GAM.mean <- summaryBy(dydt~Time+Location+Range+Treatment+combotrt,data=GAM,FUN=c(mean,standard.error),na.rm=T)
GAM.l <- split(GAM.mean,GAM.mean$combotrt)

POWER$combotrt <- as.factor(paste(POWER$Location,POWER$Range,POWER$Treatment,sep="_"))
POWER.mean <- summaryBy(RGR~Time+Location+Range+Treatment+combotrt,data=POWER,FUN=c(mean,standard.error),na.rm=T)
POWER.l <- split(POWER.mean,POWER.mean$combotrt)

POLY$combotrt <- as.factor(paste(POLY$Location,POLY$Range,POLY$Treatment,sep="_"))
names(POLY)[2] <- "Time"
POLY.mean <- summaryBy(RGR~Time+Location+Range+Treatment+combotrt,data=POLY,FUN=c(mean,standard.error),na.rm=T)
POLY.l <- split(POLY.mean,POLY.mean$combotrt)

windows(50,60);par(mfrow=c(4,2),mar=c(0.2,0.2,0.2,0.2),oma=c(7,7,3,3))
ylims=c(0,0.18);xlims=c(0,60);lwidth=2
for (i in 1:length(INTERVAL.l)){
  toplot <- INTERVAL.l[[i]]
  gam <- GAM.l[[i]]
  power <- POWER.l[[i]]
  poly <- POLY.l[[i]]
  
  #- plot interval RGR
  plotBy(RGR.mean~Time|Treatment,data=toplot,legend=F,pch=16,xlim=xlims,ylim=ylims,cex=1.5,axes=F,
         panel.first=adderrorbars(x=toplot$Time,y=toplot$RGR.mean,SE=toplot$RGR.standard.error,direction="updown"))
  legend("top",legend=INTERVAL.l[[i]]$titles[1],bty="n",cex=1.5)
  
  #- overlay GAM
  plotBy(dydt.mean~Time|Treatment,data=gam,type="l",add=T,axes=F,legend=F,lty=1,lwd=lwidth,col=c("black","red"))
  
  #- overlay POWER
  plotBy(RGR.mean~Time|Treatment,data=power,type="l",add=T,axes=F,legend=F,lty=2,lwd=lwidth,col=c("black","red"))

  #- overlay POLY
  plotBy(RGR.mean~Time|Treatment,data=poly,type="l",add=T,axes=F,legend=F,lty=3,lwd=lwidth,col=c("black","red"))
  
  
  if(i%%2==1) magaxis(side=c(1,2,3,4),labels=c(0,1,0,0),las=1,frame.plot=T,cex.axis=1.3)
  if(i%%2==0) magaxis(side=c(1,2,3,4),labels=c(0,0,0,0),las=1,frame.plot=T,cex.axis=1.3)
  if(i>=7) magaxis(side=c(1),labels=c(1),las=1,frame.plot=T,cex.axis=1.3)

}
mtext(side=1,"Time",outer=T,line=3,cex=1.3,xpd=NA)
mtext(side=2,expression(RGR~(g~g^-1~d^-1)),outer=T,line=3,cex=1.5)
legend(x=-71,y=0.85,xpd=NA,legend=c("GAM","Power","Polynomial"),bty="n",ncol=3,lty=c(1,2,3),lwd=lwidth,seg.len=3, cex=1.3)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#- PLOT

#- predicted vs. observed

#-- put observations and predictions together
combo1 <- merge(INTERVAL,GAM,by=c("Code","Time","Species","Range","Treatment","Location","Taxa"))[,c(1:7,15,23)]
names(combo1)[9] <- "predMassGam"
combo2 <- merge(combo1,POLY[,c(1,2,3)],by=c("Code","Time"))
names(combo2)[10] <- "predMassPoly"
combo <- merge(combo2,POWER[,c(1,2,11)],by=c("Code","Time"))
names(combo)[11] <- "predMassPower"

windows(60,40);par(mfrow=c(1,3))
ylims=c(0,120);xlims=ylims

#- plot Power 
plot(predMassPower~TotMass,data=combo,axes=F,ylab="",xlab="",log="xy")
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),frame.plot=T,las=1,
        xlab="Observed mass (g)",ylab="Predicted mass, Power (g)")
abline(0,1,lty=2)
lm.power <- lm(log(predMassPower)~log(TotMass),data=combo)
abline(lm.power)
r2 <- round(summary(lm.power)$adj.r.squared,3)
legend("topleft",bty="n",legend=as.expression(bquote(r^2 ~ "=" ~ 0.985)),cex=1.5)
title(main="Power")

#- plot polynomial  
plot(predMassPoly~TotMass,data=combo,axes=F,ylab="",xlab="",log="xy")
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),frame.plot=T,las=1,
        xlab="Observed mass (g)",ylab="Predicted mass, Polynomial (g)")
abline(0,1,lty=2)
lm.poly <- lm(log(predMassPoly)~log(TotMass),data=combo)
abline(lm.poly)
r2 <- round(summary(lm.poly)$adj.r.squared,3)
legend("topleft",bty="n",legend=as.expression(bquote(r^2 ~ "=" ~ 0.996)),cex=1.5)
title(main="Polynomial")

plot(predMassGam~TotMass,data=combo,axes=F,ylab="",xlab="",log="xy")
magaxis(side=c(1,2,3,4),labels=c(1,1,0,0),frame.plot=T,las=1,
        xlab="Observed mass (g)",ylab="Predicted mass, GAM (g)")
abline(0,1,lty=2)
lm.gam <- lm(log(predMassGam)~log(TotMass),data=combo)
abline(lm.gam)
r2 <- round(summary(lm.gam)$adj.r.squared,3)
legend("topleft",bty="n",legend=as.expression(bquote(r^2 ~ "=" ~ 0.997)),cex=1.5)
title(main="GAM")
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------