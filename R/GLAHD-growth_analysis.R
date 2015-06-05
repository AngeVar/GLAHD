#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-7))

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)
# 
# ###################################
# #-- self-starting power function from http://wwwuser.gwdg.de/~cscherb1/self-starting%20power%20law%20function.txt
# powermodel=function(x,a,b,c)
# {a+b*x^c}
# 
# powermodelInit=function(mCall,LHS,data){
#   xy=sortedXyData(mCall[["x"]],LHS,data)
#   lmFit1=lm(xy[,"y"]~1) #for "intercept", a
#   lmFit2=lm(log(xy[,"y"])~log(xy[,"x"])) #for b and c
#   coefs1=coef(lmFit1)
#   coefs2=coef(lmFit2)
#   a=coefs1
#   b=exp(coefs2[1])
#   c=coefs2[2]
#   value=c(a,b,c)
#   names(value)=mCall[c("a","b","c")]
#   value
# }
# 
# SSpower=selfStart(powermodel,powermodelInit,c("a","b","c"))
# 
# 
# nlm1 <- nls(TotMass~SSpower(Time,a,b,c),
#             data=crap)
###################################




# 
# #------------------------------------------------------------------------------------------------------------------------------
# #-- Is a second-order polynomial to log-transformed mass data sufficient? This may simplify subsequent analyses (ANOVA, etc)'
# 
# 
# 
# #-- plot all taxa, with taxa-specific models. The results look pretty good...
# colors <-(c("green","red"))
# dat3.l <- split(dat3,dat3$Taxa)
# windows(20,20)
# par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
# for (i in 1:length(dat3.l)){
#   toplot1 <- dat3.l[[i]]
#   toplot1$Taxa <- factor(toplot1$Taxa)
#   
#   lm1 <- lm(lnTotMass~Time*Treatment+I(Time^2)*Treatment,data=toplot1)
#   
#   newdat <- expand.grid(Taxa=levels(toplot1$Taxa),Treatment=levels(toplot1$Treatment), Time=1:60)
#   newdat$wpred <- predict(lm1,newdat,level=0)
#   
#   plotBy(lnTotMass~jitter(Time)|Treatment,data=toplot1,col=colors,legend=F,axes=F,ylim=c(-2,5),pch=20,cex=1.3)
#   plotBy(wpred ~ Time|Treatment, data=newdat, type='l',col=colors, add=TRUE, lwd=2,legend=F)
#   title(main=toplot1$Taxa[1],line=-2)
#   magaxis(side=1:4,labels=0)
#   
# }
# mtext(expression(Time~(days)),side=1,line=2,outer=T,cex=1.5)
# mtext(expression(ln~(Total~plant~mass~(g))),side=2,line=2,outer=T,cex=1.5)
# legend(x=134,y=3,legend=c("Home","Warmed"),pch=20,cex=1.5,xpd=NA,col=colors[1:2])

# 
# #-- fit linear model
# lm1 <- lm(lnTotMass~Range*Location*Treatment*Time+Range*Location*Treatment*I(Time^2),data=dat3)
# 
# #plot predictions
# newdat <- expand.grid(Range=levels(dat3$Range),Treatment=levels(dat3$Treatment),
#                       Location=levels(dat3$Location),Time=0:60)
# newdat$lnPred <- predict(lm1,newdat,level=0)
# 
# windows(20,12);par(mfrow=c(1,2))
# plotBy(lnPred~Time|Treatment,data=subset(newdat,Range=="narrow" & Location=="S"),type="l",lty=1,legend=F,lwd=2,ylim=c(-1,4))
# plotBy(lnPred~Time|Treatment,data=subset(newdat,Range=="wide" & Location=="S"),type="l",lty=2,add=T,legend=F,lwd=2)
# mtext("South",side=3)
# plotBy(lnPred~Time|Treatment,data=subset(newdat,Range=="narrow" & Location=="N"),type="l",lty=1,legend=F,lwd=2,ylim=c(-1,4))
# plotBy(lnPred~Time|Treatment,data=subset(newdat,Range=="wide" & Location=="N"),type="l",lty=2,add=T,legend=F,lwd=2)
# mtext("North",side=3)

#---------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------------------------------




# from Paine et al. 2012
###################################
# 2-parameter Power-law fit, initial mass fixed
###################################




seed.mass <- 0.4 #(45 mg, Hautier et al 2010, J. Ecology) # change this to suit!
fmla.pow2 <- as.formula(paste("~(", seed.mass, "^(1-beta) + r*x*(1-beta))^(1/(1-beta))", sep = ""))
SS.pow2   <- selfStart(fmla.pow2, initial = Init.pow2, parameters = c("r", "beta"))
xvals <- seq(1,70,by=0.5) # establish xvalues of time to use in the model fitting below



source("R/Paine.R")
##############################################################################################
##############################################################################################
#- fit growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
xlims <- c(0,70)
ylims <- c(0,100)

pdf(file="Output/Growth_power_allplants.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)

output <- data.frame(Species="",Treatment="",Location="",Taxa="",Code="",Range="",r="",beta="",AGR10g="",RGR10d="",stringsAsFactors=F)
full.output <- list()
for (i in 1:length(growth.l)){
  print(i)
  tofit <- growth.l[[i]]
  tofit <- subset(tofit,TotMass>seed.mass)
  #fit
  tmp.pow2 <- getInitial(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  fit.pow2 <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  out.pow2 <- output.pow2.gnls(fit.pow2,xvals, CI = T, M0 = seed.mass)
  
  #plot
  plot(TotMass~Time,data=tofit,pch=20,cex=3,ylab="Total plant mass (g)",xlab="Time (days)",xlim=xlims,ylim=ylims)
  lines(M~times,data=out.pow2$rates)
  mtext(text=tofit$Code[1],side=3,line=-1)
  legend("topleft",paste("r = ",round(coef(fit.pow2)[1],digits=4)),bty="n")
  legend("left",paste("beta = ",round(coef(fit.pow2)[2],digits=4)),bty="n")
  
  #compile parameters
  output[i,1:6] <- c(as.character(tofit$Species[1]),as.character(tofit$Treatment[1]),as.character(tofit$Location[1]),
                     as.character(tofit$Taxa[1]),as.character(tofit$Code[1]),as.character(tofit$Range[1]))
  output[i,"r"] <- coef(fit.pow2)[1]
  output[i,"beta"] <- coef(fit.pow2)[2]
  index_time <- findInterval(10,out.pow2$rates$times) # this was altered from "M" to "times" to get AGR and RGR at 10 days, rather than 10 grams
  index_mass <- findInterval(10,out.pow2$rates$M) # this was altered from "M" to "times" to get AGR and RGR at 10 days, rather than 10 grams
  
  output[i,"AGR10g"] <- mean(out.pow2$rates$AGR[index_mass],out.pow2$rates$AGR[index_mass+1]) # get AGR at 10 grams
  output[i,"AGR10d"] <- out.pow2$rates$AGR[index_time]           # get AGR at 10 days
  output[i,"RGR10d"] <- out.pow2$rates$RGRt[index_time]         # get RGT at 10 days
  output[i,"RGR10g"] <- mean(out.pow2$rates$RGRm[index_mass],out.pow2$rates$RGRm[index_mass+1])# get RGT at 10 grams
  
  #- grabs all output of the "rates' dataframe
  full.output[[i]] <- out.pow2$rates
  full.output[[i]]$Species <- tofit$Species[1]
  full.output[[i]]$Treatment <- tofit$Treatment[1]
  full.output[[i]]$Location <- tofit$Location[1]
  full.output[[i]]$Taxa <- tofit$Taxa[1]
  full.output[[i]]$Pot <- tofit$Pot[1]
  full.output[[i]]$Code <- tofit$Code[1]
  full.output[[i]]$Range <- tofit$Range[1]
  
  
}
dev.off()
output[,7:10] = apply(output[,7:10], 2, function(x) as.numeric(as.character(x)))
output$Species <- as.factor(output$Species)
output$Treatment <- as.factor(output$Treatment)
output$Location <- as.factor(output$Location)
output$Taxa <- as.factor(output$Taxa)
output$Code <- as.factor(output$Code)
output$Range <- as.factor(output$Range)

#- merge full.output into a dataframe
#- use this to plot response ratio of biomass to warming over time? As Miko did?
full.output.df <- do.call(rbind,full.output)

#- get the estimate of final mass to add to this output dataframe
dat.f <- subset(dat2,Date==as.Date("2015-01-05"))[,c(1,3,4,5,7,13)]
names(dat.f)[6] <- "TotMass.f"
 
output2 <- merge(output,dat.f,by=c("Species","Treatment","Location","Taxa","Code"))

#- write out a csv for Angelica to analyze
#write.csv(output2,file="Data/GLAHD_output_growth_parameters_15052015.csv",row.names=F)
output <- read.csv("Data/GLAHD_output_growth_parameters_15052015.csv")

# r.lm <- lm(r~Treatment*Location*Range,data=output)
# plot(allEffects(r.lm))
# beta.lm <- lm(beta~Treatment*Location*Range,data=output)
# plot(allEffects(beta.lm))

windows(20,12);par(mar=c(8,7,1,1))

output$Taxa <- factor(output$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)
boxplot(AGR10~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="AGR",cex.lab=2)
abline(v=16.5)
boxplot(RGR10~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="RGR",cex.lab=2)
abline(v=16.5)

##############################################################################################
##############################################################################################





##############################################################################################
##############################################################################################
#-- attempt a statisitical analysis of growth parameters.

#---
#- model for r
str(output) #- data structure look right
fm.r <- lme(log(r)~Treatment*Location*Range,random=list(~1|Species,~1|Taxa),data=output)   # R can determine the random term structure itself (most of the time!)
anova(fm.r)

#look at model diagnostics
plot(fm.r,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.r,log(r)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.r,log(r)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.r, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals

summary(fm.r) #- look at the std estiamtes of the random componenets
ranef(fm.r)   #- extract random effects. BOT and PEL had high r, while PLAT, SMIT, and LONG had lower r

#- plot effects
plot(effect("Treatment:Location:Range",fm.r)) #- There was a 2-way interaction in the S (warming increased r in wide species only), 
                                              #- and this interaction was not present in the N
plot(effect("Treatment:Location",fm.r))
plot(effect("Treatment:Range",fm.r))
#---






















##############################################################################################
##############################################################################################
# RGR decomposition
##############################################################################################
##############################################################################################
#-- merge r and beta estimates for each plant with the size dataframe, in attempt to 
#-- decompose RGR into its components
dat4 <- merge(dat3,output2,by=c("Species","Treatment","Location","Taxa","Code","Range"))
dat4$RGR <- with(dat4,r*TotMass^(beta-1)) # calculate RGR on mass basis. See Paine et al. 2012.
dat4$AGR <- with(dat4,r*(0.5^(1-beta)+r*Time*(1-beta))^(beta/(1-beta))  )
# use function to predict total plant leaf area from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
dat4$leafarea <- predict_LA(Taxa=dat4$Taxa,Treat=dat4$Treatment,Diameter=dat4$Diameter,Height=dat4$Height)
dat4$LAR <- with(dat4,leafarea/TotMass)
#dat4$ULR <- with(dat4,AGR/leafarea)
dat4$ULR <- with(dat4,RGR/LAR)

dat4.m <- summaryBy(TotMass+RGR+AGR+LAR+ULR~Treatment+Location+Range+Time,FUN=c(mean,standard.error),data=dat4)
##############################################################################################


##############################################################################################
#-- plot AGR, RGR, LAR, and ULR over mass
windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(dat4,Location=="S")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.15))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.15))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,375))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0022))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0022))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_south2.pdf")
##############################################################################################

windows(20,20);par(mfrow=c(4,2),mar=c(2,2,0,0),oma=c(5,5,4,1))
toplot <- subset(dat4,Location=="N")
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=2,"AGR",line=3)
mtext(side=3,"Narrow",line=1)
plotBy(AGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="AGR",xlab="Mass",xlim=c(0,70),ylim=c(0,5))
mtext(side=3,"Wide",line=1)

plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.3))
mtext(side=2,"RGR",line=3)
plotBy(RGR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.3))

plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,425))
mtext(side=2,"LAR",line=3)
plotBy(LAR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(50,425))

plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="narrow"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0035))
mtext(side=1,"Plant mass",line=3)
mtext(side=2,"ULR",line=3)
plotBy(ULR~TotMass|Treatment,data=subset(toplot,Range=="wide"),legend=F,pch=16,ylab="RGR",xlab="Mass",xlim=c(0,70),ylim=c(0,0.0035))
mtext(side=1,"Plant mass",line=3)
dev.copy2pdf(file="Output/RGR_decomposition_north.pdf")

##############################################################################################
##############################################################################################
##############################################################################################


##############################################################################################

dat5 <- merge(dat3,full.output.df,by=c("Species","Treatment","Location","Taxa","Code","Range"))


dat4 <- merge(dat3,output2,by=c("Species","Treatment","Location","Taxa","Code","Range"))
dat4$RGR <- with(dat4,r*TotMass^(beta-1)) # calculate RGR on mass basis. See Paine et al. 2012.
dat4$AGR <- with(dat4,r*(0.5^(1-beta)+r*Time*(1-beta))^(beta/(1-beta))  )
# use function to predict total plant leaf area from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
dat4$leafarea <- predict_LA(Taxa=dat4$Taxa,Treat=dat4$Treatment,Diameter=dat4$Diameter,Height=dat4$Height)
dat4$LAR <- with(dat4,leafarea/TotMass)
#dat4$ULR <- with(dat4,AGR/leafarea)
dat4$ULR <- with(dat4,RGR/LAR)

dat4.m <- summaryBy(TotMass+RGR+AGR+LAR+ULR~Treatment+Location+Range+Time,FUN=c(mean,standard.error),data=dat4)

full.output.df$leafarea <- predict_LA(Taxa=full.output.df$Taxa,Treat=full.output.df$Treatment,Diameter=full.output.df$Diameter,Height=full.output.df$Height)








##############################################################################################
##############################################################################################
#- compare non-linear functions
fit.pow <- nlme(TotMass ~ SS.pow2(Time, r, beta),
             fixed=list(r ~ Taxa, beta ~ 1),
             random = r ~ 1 | Code,
             start=list(fixed=c(r=rep(0.05,17),beta=rep(0.9,1))),
             data=dat3)

fit.asym <- nlme(TotMass ~ SSgompertz(Time, Asym, b2, b3),
                fixed=list(Asym  ~ Taxa, b2~1, b3 ~ 1),
                random = b3 ~ 1 | Code,
                start=list(fixed=c(Asym=rep(0.05,17),b2=0,c3=10)),
                data=dat3)
fit.gomp<- nls(TotMass ~ SSgompertz(Time, Asym, b2,b3), data = subset(dat3,Code=="BRA-21"))
##############################################################################################
##############################################################################################











#- fit them all at once in a generalized non-linear mixed effects model
fit1 <- nlme(TotMass ~ SS.pow2(Time, r, beta),
                 fixed=list(r ~ Taxa*Treatment, beta ~ 1),
                 random = r ~ 1 | Code,
                 start=list(fixed=c(r=rep(0.05,17*2),beta=rep(0.9,1))),
                 data=dat3)
fit2 <- nlme(TotMass ~ SS.pow2(Time, r, beta),
             fixed=list(r ~ Taxa*Treatment, beta ~ 1),
             random =r+beta ~ 1 | Code,
             start=list(fixed=c(r=rep(0.05,17*2),beta=rep(0.9,1))),
             data=dat3)
fit3 <- nlme(TotMass ~ SS.pow2(Time, r, beta),
             fixed=list(r ~ Location*Range*Treatment, beta ~ 1),
             random =beta ~ 1 | Code,
             start=list(fixed=c(r=rep(0.05,8),beta=rep(0.9,1))),
             data=dat3)

anova(fit1,fit2) # fit random terms for both r and beta

# Fixed effects predictions
# Make dataframe with all combinations of treatments.
newdat <- expand.grid(Taxa=levels(dat3$Taxa),Treatment=levels(dat3$Treatment), Time=7:66)
#newdat <- expand.grid(Location=levels(dat3$Location),Range=levels(dat3$Range),Treatment=levels(dat3$Treatment), Time=7:66)

newdat$wpred <- predict(fit2,newdat,level=0)

# Plot. I use jitter so not all data cover each other.
palette(topo.colors(8,alpha=0.8))
plot(TotMass~jitter(Time),data=dat3,col=as.factor(paste(Location,Range,Treatment)),main="fit1",type="n")
plotBy(wpred ~ Time|Location, data=subset(newdat,Range=="narrow" & Treatment=="Home"), type='l', add=TRUE, lwd=2,legend=F)
plotBy(wpred ~ Time|Location, data=subset(newdat,Range=="wide" & Treatment=="Home"), type='l',lty=2, add=TRUE, lwd=2,legend=F)
plotBy(wpred ~ Time|Location, data=subset(newdat,Range=="narrow" & Treatment=="Warmed"), col="red",type='l', add=TRUE, lwd=2,legend=F)
plotBy(wpred ~ Time|Location, data=subset(newdat,Range=="wide" & Treatment=="Warmed"),col="red", type='l',lty=2, add=TRUE, lwd=2,legend=F)
legend("topleft",c("Narrow-Home","Wide-Home","Narrow-Warmed","Wide-Warmed"),
       lty=c(1,2,1,2),col=c("blue","blue","red","red"))

#---------------------------------------------------------------------------
#- plot each taxa separately
newdat.l <- split(newdat,newdat$Taxa)
dat3.l <- split(dat3,dat3$Taxa)
colors <-(c("green","red"))

windows(20,20)
par(mfrow=c(5,4),mar=c(0,0,0,0),oma=c(6,6,1,1))
for (i in 1:length(dat3.l)){
  toplot1 <- dat3.l[[i]]
  toplot2 <- newdat.l[[i]]
  print(i)
  
  plotBy(TotMass~jitter(Time)|Treatment,data=toplot1,col=colors,legend=F,axes=F,ylim=c(0,100),pch=20,cex=1.3)
  plotBy(wpred ~ Time|Treatment, data=toplot2, type='l',col=colors, add=TRUE, lwd=2,legend=F)
  title(main=toplot1$Taxa[1],line=-2)
  magaxis(side=1:4,labels=0)
  
}
mtext(expression(Time~(days)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(Total~plant~mass~(g)),side=2,line=2,outer=T,cex=1.5)
legend(x=134,y=74,legend=c("Home","Warmed"),pch=20,cex=1.5,xpd=NA,col=colors[1:2])
#---------------------------------------------------------------------------















#########################
##-- attempting to replicate Fig. 3 in Paine 2012 for a subset of taxa.
### Primary task illustrated in this figure: propagating error from parameters into estimates of RGR. 
#Before propagating error, we need to check that the parameter profiles are approximately V-shaped (and that the sampling intervals of the parameters are therefore approximately multivariate normal) before proceeding. 
#To do this, we need to fit a seperate nls fit for each species. We have to use nls(), rather than gnls() or nlsList, because profile methods do not exist for those functions.
narrow.s.h <- subset(dat3,Location=="S" & Range =="narrow" & Treatment=="Home")
narrow.s.w <- subset(dat3,Location=="S" & Range =="narrow" & Treatment=="Warmed")
wide.s.h <- subset(dat3,Location=="S" & Range =="wide" & Treatment=="Home")
wide.s.w <- subset(dat3,Location=="S" & Range =="wide" & Treatment=="Warmed")



fit.nsh   <- nls(TotMass ~ SS.pow2(Time, r, beta), data = narrow.s.h)
fit.nsw   <- nls(TotMass ~ SS.pow2(Time, r, beta), data = narrow.s.w)
fit.wsh   <- nls(TotMass ~ SS.pow2(Time, r, beta), data = wide.s.h)
fit.wsw   <- nls(TotMass ~ SS.pow2(Time, r, beta), data = wide.s.w)


par(mfrow = c(4, 2), oma = c(3, 2, 3, 1))
plot(profile(fit.nsh, alphamax = 0.1))
plot(profile(fit.nsw, alphamax = 0.1))
plot(profile(fit.wsh, alphamax = 0.1))
plot(profile(fit.wsw, alphamax = 0.1))

mtext("Confidence intervals based on the profile sum of squares", side = 3, outer = TRUE)
## We see that the parofiles are nicely v-shaped

#derive output of power functions for the two treatments separately
times <- seq(min(dat3$Time),max(dat3$Time),length=101)
fit.nsh   <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = narrow.s.h)
fit.nsw   <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = narrow.s.w)
fit.wsh   <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = wide.s.h)
fit.wsw   <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = wide.s.w)

out.nsh     <- output.pow2.gnls(fit.nsh, times, CI = T, alpha = 0.05, M0=seed.mass)
out.nsw     <- output.pow2.gnls(fit.nsw, times, CI = T, alpha = 0.05, M0=seed.mass)
out.wsh     <- output.pow2.gnls(fit.wsh, times, CI = T, alpha = 0.05, M0=seed.mass)
out.wsw     <- output.pow2.gnls(fit.wsw, times, CI = T, alpha = 0.05, M0=seed.mass)



#--- plot
windows(20,20)
par(mfrow = c(4, 2), bty = "n", mar = c(4, 4, 0, 0), oma = c(0, 0, 1, 0), las = 1, cex.lab=1.3, pty = "m", tcl = 0.2, mgp = c(2.5, 0.5, 0))
COL=c("blue","red")

# (a)  Biomass over time
plot(narrow.s.h$TotMass~jitter(narrow.s.h$Time), col = COL[1], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
points(narrow.s.w$TotMass~jitter(narrow.s.w$Time), col = COL[2], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
polygon(x = c(out.nsh$rates$times, rev(out.nsh$rates$times)), y = c(out.nsh$rates$M.lo, rev(out.nsh$rates$M.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.nsw$rates$times, rev(out.nsw$rates$times)), y = c(out.nsw$rates$M.lo, rev(out.nsw$rates$M.hi)), col = alpha(COL[2],0.5), border = NA)
legend(5, 50, legend = c("Home", "Warmed"), col = COL[1:2], lty = c(1), pch = c(20), merge = T,cex=1.2)
mtext("Narrow-ranged",side=3,line=-1.5,outer=F)

plot(wide.s.h$TotMass~jitter(wide.s.h$Time), col = COL[1], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
points(wide.s.w$TotMass~jitter(wide.s.w$Time), col = COL[2], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
polygon(x = c(out.wsh$rates$times, rev(out.wsh$rates$times)), y = c(out.wsh$rates$M.lo, rev(out.wsh$rates$M.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$times, rev(out.wsw$rates$times)), y = c(out.wsw$rates$M.lo, rev(out.wsw$rates$M.hi)), col = alpha(COL[2],0.5), border = NA)
mtext("Wide-ranged",side=3,line=-1.5,outer=F)


# (b) AGR on a time basis
plot(out.nsh$rates$AGR~out.nsh$rates$times, xlab = "Time (days)", ylab = "AGR (g/day)", type = "n", xlim = c(0, 70), ylim = c(0, 2.1))
polygon(x = c(out.nsh$rates$times, rev(out.nsh$rates$times)), y = c(out.nsh$rate$AGR.lo, rev(out.nsh$rates$AGR.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$times, rev(out.wsw$rates$times)), y = c(out.wsw$rate$AGR.lo, rev(out.wsw$rates$AGR.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.nsh$rates$times, out.nsh$rates$AGR, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$times, out.wsw$rates$AGR, col = COL[2], lty = 1,lwd=2)    

plot(out.wsh$rates$AGR~out.wsh$rates$times, xlab = "Time (days)", ylab = "AGR (g/day)", type = "n", xlim = c(0, 70), ylim = c(0, 2.1))
polygon(x = c(out.wsh$rates$times, rev(out.wsh$rates$times)), y = c(out.wsh$rate$AGR.lo, rev(out.wsh$rates$AGR.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$times, rev(out.wsw$rates$times)), y = c(out.wsw$rate$AGR.lo, rev(out.wsw$rates$AGR.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.wsh$rates$times, out.wsh$rates$AGR, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$times, out.wsw$rates$AGR, col = COL[2], lty = 1,lwd=2)    



# RGR on a time basis
plot(out.nsh$rates$RGRt~out.nsh$rates$times, xlab = "Time (days)", ylab = "RGR (g/g/day)", type = "n", xlim = c(0, 70), ylim = c(0, 0.1))
polygon(x = c(out.nsh$rates$times, rev(out.nsh$rates$times)), y = c(out.nsh$rate$RGRt.lo, rev(out.nsh$rates$RGRt.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.nsw$rates$times, rev(out.nsw$rates$times)), y = c(out.nsw$rate$RGRt.lo, rev(out.nsw$rates$RGRt.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.nsh$rates$times, out.nsh$rates$RGRt, col = COL[1], lty = 1,lwd=2)    
lines(out.nsw$rates$times, out.nsw$rates$RGRt, col = COL[2], lty = 1,lwd=2)    

plot(out.wsh$rates$RGRt~out.wsh$rates$times, xlab = "Time (days)", ylab = "RGR (g/g/day)", type = "n", xlim = c(0, 70), ylim = c(0, 0.1))
polygon(x = c(out.wsh$rates$times, rev(out.wsh$rates$times)), y = c(out.wsh$rate$RGRt.lo, rev(out.wsh$rates$RGRt.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$times, rev(out.wsw$rates$times)), y = c(out.wsw$rate$RGRt.lo, rev(out.wsw$rates$RGRt.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.wsh$rates$times, out.wsh$rates$RGRt, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$times, out.wsw$rates$RGRt, col = COL[2], lty = 1,lwd=2)    


# RGR on a mass basis
plot(out.nsh$rates$M, out.nsh$rates$RGRm, xlab = "Predicted biomass (g)", ylab = "RGR (g/g/day)", type = "n", axes = F, ylim = c(0, 0.1), xlim = c(0, 50))
axis(1)
axis(2)
polygon(x = c(out.nsh$rates$M, rev(out.nsh$rates$M)), y = c(out.nsh$rates$RGRm.lo, rev(out.nsh$rates$RGRm.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.nsw$rates$M, rev(out.nsw$rates$M)), y = c(out.nsw$rates$RGRm.lo, rev(out.nsw$rates$RGRm.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.nsh$rates$M, out.nsh$rates$RGRm, col = COL[1], lty = 1,lwd=2)    
lines(out.nsw$rates$M, out.nsw$rates$RGRm, col = COL[2], lty = 1,lwd=2)    

plot(out.nsh$rates$M, out.nsh$rates$RGRm, xlab = "Predicted biomass (g)", ylab = "RGR (g/g/day)", type = "n", axes = F, ylim = c(0, 0.1), xlim = c(0, 50))
axis(1)
axis(2)
polygon(x = c(out.wsh$rates$M, rev(out.wsh$rates$M)), y = c(out.wsh$rates$RGRm.lo, rev(out.wsh$rates$RGRm.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$M, rev(out.wsw$rates$M)), y = c(out.wsw$rates$RGRm.lo, rev(out.wsw$rates$RGRm.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.wsh$rates$M, out.wsh$rates$RGRm, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$M, out.wsw$rates$RGRm, col = COL[2], lty = 1,lwd=2)    
dev.copy2pdf(file="Output/Growth_analysis_south.pdf")
  







##############################
### Figure 3               ###
### compare timing of      ###
### measuring RGRs         ###
##############################
### Primary task illustrated in this figure: propagating error from parameters into estimates of RGR. 
#Before propagating error, we need to check that the parameter profiles are approximately V-shaped (and that the sampling intervals of the parameters are therefore approximately multivariate normal) before proceeding. 
#To do this, we need to fit a seperate nls fit for each species. We have to use nls(), rather than gnls() or nlsList, because profile methods do not exist for those functions.
fit.logis.CER   <- nls(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp, subset = species == "CER")
fit.logis.GER   <- nls(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp, subset = species == "GER")
par(mfrow = c(2, 3), oma = c(3, 2, 3, 1))
plot(profile(fit.logis.CER, alphamax = 0.1))
plot(profile(fit.logis.GER, alphamax = 0.1))
mtext("Confidence intervals based on the profile sum of squares", side = 3, outer = TRUE)
mtext("Cerastium (top row) and Geranium (Bottom row) - confidence levels of 50%, 80%, 90% and 95%", side = 1, outer = TRUE)
## We see that the parofiles are nicely v-shaped
## Now we proceed to propagating error.

# Fit logistic functions for the two species using nlsList
# Look at RGR for two species nutrients, using logistic function
fit.logis4.0   <- nlsList(Y ~ SSlogis(X, Asym, xmid, scal), data = dat_asymp_spp)
out.logis4     <- output.logis.nlsList(fit.logis4.0, times = Xes_annuals_spp$X, CI = T, LOG = F, alpha = 0.05)
fit.logis.nlme <- nlme(fit.logis4.0)
# WHY is the variance and covariance among estimated parameters so much greater in nlme than in nlsList? Because it includes the among-species variance, as well as the within-species variance. thus, the vcov matrix from nlme is inappropriate for comparing among species (or other treatment groups). So, use the nlsList fit.

# compute differences in timing and magnitude of AGR
time.max.AGR.N <- out.logis4$rates[[1]]$times[out.logis4$rates[[1]]$AGR == max(out.logis4$rates[[1]]$AGR)]
time.max.AGR.Y <- out.logis4$rates[[2]]$times[out.logis4$rates[[2]]$AGR == max(out.logis4$rates[[2]]$AGR)]
time.max.AGR.Y-time.max.AGR.N # difference of 45 days
#Magnitude of peak AGR
max.AGR.N <- max(out.logis4$rates[[1]]$AGR)
max.AGR.Y <- max(out.logis4$rates[[2]]$AGR)
(max.AGR.Y-max.AGR.N)/max.AGR.N # increase of 28%


pdf("Figure 3 timing of RGR comparisons.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)
COL.CI <- "#55555530"
par(mfrow = c(2, 2), bty = "n", mar = c(4, 4, 1, 1), oma = c(0, 0, 1, 0), las = 1, family = "Helvetica", pty = "s", tcl = 0.2, mgp = c(2.5, 0.5, 0))
# Biomass over time
plot(dat_asymp_spp$X, dat_asymp_spp$Y, col = COL[7], pch = ifelse(dat_asymp_spp$species =="GER", 3, 1), xlab = "Days since sowing", ylab = "Biomass (g)", ylim = c(0, 25), xlim = c(0, 200))
for(i in 1:nrow(out.logis4$params)){
  polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$M.lo, rev(out.logis4$rates[[i]]$M.hi)), col = COL.CI, border = NA)
  lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$M, col = COL[7], lty = i)    # Logistic
}
legend(5, 23, legend = c("Geranium", "Cerastium"), col = COL[7], lty = c(2, 1), pch = c(3, 1), merge = F)
mtext("a)", adj = 0.1, line = -1)
mtext("Figure 3", adj = 0.05, cex = 1.2)

# AGR on a time basis
plot(out.logis4$rates[[2]]$times, out.logis4$rates[[2]]$AGR, xlab = "Days since sowing", ylab = "AGR (g/day)", type = "n", xlim = c(0, 200), ylim = c(0, 0.2))
for(i in 1:nrow(out.logis4$params)){
  polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$AGR.lo, rev(out.logis4$rates[[i]]$AGR.hi)), col = COL.CI, border = NA)
  lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$AGR, col = COL[7], lty = i)    # Logistic
}
mtext("b)", adj = 0.1, line = -1)

# RGR on a time basis
plot(out.logis4$rates[[2]]$times,  out.logis4$rates[[2]]$RGRt, xlab = "Days since sowing", ylab = "RGR (g/g/day)", type = "n", xlim = c(0, 200), ylim = c(0, 0.1))
for(i in 1:nrow(out.logis4$params)){
  polygon(x = c(out.logis4$rates[[i]]$times, rev(out.logis4$rates[[i]]$times)), y = c(out.logis4$rates[[i]]$RGRt.lo, rev(out.logis4$rates[[i]]$RGRt.hi)), col = COL.CI, border = NA)
  lines(out.logis4$rates[[i]]$times, out.logis4$rates[[i]]$RGRt, col = COL[7], lty = i)    # Logistic
}
mtext("c)", adj = 0.1, line = -1)

# RGR on a mass basis
plot(out.logis4$rates[[2]]$M, out.logis4$rates[[2]]$RGRm, xlab = "Predicted biomass (g)", ylab = "", type = "n", axes = F, ylim = c(0, 0.1), xlim = c(0, 11))
axis(1)
axis(2, labels = NA)
for(i in 1:nrow(out.logis4$params)){
  polygon(x = c(out.logis4$rates[[i]]$M, rev(out.logis4$rates[[i]]$M)), y = c(out.logis4$rates[[i]]$RGRm.lo, rev(out.logis4$rates[[i]]$RGRm.hi)), col = COL.CI, border = NA)
  lines(out.logis4$rates[[i]]$M, out.logis4$rates[[i]]$RGRm, col = COL[7], lty = i)    # Logistic
}
mtext("d)", adj = 0.1, line = -1)
dev.off()







#--- plot a few things for HIE science day 2015

#--- plot
windows(15,20)
par(mfrow = c(3, 1), bty = "n", mar = c(4, 4, 0, 0), oma = c(0, 0, 1, 0), las = 1, cex.lab=2.5, pty = "m", tcl = 0.2, mgp = c(2.5, 0.5, 0))
COL=c("blue","red")

# (a)  Biomass over time
plot(wide.s.h$TotMass~jitter(wide.s.h$Time), cex=1.2,col = COL[1], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
points(wide.s.w$TotMass~jitter(wide.s.w$Time), cex=1.2,col = COL[2], pch = 20, xlab = "Time (days)", ylab = "Biomass (g)", ylim = c(0, 70), xlim = c(0, 70))
polygon(x = c(out.wsh$rates$times, rev(out.wsh$rates$times)), y = c(out.wsh$rates$M.lo, rev(out.wsh$rates$M.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$times, rev(out.wsw$rates$times)), y = c(out.wsw$rates$M.lo, rev(out.wsw$rates$M.hi)), col = alpha(COL[2],0.5), border = NA)
legend(2, 70, legend = c("Home", "Warmed"), col = COL[1:2], lty = c(1), pch = c(20), merge = T,cex=1.5)

# (c) RGR vs. mass
plot(out.nsh$rates$M, out.nsh$rates$RGRm, xlab = "Biomass (g)", ylab = "RGR (g/g/day)", type = "n", axes = F, ylim = c(0, 0.1), xlim = c(0, 50))
axis(1)
axis(2)
polygon(x = c(out.wsh$rates$M, rev(out.wsh$rates$M)), y = c(out.wsh$rates$RGRm.lo, rev(out.wsh$rates$RGRm.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$M, rev(out.wsw$rates$M)), y = c(out.wsw$rates$RGRm.lo, rev(out.wsw$rates$RGRm.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.wsh$rates$M, out.wsh$rates$RGRm, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$M, out.wsw$rates$RGRm, col = COL[2], lty = 1,lwd=2)


# (b) AGR vs. mass
plot(out.nsh$rates$M, out.nsh$rates$AGR, xlab = "Biomass (g)", ylab = "AGR (g/day)", type = "n", axes = F, ylim = c(0, 2.5), xlim = c(0, 50))
axis(1)
axis(2)
polygon(x = c(out.wsh$rates$M, rev(out.wsh$rates$M)), y = c(out.wsh$rates$AGR.lo, rev(out.wsh$rates$AGR.hi)), col = alpha(COL[1],0.5), border = NA)
polygon(x = c(out.wsw$rates$M, rev(out.wsw$rates$M)), y = c(out.wsw$rates$AGR.lo, rev(out.wsw$rates$AGR.hi)), col = alpha(COL[2],0.5), border = NA)
lines(out.wsh$rates$M, out.wsh$rates$AGR, col = COL[1], lty = 1,lwd=2)    
lines(out.wsw$rates$M, out.wsw$rates$AGR, col = COL[2], lty = 1,lwd=2)










##############################################################################

#Statistical testing of growth parameters
#old analysis based on data from mass~d2h allometry with all different slopes
params<-read.csv("Data/GLAHD_output_growth_parameters_15052015.csv")
str(params)
params$Location <- factor(params$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
params$Sp_RS_EN <- as.factor(with(params,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
params$Prov_Sp_EN <- as.factor(with(params,paste(Taxa,Species)))
params$Sp_Loc_EN <- as.factor(with(params,paste(Species,Location)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize

library(nlme)
library(effects)
output2$Location <- factor(output2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
output2$Sp_RS_EN <- as.factor(with(output2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
output2$Prov_Sp_EN <- as.factor(with(output2,paste(Taxa,Species)))
output2$Sp_Loc_EN <- as.factor(with(output2,paste(Species,Location)))
#-----------------------------------------------
##model from meeting: fm1 <- lme(TotMass~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat.f) 

fm1r <- lme(r~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=params)

plot(fm1r,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r,r~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1r,r~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1r, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1r$residuals[,1])
anova(fm1r)

fm1r.l <- lme(log(r)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
plot(fm1r.l,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r.l,log(r)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1r.l,log(r)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1r.l, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1r.l$residuals[,1])
anova(fm1r.l)    

plot(allEffects(fm1r.l))      #- try to make sense of Treatment:Location:Range interaction
#------------------------------------------------

#------------------------------------------------
fm1b <- lme(beta~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1b,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1b,beta~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1b,beta~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1b, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1b$residuals[,1])
anova(fm1b)# three-way interaction significant

#------------------------------------------------

#------------------------------------------------
fm1agr10 <- lme(log(AGR10g)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2,method="REML")
#look at model diagnostics
plot(fm1agr10,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1agr10,log(AGR10g)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1agr10,log(AGR10g)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1agr10, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1agr10$residuals[,1])

anova(fm1agr10)# 3way interaction not significant? - drop?

fm2agr10 <- lme(log(AGR10g)~Treatment+Location+Range+Treatment:Location+Treatment:Range+Location:Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),
                method="ML",data=output2)
anova(fm1agr10,fm2agr10)
plot(allEffects(fm7agr10))
#------------------------------------------------

#------------------------------------------------
fm1rgr10 <- lme(RGR10d~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1rgr10,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1rgr10,RGR10d~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1rgr10,RGR10d~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1rgr10, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1rgr10$residuals[,1])
anova(fm1rgr10)
plot(allEffects(fm1rgr10))     
summary(fm1rgr10)
#------------------------------------------------

#------------------------------------------------
fm1Tmass <- lme(TotMass.f~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1Tmass,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1Tmass,TotMass.f~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1Tmass,TotMass.f~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1Tmass, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1Tmass$residuals[,1])
anova(fm1Tmass)
plot(allEffects(fm1Tmass))     
summary(fm1Tmass)

#------------------------------------------------

#------------------------------------------------

fm1agr10d <- lme(log(AGR10d)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2,method="REML")
#look at model diagnostics
plot(fm1agr10d,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1agr10d,log(AGR10d)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1agr10d,log(AGR10d)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1agr10d, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1agr10$residuals[,1])

anova(fm1agr10d)# 3way interaction significant?

#------------------------------------------------

#------------------------------------------------
fm1rgr10g <- lme(log(RGR10g)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=output2)
#look at model diagnostics
plot(fm1rgr10g,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1rgr10g,log(RGR10g)~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm1rgr10g,log(RGR10g)~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm1rgr10g, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm1rgr10g$residuals[,1])
anova(fm1rgr10g)
plot(allEffects(fm1rgr10g))     
summary(fm1rgr10g)
