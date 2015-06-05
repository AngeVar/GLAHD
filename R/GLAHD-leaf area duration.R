#------------------------------------------------------------------------------------------------------------------------------
#- Analyzes leaf area duration for the GLAHD experiment. This is the integral of total plant leaf area over time.
#------------------------------------------------------------------------------------------------------------------------------

#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")


#- read in the data, calculate the number of elapsed days
dat <- return_size_mass()
dat$Time <- as.numeric(dat$Date-min(dat$Date-7))

#- remove data with fewer than 9 observations through time (the maximum possible)
obs <- unname(table(dat$Code)) # get the frequency of observations for each pot
names <- names(table(dat$Code))# get the associated name of each pot
keeps <- names[which(obs>=9)]    # return a vector of pot names with more than n observations
dat2 <- subset(dat,Code %in% keeps) # subset dataframe
dat2$Code <- factor(dat2$Code)




#----------------------------------
#-- average across groups and plot
#- average across four groups (N vs. S, Wide vs. Narrow)
dat.m <- summaryBy(leafArea+TotMass~Time+Date+Location+Range+Treatment,FUN=c(mean,standard.error),data=dat2)
dat.m$combo <- factor(paste(dat.m$Location,dat.m$Range,sep="_"),levels=c("S_wide","S_narrow","N_wide","N_narrow"))


#- make 4-panel plot of leaf area over time 
dat.m.l <- split(dat.m,dat.m$combo)
windows(20,20);par(mfrow=c(2,2),mar=c(4,5,1,1),oma=c(4,3,1,1))
labels <- c("South Wide","South Narrow","North Wide","North Narrow")
for (i in 1:length(dat.m.l)){
  toplot <- dat.m.l[[i]]
  
  plotBy(leafArea.mean~Time|Treatment,type="b",pch=16,legend=F,data=toplot,ylim=c(0,6000),cex=1.4,cex.lab=1.3,axes=F,
         ylab="",xlab="",
         panel.first=adderrorbars(x=toplot$Time,y=toplot$leafArea.mean,SE=toplot$leafArea.standard.error,direction="updown"))
  if (i==1) legend("topleft",legend=c("Home","Warmed"),pch=16,col=c("blue","red"),cex=1.4)
  title(main=labels[i],line=0.2,xpd=F)
  magaxis(side=c(1,2,4),labels=c(1,1,0),las=1,box=T)
}
title(ylab=expression(Total~leaf~area~(cm^2)),outer=T,line=-1,cex.lab=2)
title(xlab=expression(Time~(days)),outer=T,line=1,cex.lab=2)
#----------------------------------




#----------------------------------
#- loop over each pot, fit a polynomial function to the leaf area vs. time relationship, and integrate it. Integrated leaf area over time is calle LAD.i by me.

polyfn <- function(x,params){
  y <- params[1]+params[2]*x+params[3]*x^2+params[4]*x^3+params[5]*x^4
  return(y)
}

dat.l <- split(dat2,dat2$Code)
integrations <- list()
LAD <- c()
Code <- c()
for (i in 1:length(dat.l)){
  lm1 <- lm(leafArea~Time+I(Time^2)+I(Time^3)+I(Time^4),data=dat.l[[i]])  # fit fifth order polynomial to data
#  newdata <- data.frame(Time=seq(min(dat2$Time),max(dat2$Time),length=101))
#   newdata$pred <- polyfn(x=newdata$Time,params=unname(coef(lm1)))
#   plot(leafArea~Time,data=dat.l[[i]])
#   lines(pred~Time,data=newdata)
  
  LAD[i] <- integrate(f=polyfn,lower=7,upper=66,params=unname(coef(lm1)))[[1]] # integrate the 5th order poly from 7 to 66 days to get LAD
  Code[i] <- as.character(dat.l[[i]]$Code[1])
#LAD <- integrations[[i]][[i]]
}
df <- data.frame(LAD.i=LAD,Code) # put LAD.i and Code together in a dataframe, for easier merging later
#----------------------------------




#----------------------------------
# statistical anlaysis of LAD

#- calculate LAD as the sum of leafArea measured over time for each pot. This is a poor-man's integral.
dat.lad <- summaryBy(leafArea~Species+Treatment+Location+Taxa+Pot+Code+Range,FUN=sum,keep.names=T,data=dat2)
names(dat.lad)[which(names(dat.lad)=="leafArea")] <- "LAD" # rename the sum of leaf area over time to be leaf area duration
dat.lad <- merge(dat.lad,df,by="Code") # get LAD.i (the real integrated leaf area over time)
#- looks like an interaction between warming and taxa (or warming x location x range?)
windows(20,12);par(mar=c(8,7,1,1))
boxplot(LAD.i~Treatment+Taxa,data=dat.lad,las=2,ylab=c("Leaf area duration"),col=c("grey","red"))



#------
#- fit the statistical model with lmer in lme4 package
library(arm)
library(lme4)
library(LMERConvenienceFunctions)
library(lmerTest)
#fit null model, with no fixed effects, just random effects. Then fit full model
m0 <- lmer(LAD.i~1+(1|Taxa),data=dat.lad)#null model
mF <- lmer(LAD.i~Treatment*Location*Range+(1|Taxa),data=dat.lad,REML=T)#full model
mcp.fnc(mF) # resid-pred plot looks good. qq-plot looks okay, but could be better

#- remove some outliers and try again
dat.lad2 <- romr.fnc(mF, dat.lad, trim = 2.5)
dat.lad2$n.removed
dat.lad2$percent.removed
dat.lad2<-dat.lad2$data

mF2 <- lmer(LAD.i~Treatment*Location*Range+(1|Taxa),data=dat.lad2,REML=T)#full model
mcp.fnc(mF2) # resid-pred plot looks good. qq-plot looks okay, but could be better

anova(mF2, ddf = "Kenward-Roger") #- three way interaction between Treatment, Range, and Locaiton is highly significant with high F-value
anova(mF2, ddf = "lme4") #- as per Doug Bates' theoretical concerns, this doesn't calcualte ddf, and thus can't get p-values etc.


mF.back <- step(mF2) # backwards elimation of non-significant effects. print object and note that all terms were kept.
plot(mF.back) # cool plot, but actually not meaningful
#get marginal and conditional r2 values
source("W:/WorkingData/GHS39/GLAHD/Share/R/rsquared_glmm.R")
r.squared(mF2) # marginal r2 (fixed effects only) is 0.51, conditional r2 (fixed and random effects) is 0.87


#-- attempting to do post-hoc tests, but I don't really understand them. I don't know what these mean.
mF2.ph <- mcposthoc.fnc(model = mF2, var = list(ph1 = c("Treatment")))
smry <- summary(object = mF2.ph, term = "TreatmentWarmed:LocationS:Rangewide")
smry <- summary(object = mF2.ph, term = "TreatmentHome:LocationS:Rangewide")

### the take home message is that warming increased leaf area durating for wide taxa in the south relative to narrow taxa in the south
### this interaction was not observed in the north (in fact, it was reversed)

#----------------------------------