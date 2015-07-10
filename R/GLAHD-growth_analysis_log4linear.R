#-------------------------------------------------------------------------------------
#- This script provides an alternative curve fit for growth analysis.
#- This works quickly and simply because it uses a LINEAR model rather than a nonlinear model.
#- This assumes that a plot of log(Mass)~Time can be adequately fit with a polynomial of some order (e.g., 2nd, 3rd).
#- This assumption seems to be correct in this case, for our data.
#- I don't know whether it's "better" than a power model, but it is much more convenient!
#-------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")

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


#-------------------------------------------------------------------------------------
#-- function to interpret log-polynomial curve fits that are FOURTH ORDER
output.log_4order <- function(X, Y, params, times,Code){
  a <- params[1]; b <- params[2]; c <- params[3]; d <- params[4]; e<- params[5]
  #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
  fitted <- exp(a + b*X + c*X^2 + d*X^3 +e*X^4)
  resid  <- Y - fitted
  data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
  #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
  eq <- bquote(a+bx+cx^2+dx^3+ex^4)
  mss  <- sum((fitted - mean(fitted))^2)
  rss  <- sum(resid^2)
  R2   <- mss/(mss + rss)
  rmse <- sqrt(rss)
  N <- length(X)
  logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
  AIC  <- -2 * logLik  + 2 * 3 # four parameters
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  #temp <- a + b*X + c*X^2 
  rates = data.frame(
    times = times,
    M    =  exp(a+b*times+c*times^2+d*times^3+e*times^4))                  # from Hunt- Plant Growth Curves
  rates$AGR  =  with(rates,M*(b+2*c*times+3*d*times^2+4*e*times^3))            # from Hunt- Plant Growth Curves
  rates$RGR <- rates$AGR/rates$M
  rates$Code <- Code
  out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
  return(out)
}
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#- Calculate interval-based RGR for each plant
dat3 <- dat3[with(dat3,order(Code,Date)),] #- make sure the observations for each plant are in order
growth.l <- split(dat3,dat3$Code)          #- split into a list for each plant
for (i in 1:length(growth.l)){
  #- use diff() to calculate RGR. Note that I've added a trailing "NA", or else the vector would be too short to fit
  #-   the data frame.
  growth.l[[i]]$RGR <- c(diff(growth.l[[i]]$lnTotMass),NA)/c(unname(diff(growth.l[[i]]$Date)),NA) 
  growth.l[[i]]$dDiameter <- c(diff(growth.l[[i]]$Diameter),NA) 
  
}
dat4 <- do.call(rbind,growth.l)
windows(30,40)
plotBy(RGR~jitter(as.numeric(Date))|Code,type="b",data=dat4,legend=F);abline(h=0) 
smoothScatter(x=dat4$Time,y=dat4$RGR,ylab="RGR",xlab="Time since treatment began (days)",cex.lab=1.5,
              xlim=c(0,60));abline(h=0)
#- it looks like a third-order polynomial may fit the data best, which would produce
#   and RGR that changes over time according to a second-order polynomial


#-------------------------------------------------------------------------------------
#- fit log-4polynomial growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
log4fits <- output <- data.out <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  log4fits[[i]] <- lm(log(TotMass)~Time+I(Time^2)+I(Time^3)+I(Time^4),data=tofit) #- fit log-polynomial
  
  #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
  output[[i]] <- output.log_4order(X=tofit$Time,Y=tofit$TotMass,
                                   params=unname(coef(log4fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$rates
  data.out[[i]] <- output.log_4order(X=tofit$Time,Y=tofit$TotMass,
                                     params=unname(coef(log4fits[[i]])),times=seq(1:60),Code=tofit$Code[1])$data
  data.out[[i]]$Code <- tofit$Code[1]
}
#- put rates dataframe together. This has the timecourse of mass, AGR, and RGR for each plant
rates.df <-do.call(rbind,output)
data.df <- do.call(rbind,data.out)
data.df$fit_type <- "log-4polynomial"
params <- as.data.frame(do.call(rbind,lapply(log4fits,coef))) #- get polynomial parameters. Perhaps not that useful alone
params$Code <- unique(rates.df$Code) # add the code variable to the parameter estimate
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#- plot the development for each treatment. Note how effect on RGR is particularly time-dependant
rates.df2 <- merge(rates.df,subset(dat2,Date==as.Date("2014-11-17")),by="Code")[,c(1:6,8:10,17,19)]#,16)] # merge with dat2 to get Treatment, Location, etc

rates.trt <- summaryBy(M+RGR+AGR~times+Treatment+Location+Range,data=rates.df2,FUN=mean,keep.names=T)
rates.trt$combotrt <- as.factor(paste(rates.trt$Location,rates.trt$Range,rates.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(RGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)
plotBy(AGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
plotBy(M~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
plotBy(M~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#- plot RGR via the segment method for each taxa-combination. Overlay model outputs.
dat4$combotrt <- as.factor(paste(dat4$Location,dat4$Range,dat4$Treatment,sep="_"))
dat4.l <- split(dat4,dat4$combotrt)
names(rates.trt)[1] <- "Time"
rates.trt.l <- split(rates.trt,rates.trt$combotrt)
windows(40,40);par(mfrow=c(4,2),mar=c(0.2,0.2,0.2,0.2),oma=c(7,7,3,3))
ylims=c(0,0.3)
for (i in 1:length(dat4.l)){
  smoothScatter(x=dat4.l[[i]]$Time,y=dat4.l[[i]]$RGR,ylab="",xlab="",ylim=ylims,axes=F,
                xlim=c(0,60));abline(h=0)
  lines(rates.trt.l[[i]]$RGR~rates.trt.l[[i]]$Time,lwd=3)
  legend("top",legend=dat4.l[[i]]$combotrt[1],bty="n")
  if(i%%2==1) axis(2,las=1)
  if(i%%2==0) axis(4,las=1)
  if(i>6) axis(1)
}
mtext(side=1,"Time (days after treatment began)",outer=T,line=4,cex=1.5,xpd=NA)
mtext(side=2,"RGR",outer=T,line=3,cex=1.5)
#----------------------- --------------------------------------------------------------



#-------------------------------------------------------------------------------------
#- plot mass and RGR over time for each taxa-treatment combination, overlay data and model output
dat4$RGR_time <- dat4$Time+5
dat4.l.taxa <- split(dat4,dat4$Taxa)
rates.taxa <- summaryBy(M+RGR+AGR~times+Treatment+Location+Range+Taxa,data=rates.df2,FUN=mean,keep.names=T)
rates.taxa.l <- split(rates.taxa,rates.taxa$Taxa)
rates.l <- split(rates.df2,rates.df2$Taxa)

# loop over each taxa, plot a set of four figures, export to a pdf
pdf(file="C:/Repos/GLAHD/Output/RGR_M_Taxa_filtered_log4.pdf")

layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
       widths=rep(1,4), heights=rep(1,4))
for (i in 1:length(dat4.l.taxa)){
  ylims=c(0,max(dat4.l.taxa[[i]]$TotMass+5))
  #- plot home mass over time
  plot(TotMass~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
  plotBy(M~times|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(M~times,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
  
  #- plot warmed mass over time
  plot(TotMass~Time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Warmed",sep='_'))
  plotBy(M~times|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(M~times,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
  
  #- plot home RGR over time
  ylims=c(0,max(dat4.l.taxa[[i]]$RGR+0.1,na.rm=T))
  
  plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Home"),ylim=ylims,col="blue",xlim=c(0,60),
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Home")$Taxa[1],"Home",sep='_'))
  plotBy(RGR~times|Code,data=subset(rates.l[[i]],Treatment=="Home"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(RGR~times,data=subset(rates.taxa.l[[i]],Treatment=="Home"),lwd=3)
  
  #- plot warmed RGR over time
  plot(RGR~RGR_time,data=subset(dat4.l.taxa[[i]],Treatment=="Warmed"),ylim=ylims,col="red",xlim=c(0,60),
       main=paste(subset(dat4.l.taxa[[i]],Treatment=="Warmed")$Taxa[1],"Home",sep='_'))
  plotBy(RGR~times|Code,data=subset(rates.l[[i]],Treatment=="Warmed"),lwd=0.5,add=T,type="l",legend=F,col="grey")
  lines(RGR~times,data=subset(rates.taxa.l[[i]],Treatment=="Warmed"),lwd=3)
  
}
dev.off()

#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#- Do analysis suggested by Mark, where interval-based RGR estimates are compared to model estimates
#-  The basic question: is a 4th order polynomial model sufficient to describe the data,
#-      particularly the peaked RGR response to warming during the second growth interval?


#- the interval-based RGR estimate for the segment in question
rgr.seg2 <- subset(dat4,Date==as.Date("2014-11-17"))
names(rgr.seg2)[18] <- "RGR.seg"

#- get the modeled RGR for each day of the growth interval in question for each plant
#remember that time "7" is 2014-11-07
table(dat2$Date) # second growth interval is from 2014-11-17 to 2014-11-26, which is "times" 11 to 20
rgr.mod.seg2 <- subset(rates.df2,times>=18 & times <=27)
rgr.mod.seg2.plant <- summaryBy(M+AGR+RGR~Code+Species+Treatment+Location+Taxa+Range,data=rgr.mod.seg2,keep.names=T)

#- put them together
rgr2 <- merge(rgr.seg2,rgr.mod.seg2.plant,by=c("Code","Species","Treatment","Location","Taxa","Range"))

#- plot both estimates of RGR for each taxa

#plot RGR seg, then RGR modeled
ylims=c(0,0.2)
xlims=c(0,35)
colors <- c("blue","red")
windows(10,10);par(mfrow=c(2,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)

#RGR segmented
boxplot(RGR.seg~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[seg]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=F,las=2,cex.axis=1.5)
abline(v=16.4)

#RGR modeled
boxplot(RGR~Treatment+Taxa,data=rgr2,ylim=ylims,xlim=xlims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[mod]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(rgr2$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)


#- by-plot
windows(20,20)
plotBy(RGR~RGR.seg|Location,data=rgr2,xlab="RGR-seg",ylab="RGR-mod",xlim=c(0,0.22),ylim=c(0,0.22),cex.lab=2)
abline(0,1)

#- compare statistical analyes
library(nlme)
rgr2$Sp_RS_EN <- as.factor(with(rgr2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rgr2$Prov_Sp_EN <- as.factor(with(rgr2,paste(Taxa,Species)))
fm.rgr.mod <- lme(RGR~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
plot(fm.rgr.mod,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.rgr.mod,RGR~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.rgr.mod,RGR~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.rgr.mod, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.rgr.mod$residuals[,1])
anova(fm.rgr.mod)  

fm.rgr.seg <- lme(RGR.seg~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rgr2)
plot(fm.rgr.seg,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.rgr.seg,RGR.seg~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.rgr.seg,RGR.seg~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.rgr.seg, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.rgr.seg$residuals[,1])
anova(fm.rgr.seg)  

#-------------------------------------------------------------------------------------------------------

#- plot the first four segements
rgr.seg14 <- subset(dat4,Date <=as.Date("2014-12-01"))
names(rgr.seg14)[17] <- "RGR.seg"
windows(20,20)
plotBy(RGR.seg~Date|Code,data=rgr.seg14,type="b",legend=F)

#-------------------------------------------------------------------------------------------------------
#- compare modeled and segmented data again, but for the second and third growth intervals
#- the interval-based RGR estimate for the segment in question
rgr.seg23 <- subset(dat4,Date>=as.Date("2014-11-17") & Date <=as.Date("2014-12-01"))
names(rgr.seg23)[18] <- "RGR.seg"

#- get the modeled RGR for each day of the growth interval in question for each plant
#remember to adjust time, time "7" is 2014-11-07 if first date is labelled 7
table(dat2$Date) # second growth interval is from 2014-11-17 to 2014-11-26, which is "times" 11 to 20
rgr.mod.seg23 <- subset(rates.df2,times >=18 & times <= 32)
rgr.mod.seg23.plant <- summaryBy(M+AGR+RGR~Code+Species+Treatment+Location+Taxa+Range,data=rgr.mod.seg23,keep.names=T)

#- put them together
rgr23 <- merge(rgr.seg23,rgr.mod.seg23.plant,by=c("Code","Species","Treatment","Location","Taxa","Range"))

#plot RGR seg, then RGR modeled
ylims=c(0,0.2)
colors <- c("blue","red")
windows(20,20);par(mfrow=c(2,1),mar=c(0.5,6,0.5,3),oma=c(6,0,0,0),cex.axis=1.2)

#RGR segmented
boxplot(RGR.seg~Treatment+Taxa,data=rgr23,ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[seg]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=F,las=2,cex.axis=1.5)
abline(v=16.4)

#RGR modeled
boxplot(RGR~Treatment+Taxa,data=rgr23,ylim=ylims,
        axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(RGR[mod]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(rgr12$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)


#- by-plot
windows(20,20);par(mar=c(6,6,0,0))
plotBy(RGR~RGR.seg|Location,data=rgr23,xlab="RGR-seg",ylab="RGR-mod",xlim=c(0,0.22),ylim=c(0,0.22),cex.lab=2)
abline(0,1)

#-------------------------------------------------------------------------------------------------------
rates.df2 <- merge(rates.df,subset(dat2,Date==as.Date("2014-11-17")),by="Code")[,c(1:6,8:10,17,19)]#,16)] # merge with dat2 to get Treatment, Location, etc
rates.df2$LAR <- rates.df2$leafArea/rates.df2$M
rates.trt <- summaryBy(LAR+M+RGR+AGR~times+Treatment+Location+Range,data=rates.df2,FUN=mean,keep.names=T)
rates.trt$combotrt <- as.factor(paste(rates.trt$Location,rates.trt$Range,rates.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(1,1))
plotBy(RGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)
plotBy(AGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
plotBy(M~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft",pch=15)
plotBy(LAR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright",pch=15)

plotBy(log(LAR)~RGR|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="bottomright",pch=15)
##-------------------------------------------------------------------------------------------------------
#Biomass enhancement ratio (mean of all data)
ber<- summaryBy(M~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]/bio[,5]
BER<- bio[,c(6:11)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

#Relative enhancement of biomass
windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,2), xlim=c(0,67))
mtext(text=expression(M[W]:M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,2),yaxt='n', xlim=c(0,67))
abline(h=1)
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Biomass Enhancement Ratio",side=3,outer=T,cex=1,adj=0.5,line=1)

#Absolute enhancement of biomass
ber<- summaryBy(M~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,10]-bio[,5]
BER<- bio[,c(6:11)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)


windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-5,20), xlim=c(0,67))
mtext(text=expression(M[W]-M[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-5,20),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of M",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

#Relative enhancement of RGR
ber<- summaryBy(M+RGR~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,12]/bio[,6]
BER<- bio[,c(1:6,13)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,2), xlim=c(0,67))
mtext(text=expression(RGR[W]:RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,2),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=1)

#Absolute enhancement of RGR
ber<- summaryBy(M+RGR~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,12]-bio[,6]
BER<- bio[,c(1:6,13)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-0.02,0.04), xlim=c(0,67))
mtext(text=expression(RGR[W]-RGR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-0.02,0.04),yaxt='n', xlim=c(0,67))
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of RGR",side=3,outer=T,cex=1,adj=0.5,line=1)
abline(h=0)

#Relative enhancement of LAR
ber<- summaryBy(M+LAR~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,12]/bio[,6]
BER<- bio[,c(1:6,13)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(0.5,1.5), xlim=c(0,67))
mtext(text=expression(LAR[W]:LAR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=1)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(0.5,1.5),yaxt='n', xlim=c(0,67))
abline(h=1)
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Relative Enhancement of LAR",side=3,outer=T,cex=1,adj=0.5,line=1)

#Absolute enhancement of LAR
ber<- summaryBy(M+LAR~times+Range+Location+Treatment, data=rates.df2, FUN=mean, keep.names=T) 
berh <- subset(ber, Treatment == "Home")
berw <- subset (ber, Treatment == "Warmed")  
bio <- cbind(berh,berw)
bio$ber<- bio[,12]-bio[,6]
BER<- bio[,c(1:6,13)]
bersum<- summaryBy(ber~Range+Location+times, data=BER, FUN=mean, keep.names=T)

windows(10,6);par(mfrow=c(1,2),mar=c(2,0,2,0),oma=c(5,9,3,5),cex.axis=1)

plotBy(ber~times|Range,data=subset(BER, Location =="S"),col=c("black","red","blue","orange"),
       legend=F,type="o", main="South", ylim=c(-25,50), xlim=c(0,67))
mtext(text=expression(LAR[W]-LAR[H]),side=2,outer=T,cex=1,adj=0.5,line=3)
abline(h=0)
plotBy(ber~times|Range,data=subset(BER,Location =="N"),col=c("black","red","blue","orange"),
       legendwhere="topright",type="o", main = "North", ylim=c(-25,50),yaxt='n', xlim=c(0,67))
abline(h=0)
mtext(text="Time",side=1,outer=T,cex=1,adj=0.5,line=1)
mtext(text="Absolute Enhancement of LAR",side=3,outer=T,cex=1,adj=0.5,line=1)
##---------------------------------------------------------------------------------------------------



#not modified yet
#Extract data from day 17 (the second harvest) i.e. Day 14 of the experiment

dat4<- subset(DEfits.all,Date==as.Date("2014-11-17"))

dat4$Location <- factor(dat4$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
dat4$Sp_RS_EN <- as.factor(with(dat4,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
dat4$Prov_Sp_EN <- as.factor(with(dat4,paste(Taxa,Species)))
dat4$Sp_Loc_EN <- as.factor(with(dat4,paste(Species,Location)))

fm1r.l <- lme(log(r)~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=dat4)
plot(fm1r.l,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1r.l,log(r)~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1r.l,log(r)~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1r.l, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1r.l$residuals[,1])
anova(fm1r.l)                 #3-way interaction almost there

plot(allEffects(fm1r.l))      #- try to make sense of Treatment:Location:Range interaction


