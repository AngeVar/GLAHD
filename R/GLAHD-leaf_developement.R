#------------------------------------------------------------------------------------------------------------------------------
#- Reads in and processes the leaf growth dataset based on leaf area measured by pictures and imageJ
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")

#- read data, do simple manipulation
la <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/Leaf photos/GLAHD-LeafArea.csv")
la$Date <- as.Date(la$Date,format="%d/%m/%Y")

la$Taxa <- unlist(strsplit(x=as.character(la$Code),split="-"))[seq(from=1,to=nrow(la)*2,by=2)]
la$Taxa <- factor(la$Taxa,levels=c("BTER","ACAM","BOT","LONG","CTER","DCAM","BRA","PLAT"))
la$Pot <- as.numeric(unlist(strsplit(x=as.character(la$Code),split="-"))[seq(from=2,to=nrow(la)*2,by=2)])
la$Treat <- as.factor(ifelse(la$Pot < 20, "Home","Warmed"))
la$Time <- as.numeric(la$Date-min(la$Date))
                   
linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),
                     Range= c("wide","wide",rep("narrow",6)))


allom <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom2 <- unique(subset(allom,select=c("Location","Species","Taxa")))

la.1 <- merge(la,allom2,by=c("Taxa"))
la.2 <- merge(la.1,linkdf,by=c("Species"))
la <- la.2
la$Combo <- factor(paste(la$Taxa,la$Treat,sep="-"),levels=c("BTER-Home","BTER-Warmed","ACAM-Home","ACAM-Warmed","BOT-Home","BOT-Warmed","LONG-Home","LONG-Warmed"
                                                            ,"CTER-Home","CTER-Warmed","DCAM-Home","DCAM-Warmed","BRA-Home","BRA-Warmed","PLAT-Home","PLAT-Warmed"))
la$Combo2 <- factor(paste(la$Location,la$Range,la$Treat,sep="-"))                

#- plot
plotBy(Area~Date|Code,data=la,legend=F,type="b")
with(la,plot(Area~Date,col=Treat))

#------------------------------------------------------------------
#- average across species, treatments, and dates, then plot
la.m <- summaryBy(Area~Date+Treat+Taxa,data=la,FUN=mean,keep.names=T)
la.m.l <- split(la.m,la.m$Taxa)
la.l <- split(la,la$Taxa)

windows(12,12)
par(mfrow=c(4,2),mar=c(2,4,1,1))
xlims <- c(as.Date("2014-12-20"),as.Date("2015-1-15"))
for (i in 1:length(la.l)){
  dat <- la.l[[i]]
  ylims <- range(dat$Area,na.rm=T)
  plotBy(Area~Date|Code,legend=F,data=subset(dat,Treat=="Home"),axes=F,col="blue",type="b",pch=15,ylim=ylims,xlim=xlims)
  plotBy(Area~Date|Code,legend=F,add=T,data=subset(dat,Treat=="Warmed"),col="red",type="b",pch=15)
  axis(2);box()
  mtext(side=2,dat$Taxa[1],line=-4,las=1)
  axis.Date(side=1,at=seq.Date(from=as.Date("2014-12-20"),to=as.Date("2015-01-15"),by="week"))
}
#------------------------------------------------------------------




#------------------------------------------------------------------
#- try fitting using the Richards function
library(minpack.lm)
la.l <- split(la,la$Combo)

richards <- function(x,a,b,c,d){
  y <- a*(1+exp(b-c*x))^(-1/d)
  return(y)
}

r.out <- richards(x=xvals,a=25,b=3,c=0.5,d=2)
lines(r.out~xvals)

fit1 <- nlsLM(Area~richards(Time,a,b,c,d),data=la.l[[1]],
    start=list(a=25,b=3,c=0.5,d=2))
plot(la.l[[1]]$Area~la.l[[1]]$Time)
newdata <- data.frame(Time=xvals,newdata=newdata)
newdata$Area <- predict(fit1,newdata=newdata)
lines(Area~Time,data=newdata)

#------------------------------------------------------------------




#-------------------------------------------------------------------
#- Try Gompertz exponential. Code from Paine et al. nonlinear fitting paper
la.all <- la[complete.cases(la),]
la.l <- split(la.all,la.all$Combo)
tmp.gomp <- list()
fit.gomp <- list()
out.gomp <- list()
xvals <- seq(min(la$Time),max(la$Time),length=101)
for (i in 1:length(la.l)){
  print(i)
  tmp.gomp[[i]] <- getInitial(Area ~ SSgompertz(Time, Asym, b2, b3), data = la.l[[i]])
  fit.gomp[[i]] <- gnls(Area ~ SSgompertz(Time, Asym, b2, b3), data = la.l[[i]], weights= varExp(form = ~ fitted(.)))
  out.gomp[[i]] <- output.gomp.gnls(fit.gomp[[i]], xvals, CI = T)
}

windows(12,12)
par(mfrow=c(3,1),mar=c(4,4,1,1))
ivalues <- seq(1,15,by=2)
for(a in 1:length(ivalues)){
  i <- ivalues[a]
  toplot <- rbind(la.l[[i]],la.l[[i+1]])
  plotBy(Area~Time|Treat,data=toplot, xlab = "", ylab = "Area (cm2)",type="p",legend=F,pch=16,cex=0.5)
  title(main=toplot$Taxa[1],line=-1)
  
  
  #-- plot biomass
  #home
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$M,  col = "black")     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$M.hi,  col = "black",lty=2)     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$M.lo,  col = "black",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i]]$rates$times, rev(out.gomp[[i]]$rates$times)), y = c(out.gomp[[i]]$rates$M.lo, rev(out.gomp[[i]]$rates$M.hi)), col = alpha("grey",0.5), border = NA)
  
  #warmed
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$M,  col = "red")     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$M.hi,  col = "red",lty=2)     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$M.lo,  col = "red",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i+1]]$rates$times, rev(out.gomp[[i+1]]$rates$times)), y = c(out.gomp[[i+1]]$rates$M.lo, rev(out.gomp[[i+1]]$rates$M.hi)), col = alpha("red",0.5), border = NA)
  #---
  
  #-- plot AGR
  plot(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$AGR,type="n",ylim=c(0,max(out.gomp[[i+1]]$rates$AGR.hi+1)),ylab="AGR",xlab="")
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$AGR,  col = "black")     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$AGR.lo,  col = "black",lty=2)     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$AGR.hi,  col = "black",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i]]$rates$times, rev(out.gomp[[i]]$rates$times)), y = c(out.gomp[[i]]$rates$AGR.lo, rev(out.gomp[[i]]$rates$AGR.hi)), col = alpha("grey",0.5), border = NA)

  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$AGR,  col = "red")     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$AGR.lo,  col = "red",lty=2)     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$AGR.hi,  col = "red",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i+1]]$rates$times, rev(out.gomp[[i+1]]$rates$times)), y = c(out.gomp[[i+1]]$rates$AGR.lo, rev(out.gomp[[i+1]]$rates$AGR.hi)), col = alpha("red",0.5), border = NA)
  #-- 
  
  
  #-- plot RGR
  plot(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$RGRt,type="n",ylim=c(0,max(out.gomp[[i+1]]$rates$RGRt.hi)),ylab="RGR",xlab="")
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$RGRt,  col = "black")     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$RGRt.lo,  col = "black",lty=2)     # gompertz
  lines(out.gomp[[i]]$rates$times,  out.gomp[[i]]$rates$RGRt.hi,  col = "black",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i]]$rates$times, rev(out.gomp[[i]]$rates$times)), y = c(out.gomp[[i]]$rates$RGRt.lo, rev(out.gomp[[i]]$rates$RGRt.hi)), col = alpha("grey",0.5), border = NA)
  
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$RGRt,  col = "red")     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$RGRt.lo,  col = "red",lty=2)     # gompertz
  lines(out.gomp[[i+1]]$rates$times,  out.gomp[[i+1]]$rates$RGRt.hi,  col = "red",lty=2)     # gompertz
  polygon(x = c(out.gomp[[i+1]]$rates$times, rev(out.gomp[[i+1]]$rates$times)), y = c(out.gomp[[i+1]]$rates$RGRt.lo, rev(out.gomp[[i+1]]$rates$RGRt.hi)), col = alpha("red",0.5), border = NA)
  #-- 
  
  
  
  
}










#------------------------------------------------------------------
# playing around to find a function to use
nlsfit <- nls(Area ~  SSlogis(Time, Asym, xmid, scal),
                data=subset(la,Code=="PLAT-26"))
tofit <- subset(la,Code=="BTER-6")
nlsfit <- nls(Area ~  SSlogis(Time, Asym, xmid, scal),
              data=tofit)


tmp.gomp<- getInitial(Area ~ SSgompertz(Time, Asym, b2, b3), data=tofit)
fit.gomp <- gnls(Area ~ SSgompertz(Time, Asym, b2, b3), data = tofit, weights= varExp(form = ~ fitted(.)))
out.gomp <- output.gomp.gnls(fit.gomp, xvals, CI = T)
# Make dataframe with X variable that we wish to predict Y values for.
# Make sure it has the same name as in the dataframe we used to fit the model!
newdat <- data.frame(Time=seq(0,21, length=101))


# Predict from the fitted model for the new dataframe.
newdat$Areapred <- predict(nlsfit, newdata=newdat)


# Add a line.
with(subset(la,Code=="BTER-6"), plot(Time, Area, xlim=c(0, 30), ylim=c(0,100)))
with(newdat, lines(Time, Areapred, col="blue"))
lines(out.gomp$rates$times,  out.gomp$rates$M,  col = "black")     # gompertz
#------------------------------------------------------------------













#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#- fit each leaf separately, analyze parameter outputs
la2 <- subset(la,Code!="BTER-25")
la2$Code <- factor(la2$Code)
la.l <- split(la2,la2$Code)
fit <- coefs <- list()

#make a big pdf with all of the fits for each curve.
pdf(file="./output/leaf_dev_fits.pdf")
for (i in 1:length(la.l)){
  dat <- la.l[[i]]
  fit[[i]] <- nls(Area ~  SSlogis(Time, Asym, xmid, scal),data=dat)
  coefs[[i]] <- coef(fit[[i]])
  #coefs[[i]]$Code <- as.character(dat$Code[[1]])
  
  newdat <- expand.grid(Code=dat$Code[1], Time=0:21)
  newdat$Areapred <- predict(fit[[i]],newdat,level=0)
  
  plot(Area~Time,data=dat)
  lines(Areapred~Time,data=newdat)
  title(main=dat$Code[1])
}
dev.off()

fits <- as.data.frame(do.call(rbind,coefs))
fits$Code <- fits$origin <- fits$range <- NA
for (i in 1:nrow(fits)){
  fits$Code[i] <- as.character(la.l[[i]]$Code[1])
}
fits$Taxa <- unlist(strsplit(x=as.character(fits$Code),split="-"))[seq(from=1,to=nrow(fits)*2,by=2)]
fits$Taxa <- factor(fits$Taxa,levels=c("BTER","ACAM","BOT","LONG","CTER","DCAM","BRA","PLAT"))
fits$Pot <- as.numeric(unlist(strsplit(x=as.character(fits$Code),split="-"))[seq(from=2,to=nrow(fits)*2,by=2)])
fits$Treat <- as.factor(ifelse(fits$Pot < 20, "Home","Warmed"))
narrow <- c("BOT","LONG","BRA","PLAT")
wide <- c("BTER","ACAM","CTER","DCAM")
north <- c("CTER","DCAM","BRA","PLAT")
south <- c("BTER","ACAM","BOT","LONG")
for (i in 1:nrow(fits)){
  fits$Code[i] <- as.character(la.l[[i]]$Code[1])
  fits$range[i] <- ifelse(fits$Taxa[i] %in% narrow,"Narrow","Wide")
  fits$origin[i] <- ifelse(fits$Taxa[i] %in% north,"North","South")
  
}
fits$range <- as.factor(fits$range)
fits$origin <- as.factor(fits$origin)
# warming reduces the asymptote overall, but this varies strongly by taxa
lm.asym <- lm(Asym~Taxa*Treat,data=fits)


library(visreg)
anova(lm.asym)
plot(allEffects(lm.asym))
effect("Treat",lm.asym)
effect(Treat:Taxa,lm.asym)
windows(20,12)
visreg(lm.asym, "Taxa", by="Treat", overlay=TRUE,
       points=list(col=c("blue","red")),
       fill=list(col=c("blue","red")),
       line=list(col=c("blue","red")))

# warming reduces xmid, interaction with taxa is marginally significant
lm.xmid <- lm(xmid~Taxa*Treat,data=fits)
anova(lm.xmid)
plot(allEffects(lm.xmid))
effect("Treat",lm.xmid)
visreg(lm.xmid, "Taxa", by="Treat", overlay=TRUE,
       points=list(col=c("blue","red")),
       fill=list(col=c("blue","red")),
       line=list(col=c("blue","red")))

# warming reduces scale, no interaction with taxa
lm.scal <- lm(scal~Taxa*Treat,data=fits)
lm.scal2 <- lm(scal~range*origin*Treat,data=fits) #warming only reduces scale in the south
anova(lm.scal)
plot(allEffects(lm.scal))
effect("Treat",lm.scal)
visreg(lm.scal, "Taxa", by="Treat", overlay=TRUE,
       points=list(col=c("blue","red")),
       fill=list(col=c("blue","red")),
       line=list(col=c("blue","red")))


#- the big effect of warming is to reduce the scale parameter, particularly in the south
#- this means that the rate of leaf development was increased by warming in the south
windows(12,12)
Asym <- c(100,100)
xmid <- c(9,8)
scal <- c(3,2) 
xvals <- seq(0,21,length=101)
pred1 <- Asym[1]/(1+exp((xmid[1]-xvals)/scal[1]))
pred2 <- Asym[2]/(1+exp((xmid[2]-xvals)/scal[2]))
plot(pred1~xvals,type="l",col="black",cex.lab=1.5,ylab="Leaf area",xlab="Time in days")
lines(pred2~xvals,type="l",col='red')
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------













#------------------------------------------------------------------
#-- fit non-linear mixed effects model
la.all <- la[complete.cases(la),]
fitnlme1 <- nlme(Area ~ SSlogis(Time, Asym, xmid, scal),
                 fixed=list(Asym ~ Taxa*Treat, xmid ~ Taxa, scal ~ Taxa),
                 random = Asym ~ 1 | Code,
                 start=list(fixed=c(Asym=rep(60,16),xmid=rep(4,8),scal=rep(2,8))),
                 data=la.all)


# Fixed effects predictions
# Make dataframe with all combinations of Diet and Time.
newdat <- expand.grid(Taxa=levels(la$Taxa), Time=0:21,Treat=levels(la$Treat))
newdat$Areapred <- predict(fitnlme1,newdat,level=0)

# Plot. I use jitter so not all data cover each other.
windows(12,12)
palette(c("blue","red","lightgreen","darkorange","violet","yellow","darkgreen","black"))
with(subset(la,Treat=="Home"), plot(jitter(Time), Area, pch=19, col=Taxa, main="fitnlme1",xlab="Days"))
with(subset(la,Treat=="Warmed"), points(jitter(Time), Area, pch=1, col=Taxa, main="fitnlme1",xlab="Days"))

plotBy(Areapred ~ Time|Taxa, data=subset(newdat,Treat=="Home"), type='l', add=TRUE, lwd=2,legend=F)
plotBy(Areapred ~ Time|Taxa, data=subset(newdat,Treat=="Warmed"),lty=2, type='l', add=TRUE, lwd=2,legend=F)
legend("topleft",legend=levels(la$Taxa),col=as.factor(levels(la$Taxa)),pch=19)
