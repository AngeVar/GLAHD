#supplementary figure allometry
allom <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
allom$logd2h <- log10(allom$d2h)
allom$logTM <- log10(allom$Totmass)
allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                ifelse(allom$Pot>=40,"Pre","Warmed")))
# allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
#                                          "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))

allom$Taxa <- factor(allom$Taxa,levels=c("BRA","CTER","DTER","PEL","ETER","DCAM","PLAT",
                                         "ECAM","FCAM","BOT","ATER","BTER","LONG","ACAM","BCAM","SMIT","CCAM"))
allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

allom.l <-split(allom,allom$Taxa)

combostrop<- c("BRA",NA,"CTER","DTER","PEL",NA,"ETER","DCAM","PLAT",NA,
               "ECAM","FCAM")
combostemp <- c("BOT",NA,"ATER","BTER","LONG",NA,"ACAM","BCAM","SMIT",NA,"CCAM")
combos<-c(combostrop,rep(NA,4),combostemp)
figIDs<- c("a",NA,"g","h","b", NA, "i","j","c",NA, "k","l",NA,NA,NA,NA,"d",NA,"m","n","e",NA,"o","p","f",NA,"q")


ylims <- c(-1.5,3.5)
xlims <- c(-1.5,3)
windows(8.27,11.69)
par(mfrow=c(8,8),mar=c(0,0,0,0),oma=c(8,8,4,4))
layout(matrix(c(1:28), nrow=7, ncol=4,byrow=T),
       heights=c(1,1,1,0.3,1,1,1),
       widths=c(1,0.3,1,1,1))


for (i in 1:length(combos)){
  dat2 <- subset(allom,Taxa==as.character(combos[i]))
  plotBy(logTM~logd2h|Treat, data= dat2,col=c("black","darkgrey","red"),type="p",pch=19,cex=1.4,
         xlim=xlims,ylim=ylims,axes=FALSE,xlab="",ylab="", legend=F)
  if(nrow(dat2)>=1){
    legend("topright",legend=figIDs[i], bty='n', cex=1.5)}
  par(new=T)
  if(nrow(dat2)>=1){
    predline(lm(logTM~logd2h,data=dat2),col=alpha("black",alpha=0.4), lwd=1)}
  #first plot
  ifelse(dat2$Taxa %in% c("BRA","PEL","BOT","LONG"),
         magaxis(side=c(1,2),labels=c(0,1),frame.plot=T,las=1,cex.axis=1.2),
         ifelse(dat2$Taxa %in% c("CCAM","BCAM","ECAM","FCAM"),
                magaxis(side=c(1,2),labels=c(1,0),frame.plot=T,las=1,cex.axis=1.2),
                ifelse(dat2$Taxa %in% c("SMIT","PLAT"),
                       magaxis(side=c(1,2),labels=c(1,1),frame.plot=T,las=1,cex.axis=1.2),
                       magaxis(side=c(1,2),labels=c(0,0),frame.plot=T,las=1,cex.axis=1.2))))
  
  legend("topleft",legend=dat2$Taxa[1],bty="n",cex=1.5)
  
}
mtext(expression(log[10]~(diameter^2~"*"~height)),side=2,line=3.5,outer=T,cex=1.4)
mtext(expression(log[10]~(Total~mass~(g))),side=1,line=3.5,outer=T,cex=1.3)
mtext("Narrow",side=3,line=1,outer=T,cex=1.2, adj=0.08)
mtext("Wide",side=3,line=1,outer=T,cex=1.2, adj=0.72)
text(9.2,y=23,labels="Tropical", xpd=NA, srt=-90, pos=2, cex=1.9)
text(9.2,y=4.5,labels="Temperate", xpd=NA, srt=-90, pos=2, cex=1.9)

legend(x=3.2,y=2.4,legend=c(expression(Warmed~(+3.5~degree~C)),"Home", "Pre-treatment"),cex=1.4,
       xpd=NA,fill=c("red","black","darkgrey"), border="black",
       bty="n")


# 
# #analysis
# library(smatr)
# slope.com(y=allom.2$logTM,x=allom.2$logd2h,groups=allom.2$Treat) # returns a p-value of 0.3
# fit1<- sma(logTM~logd2h*Taxa, data=allom.2)
# 
# plot(fit1, which='residual') 
# plot(fit1, which="qq")
# 
# library(smatr)
# fit1 <- slope.com(y=allom.2$logTM,x=allom.2$logd2h,groups=allom.2$Taxa) # slopes do not differ (p = 0.1)
# fit2 <- elev.com(y=allom.2$logTM,x=allom.2$logd2h,groups=allom.2$Taxa)  # intercepts DO differ (p < 0.0001)
# sma(logTM~logd2h+Taxa, data=allom.2)
# 
# 
# ancova.full <- lm(logTM~logd2h*Taxa*Treat,data=allom.2) # most higher order terms not significant
# ancova.2 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=allom.2) # drop 3-way interaction
# ancova.3 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat,data=allom.2)#drop Taxa:Treat
# ancova.4 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Treat,data=allom.2)#drop logd2h:Taxa
# 
# #seems like allometires should have taxa specific intercept and treatment specific slopes?
# 
# 
# ancova.5 <- lm(logTM~logd2h+Taxa+logd2h:Treat,data=allom.2) # significance of logd2h:Treat interaction goes away after dropping non-significant terms
# ancova.6 <- lm(logTM~logd2h+Taxa,data=allom.2)
# ancova.y <- lm(logTM~logd2h*Taxa+Taxa,data=allom.2)
# anova(ancova.full)
# anova(ancova.full,ancova.6) # full ancova has slightly lower AIC score, but is not significantly different than the full ancova.
# plot(ancova.4) # assumptions are met pretty well. It's not perfect, but it's good.
# anova(ancova.6)
# 
# 
# layout(matrix(c(1:28), nrow=7, ncol=4,byrow=T),
# heights=c(1,1,1,0.3,1,1,1),
# widths=c(1,0.3,1,1,1))
#  
