#supplementary figure allometry

allom <- read.csv("Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
allom$logd2h <- log10(allom$d2h)
allom$logTM <- log10(allom$Totmass)
allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                ifelse(allom$Pot>=40,"Pre","Warmed")))
allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

#- split into list across all taxa for plotting
allom.2 <- subset(allom,Treat!="Pre")
allom.2$Treat <- factor(allom.2$Treat)
allom.l <- split(allom.2,allom.2$Taxa)

windows(30,30)
par(mfrow=c(5,4), mar=c(2,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(-0.5,2.1)
xlims <- c(-.5,3)
palette(c("black","red"))
for (i in 1:length(allom.l)){
  toplot <- allom.l[[i]]
  plotBy(logTM~logd2h|Treat,data=toplot,type="p",pch=19,xlim=xlims,ylim=ylims,axes=F,
         ylab="H",xlab="",legend=F, cex=2)
  magaxis(side=c(1:4),labels=c(1,1,0,0))
  mtext(text=toplot$Taxa,side=1,line=-1.5)
  
  fit.h <- lm(logTM~logd2h,data=subset(toplot,Treat=="Home"))
  predline(fit.h,col=alpha("black",alpha=0.4), lwd=1)
  
  fit.w <- lm(logTM~logd2h,data=subset(toplot,Treat=="Warmed"))
  predline(fit.w,col=alpha("red",alpha=0.4),lwd=1)
  
}
mtext(expression(log[10]~(diameter^2~"*"~height)),side=1,line=2,outer=T,cex=1.5)
mtext(expression(log[10]~(Total~mass~(g))),side=2,line=2,outer=T,cex=1.5)
legend(x=4,y=2,legend=c(expression(Warmed~(+3.5~degree~C)),"Home"),pch=22,pt.cex=2,pt.bg=c(alpha("red",0.6),alpha("black",0.6)),
       bty="n",cex=1.5, xpd=T)


#analysis
ancova.full <- lm(logTM~logd2h*Taxa*Treat,data=allom.2) # most higher order terms not significant
ancova.2 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat+Taxa:Treat,data=allom.2) # drop 3-way interaction
ancova.3 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Taxa+logd2h:Treat,data=allom.2)#drop Taxa:Treat
ancova.4 <- lm(logTM~logd2h+Taxa+Treat+logd2h:Treat,data=allom.2)#drop logd2h:Taxa

#seems like allometires should have taxa specific intercept and treatment specific slopes?


ancova.5 <- lm(logTM~logd2h+Taxa+logd2h:Treat,data=allom.2) # significance of logd2h:Treat interaction goes away after dropping non-significant terms
ancova.6 <- lm(logTM~logd2h+Taxa,data=allom.2)
ancova.y <- lm(logTM~logd2h*Taxa+Taxa,data=allom.2)
anova(ancova.full)
anova(ancova.full,ancova.6) # full ancova has slightly lower AIC score, but is not significantly different than the full ancova.
plot(ancova.4) # assumptions are met pretty well. It's not perfect, but it's good.
anova(ancova.6)
