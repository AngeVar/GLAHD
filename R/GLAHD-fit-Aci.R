#- to install the newest version of plantecophys
#library(devtools)
#install_bitbucket("remkoduursma/plantecophys")

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
p.flag <- T # use parallel processing?
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
acifits <- merge(fits.params,growth2,by=c("Code","Taxa"))
acifits$JtoV <- with(acifits,Jmax/Vcmax)
acifits$Location <- factor(acifits$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
acifits$Sp_RS_EN <- as.factor(with(acifits,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
acifits$Prov_Sp_EN <- as.factor(with(acifits,paste(Taxa,Species)))


#- does Vcmax or Jmax change with plant size? note that there will be a huge N vs. S effect in Vcmax
pairs(acifits[,c("Vcmax","Jmax","d2h","TotMass","leafArea")])



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






#----------------------------------------------------------------------------------
#-- plot Vcmax and Jmax for each taxa
windows(16,14);par(mfrow=c(3,1),mar=c(6,5,2,3),cex.lab=1.3,cex.axis=1.3,las=1)
#colors <- c(rep("white",5),rep("grey",3),rep("white",6),rep("grey",3))
colors <- c("blue","red")

#plot Vcmax
ylims=c(10,35)
boxplot(Vcmax~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(V["c,max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(acifits$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)

#plot Jmax
boxplot(Jmax~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
title(ylab=expression(J["max"]),cex.lab=2,line=2.5)
axis(side=1,at=seq(from=1.5,to=34.5,by=2),labels=levels(acifits$Taxa),las=2,cex.axis=1.5)
abline(v=16.4)

#plot Jmax:Vcmax
boxplot(JtoV~Treat+Taxa,data=acifits,axes=F,las=2,col=colors)
magaxis(c(2,3,4),labels=c(1,0,1),box=T,las=1)
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
#- create ANOVA table of results, embed that table into a word document. 
#- the table needs cleaning in word to embed the degrees of freedom and label the columns, etc.
vcmax.o <- anova(fm.vcmax)

extract.lme <- function(model){
  mod.o <- anova(model)
  
  ndf <- mod.o[2:nrow(mod.o),1]
  ddf <- mod.o[2:nrow(mod.o),2]
  df <- (paste(ndf,ddf,sep=", "))
  Fvalue <-round(mod.o[2:nrow(mod.o),3],1)
  Pvalue <-round(mod.o[2:nrow(mod.o),4],3)
  Pvalue[which(Pvalue==0)] <- "<0.001"
  return(data.frame("df"=df, "F" = Fvalue,"P" = Pvalue))
}


#- make a big table
gxtable <- do.call(cbind,lapply(list(fm.vcmax,fm.jmax,fm.jtov,fm.Asat,fm.Rmass),FUN=extract.lme))
row.names(gxtable) <- row.names(anova(fm.vcmax))[2:8]

library(R2wd)
wdGet()
wdTable(gxtable,autoformat=2)
#----------------------------------------------------------------------------------

