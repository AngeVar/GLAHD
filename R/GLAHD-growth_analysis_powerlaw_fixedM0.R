#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment
#------------------------------------------------------------------------------------------------------------------------------
source("R/loadLibraries.R")

#get taxa specific initial total mass from the pre-treatment harvest
allom <- read.csv("C:/Repos/GLAHD/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
pre <- subset(allom, Treatment == "Pre-treatment")
seed_mass<- summaryBy(Totmass~Taxa, FUN=min, data=pre, var.names = "seedmass", keep.names=T) # tried the minimum instead of the mean?

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-7)) #finds first date and labels it as Time 7 i.e. 07112014 is Day 7

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)

#dat4<- merge(dat3,seed_mass, by="Taxa")#add taxa specific initial mass to data



source("R/Paine.R")
##############################################################################################
##############################################################################################
#- fit growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
xlims <- c(0,70)
ylims <- c(0,100)

pdf(file="C:/Repos/GLAHD/Output/Growth_power_allplants.pdf", paper = "a4", width = 14/2.54, height = 14/2.54)

output <- data.frame(Species="",Treatment="",Location="",Taxa="",Code="",Range="",r="",beta="",AGR10g="",AGR10d="",RGR10d="",RGR10g="",stringsAsFactors=F)
full.output <- list()
for (i in 1:length(growth.l)){
  print(i)
  tofit <- growth.l[[i]]
  #seedmass <- seed_mass[which(tofit$Taxa[1] == seed_mass$Taxa),2] # get the right seedmass out of the dataframe above
  seedmass <- ifelse(min(tofit$TotMass-0.1)>0,min(tofit$TotMass-0.1),0.1) # define the seedmass for each plant, but force it to be positive
  #tofit <- subset(tofit,TotMass>seed.mass)
  fmla.pow2 <- as.formula(paste("~(", seedmass, "^(1-beta) + r*x*(1-beta))^(1/(1-beta))", sep = ""))
  SS.pow2   <- selfStart(fmla.pow2, initial = Init.pow2, parameters = c("r", "beta"))
  xvals <- seq(1,70,by=0.5) # establish xvalues of time to use in the model fitting below
  
  #fit
  tmp.pow2 <- getInitial(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  fit.pow2 <- gnls(TotMass ~ SS.pow2(Time, r, beta), data = tofit)
  out.pow2 <- output.pow2.gnls(fit.pow2,xvals, CI = T, M0 = seedmass) # i'm not sure this works. r and beta should be good, AGR and RGR are questionable.
  
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
output[,7:12] = apply(output[,7:12], 2, function(x) as.numeric(as.character(x)))
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



#Some plots
windows(20,12);par(mar=c(8,7,1,1))

output2$Taxa <- factor(output2$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=output2,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=output,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)

