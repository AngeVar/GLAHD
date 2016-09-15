#-------------------------------------------------------------------------------------
#- This script provides an alternative curve fit for growth analysis via a GAM.
#- 
#-------------------------------------------------------------------------------------


#- load libraries from script
source("R_cleaned/GLAHD_LoadLibraries.R")
source("R_cleaned/fitGAM/derivSimulCI.R")
source("R_cleaned/fitGAM/plotCIdate.R")
source("R_cleaned/fitGAM/smoothplot.R")

# for GAM
Library(mgcv)



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
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
#   smoothplot(Time, lnTotMass, data=tofit, kgam=kgam)
#   title(main=tofit$Code[1])
  
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
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
