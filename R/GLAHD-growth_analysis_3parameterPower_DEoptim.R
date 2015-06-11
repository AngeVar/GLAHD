#------------------------------------------------------------------------------------------------------------------------------
# This script analyzes the growth data from the GLAHD glasshouse experiment.
#  - I don't believe the parameter estimates produced from this script.
#     plants that clearly have an increased rgr don't here, as the model
#     just fits a higher initial mass. Ugh. I think we need to fix the initial mass.
#  - NOTE- you can avoid the long model fitting process. Just skip to the bottom and read in the csv file.
#
#  - I'm attempting to fit the three-parameter power law function, but
#     parameter equifinality is a problem (some curves fit with a high
#     initial mass and a low r, which I think is wrong). One approach is
#     to fix the initial mass parameter at a common and low value, which 
#     works well but feels somewhat contrived- how do you choose a "correct" 
#     initial mass?
#  - This script takes a different approach and fits the full 3-parameter model,
#     but uses an evolutionary algorithm to fully explore the parameter space and 
#     (hopefully) avoid local minima to find the true global minimum (find the "true"
#     parameter values). Unfortunately this takes a lot of computing power, and is slow.
#  - In doing this, I noticed that the power function has a discontinuity at values of 
#    beta just above 1.0 (say, 1.001), where the predicted mass values are undefined.
#    This means that the evolutionary algorithm must be restricted out of that space, 
#    so I set an upper limit to the allowed beta values at 0.999 . This prevents our plants
#    from growing at rates that exceed the true exponential, which is fine.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")
library(DEoptim)

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-7)) #finds first date and labels it as Time 7 i.e. 07112014 is Day 7

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)




#--------------------------
#-- an example?
#- get an example curve to play with
junk <- subset(dat3,Code=="DTER-21")

#- function to fit the 3-parameter power law. From Paine et al. 2012.
fit.power3 <- function(p, Time, TotMass){
  M0   <- exp(p[1])
	beta <- exp(p[2])
	r    <- exp(p[3])
	y.pred <- (M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta))
	RSS    <- sum((TotMass-y.pred)*(TotMass-y.pred))/length(Time)
	return(RSS)
}

#- function to plot the result
fit.power.pred <- function(p, Time, TotMass){
  M0   <- exp(p[1])
  beta <- exp(p[2])
  r    <- exp(p[3])
  y.pred <- (M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta))
  RSS    <- sum((TotMass-y.pred)*(TotMass-y.pred))/length(Time)
  return(y.pred)
}
# # Provide the routine with some starting guesses for parameter values
# p <- log(c(0.08, 0.9, 0.4))
# 
# #- fit via optim(). I don't think this works very well.
# tmp <- optim(p, fit.power3, method="Nelder-Mead", Time=junk$Time, TotMass=junk$TotMass,
#              control=list(maxit=20000,trace=0))
# # Output parameter values
# best.pow.optim <- exp(tmp$par)

#- define upper and lower bounds and initial population estimate based on the output of optim()

# n <- length(lower) * 10
# Npop <- data.frame()
# Npop[1:n,1] <- rnorm(n=n,mean=best.pow.optim[1],sd=best.pow.optim[1]/5)
# Npop[1:n,2] <- rnorm(n=n,mean=best.pow.optim[2],sd=best.pow.optim[2]/5)
# Npop[1:n,3] <- rnorm(n=n,mean=best.pow.optim[3],sd=best.pow.optim[3]/5)
# Npop <- log(as.matrix(Npop))



set.seed(1234) # set seed for repeatablility
DEoptim.control <- list(VTR = -Inf, strategy = 2,itermax=500,trace=TRUE,CR=0.9)
out <- DEoptim(fit.power3,lower,upper,Time=junk$Time,TotMass=junk$TotMass)
out.par <- unname(out$optim$bestmem)



plot(junk$TotMass~junk$Time)
points(fit.power.pred(p=out.par,Time=junk$Time,TotMass=junk$TotMass)~junk$Time,col="red")



#- well that worked, so let's apply it to all our curves
#--------------------------







# #--- attempt to use parallel processing to speed up thecurve fitting process.
# THis is much faster, but I can't make sense of the output, as the plants
#  aren't fit in order. So I have no idea which parameters correspond to which
#  plant.
# p.flag <- T # use parallel processing?
# if(p.flag==T){
#   library(doParallel)
#   cl <- makeCluster(3)
#   registerDoParallel(cl)
# }
# 
# #- fit the curves via parallel processing, returns a list of output
# dat.l <- split(dat3,dat3$Code)
# output <- out.par <- list()
# DEoptim.control <- list(VTR = -Inf, strategy = 2,NP=300,itermax=500,trace=100,CR=0.9)
# output <- foreach(i=1:length(dat.l),.packages="DEoptim") %dopar% 
#   DEoptim(fit.power3,lower,upper,Time=dat.l[[i]]$Time,TotMass=dat.l[[i]]$TotMass)





#----------------------------------------------------------------------------------------------------------------
#- Loop over each plant and estimate the parameters by an evolutionary algorithm. This is slow (5-10 min?).
#- For testing purposes, reduce NP in DE.optim.control to 30. This gives approximately the same answers but is
#-   much much faster. Judging from my system monitor, DEoptim must use parallel processing under the hood.
#----------------------------------------------------------------------------------------------------------------

DEoptim.control <- list(VTR = -Inf, strategy = 2,NP=300,itermax=500,trace=100,CR=0.9)
lower <- log(c(0.001,0.5,0.001))
upper <- log(c(0.5,0.999,0.5))
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
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------




#---- start from here

#- merge it all the treatment code information from dat3
DEfits <- merge(params,subset(dat3,Date==as.Date("2014-11-17")),by="Code")

#- write out a csv for future work, as teh evolutionary algorithm takes a long time.
#write.csv(DEfits,row.names=F,file="C:/Repos/GLAHD/Output/DE_power_fits.csv")
DEfits <- read.csv("C:/Repos/GLAHD/Output/DE_power_fits.csv")

#Some plots. These make no sense. Initial mass is all over the place and explains why r doens't respond to treatment
windows(20,12);par(mar=c(8,7,1,1))

output$Taxa <- factor(output$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))#boxplot(r~Treatment*Range*Location,data=output,col=c("blue","red"),las=2,ylab="r",las=1)
boxplot(r~Treatment*Taxa,data=DEfits,col=c("blue","red"),las=2,ylab="r",cex.lab=2)
abline(v=16.5)
boxplot(beta~Treatment*Taxa,data=DEfits,col=c("blue","red"),las=2,ylab="beta",cex.lab=2)
abline(v=16.5)
boxplot(M0~Treatment*Taxa,data=DEfits,col=c("blue","red"),las=2,ylab="initial mass",cex.lab=2)

#- ugh. This is bad.
plotBy(r~M0|Treatment,data=DEfits,legendwhere="topright")
plotBy(r~beta|Treatment,data=DEfits,legendwhere="topright")
