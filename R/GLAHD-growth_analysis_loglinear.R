#-------------------------------------------------------------------------------------
#- This script provides an alternative curve fit for growth analysis.
#- This works quickly and simply because it uses a LINEAR model rather than a nonlinear model.
#- This assumes that a plot of log(Mass)~Time can be adequately fit with a polynomial.
#- This assumption seems to be correct in this case, for our data.
#- I don't know whether it's "better" than a power model, but it is much more convenient!
#-------------------------------------------------------------------------------------


#- load libraries from script
source("R/loadLibraries.R")

#- read in the data, do a few conversions
dat2 <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-7)) #finds first date and labels it as Time 7 i.e. 07112014 is Day 7

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)) # get the frequency of observations for each pot
names <- names(table(dat2$Code))# get the associated name of each pot
keeps <- names[which(obs>6)]    # return a vector of pot names with more than n observations
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)


# 
# #-------------------------------------------------------------------------------------
# #-- function to interpret log-polynomial curve fits
# output.log_lin <- function(X, Y, params, times,Code){
#   a <- params[1]; b <- params[2]; c <- params[3]
#   #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
#   fitted <- exp(a + b*X + c*X^2)
#   resid  <- Y - fitted
#   data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
#   #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
#   eq <- bquote(a+bx+cx^2)
#   mss  <- sum((fitted - mean(fitted))^2)
#   rss  <- sum(resid^2)
#   R2   <- mss/(mss + rss)
#   rmse <- sqrt(rss)
#   N <- length(X)
#   logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
#   AIC  <- -2 * logLik  + 2 * 3 # three parameters
#   summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
#   temp <- a + b*X + c*X^2 # fix from here down
#   rates = data.frame(
#     times = times,
#     M    =  exp(a+b*times+c*times^2))                  # from Hunt- Plant Growth Curves
#   rates$AGR  =  with(rates,M*(b+2*c*times))            # from Hunt- Plant Growth Curves
#   rates$RGR <- rates$AGR/rates$M
#   rates$Code <- Code
#   out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
#   return(out)
# }
# #-------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------
#-- function to interpret log-polynomial curve fits that are THIRD ORDER
output.log_3order <- function(X, Y, params, times,Code){
  a <- params[1]; b <- params[2]; c <- params[3]; d <- params[4]
  #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
  fitted <- exp(a + b*X + c*X^2 + d*X^3)
  resid  <- Y - fitted
  data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
  #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
  eq <- bquote(a+bx+cx^2+dx^3)
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
    M    =  exp(a+b*times+c*times^2+d*times^3))                  # from Hunt- Plant Growth Curves
  rates$AGR  =  with(rates,M*(b+2*c*times+3*d*times^2))            # from Hunt- Plant Growth Curves
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
}
dat4 <- do.call(rbind,growth.l)
plotBy(RGR~jitter(as.numeric(Date))|Code,type="b",data=dat4,legend=F);abline(h=0) 
#- it looks like a third-order polynomial may fit the data best, which would produce
#   and RGR that changes over time according to a second-order polynomial

#-------------------------------------------------------------------------------------


# 
# #-------------------------------------------------------------------------------------
# #- fit log-polynomial growth model to each plant one at a time, extract parameters, and analyze them.
# growth.l <- split(dat3,dat3$Code)
# loglinfits <- output <- data.out <- list()
# for(i in 1:length(growth.l)){
#   tofit <- growth.l[[i]]
#   loglinfits[[i]] <- lm(log(TotMass)~Time+I(Time^2),data=tofit) #- fit log-polynomial
#   
#   #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
#   output[[i]] <- output.log_lin(X=tofit$Time,Y=tofit$TotMass,
#                               params=unname(coef(loglinfits[[i]])),times=seq(1:70),Code=tofit$Code[1])$rates
#   data.out[[i]] <- output.log_lin(X=tofit$Time,Y=tofit$TotMass,
#                                 params=unname(coef(loglinfits[[i]])),times=seq(1:70),Code=tofit$Code[1])$data
#   data.out[[i]]$Code <- tofit$Code[1]
# }
# #- put rates dataframe together. This has the timecourse of mass, AGR, and RGR for each plant
# rates.df <-do.call(rbind,output)
# data.df <- do.call(rbind,data.out)
# data.df$fit_type <- "log-polynomial"
# params <- as.data.frame(do.call(rbind,lapply(loglinfits,coef))) #- get polynomial parameters. Perhaps not that useful alone
# params$Code <- unique(rates.df$Code) # add the code variable to the parameter estimate
# #-------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------
#- fit log-3polynomial growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
log3fits <- output <- data.out <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  log3fits[[i]] <- lm(log(TotMass)~Time+I(Time^2)+I(Time^3),data=tofit) #- fit log-polynomial
  
  #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
  output[[i]] <- output.log_3order(X=tofit$Time,Y=tofit$TotMass,
                                params=unname(coef(log3fits[[i]])),times=seq(1:70),Code=tofit$Code[1])$rates
  data.out[[i]] <- output.log_3order(X=tofit$Time,Y=tofit$TotMass,
                                  params=unname(coef(log3fits[[i]])),times=seq(1:70),Code=tofit$Code[1])$data
  data.out[[i]]$Code <- tofit$Code[1]
}
#- put rates dataframe together. This has the timecourse of mass, AGR, and RGR for each plant
rates.df <-do.call(rbind,output)
data.df <- do.call(rbind,data.out)
data.df$fit_type <- "log-3polynomial"
params <- as.data.frame(do.call(rbind,lapply(log3fits,coef))) #- get polynomial parameters. Perhaps not that useful alone
params$Code <- unique(rates.df$Code) # add the code variable to the parameter estimate
#-------------------------------------------------------------------------------------



#-------------------------------------------------------------------------------------
#- plot the development for each treatment. Note how effect on RGR is particularly time-dependant
rates.df2 <- merge(rates.df,subset(dat2,Date==as.Date("2014-11-17")),by="Code")[,c(1:6,8:10,16)] # merge with dat2 to get Treatment, Location, etc

rates.trt <- summaryBy(M+RGR+AGR~times+Treatment+Location+Range,data=rates.df2,FUN=mean,keep.names=T)
rates.trt$combotrt <- as.factor(paste(rates.trt$Location,rates.trt$Range,rates.trt$Treatment,sep="_"))

windows(40,30);par(mfrow=c(3,1))
plotBy(RGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright")
plotBy(AGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft")
plotBy(M~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft")
#-------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------
#- compare fits of the 3-parameter power function and the log-polynomial function.


#- read in the 3-parameter power law fits as a comparison
powfits <- read.csv("C:/Repos/GLAHD/R/DE_power_fits.csv")
powfits <- powfits[,1:4]
powfits2 <- merge(powfits,dat3,by="Code")
head(powfits2)
powfits2$predmass <- with(powfits2,(M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta)))

windows(30,40);par(mfrow=c(2,1),mar=c(5,6,1,1))
#- linear scale
plot(predmass~TotMass,data=powfits2,pch=3,col="red",cex.lab=1.5,
     xlab="Observed mass (g)",ylab="Predicted mass (g)")
points(fitted~observed,data=data.df,pch=3,col="black")
abline(0,1)
legend("topleft",legend=c("log-3poly","3-power"),pch=3,col=c("black","red"))
#- log scale
plot(predmass~TotMass,data=powfits2,pch=3,col="red",cex.lab=1.5,log="xy",
     xlab="Observed mass (g)",ylab="Predicted mass (g)")
points(fitted~observed,data=data.df,pch=3,col="black")
abline(0,1)
legend("topleft",legend=c("log-3poly","3-power"),pch=3,col=c("black","red"))

summary(lm(predmass~TotMass,data=powfits2)) #- note that the power function still has a higher r2 and fits better at higher mass values
summary(lm(fitted~observed,data=data.df))#- note that the 3rd order poly's r2 is still very good but lower than the powers. 
                                         #- it fits better at low values of mass, which is expected given that it is fit to log-transformed values.
#-------------------------------------------------------------------------------------