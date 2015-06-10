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



#-------------------------------------------------------------------------------------
#-- function to interpret log-polynomial curve fits
output.log_lin <- function(X, Y, params, times,Code){
  a <- params[1]; b <- params[2]; c <- params[3]
  #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
  fitted <- exp(a + b*X + c*X^2)
  resid  <- Y - fitted
  data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
  #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
  eq <- bquote(a+bx+cx^2)
  mss  <- sum((fitted - mean(fitted))^2)
  rss  <- sum(resid^2)
  R2   <- mss/(mss + rss)
  rmse <- sqrt(rss)
  N <- length(X)
  logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
  AIC  <- -2 * logLik  + 2 * 3 # three parameters
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  temp <- a + b*X + c*X^2 # fix from here down
  rates = data.frame(
    times = times,
    M    =  exp(a+b*times+c*times^2))                  # from Hunt- Plant Growth Curves
  rates$AGR  =  with(rates,M*(b+2*c*times))            # from Hunt- Plant Growth Curves
  rates$RGR <- rates$AGR/rates$M
  rates$Code <- Code
  out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
  return(out)
}
#-------------------------------------------------------------------------------------





#-------------------------------------------------------------------------------------
#- fit log-polynomial growth model to each plant one at a time, extract parameters, and analyze them.
growth.l <- split(dat3,dat3$Code)
loglinfits <- output <- list()
for(i in 1:length(growth.l)){
  tofit <- growth.l[[i]]
  loglinfits[[i]] <- lm(log(TotMass)~Time+I(Time^2),data=tofit) #- fit log-polynomial
  
  #- extract the output using a new function which calculates AGR and RGR from the polynomial fit parameters
  output[[i]] <- output.log_lin(X=tofit$Time,Y=tofit$TotMass,
                              params=unname(coef(loglinfits[[i]])),times=seq(1:70),Code=tofit$Code[1])$rates
}
#- put rates dataframe together. This has the timecourse of mass, AGR, and RGR for each plant
rates.df <-do.call(rbind,output)
params <- as.data.frame(do.call(rbind,lapply(loglinfits,coef))) #- get polynomial parameters. Perhaps not that useful alone
params$Code <- unique(rates.df$Code) # add the code variable to the parameter estimate
#-------------------------------------------------------------------------------------




#-------------------------------------------------------------------------------------
#- plot the development for each treatment. Note how effect on RGR is particularly time-dependant
rates.df2 <- merge(rates.df,subset(dat2,Date==as.Date("2014-11-17")),by="Code")[,c(1:6,8:10,16)] # merge with dat2 to get Treatment, Location, etc

rates.trt <- summaryBy(M+RGR+AGR~times+Treatment+Location+Range,data=rates.df2,FUN=mean,keep.names=T)
rates.trt$combotrt <- as.factor(paste(rates.trt$Location,rates.trt$Range,rates.trt$Treatment,sep="_"))
plotBy(RGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topright")
plotBy(AGR~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft")
plotBy(M~times|combotrt,data=rates.trt,col=c("red","black","blue","green","orange","cyan","grey","yellow"),
       legendwhere="topleft")
#-------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------
#-- pull out RGR at 10 days and analyze
rates_10d <- subset(rates.df2,times==10)
rates_10d$Sp_RS_EN <- as.factor(with(rates_10d,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
rates_10d$Prov_Sp_EN <- as.factor(with(rates_10d,paste(Taxa,Species)))



#- analyze RGR
fm.rgr <- lme(RGR~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=rates_10d)
#look at model diagnostics
plot(fm.rgr,resid(.,type="p")~fitted(.) | Treatment,abline=0)   #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm.rgr,RGR~fitted(.)|Species,abline=c(0,1))            #predicted vs. fitted for each species
plot(fm.rgr,RGR~fitted(.),abline=c(0,1))                    #overall predicted vs. fitted
qqnorm(fm.rgr, ~ resid(., type = "p"), abline = c(0, 1))       #qqplot to assess normality of residuals
hist(fm.rgr$residuals[,1])
anova(fm.rgr) #- our three way interaction is now marginal. (it comes back if you choose 1 day rather than 10!)
plot(allEffects(fm.rgr))     
summary(fm.rgr)
#-------------------------------------------------------------------------------------