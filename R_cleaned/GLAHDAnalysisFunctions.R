#AnalysisFunctions

# function to predict total plant mass from information on taxa, 
# treatment, height, and diameter. Takes a series of vectors, returns a vector

predict_mass <- function(Taxa,Treat,Diameter,Height,model_flag="complex"){
  
  #- read in the data, do a few conversions
  allom <- read.csv("Data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
  allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
  allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  allom$logd2h <- log10(allom$d2h)
  allom$logTM <- log10(allom$Totmass)
  allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                  ifelse(allom$Pot>=40,"Pre","Warmed")))
  allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  
  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  allom <- merge(allom,linkdf,by="Species")
  allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
  
  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  lm.taxa_complex <- lm(logTM~logd2h*Taxa,data=allom) #taxa specfic intercepts and slopes
  lm.taxa_simple  <- lm(logTM~logd2h+Taxa,data=allom) #taxa specific intercepts, common slope
  # -----------------------------------
  
  newdat.s <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat.s$logd2h <- with(newdat.s,log10((Diameter/10)^2*(Height)))
  newdat.s$Totmass <- 10^predict(lm.taxa_simple,newdat.s)
  
  
  newdat.c <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat.c$logd2h <- with(newdat.c,log10((Diameter/10)^2*(Height)))
  newdat.c$Totmass <- 10^predict(lm.taxa_complex,newdat.c)
  if(model_flag=="complex") {return(newdat.c$Totmass)} # simple or complex model? Should taxa have different slopes?
  if(model_flag=="simple") {return(newdat.s$Totmass)} # simple or complex model? Should taxa have different slopes?
  
}

# function to predict total plant leaf area from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
predict_LA <- function(Taxa,Treat,Diameter,Height){
  
  #- read in the data, do a few conversions
  allom <- read.csv("data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
  allom$Totmass <- base::rowSums(allom[,11:13]) #total mass is the sum of leaf, stem, and root mass
  allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  allom$logd2h <- log10(allom$d2h)
  allom$logLA <- log10(allom$Leafarea)
  allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                  ifelse(allom$Pot>=40,"Pre","Warmed")))
  allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  
  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  allom <- merge(allom,linkdf,by="Species")
  allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
  
  
  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  allom.2 <- subset(allom,Treatment!="Pre-treatment")
  #get mean of pre-treatment trees and add to both treatments
  pretre<- subset(allom, Treat== "Pre")
  sumpre1<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T);sumpre2<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T)
  sumpre1$Treat<- "Home";sumpre2$Treat<- "Warmed"; sumpre3<- droplevels(rbind(sumpre1,sumpre2))
  sumpre3$Date<-as.Date("2014-11-06");sumpre3$Treatment<-"Pre-treatment";sumpre3$Location<-c(rep("S",8),rep("N",9), rep("S",8),rep("N",9));sumpre3$Notes<- " "
  sumpre3$Species <-c(rep(c("TER","TER","CAM","CAM","CAM","BOT","LONG","SMIT","TER","TER","TER","CAM","CAM","CAM","BRA","PEL","PLAT"),2))
  sumpre3$Code <- paste(sumpre3$Taxa,sumpre3$Pot, sep="-");sumpre3$Range<-c(rep(c(rep("wide",5),rep("narrow",3),rep("wide",6),rep("narrow",3)),2))
  allom.3<- droplevels(rbind.fill(sumpre3,allom.2))
  allom.3$Treat<-as.factor(allom.3$Treat)
  
  #   allom.2$Treatment <- factor(allom.2$Treatment)
  #   ancova.full <- lm(logLA~logd2h*Taxa*Treat,data=allom.2) # 3-way term significant, can't be dropped.
  #                                                           # this means that the warming effect on the slope is taxa-specific
  ancova.full <- lm(logLA~logd2h*Taxa*Treat,data=allom.3)     # fit model including pre-treatment mean added to both treatments
  #-----------------------------------
  
  newdat <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat$logd2h <- with(newdat,log10((Diameter/10)^2*(Height)))
  newdat$leafArea <- 10^predict(ancova.full,newdat)
  return(newdat$leafArea)
}

# function to predict total plant leaf mass from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
predict_leaf_mass <- function(Taxa,Treat,Diameter,Height,model_flag="simple"){

  #- read in the data, do a few conversions
  allom <- read.csv("data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
  allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  allom$logd2h <- log10(allom$d2h)
  allom$logLM <- log10(allom$Leafmass)
  allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                  ifelse(allom$Pot>=40,"Pre","Warmed")))
  allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))

  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  allom <- merge(allom,linkdf,by="Species")
  allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  lm.taxa_complex <- lm(logLM~logd2h*Taxa,data=allom)
  lm.taxa_simple  <- lm(logLM~logd2h+Taxa,data=allom)
  # -----------------------------------

  newdat.s <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat.s$logd2h <- with(newdat.s,log10((Diameter/10)^2*(Height)))
  newdat.s$leafMass <- 10^predict(lm.taxa_simple,newdat.s)

  newdat.c <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat.c$logd2h <- with(newdat.c,log10((Diameter/10)^2*(Height)))
  newdat.c$leafMass <- 10^predict(lm.taxa_complex,newdat.c)
  if(model_flag=="complex") {return(newdat.c$leafMass)} # simple or complex model? Should taxa have different slopes?
  if(model_flag=="simple") {return(newdat.s$leafMass)} # simple or complex model? Should taxa have different slopes?

}

# function to predict total plant stem mass from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
predict_stem_mass <- function(Taxa,Treat,Diameter,Height){

  #- read in the data, do a few conversions
  allom <- read.csv("data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
  allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  allom$logd2h <- log10(allom$d2h)
  allom$logSM <- log10(allom$Stemmass)
  allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                  ifelse(allom$Pot>=40,"Pre","Warmed")))
  allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))

  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  allom <- merge(allom,linkdf,by="Species")
  allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  #-- allometry ANCOVA for all taxa
  allom.2 <- subset(allom,Treatment!="Pre-treatment")
  #get mean of pre-treatment trees and add to both treatments
  pretre<- subset(allom, Treat== "Pre")
  sumpre1<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T);sumpre2<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T)
  sumpre1$Treat<- "Home";sumpre2$Treat<- "Warmed"; sumpre3<- droplevels(rbind(sumpre1,sumpre2))
  sumpre3$Date<-as.Date("2014-11-06");sumpre3$Treatment<-"Pre-treatment";sumpre3$Location<-c(rep("S",8),rep("N",9), rep("S",8),rep("N",9));sumpre3$Notes<- " "
  sumpre3$Species <-c(rep(c("TER","TER","CAM","CAM","CAM","BOT","LONG","SMIT","TER","TER","TER","CAM","CAM","CAM","BRA","PEL","PLAT"),2))
  sumpre3$Code <- paste(sumpre3$Taxa,sumpre3$Pot, sep="-");sumpre3$Range<-c(rep(c(rep("wide",5),rep("narrow",3),rep("wide",6),rep("narrow",3)),2))
  allom.3<- droplevels(rbind.fill(sumpre3,allom.2))
  allom.3$Treat<-as.factor(allom.3$Treat)

  lm.taxa <- lm(logSM~logd2h+Taxa+Treat+logd2h:Treat+Taxa:Treat,data=allom.3)
  # -----------------------------------

  newdat <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat$logd2h <- with(newdat,log10((Diameter/10)^2*(Height)))
  newdat$stemMass <- 10^predict(lm.taxa,newdat)

  return(newdat$stemMass)

}

# function to predict total plant root mass from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
predict_root_mass <- function(Taxa,Treat,Diameter,Height){

  #- read in the data, do a few conversions
  allom <- read.csv("data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
  allom$d2h <- with(allom,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  allom$logd2h <- log10(allom$d2h)
  allom$logRM <- log10(allom$Rootmass)
  allom$Treat <- as.factor(ifelse(allom$Pot < 20, "Home",
                                  ifelse(allom$Pot>=40,"Pre","Warmed")))
  allom$Taxa <- factor(allom$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                           "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))

  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  allom <- merge(allom,linkdf,by="Species")
  allom[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA

  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  allom.2 <- subset(allom,Treatment!="Pre-treatment")
  #get mean of pre-treatment trees and add to both treatments
  pretre<- subset(allom, Treat== "Pre")
  sumpre1<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T);sumpre2<- summaryBy(.~Taxa, data=pretre, FUN=mean, keep.names=T)
  sumpre1$Treat<- "Home";sumpre2$Treat<- "Warmed"; sumpre3<- droplevels(rbind(sumpre1,sumpre2))
  sumpre3$Date<-as.Date("2014-11-06");sumpre3$Treatment<-"Pre-treatment";sumpre3$Location<-c(rep("S",8),rep("N",9), rep("S",8),rep("N",9));sumpre3$Notes<- " "
  sumpre3$Species <-c(rep(c("TER","TER","CAM","CAM","CAM","BOT","LONG","SMIT","TER","TER","TER","CAM","CAM","CAM","BRA","PEL","PLAT"),2))
  sumpre3$Code <- paste(sumpre3$Taxa,sumpre3$Pot, sep="-");sumpre3$Range<-c(rep(c(rep("wide",5),rep("narrow",3),rep("wide",6),rep("narrow",3)),2))
  allom.3<- droplevels(rbind.fill(sumpre3,allom.2))
  allom.3$Treat<-as.factor(allom.3$Treat)

  lm.taxa <- lm(logRM~logd2h*Taxa+Treat,data=allom.3)
  # -----------------------------------

  newdat <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat$logd2h <- with(newdat,log10((Diameter/10)^2*(Height)))
  newdat$rootMass <- 10^predict(lm.taxa,newdat)

  return(newdat$rootMass)
}

#- read and process the height and diameter dataset. Estimate total mass based on size by calling predict_mass()
return_size_mass <- function(model_flag="simple"){
  
  #- read in the data, do a few conversions
  
  dat <- read.csv("data/GHS39_GLAHD_MAIN_HEIGHTDIAMETER_20141106-20150105_L2.1.csv")
  dat$Date <- as.Date(dat$Date)
  dat$d2h <- with(dat,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  
  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  dat <- merge(dat,linkdf,by="Species")
  
  # get rid of the initial harvest and periodic harvest data
  dat2 <- subset(dat,Date!=as.Date("2014-11-06") & Date != as.Date("2014-11-20")& ToFit != 0)
  dat2$Treatment <- factor(dat2$Treatment)
  dat2$Code <- factor(dat2$Code)
  dat2$Taxa <- factor(dat2$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  
  
  dat2 <- dat2[with(dat2,order(Taxa)),]
  
  dat2$TotMass <- predict_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height,model_flag=model_flag)
  dat2$leafArea <- predict_LA(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  
  return(dat2)
}


#- read and process the height and diameter dataset. Estimate total mass based on size by calling predict_mass()
return_size_mass_all <- function(model_flag="simple"){

  #- read in the data, do a few conversions
  dat <- read.csv("data/GHS39_GLAHD_MAIN_HEIGHTDIAMETER_20141106-20150105_L2.1.csv")
  dat$Date <- as.Date(dat$Date)
  dat$d2h <- with(dat,(Diameter/10)^2*(Height)) #calculate d2h in cm3

  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  dat <- merge(dat,linkdf,by="Species")

  # get rid of the initial harvest and periodic harvest data
  dat2 <- subset(dat,Date!=as.Date("2014-11-06") & Date != as.Date("2014-11-20")& ToFit != 0)
  dat2$Treatment <- factor(dat2$Treatment)
  dat2$Code <- factor(dat2$Code)
  dat2$Taxa <- factor(dat2$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))


  dat2 <- dat2[with(dat2,order(Taxa)),]

  dat2$TotMass <- predict_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height,model_flag=model_flag)
  dat2$leafArea <- predict_LA(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  dat2$leafMass <- predict_leaf_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  dat2$stemMass <- predict_stem_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  dat2$rootMass <- predict_root_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  return(dat2)
}

# #fit 3 parameter power function
# fit.power3 <- function(p, Time, TotMass){
#   M0   <- exp(p[1])
#   beta <- exp(p[2])
#   r    <- exp(p[3])
#   y.pred <- (M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta))
#   RSS    <- sum((TotMass-y.pred)*(TotMass-y.pred)/y.pred)/length(Time)#edited by JED 11 June 2015
#   return(RSS)
# }
# #plot results from 3 parameter power function
# fit.power.pred <- function(p, Time, TotMass){
#   M0   <- exp(p[1])
#   beta <- exp(p[2])
#   r    <- exp(p[3])
#   y.pred <- (M0^(1-beta) + r*Time*(1-beta))^(1/(1-beta))
#   return(y.pred)
# }
# 
# # from Paine et al. 2012
# # Self-start function for 2-parameter power-law fit (initial mass fixed by used). This self-start for this model is not provided by Pinhero and Bates (2000), so I wrote one. 
# Init.pow2 <- function(mCall, LHS, data){
#   xy <- sortedXyData(mCall[["x"]], LHS, data)
#   if(nrow(xy) < 4) {stop("Too few distinct x values to fit a power-law")}
#   r    <-coef(lm(log(y) ~ x, xy))[2]    # Use the slope from a log fit to the data as an initial guess for r
#   beta <- 0.9                             # give initial guess of beta as 0.9. don't use beta = 1, as it's undefined (1/0)
#   value <- c(r, beta) 
#   names(value) <- mCall[c("r", "beta")]
#   return(value)
# }
# 
# 
# output.pow2.gnls <- function(fit, times, CI = F, alpha = 0.05, M0){
#   params <- c(coef(fit), M0 = M0)
#   r = params[1]; beta = params[2]
#   fitted <- fit$fitted
#   resid  <- fit$residual
#   data <- data.frame(fitted = fitted, resid = resid)
#   eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
#   mss  <- sum((fitted - mean(fitted))^2)
#   rss  <- sum(resid^2)
#   R2   <- mss/(mss + rss)
#   rmse <- sqrt(rss)
#   AIC <- AIC(fit)
#   summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
#   temp <- M0^(1-beta) + r*times*(1-beta)
#   rates = data.frame(
#     times = times,
#     M    =     temp^(1/(1-beta)),
#     AGR  = r * temp^(beta/(1-beta)))
#   rates$RGRt <- rates$AGR/rates$M
#   rates$RGRm  <-  r * rates$M^(beta-1)
#   if(CI == T){
#     cov <- fit$varBeta
#     cov <- cbind(rbind(cov, 0), 0) # variance-covariance for M0 is 0 - it's fixed. 
#     x <- data.frame(rmvnorm(n=1000, mean=params, sigma=cov))
#     M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
#     for(i in 1:nrow(x)){
#       x.i  <- x[i,]
#       temp.i   <- x.i$M0^(1-x.i$beta) + x.i$r * times *(1-x.i$beta)
#       M[i,]    <-         temp.i^(1/(1-x.i$beta))
#       AGR[i,]  <- x.i$r * temp.i^(x.i$beta/(1-x.i$beta))
#       RGRt[i,] <- AGR[i,]/M[i,]
#       RGRm[i,] <- x.i$r * M[i,]^(x.i$beta - 1)
#     }
#     CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
#     out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
#   } else {
#     out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
#   }
#   return(out)
# }
# 
# # this function returns confidence envelopes around growth trajectories, and growth rates. 
# summarizer <- function(dat, alpha){
#   n <- length(dat)
#   quantiles <- c(alpha/2, 1-(alpha/2))
#   CIs <- data.frame(matrix(NA, ncol(dat[[1]]), n*2))
#   names(CIs) <- paste(rep(names(dat), each = 2), c("lo", "hi"), sep = ".")
#   for(i in 1:n){
#     CIs[,(2*i-1):(2*i)] <- t(apply(dat[[i]],    2, quantile, quantiles, na.rm = T))
#   }
#   return(CIs)
# }



#- this function reads in and returns the Rdark spot measurements. Returns a dataframe.
getRdark <- function(){

  #- read data
  gx <- read.csv(file="data/GHS39_GLAHD_MAIN_GX-RDARK_20141215_L2.csv")

  #- get the first bit of the code (the taxa)
  gx$Taxa <- unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=1,to=nrow(gx)*2,by=2)]
  gx$Taxa <- factor(gx$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  gx$Pot <- as.numeric(unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=2,to=nrow(gx)*2,by=2)])
  gx$Treatment <- as.factor(ifelse(gx$Pot < 20, "Home","Warmed"))

  #- average over sub-replicate logs
  gx1 <- summaryBy(.~Code+Taxa+Treatment+Pot,data=gx,keep.names=T)
  gx1$R <- gx1$Photo*-1


  #---------------------------
  #- get the leaf mass data
  lm <- read.csv(file="data/GHS39_GLAHD_MAIN_BIOMASS-RDLEAVES_20141222_L2.csv")

  gx2 <- merge(gx1,lm, by=c("Code","Taxa","Pot"))
  gx2$Rmass <- with(gx2,R*Area/10000/Leafmass*1000*1000)
  gx2$SLA <- with(gx2,(Area/(Leafmass/1000)))
  #---------------------------


  #---------------------------
  #- prepare data for statistical analysis
  #- get the growth data, mostly just for the treatment codes
  growth <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
  growth2 <- summaryBy(d2h+TotMass+leafArea~Species+Treatment+Location+Taxa+Code+Range,keep.names=T,data=subset(growth,Date >= as.Date("2014-12-8") & Date <=as.Date("2014-12-20")))

  #- merge size totalmass and leafarea data into dataframe with Rdark values
  gx3 <- merge(gx2,growth2,by=c("Code","Taxa","Treatment"))


  gx3$Location <- factor(gx3$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
  gx3$Sp_RS_EN <- as.factor(with(gx3,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
  gx3$Prov_Sp_EN <- as.factor(with(gx3,paste(Taxa,Species)))
  return(gx3)
}


#- this function reads and processes the Asat spot measurements. Returns a dataframe.
getAsat <- function(){


  #- read data
  gx <- read.csv(file="data/GHS39_GLAHD_MAIN_GX-ASAT_20141204_L2.csv")

  #- get the first bit of the code (the taxa)
  gx$Taxa <- unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=1,to=2071,by=2)]
  gx$Taxa <- factor(gx$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  gx$Pot <- as.numeric(unlist(strsplit(x=as.character(gx$Code),split="-"))[seq(from=2,to=2072,by=2)])
  gx$Treat <- as.factor(ifelse(gx$Pot < 20, "Home","Warmed"))
  gx$Location <- as.factor(ifelse(gx$Taxa == "ACAM" | gx$Taxa == "BCAM"| gx$Taxa == "CCAM"|
                                    gx$Taxa == "ATER"| gx$Taxa == "BTER"| gx$Taxa == "LONG"|
                                    gx$Taxa == "SMIT" | gx$Taxa == "BOT", "S","N"))
  gx$Range <- as.factor(ifelse(gx$Taxa == "LONG"| gx$Taxa == "SMIT" |gx$Taxa == "BOT"| gx$Taxa == "PEL"
                               | gx$Taxa == "PLAT"|gx$Taxa =="BRA","Narrow","Wide"))
  #- average over sub-replicate logs
  gx2 <- summaryBy(.~Code+Taxa+Treat+Location+Range,data=gx,keep.names=T)





  #---------------------------
  #- prepare data for statistical analysis
  #- get the growth data, mostly just for the treatment codes
  growth <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
  growth2 <- summaryBy(d2h+TotMass+leafArea~Species+Treatment+Location+Taxa+Code+Range,keep.names=T,data=subset(growth,Date >= as.Date("2014-12-8") & Date <=as.Date("2014-12-20")))

  #- merge size totalmass and leafarea data into dataframe with aci values
  gx2$Range <- as.factor(tolower(gx2$Range))
  gx3 <- merge(gx2,growth2,by=c("Code","Taxa","Location","Range"))


  gx3$Location <- factor(gx3$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
  gx3$Sp_RS_EN <- as.factor(with(gx3,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
  gx3$Prov_Sp_EN <- as.factor(with(gx3,paste(Taxa,Species)))

  return(gx3)
}

# 
# 
# 
# #-------------------------------------------------------------------------------------
# #-- function to interpret log-polynomial curve fits that are THIRD ORDER
# output.log_3order <- function(X, Y, params, times,Code){
#   a <- params[1]; b <- params[2]; c <- params[3]; d <- params[4]
#   #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
#   fitted <- exp(a + b*X + c*X^2 + d*X^3)
#   resid  <- Y - fitted
#   data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
#   #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
#   eq <- bquote(a+bx+cx^2+dx^3)
#   mss  <- sum((fitted - mean(fitted))^2)
#   rss  <- sum(resid^2)
#   R2   <- mss/(mss + rss)
#   rmse <- sqrt(rss)
#   N <- length(X)
#   logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
#   AIC  <- -2 * logLik  + 2 * 3 # four parameters
#   summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
#   #temp <- a + b*X + c*X^2 
#   rates = data.frame(
#     times = times,
#     M    =  exp(a+b*times+c*times^2+d*times^3))                  # from Hunt- Plant Growth Curves
#   rates$AGR  =  with(rates,M*(b+2*c*times+3*d*times^2))            # from Hunt- Plant Growth Curves
#   rates$RGR <- rates$AGR/rates$M
#   rates$Code <- Code
#   out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
#   return(out)
# }
# #-------------------------------------------------------------------------------------
# 
# 
# 
# 
# #-------------------------------------------------------------------------------------
# #-- function to interpret log-polynomial curve fits that are FOURTH ORDER
# output.log_4order <- function(X, Y, params, times,Code){
#   a <- params[1]; b <- params[2]; c <- params[3]; d <- params[4]; e<- params[5]
#   #fitted <- (M0^(1-beta)   + r  * X *(1-beta))^(1/(1-beta))
#   fitted <- exp(a + b*X + c*X^2 + d*X^3 +e*X^4)
#   resid  <- Y - fitted
#   data <- data.frame(time=X,observed = Y, fitted = fitted, resid = resid)
#   #eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
#   eq <- bquote(a+bx+cx^2+dx^3+ex^4)
#   mss  <- sum((fitted - mean(fitted))^2)
#   rss  <- sum(resid^2)
#   R2   <- mss/(mss + rss)
#   rmse <- sqrt(rss)
#   N <- length(X)
#   logLik <- -N * (log(2 * pi) + 1 - log(N) + log(sum(resid^2)))/2
#   AIC  <- -2 * logLik  + 2 * 3 # four parameters
#   summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
#   #temp <- a + b*X + c*X^2 
#   rates = data.frame(
#     times = times,
#     M    =  exp(a+b*times+c*times^2+d*times^3+e*times^4))                  # from Hunt- Plant Growth Curves
#   rates$AGR  =  with(rates,M*(b+2*c*times+3*d*times^2+4*e*times^3))            # from Hunt- Plant Growth Curves
#   rates$RGR <- rates$AGR/rates$M
#   rates$Code <- Code
#   out <- list(params = params[-4], summary = summary, equation = eq, data = data, rates = rates)
#   return(out)
# }
# #-------------------------------------------------------------------------------------
# 
