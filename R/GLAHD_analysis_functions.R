
# function to predict total plant mass from information on taxa, treatment, height, and diameter. Takes a series of vectors, returns a vector
predict_mass <- function(Taxa,Treat,Diameter,Height,model_flag="simple"){
  
  #- read in the data, do a few conversions
  allom <- read.csv("C:/Repos/GLAHD/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
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
  allom$LMF <- with(allom,Leafmass/Totmass)
  
  #-----------------------------------
  #-- allometry ANCOVA for all taxa
  lm.taxa_complex <- lm(logTM~logd2h*Taxa,data=allom)
  lm.taxa_simple  <- lm(logTM~logd2h+Taxa,data=allom)
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
  allom <- read.csv("C:/Repos/GLAHD/Data/Harvests/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
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
  allom.2$Treatment <- factor(allom.2$Treatment)
  ancova.full <- lm(logLA~logd2h*Taxa*Treat,data=allom.2) # 3-way term significant, can't be dropped.
                                                          # this means that the warming effect on the slope is taxa-specific
  
  #-----------------------------------
  
  newdat <- data.frame(Taxa=Taxa,Treat=Treat,Diameter=Diameter,Height=Height)
  newdat$logd2h <- with(newdat,log10((Diameter/10)^2*(Height)))
  newdat$leafArea <- 10^predict(ancova.full,newdat)
  return(newdat$leafArea)
}


#- read and process the height and diameter dataset. Estimate total mass based on size by calling predict_mass()
return_size_mass <- function(model_flag="simple"){
  
  #- read in the data, do a few conversions
  dat <- read.csv("C:/Repos/GLAHD/Data/HeightDiam/GHS39_GLAHD_MAIN_HEIGHT&DIAMETER_20141106-20150105_L1.csv")
  dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
  dat$d2h <- with(dat,(Diameter/10)^2*(Height)) #calculate d2h in cm3
  
  linkdf <- data.frame(Species = c("CAM","TER","BOT","BRA","LONG","PEL","PLAT","SMIT"),Range= c("wide","wide",rep("narrow",6)))
  dat <- merge(dat,linkdf,by="Species")
  
  # get rid of the initial harvest and periodic harvest data
  dat2 <- subset(dat,Date!=as.Date("2014-11-06") & Date != as.Date("2014-11-20"))
  dat2$Treatment <- factor(dat2$Treatment)
  dat2$Code <- factor(dat2$Code)
  dat2$Taxa <- factor(dat2$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                         "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
  
                 
  dat2 <- dat2[with(dat2,order(Taxa)),]
  
  dat2$TotMass <- predict_mass(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height,model_flag=model_flag)
  dat2$leafArea <- predict_LA(Taxa=dat2$Taxa,Treat=dat2$Treatment,Diameter=dat2$Diameter,Height=dat2$Height)
  
  return(dat2)
}



# from Paine et al. 2012
# Self-start function for 2-parameter power-law fit (initial mass fixed by used). This self-start for this model is not provided by Pinhero and Bates (2000), so I wrote one. 
Init.pow2 <- function(mCall, LHS, data){
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4) {stop("Too few distinct x values to fit a power-law")}
  r    <-coef(lm(log(y) ~ x, xy))[2]    # Use the slope from a log fit to the data as an initial guess for r
  beta <- 0.9                             # give initial guess of beta as 0.9. don't use beta = 1, as it's undefined (1/0)
  value <- c(r, beta) 
  names(value) <- mCall[c("r", "beta")]
  return(value)
}


output.pow2.gnls <- function(fit, times, CI = F, alpha = 0.05, M0){
  params <- c(coef(fit), M0 = M0)
  r = params[1]; beta = params[2]
  fitted <- fit$fitted
  resid  <- fit$residual
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste((.(round(r * (1-beta), 3))*t)^.(round(1/(1-beta), 2))))
  mss  <- sum((fitted - mean(fitted))^2)
  rss  <- sum(resid^2)
  R2   <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  temp <- M0^(1-beta) + r*times*(1-beta)
  rates = data.frame(
    times = times,
    M    =     temp^(1/(1-beta)),
    AGR  = r * temp^(beta/(1-beta)))
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm  <-  r * rates$M^(beta-1)
  if(CI == T){
    cov <- fit$varBeta
    cov <- cbind(rbind(cov, 0), 0) # variance-covariance for M0 is 0 - it's fixed. 
    x <- data.frame(rmvnorm(n=1000, mean=params, sigma=cov))
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      temp.i   <- x.i$M0^(1-x.i$beta) + x.i$r * times *(1-x.i$beta)
      M[i,]    <-         temp.i^(1/(1-x.i$beta))
      AGR[i,]  <- x.i$r * temp.i^(x.i$beta/(1-x.i$beta))
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,] <- x.i$r * M[i,]^(x.i$beta - 1)
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}

# this function returns confidence envelopes around growth trajectories, and growth rates. 
summarizer <- function(dat, alpha){
  n <- length(dat)
  quantiles <- c(alpha/2, 1-(alpha/2))
  CIs <- data.frame(matrix(NA, ncol(dat[[1]]), n*2))
  names(CIs) <- paste(rep(names(dat), each = 2), c("lo", "hi"), sep = ".")
  for(i in 1:n){
    CIs[,(2*i-1):(2*i)] <- t(apply(dat[[i]],    2, quantile, quantiles, na.rm = T))
  }
  return(CIs)
}



#- this function reads in and returns the Rdark spot measurements. Returns a dataframe.
getRdark <- function(){
  
  #- read data
  gx <- read.csv(file="Data/GasEx/Rdark/GLAHD-Rdark-151214-compiled.csv")
  
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
  lm <- read.csv(file="Data/GasEx/Rdark/Dry Weight_22Dec2014_CRamig.csv")
  
  gx1$Code2 <- paste(gx1$Taxa,gx1$Pot,sep="")
  
  gx2 <- merge(gx1,lm,by.x="Code2",by.y="Code")
  gx2$Rmass <- with(gx2,R*Area/10000/leafDW*1000*1000)
  gx2$SLA <- with(gx2,(Area/(leafDW/1000))) 
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
  gx <- read.csv(file="C:/Repos/GLAHD/Data/GasEx/Asat/GHS39_GLAHD_MAIN_Asat_04122014_L1.csv")
  
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
