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
