###------------------------------------------------------------------------------------------------------------------
#output from GAM model:"Code" "Time" "dydt" "predMass" "Species" "Treatment" "Location" "Taxa" "Range" "AGR" 
source("R/GLAHD_gamfits.R") 
library(effects)
###------------------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------------------
#output from allometries: "Species" "Date" "Treatment" "Location" "Taxa" "Pot" "Code" Height" "Diameter" "Comment" "ToFit" "d2h" "Range" "TotMass" "leafArea" "leafMass" "stemMass" "rootMass" "Time" 
dat2 <- return_size_mass_all(model_flag="complex") # use common slope allometry ("simple") or taxa-specific slope ("complex")
dat2$Time <- as.numeric(dat2$Date-(min(dat2$Date)-1)) #finds first date and labels it as Time 1 i.e. 07112014 is Day 1

#- remove data with fewer than 6 observations through time
obs <- unname(table(dat2$Code)); names <- names(table(dat2$Code));keeps <- names[which(obs>6)]
dat3 <- subset(dat2,Code %in% keeps) # subset dataframe
dat3$Code <- factor(dat3$Code)
dat3$lnTotMass <- log(dat3$TotMass)

#- Calculate some new variables
dat3$LMF<- with(dat3, leafMass/TotMass)             #Leaf Mass Fraction
dat3$SMF<- with(dat3, stemMass/TotMass)             #Stem Mass Fraction
dat3$RMF<- with(dat3, rootMass/TotMass)             #Root Mass Fraction
dat3$LMA<- with(dat3, leafMass/leafArea)            #Leaf Mass per Area
dat3$SLA<- with(dat3, leafArea/leafMass)            #Specific leaf area
dat3$LAR<- with(dat3, leafArea/TotMass)             #Leaf Area Ratio
dat3$RSR<- with(dat3, rootMass/(stemMass+leafMass)) #Root:Shoot Ratio

rate <- merge(gamfits2,dat3,by=c("Code","Time","Species","Location","Treatment","Taxa","Range"))
rate$ULR <- with(rate,dydt/LAR)
###------------------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------------------
#Harvested data:"Species Date Treatment Location Taxa Pot Code Height Diameter Comment ToFit d2h Range TotMass leafArea leafMass stemMass rootMass Time

data <- read.csv("data/GHS39_GLAHD_MAIN_BIOMASS_20141106-20150116_L1.csv")
data$Date <- as.Date(data$Date)
data$Totmass <- base::rowSums(data[,11:13]) #total mass is the sum of leaf, stem, and root mass
data$d2h <- with(data,(Diameter/10)^2*(Height)) #calculate d2h in cm3
data$Treat <- as.factor(ifelse(data$Pot < 20, "Home",
                               ifelse(data$Pot>=40,"Pre","Warmed")))
data$Taxa <- factor(data$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                       "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
data[24,] <- NA # get rid of BOT-45, which has crazy high LAR and SLA
data$Range<- as.factor(ifelse(data$Species == "TER"|data$Species == "CAM", "wide","narrow"))
data2<-droplevels(subset(data,Treat!="Pre"))

###------------------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------------------
#Gasexchange data

Rdark<-getRdark()                                  
Asat<- getAsat()                                   
source("R_cleaned/GLAHD_ACIfit.R") #this takes a while     

SLARd<- Rdark[,c("Code","SLA")]                    #extract the codes we have leaf SLA for
acifits<-merge(acifit,SLARd, by="Code")            #give SLA to the codes that we have Jmax and Vcmax data for
Asatm<-merge(Asat,SLARd, by="Code")                #give SLA to the codes that we have Asat data for
###------------------------------------------------------------------------------------------------------------------


###------------------------------------------------------------------------------------------------------------------


#add nesting for mixed model
gamfits2$Location <- factor(gamfits2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2$Sp_RS_EN <- as.factor(with(gamfits2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2$Prov_Sp_EN <- as.factor(with(gamfits2,paste(Taxa,Species)))
gamfits2$Sp_Loc_EN <- as.factor(with(gamfits2,paste(Species,Location)))

rate$Location <- factor(rate$Location,levels=c("S","N")) 
rate$Sp_RS_EN <- as.factor(with(rate,paste(Species,Range)))   
rate$Prov_Sp_EN <- as.factor(with(rate,paste(Taxa,Species)))
rate$Sp_Loc_EN <- as.factor(with(rate,paste(Species,Location)))

data2$Location <- factor(data2$Location,levels=c("S","N")) 
data2$Sp_RS_EN <- as.factor(with(data2,paste(Species,Range)))   
data2$Prov_Sp_EN <- as.factor(with(data2,paste(Taxa,Species)))
data2$Sp_Loc_EN <- as.factor(with(data2,paste(Species,Location)))

###########################################################

#gamfits2 has data interpolated for every day

#rate has allometrically estimated data only for the days on which we measured Height and Diameter: 
#Day 1 11 20 25 32 39 46 53 60 plus dydt (RGR), AGR and predmass from the GAM 

#data2 has harvested, actually measured data
#Rdark has Rmass and SLA for gax exchange leaves
#Asat has Asat at growth T
#acifit has Vcmax, Jmax and JtoV
#acifits has Vcmax and Jmax on a mass basis

###########################################################
