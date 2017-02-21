###------------------------------------------------------------------------------------------------------------------
#output from GAM model:"Code" "Time" "dydt" "predMass" "Species" "Treatment" "Location" "Taxa" "Range" "AGR" 
source("R_cleaned/GLAHD_gamfits.R") 
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

dat <- read.csv("Data/GHS39_GLAHD_MAIN_GX-RVT_20150114-20150120_L2.csv")
dat$Date <- as.Date(dat$Date,format="%d/%m/%Y")
dat$R_area <- dat$Photo*-1

#- read in the leaf mass data, calculate R per unit leaf mass
leafmass <- read.csv("Data/GHS39_GLAHD_MAIN_BIOMASS-RVTLEAVES_20150120_L2.csv")
r <- merge(dat,leafmass,by="Code")
r$Taxa <- as.factor(sapply(1:nrow(r),function(i)strsplit(as.character(r$Code[i]),split="-")[[1]][1])) # get the "Taxa" factor out of Code. Note this could have been a loop

r$R_mass <- r$R_area/10000*r$Area/r$Leafmass*1000 #note mass was in g
rt <- subset(r,Taxa!="SMIT" & Taxa !="ACAM")
rt$Taxa <- factor(rt$Taxa,levels=c("BOT","BTER","BRA","CTER"))
rt$lnRmass <- log(rt$R_mass)
rt <- subset(rt,Code!="CTER-4" & Code!="CTER-22" & Code!="CTER-27") #CTER4 looks weird... peaks much much lower than others. CTER-22 also looks bad
rt$Code <- factor(rt$Code)

q10mods <- list() #simple q10
arrmods <- list() #simple arrhenius equation
mtmods <- list()  #modified Q10 to allow Q10 to decline with measurement temperature. Mark's preferred model?
polymods <- list() #polynomial equation popularized by Owen Atkin. Basically an arrhenius where the activation energy can decline with temperature


xvals <- seq(13,45,length=101)
rt.l <- split(rt,rt$Code)
for (i in 1:length(rt.l)){
  
  dat <- subset(rt.l[[i]],Tleaf<40)
  print(i)
  #fit the models
  q10mods[[i]] <- nls(R_mass~ Rref*Q10^((Tleaf-15)/10),data=dat,start=list(Rref=10,Q10=2))
  arrmods[[i]] <- nls(R_mass~Rref*exp((Ea/(0.008314*(273.15+15))*(1-((273.15+15))/(273.15+Tleaf)))),data=dat,start=list(Rref=10,Ea=45))
  mtmods[[i]] <- nls(R_mass~Rref*((x-y*(Tleaf-15)/2))^((Tleaf-15)/10),data=dat,start=list(x=3,y=0.03,Rref=10))
  polymods[[i]] <- lm(lnRmass~Tleaf+I(Tleaf^2),data=dat)
}


#get coefs
q10.params <- as.data.frame(do.call(rbind,lapply(q10mods,coef)))
arr.params <- as.data.frame(do.call(rbind,lapply(arrmods,coef)))
mt.params <- as.data.frame(do.call(rbind,lapply(mtmods,coef)))
poly.params <- as.data.frame(do.call(rbind,lapply(polymods,coef)))

#rename parameter columns
names(q10.params) <- paste("Q10.",names(q10.params),sep="")
names(arr.params) <- paste("Arr.",names(arr.params),sep="")
names(mt.params) <- paste("mt.",names(mt.params),sep="")
names(poly.params) <- c("poly.a","poly.b","poly.c")

#combine fitted parameters, add chambers and treatments
params <- cbind(q10.params,arr.params,mt.params,poly.params)
params$Code  <- as.factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Code[1]))) # get the Code. Note this could have been a loop
params$Treatment  <- as.factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Treatment[1]))) # get the Code. Note this could have been a loop
params$Taxa  <- factor(sapply(1:length(rt.l),function(i)as.character(rt.l[[i]]$Taxa[1])),levels=c("BTER","BOT","CTER","BRA")) # get the Code. Note this could have been a loop

#get AIC values
q10.AIC <- abs(as.data.frame(do.call(rbind,lapply(q10mods,AIC))))
arr.AIC <- abs(as.data.frame(do.call(rbind,lapply(arrmods,AIC))))
mt.AIC <- abs(as.data.frame(do.call(rbind,lapply(mtmods,AIC))))
poly.AIC <- abs(as.data.frame(do.call(rbind,lapply(polymods,AIC))))

#rename AIC columns
names(q10.AIC) <- "Q10.AIC"
names(arr.AIC) <- "arr.AIC"
names(mt.AIC) <- "mt.AIC"
names(poly.AIC) <- "poly.AIC"

fitted <- cbind(params,q10.AIC,arr.AIC,mt.AIC,poly.AIC)
fitted <- fitted[,11:17] # Note that Mark's modified Q10 tends to have the lowest AIC for each curve, but not always

# mark's modified Q10 has the lowest overall AIC.
#summaryBy(Q10.AIC+arr.AIC+mt.AIC+poly.AIC~1,FUN=sum,data=fitted)

# parameter estimates
params.m <- summaryBy(mt.x + mt.y + mt.Rref~Taxa+Treatment,FUN=c(mean,standard.error),data=params)
params$R25 <- with(params,mt.Rref*((mt.x-mt.y*(25-15)/2))^((25-15)/10))
params$Q10_25 <- with(params,(mt.x-mt.y*(25-15)/2))
params$R18.5 <- with(params,mt.Rref*((mt.x-mt.y*(18.5-15)/2))^((18.5-15)/10))
params$Q10_18.5 <- with(params,(mt.x-mt.y*(18.5-15)/2))
params$R28.5 <- with(params,mt.Rref*((mt.x-mt.y*(28.5-15)/2))^((28.5-15)/10))
params$Q10_28.5 <- with(params,(mt.x-mt.y*(28.5-15)/2))


tvals <- seq(from=13,to=40,by=1)
params.list <- split(params,params$Code)
predR <- list()
for(i in 1:length(params.list)){
  predR[[i]] <- data.frame("Code"=params.list[[i]]$Code[1],"Treatment"=params.list[[i]]$Treatment[1],
                           "Taxa"=params.list[[i]]$Taxa[1],"T"=tvals,
                           "R"=with(params.list[[i]],exp(poly.a+poly.b*tvals+poly.c*tvals^2)),
                           "Q10"=with(params.list[[i]],exp(10*(poly.b+2*poly.c*tvals))))
}
predR.df <- do.call(rbind,predR)
predR.df$logR <- log(predR.df$R)



###------------------------------------------------------------------------------------------------------------------


#add nesting for mixed model
gamfits2$Location <- factor(gamfits2$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2$Sp_RS_EN <- as.factor(with(gamfits2,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2$Prov_Sp_EN <- as.factor(with(gamfits2,paste(Taxa,Species)))
gamfits2$Sp_Loc_EN <- as.factor(with(gamfits2,paste(Species,Location)))

gamfits2mass$Location <- factor(gamfits2mass$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
gamfits2mass$Sp_RS_EN <- as.factor(with(gamfits2mass,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
gamfits2mass$Prov_Sp_EN <- as.factor(with(gamfits2mass,paste(Taxa,Species)))
gamfits2mass$Sp_Loc_EN <- as.factor(with(gamfits2mass,paste(Species,Location)))

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
