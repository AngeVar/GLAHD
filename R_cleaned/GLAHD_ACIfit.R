#-------------------------------------------------------------------------------------
#- This script provides ACi curve fit for GLAHD
#- 
#-------------------------------------------------------------------------------------

#- to install the newest version of plantecophys
# library(devtools)
# install_bitbucket("remkoduursma/plantecophys")


#- read data from csv files. These are rawish- they are the files right off the machine, but manually processed to remove remarks
#- and code which points to exclude from the A:Ci curve fitting.
aci1 <- read.csv(file="data/GHS39_GLAHD_MAIN_GX-ACI_20141208_L2.csv")
aci2 <- read.csv(file="data/GHS39_GLAHD_MAIN_GX-ACI_20141209_L2.csv")
aci3 <- read.csv(file="data/GHS39_GLAHD_MAIN_GX-ACI_20141210_L2.csv")

aci <- rbind(aci1,aci2,aci3)

#- get the first bit of the code (the taxa)
aci$Taxa <- unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=1,to=nrow(aci)*2,by=2)]
aci$Taxa <- factor(aci$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
aci$Pot <- as.numeric(unlist(strsplit(x=as.character(aci$Code),split="-"))[seq(from=2,to=nrow(aci)*2,by=2)])
aci$Treat <- ifelse(aci$Pot < 20, "Home","Warmed")

#- sort by code
aci <- aci[with(aci,order(Code)),]

#- extract the data to fit (exclude lines that I've manually decided to remove from the curve fitting)
aci.fit <- subset(aci,toFit==1)
aci.fit$PPFD<-aci.fit$PARi



#--- attempt to use parallel processing to speed up the A:Ci curve fitting process.
p.flag <- F # use parallel processing?
if(p.flag==T){
  library(doParallel)
  cl <- makeCluster(2)
  registerDoParallel(cl)
}

if(p.flag==T){
  starttime <- Sys.time()
  fits.p <- list()
  aci.fit.list <- split(aci.fit,aci.fit$Code)
  crap <- foreach(i=1:length(aci.fit.list)) %dopar% 
    plantecophys::fitaci(aci.fit.list[[i]],PPFD="PARi",Tcorrect=F,citransition=450)
  
  Sys.time()-starttime # took 55 seconds to fit all curves, rather than 1.6 mintues on 1 core.
  fits.params <- as.data.frame(do.call(rbind,lapply(crap,FUN=coef)))
  fits.params$Code <- levels(aci.fit$Code)
}#-----


#- fit the aci curves for each plant
if(p.flag==F){
  starttime <- Sys.time()
  fits <- fitacis(aci.fit,group="Code",fitmethod="bilinear",Tcorrect=F,citransition=450)
  fits.params <- coef(fits) #extract the Vcmax and Jmax parameters
  Sys.time()-starttime
  
}

#- assign other factor variables based on the "Code".
fits.params$Taxa <- unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=1,to=nrow(fits.params)*2,by=2)]
fits.params$Taxa <- factor(fits.params$Taxa,levels=c("ATER","BTER","ACAM","BCAM","CCAM","BOT","LONG","SMIT",
                                                     "CTER","DTER","ETER","DCAM","ECAM","FCAM","BRA","PEL","PLAT"))
fits.params$Pot <- as.numeric(unlist(strsplit(x=as.character(fits.params$Code),split="-"))[seq(from=2,to=nrow(fits.params)*2,by=2)])
fits.params$Treat <- ifelse(fits.params$Pot < 20, "Home","Warmed")


#----------------------------------------------------------------------------------
#- process Aci fits for statistical analysis

#- get the growth data, mostly just for the treatment codes
growth <- return_size_mass(model_flag="simple") # use common slope allometry ("simple") or taxa-specific slope ("complex")
growth2 <- summaryBy(d2h+TotMass+leafArea~Species+Treatment+Location+Taxa+Code+Range,keep.names=T,data=subset(growth,Date >= as.Date("2014-12-8") & Date <=as.Date("2014-12-20")))

#- merge size totalmass and leafarea data into dataframe with aci values
acifit <- merge(fits.params,growth2,by=c("Code","Taxa"))
acifit$JtoV <- with(acifit,Jmax/Vcmax)
acifit$Location <- factor(acifit$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
acifit$Sp_RS_EN <- as.factor(with(acifit,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
acifit$Prov_Sp_EN <- as.factor(with(acifit,paste(Taxa,Species)))

