#------------------------------------------------------------------------------------------------------------------------------
#- Reads in and processes the leaf gas exchange data from the heatwave
#- specifically quantifies non-stomatal limitation to photosynthesis
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")


#- read data
gx <- read.csv(file="W:/WorkingData/GHS39/GLAHD/Share/Data/Heatwave/GLADH heatwave combined cleaned gasex wp data.csv")
gx$Ca <- 400
gx$Ci <- gx$cica*gx$Ca
gx$date <- as.Date(gx$date,format="%d/%m/%Y")

#- define function to return Vcmax based on a single point A:Ci curve as in Zhou et al. 2013 Tree Phys.
returnVcmax <- function(A,Ci,Tleaf){
  Tk <- Tleaf + 274.15
  R <- 8.314                                                       #universal gas constant in J mol-1 K-1
  gammastar <- 42.75*exp( (37830*(Tk-298))/(298*R*Tk) )            #calculate gamma star as in bernacchi et al. 2001
  Ko <- exp(20.3-36.38/(R/1000*Tk))                                     #calculate Ko as in Bernacchi et al. 2001
  Kc <- exp(38.05-79.43/(R/1000*Tk))                                    #calculate Kc as in Barnacchi et al. 2001
  Km <- Kc*(1+210/Ko)                                              # calculate Km. Assumes [O2] = 210 mmol mol-1
  
  Vcmax <- A*(Ci+Km)/(Ci-gammastar)                                #calculate apparent Vcmax. This assumes Rd is zero (or negligible)
  return(Vcmax)
}




#- get apparent Vcmax for each observation
gx$Vcmax_a <- returnVcmax(A=gx$photo,Ci=gx$Ci,Tleaf=gx$tleaf)

#- get the "maximum" Vcmax as the 95th percentile of the apparent Vcmax of the well-watered treatments
gx.list <- split(gx,gx$species)
species <- c()
Vcmax_max <- c()
Tleaf_at_Vcmax <- c()
for(i in 1:length(gx.list)){
  dat <- subset(gx.list[[i]],date==as.Date("2015-01-28"))
  species[i] <- as.character(dat$species[1])
  Vcmax_max[i] <- unname(quantile(dat$Vcmax_a,probs=0.5))
  Tleaf_at_Vcmax[i] <- mean(dat$tleaf)
  
}
df2 <- data.frame(species=species,Vcmax_max=Vcmax_max,Tleaf_at_Vcmax = Tleaf_at_Vcmax)

gx2 <- merge(gx,df2,by=c("species"))
gx3 <- subset(gx2,photo>-10)

gx3$Photo_a <- Photosyn(VPD=gx3$vpdl,Ca=gx3$Ca,PPFD=1800,Tleaf=gx3$tleaf,Ci=gx3$Ci,Vcmax=gx3$Vcmax_max,Tcorrect=T)$ALEAF
gx3$NSL <- with(gx3,photo/Photo_a)
gx3.5 <- gx3 #subset(gx3,NSL> -5 & NSL < 2)
gx4 <- summaryBy(NSL+Vcmax_a+photo+tleaf~species+date,data=gx3.5,FUN=c(mean,standard.error))




windows(20,35);par(mfrow=c(2,1))
plotBy(Vcmax_a.mean~date|species,data=gx4,pch=15,type="b",legendwhere="topright",ylab="Apparent Vcmax at leaf T",ylim=c(-50,300),
       panel.first=adderrorbars(x=gx4$date,y=gx4$Vcmax_a.mean,SE=gx4$Vcmax_a.standard.error,direction="updown"))
abline(h=0)
plotBy(tleaf.mean~date|species,data=gx4,pch=15,type="b",legend=F,ylab="Leaf temperature (deg C)")

# 
# #adjust to make initial value zero
# bot1 <- mean(subset(gx4,species=="bot" & date==as.Date("2015-01-28"))$NSL.mean)
# cam1 <- mean(subset(gx4,species=="cam" & date==as.Date("2015-01-28"))$NSL.mean)
# smit1 <- mean(subset(gx4,species=="smit" & date==as.Date("2015-01-28"))$NSL.mean)
# ter1 <- mean(subset(gx4,species=="ter" & date==as.Date("2015-01-28"))$NSL.mean)
# 
# initials <- summaryBy(NSL.mean~species,data=subset(gx4,date==as.Date("2015-01-28")),FUN=mean,keep.names=F)
# gx5 <- merge(gx4,initials,by="species")
# gx5$NSL <- with(gx5,NSL.mean-(NSL.mean.mean-1))

#- this doesn't really work, given that Vcmax is so far above the optimum
plotBy(NSL.mean~date|species,data=gx4,pch=15,legendwhere="bottomright",type="b",
       panel.first=adderrorbars(x=gx4$date,y=gx4$NSL.mean,SE=gx4$NSL.standard.error,direction="updown"))
