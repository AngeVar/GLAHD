#------------------------------------------------------------------------------------------------------------------------------
# This script looks at the measurements of node number in the GLAHD experiment
#------------------------------------------------------------------------------------------------------------------------------



#- load libraries from script
source("R/loadLibraries.R")

#- read in the data
nodes <- read.csv("Data/Nodes/GHS39_GLAHD_MAIN_NODE_NUMBER_20150112_L1.csv")
nodes$Taxa <- factor(nodes$Taxa,levels=c("BTER","ACAM","BOT","LONG",
                                         "CTER","DCAM","BRA","PLAT"))
nodes$Range <- factor(ifelse(nodes$Species == "TER"|nodes$Species == "CAM", "wide","narrow"))

#- plot node number
#- it looks like warming increased node # in the south but not the north
windows(12,12);par(mar=c(8,4,1,1))
boxplot(Nodes~Treatment*Taxa,data=nodes,col=c("blue","red"),las=2,ylab="Nodes (#)",ylim=c(0,70))
abline(v=8.5)
text(x=4,y=70,"South")
text(x=12,y=70,"North")
legend("bottomright",legend=c("Home","Warmed"),fill=c("blue","red"))

nodes$Location <- factor(nodes$Location,levels=c("S","N")) # relevel Location so that "S" is the first level and "N" is the second
nodes$Sp_RS_EN <- as.factor(with(nodes,paste(Species,Range)))   # use "explicit nesting" to create error terms of species:rangesize and prov:species:rangesize
nodes$Prov_Sp_EN <- as.factor(with(nodes,paste(Taxa,Species)))


fm1n <- lme(Nodes~Treatment*Location*Range,random=list(~1|Sp_RS_EN,~1|Prov_Sp_EN),data=nodes)#, method="ML")
plot(fm1n,resid(.,type="p")~fitted(.) | Treatment,abline=0)     #resid vs. fitted for each treatment. Is variance approximately constant?
plot(fm1n,Nodes~fitted(.)|Species,abline=c(0,1))               #predicted vs. fitted for each species
plot(fm1n,Nodes~fitted(.),abline=c(0,1))                       #overall predicted vs. fitted
qqnorm(fm1n, ~ resid(., type = "p"), abline = c(0, 1))          #qqplot to assess normality of residuals
hist(fm1n$residuals[,1])
anova(fm1n)    
