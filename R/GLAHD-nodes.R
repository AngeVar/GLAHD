#------------------------------------------------------------------------------------------------------------------------------
# This script looks at the measurements of node number in the GLAHD experiment
#------------------------------------------------------------------------------------------------------------------------------



#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")

#- read in the data
nodes <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/Nodes/GHS39_GLAHD_MAIN_NODE_NUMBER_20150112_L1.csv")
nodes$Taxa <- factor(nodes$Taxa,levels=c("BTER","ACAM","BOT","LONG",
                                         "CTER","DCAM","BRA","PLAT"))

#- plot node number
#- it looks like warming increased node # in the south but not the north
windows(12,12);par(mar=c(8,4,1,1))
boxplot(Nodes~Treatment*Taxa,data=nodes,col=c("blue","red"),las=2,ylab="Nodes (#)",ylim=c(0,70))
abline(v=8.5)
text(x=4,y=70,"South")
text(x=12,y=70,"North")
legend("bottomright",legend=c("Home","Warmed"),fill=c("blue","red"))
