#------------------------------------------------------------------------------------------------------------------------------
# This script reads and plots the height and diameter data for the GLAsshouse Heating and Distribution (GLAHD) project.
#------------------------------------------------------------------------------------------------------------------------------


#- load libraries from script
source("W:/WorkingData/GHS39/GLAHD/Share/R/loadLibraries.R")

#- read in the data, do a few conversions
dat2 <- return_size_mass()

#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
# Process and plot all taxa individually
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

                     
#- average across taxa and treatments
dat.m <- summaryBy(Height+Diameter+d2h+TotMass~Taxa+Location+Treatment+Date+Range,data=dat2,FUN=c(mean,standard.error))
dat.m$taxa_loc <- as.factor(paste(dat.m$Location,dat.m$Taxa,sep="-"))
dat.m$taxa_loc <- factor(dat.m$taxa_loc,levels=c("S-ATER","S-BTER","S-ACAM","S-BCAM","S-CCAM","S-BOT","S-LONG","S-SMIT",
                                       "N-CTER","N-DTER","N-ETER","N-DCAM","N-ECAM","N-FCAM","N-BRA","N-PEL","N-PLAT"))
#- split into list across all taxa for plotting
dat.l <- split(dat.m,dat.m$taxa_loc)

#------------------------------------------------------------------------------
#- plot Total mass for each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,90)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
palette(c("black","red"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(TotMass.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$TotMass.mean,
                                  SE=toplot$TotMass.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=3,line=-1.5)
  if (i==length(dat.l)){
    mtext(text="Total plant mass (g)",side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/TotalMass_all_taxa.pdf")

#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#- plot heights for each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,200)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
palette(c("black","red"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(Height.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$Height.mean,
               SE=toplot$Height.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(dat.l)){
    mtext(text="Height (cm)",side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Heights_all_taxa.pdf")

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#- plot diameters for each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,17)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(Diameter.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$Diameter.mean,
                                  SE=toplot$Diameter.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=1,line=-1.5)
  if (i==length(dat.l)){
    mtext(text="Diameter (mm)",side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Diameters_all_taxa.pdf")

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#- plot d2h for each taxa on a separate panel. Results in a huge figure
windows(30,30)
par(mfrow=c(5,4), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,320)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(d2h.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$d2h.mean,
                                  SE=toplot$d2h.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Taxa,sep="-"),side=3,line=-1.5)
  if (i==length(dat.l)){
    mtext(text=expression(d^2*h~(cm^3)),side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/d2h_all_taxa.pdf")

#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
# end
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
# Process and plot narrow vs. wide in N vs. S

#- average across location (N vs. S), treatment, and range (wide vs. narrow)
dat.m <- summaryBy(Height+Diameter+TotMass+d2h~Location+Treatment+Date+Range,data=dat2,FUN=c(mean,standard.error))

#- split into list across all taxa for plotting
dat.l <- split(dat.m,as.factor(paste(dat.m$Location,dat.m$Range,sep="-")))

#- plot height
windows(30,30)
par(mfrow=c(2,2), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,200)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(Height.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$Height.mean,
                                  SE=toplot$Height.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Range,sep="-"),side=1,line=-1.5)
  if (i==length(dat.l)){
    mtext(text=expression(Height~(cm)),side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Height_groups.pdf")


#- plot diameter
windows(30,30)
par(mfrow=c(2,2), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,15)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(Diameter.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$Diameter.mean,
                                  SE=toplot$Diameter.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Range,sep="-"),side=1,line=-1.5)
  if (i==length(dat.l)){
    mtext(text=expression(Diameter~(mm)),side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Diameter_groups.pdf")



#- plot d2h
windows(30,30)
par(mfrow=c(2,2), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,275)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(d2h.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$d2h.mean,
                                  SE=toplot$d2h.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Range,sep="-"),side=3,line=-1.5)
  if (i==length(dat.l)){
    mtext(text=expression(d^2*h~(cm^3)),side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
    
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/d2h_groups.pdf")

#- plot total mass

windows(30,30)
par(mfrow=c(2,2), mar=c(0.3,2,0.3,0.8), oma=c(5,6,2,2.5))
ylims <- c(0,90)
xlims <- c(as.Date("2014-11-1"),as.Date("2015-1-10"))
palette(c("black","red"))
for (i in 1:length(dat.l)){
  toplot <- dat.l[[i]]
  plotBy(TotMass.mean~Date|Treatment,data=toplot,type="b",pch=15,xlim=xlims,ylim=ylims,
         ylab="H",xlab="",legend=F,
         panel.first=adderrorbars(x=toplot$Date,y=toplot$TotMass.mean,
                                  SE=toplot$TotMass.standard.error,direction="updown"))
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="month"),labels=F)
  mtext(text=paste(toplot$Location,toplot$Range,sep="-"),side=3,line=-1.5)
  if (i==length(dat.l)){
    mtext(text="Total plant mass (g)",side=2,outer=T,line=2,cex=2)
    mtext(text="Date",side=1,outer=T,line=2,cex=2)
  }
}
dev.copy2pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/TM_groups.pdf")
#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
#- plot height and diameter over time for each individual plant, giving each plant a page on a pdf. Results in a huge pdf.

pdf(file="W:/WorkingData/GHS39/GLAHD/Share/Output/Height_diam_eachpot.pdf")


par(mfrow=c(2,1),mar=c(2,4,1,1),oma=c(4,1,2,1))
xlims <- c(min(dat2$Date)-2,max(dat2$Date)+2)

dat.l.all <- split(dat2,dat2$Code)
for (i in 1:length(dat.l.all)){
  plotBy(Height~Date,data=dat.l.all[[i]],type="b",xlim=xlims,xlab="",legend=F)
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="week"),labels=F)
  plotBy(Diameter~Date,data=dat.l.all[[i]],type="b",xlim=xlims,legend=F)
  axis.Date(side=1,at=seq(from=xlims[1],to=xlims[2],by="week"),labels=T)
  
  title(main=dat.l.all[[i]]$Code[1],outer=T)

}

dev.off()

#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------
