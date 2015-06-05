#This script reads absorbance data from photospec and calculates TNC#

library(plyr)
library(HH)
library(doBy)

#--------------------------------------------------------------------------------------

#First run was only standard material

stdrun<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_intstd_5May2015.csv", header = F)
stdrun <- rename(stdrun[-c(1:30,43:45,160:162),c(1:3)], c("V1"="ID", "V2"="Class","V3"="Abs490"))
stdrun$ID <- as.numeric(as.character(stdrun$ID));stdrun$Abs490 <- as.numeric(as.character(stdrun$Abs490));stdrun$Class <- droplevels(stdrun$Class)
stdrun$rep<- as.factor(rep(c(1,2,3)))
#plot calibration curve
data <- subset(stdrun, Class=="c");plot(data$Abs490,data$ID)
cal <-lm(ID~Abs490, data=data);abline(cal)
data$curve <- as.factor(c(rep(1,18),rep(2,18),rep(3,18)));ancova(ID ~ Abs490 + curve, data=data)

#Curves are not different - use slope off all datapoints and calculate glucose concentration (g/ml)
stdrun$Vol <- as.numeric(ifelse(stdrun$Class == "su","6.5", "5"))
stdrun$glugml <- as.numeric(cal$coef[2])*stdrun$Abs490+as.numeric(cal$coef[1])

#calculate the amount of glucose in sample (mg/g)
stdmass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/refmass.csv")
stdmass <- rename(stdmass, c("ID"="Code"));stdmass <- rename(stdmass, c("Order"="ID"))
stdmass<- subset(stdmass, is.na(date))
std <- merge(stdrun, stdmass, by="ID")
std$glumgg <- (std$glugml*5*std$Vol)/std$Mass
std$prcglu<- std$glugml*std$Vol*5*0.1/std$Mass

#correct for contamination by subtracting the average absorbance of the blanks
blanks <- subset(std, ID < 0)
bstd<-summaryBy(Abs490+glugml~ID+Class, data=blanks, FUN=c(mean, sd)); bstdst<-subset(bstd,Class=="st")

std$corrabs <- ifelse(std$Class == "st", std$Abs490-0.07365, std$Abs490)
std$corrglugml <- as.numeric(cal$coef[2])*std$corrabs+as.numeric(cal$coef[1])
std$corrglumgg <- (std$corrglugml*5*std$Vol)/std$Mass
std$corrprcglu<- std$corrglugml*std$Vol*5*0.1/std$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv<- std[std$ID>0,c(1,2,4,6,7,12,13,16,17)]
TNCsu<- csv[csv$Class == "su",];TNCst<- csv[csv$Class == "st",]
TNC <- merge(TNCsu,TNCst, by=c("ID","rep"));
TNC$TNC <- TNC$corrprcglu.x+TNC$corrprcglu.y;
tnc<-TNC[,c(1,2,17)]
stdcsv<- merge(csv,tnc, by=c("ID","rep"))

#---------------------------------------------------------------------------------------

#Round 1: 20150507

r1<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round1_7May2015.csv", header = F)
r1 <- rename(r1[-c(1:36,115:120),c(1:3)], c("V1"="ID", "V2"="Class","V3"="Abs490"))
r1$ID <- as.numeric(as.character(r1$ID));r1$Abs490 <- as.numeric(as.character(r1$Abs490));r1$Class <- droplevels(r1$Class)
r1$rep <- as.factor(rep(c(1,2,3)))

#plot calibration curve
data1 <- subset(r1, Class=="c" &ID !=100); plot(data1$Abs490,data1$ID)
calcurve <-lm(ID~Abs490, data=data1);abline(calcurve)
data1$curve <- as.factor(c(rep(1,15),rep(2,15),rep(3,15)));ancova(ID ~ Abs490 + curve, data=data1)

#Curves are not different - use slope off all datapoints and calculate glucose concentration (g/ml)
r1$Vol <- as.numeric(ifelse(r1$Class == "su","6.75", "5"))
r1$glugml <- as.numeric(calcurve$coef[2])*r1$Abs490+as.numeric(calcurve$coef[1])

#calculate the amount of glucose in sample (mg/g)
mass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/TNCMass.csv")
mass <- rename(mass[c(1:144),c(1,5,6)], c("ID"="Code"));mass <- rename(mass, c("Order"="ID"))
run1 <- subset(r1, Class== "su"| Class== "st"); round1 <- merge(run1, mass, by="ID")
round1$glumgg <- (round1$glugml*5*round1$Vol)/round1$Mass
round1$prcglu<- (round1$glugml*round1$Vol*5*0.1)/round1$Mass

#correct for contamination by subtracting the average absorbance of the blanksmple
blanks <- subset(run1, ID < 0)
#remove -2.2 - twice as high as the others reps.
blanks<- subset(blanks,as.integer(rownames(blanks))!=233)
summaryBy(Abs490+glugml~Class, data=blanks, FUN=c(mean, sd))

round1$corrabs <- ifelse(round1$Class == "st", round1$Abs490-0.064060000, round1$Abs490 -0.004733333)
round1$corrglugml <- as.numeric(calcurve$coef[2])*round1$corrabs+as.numeric(calcurve$coef[1])
round1$corrglumgg <- (round1$corrglugml*5*round1$Vol)/round1$Mass
round1$corrprcglu<- round1$corrglugml*round1$Vol*5*0.1/round1$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv1<- round1[round1$ID>0,c(1,2,4,6,7,9,10,13,14)]
TNCsu1<- csv1[csv1$Class == "su",];TNCst1<- csv1[csv1$Class == "st",]
TNC1 <- merge(TNCsu1,TNCst1, by=c("ID","rep"));
TNC1$TNC <- TNC1$corrprcglu.x+TNC1$corrprcglu.y;
tnc1<-TNC1[,c(1,2,17)]
r1csv<- merge(csv1,tnc1, by=c("ID","rep"))


#------------------------------------------------------------------------------

#Round 2: 20150513

r2<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round2_sugars+starch_13May2015.csv", header = F);r2[90,1]<-"34"
r2<- rename(r2[-c(1:36,55:57,154:156,268:270),c(1:3)],c("V1"="ID", "V2"="Class","V3"="Abs490"))
r2$Class <- c(rep("c",18),rep("su",78), rep("c",36), rep("st",75),rep("c",18))
r2$Class <- as.factor(r2$Class); r2$ID <- as.numeric(as.character(r2$ID)); r2$Abs490 <- as.numeric(as.character(r2$Abs490))
r2$rep <- as.factor(rep(c(1,2,3)))

#plot calibration curve
data2 <- subset(r2, Class=="c"); plot(data2$Abs490,data2$ID)
calcurve2 <-lm(ID~Abs490, data=data2);abline(calcurve2)
data2$curve <- as.factor(c(rep(1,18),rep(2,18),rep(3,18),rep(4,18)));ancova(ID ~ Abs490 + curve, data=data2) 
ancova(ID ~ Abs490 + curve, data=droplevels(subset(data2,curve=="3"|curve=="4")))

# Curves are different - use slope off curve 3and4 to calculate glucose concentration (g/ml)
curve<- subset(data2,curve ==3|curve==4)
plot(curve$Abs490,curve$ID);calcurve2 <-lm(ID~Abs490, data=curve);abline(calcurve2)
r2$Vol <- as.numeric(ifelse(r2$Class == "su","6.75", "5"))
r2$glugml <- as.numeric(calcurve2$coef[2])*r2$Abs490+as.numeric(calcurve2$coef[1])

#calculate the amount of glucose in sample (mg/g)
mass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/TNCMass.csv")
mass <- rename(mass[c(1:144),c(1,5,6)], c("ID"="Code"));mass <- rename(mass, c("Order"="ID"))
run2 <- subset(r2, Class== "su"| Class== "st"); round2 <- merge(run2, mass, by="ID")
round2$glumgg <- (round2$glugml*5*round2$Vol)/round2$Mass
round2$prcglu<- (round2$glugml*round2$Vol*5*0.1)/round2$Mass

#correct for contamination by subtracting the average absorbance of the blanksmple
blanks <- subset(run2, ID < 0)
summaryBy(Abs490+glugml~Class, data=blanks, FUN=c(mean, sd))

round2$corrabs <- ifelse(round2$Class == "st", round2$Abs490-0.09566667, round2$Abs490 -0.01221667)
round2$corrglugml <- as.numeric(calcurve$coef[2])*round2$corrabs+as.numeric(calcurve$coef[1])
round2$corrglumgg <- (round2$corrglugml*5*round2$Vol)/round2$Mass
round2$corrprcglu<- round2$corrglugml*round2$Vol*5*0.1/round2$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv2<- round2[round2$ID>0,c(1,2,4,6,7,9,10,13,14)]
TNCsu2<- csv2[csv2$Class == "su",];TNCst2<- csv2[csv2$Class == "st",]
TNC2 <- merge(TNCsu2,TNCst2, by=c("ID","rep"));
TNC2$TNC <- TNC2$corrprcglu.x+TNC2$corrprcglu.y;
tnc2<-TNC2[,c(1,2,17)]
r2csv<- merge(csv2,tnc2, by=c("ID","rep"))

#------------------------------------------------------------------------------

#Round 3: 20150519

r3st<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round3+4_starch_19May2015.csv", header = F)
r3su<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round3+4_sugars_19May2015.csv", header = F)
r3st<- r3st[-c(1:33,202:207,223:225),c(1:3)]
r3su<- r3su[-c(1:33,211:213),c(1:3)]
r3st<- rename(r3st, c("V1"="ID", "V2"="Class","V3"="Abs490"))
r3su <- rename(r3su, c("V1"="ID", "V2"="Class","V3"="Abs490"))
r3st$Class <- c(rep("c",18),rep("st",141), rep("c",18),rep("oldc",6))
r3su$Class <- c(rep("c",18),rep("su",141), rep("c",18))
r3 <- rbind(r3st,r3su)
r3$Class <- as.factor(r3$Class); r3$ID <- as.numeric(as.character(r3$ID)); r3$Abs490 <- as.numeric(as.character(r3$Abs490))
r3$rep <- as.factor(rep(c(1,2,3)))

#plot calibration curve
data3 <- subset(r3, Class=="c"); plot(data3$Abs490,data3$ID)
calcurve3 <-lm(ID~Abs490, data=data3);abline(calcurve3)
data3$curve <- as.factor(c(rep(1,18),rep(2,18),rep(3,18),rep(4,18)));ancova(ID ~ Abs490 + curve, data=data3) 

#Curves are not different - use slope off all datapoints and calculate glucose concentration (g/ml)
r3$Vol <- as.numeric(ifelse(r3$Class == "su","6.75", "5"))
r3$glugml <- as.numeric(calcurve3$coef[2])*r3$Abs490+as.numeric(calcurve3$coef[1])

#calculate the amount of glucose in sample (mg/g)
mass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/TNCMass.csv")
mass <- rename(mass[c(1:144),c(1,5,6)], c("ID"="Code"));mass <- rename(mass, c("Order"="ID"))
run3 <- subset(r3, Class== "su"| Class== "st"); round3 <- merge(run3, mass, by="ID")
round3$glumgg <- (round3$glugml*5*round3$Vol)/round3$Mass
round3$prcglu<- (round3$glugml*round3$Vol*5*0.1)/round3$Mass

#correct for contamination by subtracting the average absorbance of the blank sample
blanks <- subset(run3, ID < 0)
summaryBy(Abs490+glugml~Class, data=blanks, FUN=c(mean, sd))

round3$corrabs <- ifelse(round3$Class == "st", round3$Abs490-0.084966667, round3$Abs490 -0.003133333)
round3$corrglugml <- as.numeric(calcurve$coef[2])*round3$corrabs+as.numeric(calcurve$coef[1])
round3$corrglumgg <- (round3$corrglugml*5*round3$Vol)/round3$Mass
round3$corrprcglu<- round3$corrglugml*round3$Vol*5*0.1/round3$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv3<- round3[round3$ID>0,c(1,2,4,6,7,9,10,13,14)]
TNCsu3<- csv3[csv3$Class == "su",];TNCst3<- csv3[csv3$Class == "st",]
TNC3 <- merge(TNCsu3,TNCst3, by=c("ID","rep"));
TNC3$TNC <- TNC3$corrprcglu.x+TNC3$corrprcglu.y;
tnc3<-TNC3[,c(1,2,17)]
r3csv<- merge(csv3,tnc3, by=c("ID","rep"))


#---------------------------------------------------------------------------------

#Round 4: 20150521

r4su<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round5_sugars_21May2015.csv", header = F)
r4st<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_round5_starch_21May2015.csv", header = F)
r4su<- r4su[-c(1:33,94:99,208:213,223:255),c(1:3)]
r4st<- r4st[-c(1:33,124:129,235:237),c(1:3)]
r4su <- rename(r4su, c("V1"="ID", "V2"="Class","V3"="Abs490"))
r4st<- rename(r4st, c("V1"="ID", "V2"="Class","V3"="Abs490"))
r4su$Class <- c(rep("c",18),rep("su",141), rep("c",18))
r4st$Class <- c(rep("c",18),rep("st",141), rep("c",18),rep("oldc",18))
r4 <- rbind(r4su,r4st);r4$ID<-as.character(r4$ID)
r4$ID<-ifelse(r4$ID=="-8","62",r4$ID)
r4$ID<-ifelse(r4$ID=="5h100","100",ifelse(r4$ID=="5h80","80",ifelse(r4$ID=="5h60","60",ifelse(r4$ID=="5h40","40",ifelse(r4$ID=="5h20","20",ifelse(r4$ID=="5h0","0",r4$ID))))))
r4$Class <- as.factor(r4$Class); r4$ID <- as.numeric(r4$ID); r4$Abs490 <- as.numeric(as.character(r4$Abs490))
r4$rep <- as.factor(rep(c(1,2,3)));

#plot calibration curve
data4 <- subset(r4, Class=="c"); plot(data4$Abs490,data4$ID)
calcurve4 <-lm(ID~Abs490, data=data4);abline(calcurve4)
data4$curve <- as.factor(c(rep(1,18),rep(2,18),rep(3,18),rep(4,18)));ancova(ID ~ Abs490 + curve, data=data4) 
#data4$curve <- as.factor(c(rep(1,36),rep(2,36)));ancova(ID ~ Abs490 + curve, data=data4)
# ancova(ID ~ Abs490 + curve, data=droplevels(subset(data4, curve == "1"|curve =="2")))
# ancova(ID ~ Abs490 + curve, data=droplevels(subset(data4, curve == "3"|curve =="4")))

#Curves are different - use slope 1and2 to calculate glucose concentration (g/ml) in sugar and slope 3and4 for starch
calcurve4.1<-lm(ID~Abs490, data=droplevels(subset(data4, curve == "1"|curve =="2")))
calcurve4.2<-lm(ID~Abs490, data=droplevels(subset(data4, curve == "3"|curve =="4")))
r4$Vol <- as.numeric(ifelse(r4$Class == "su","6.75", "5"))
r4$glugml <- ifelse(r4$Class == "su",as.numeric(calcurve4.1$coef[2])*r4$Abs490+as.numeric(calcurve4.1$coef[1]),
                     as.numeric(calcurve4.2$coef[2])*r4$Abs490+as.numeric(calcurve4.2$coef[1]))

#calculate the amount of glucose in sample (mg/g)
mass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/TNCMass.csv")
mass <- rename(mass[c(1:144),c(1,5,6)], c("ID"="Code"));mass <- rename(mass, c("Order"="ID"))
run4 <- subset(r4, Class== "su"| Class== "st"); round4 <- merge(run4, mass, by="ID")
round4$glumgg <- (round4$glugml*5*round4$Vol)/round4$Mass
round4$prcglu<- (round4$glugml*round4$Vol*5*0.1)/round4$Mass

#correct for contamination by subtracting the average absorbance of the blank sample
blanks <- subset(run4, ID < 0)
#remove -6.3 - six times as high as the others reps.
blanks2<- subset(blanks,as.integer(rownames(blanks))!=198)
summaryBy(Abs490+glugml~Class, data=blanks2, FUN=c(mean, sd))

round4$corrabs <- ifelse(round4$Class == "st", round4$Abs490-0.06678333, round4$Abs490 -0.01030000)
round4$corrglugml <- as.numeric(calcurve$coef[2])*round4$corrabs+as.numeric(calcurve$coef[1])
round4$corrglumgg <- (round4$corrglugml*5*round4$Vol)/round4$Mass
round4$corrprcglu<- round4$corrglugml*round4$Vol*5*0.1/round4$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv4<- round4[round4$ID>0,c(1,2,4,6,7,9,10,13,14)]
TNCsu4<- csv4[csv4$Class == "su",];TNCst4<- csv4[csv4$Class == "st",]
TNC4 <- merge(TNCsu4,TNCst4, by=c("ID","rep"));
TNC4$TNC <- TNC4$corrprcglu.x+TNC4$corrprcglu.y;
tnc4<-TNC4[,c(1,2,17)]
r4csv<- merge(csv4,tnc4, by=c("ID","rep"))

#------------------------------------------------------------------------------
#Practice run: 20150602

stdrun2<-read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/RAW/TNC_intstd_2june2015.csv", header = F)
stdrun2 <- rename(stdrun2[-c(1:27,105,154,174,225),c(1:3)], c("V1"="ID", "V2"="Class","V3"="Abs490"))
stdrun2$ID <- as.numeric(as.character(stdrun2$ID));stdrun2$Abs490 <- as.numeric(as.character(stdrun2$Abs490));stdrun2$Class <- droplevels(stdrun2$Class)
stdrun2$rep<- as.factor(c(rep(c(1,2,3),25),1,2, rep(c(1,2,3),16),1,2,rep(c(1,2,3),5),1,2,rep(c(1,2,3),16),1,2,rep(c(1,2,3),18)))
stdrun2$Class <- as.factor(c(rep("c",18),rep("su",107), rep("st",105),rep("c",18)))

#plot calibration curve
data0 <- subset(stdrun2, Class=="c");plot(data0$Abs490,data0$ID)
cal <-lm(ID~Abs490, data0=data0);abline(cal)
data0$curve <- as.factor(c(rep(1,18),rep(2,18)));ancova(ID ~ Abs490 + curve, data=data0)

#Curves are not different - use slope off all datapoints and calculate glucose concentration (g/ml)
stdrun2$Vol <- as.numeric(ifelse(stdrun2$Class == "su","6.5", "5"))
stdrun2$glugml <- as.numeric(cal$coef[2])*stdrun2$Abs490+as.numeric(cal$coef[1])

#calculate the amount of glucose in sample (mg/g)
stdmass2 <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/refmass.csv")
stdmass2 <- rename(stdmass2, c("ID"="Code"));stdmass2 <- rename(stdmass2, c("Order"="ID"))
mass2 <- stdmass2[!is.na(stdmass2$date),]
stdrun2<- stdrun2[stdrun2$Class!="c",]
std2 <- merge(stdrun2, mass2, by="ID")
std2$glumgg <- (std2$glugml*5*std2$Vol)/std2$Mass
std2$prcglu<- std2$glugml*std2$Vol*5*0.1/std2$Mass

#correct for contamination by subtracting the average absorbance of the blanks
std2$ID <- as.factor(std2$ID)
blanks0 <- std2[std2$ID == 1|std2$ID ==12|std2$ID ==13|std2$ID ==24|std2$ID ==25|std2$ID ==36,]
bstd2<-summaryBy(Abs490+glugml~ID+Class, data=blanks0, FUN=c(mean, sd))
blanks0[order(blanks0$Class),]
#Starch Blanks 12 and 13 unusually high, mean of remaining blanks used for correction
blanks00<- blanks0[blanks0$ID!= 12 & blanks0$ID!= 13,]
summaryBy(Abs490+glugml~Class, data=blanks00, FUN=c(mean, sd))

std2$corrabs <- ifelse(std2$Class == "st", std2$Abs490-0.09088000, std2$Abs490-0.01321667)
std2$corrglugml <- as.numeric(cal$coef[2])*std2$corrabs+as.numeric(cal$coef[1])
std2$corrglumgg <- (std2$corrglugml*5*std2$Vol)/std2$Mass
std2$corrprcglu<- std2$corrglugml*std2$Vol*5*0.1/std2$Mass

#print glucose concentration (mg/g) and (%) and TNC (%)
csv0<- std2[,c(1,2,4,6,7,12,13,16,17)]
TNCsu<- csv0[csv0$Class == "su",];TNCst<- csv0[csv0$Class == "st",]
TNC <- merge(TNCsu,TNCst, by=c("ID","rep"));
TNC$TNC <- TNC$corrprcglu.x+TNC$corrprcglu.y;
tnc<-TNC[,c(1,2,17)]
stdcsv2<- merge(csv0,tnc, by=c("ID","rep"))
stdcsv2 <- droplevels(subset(stdcsv2,ID != 1 & ID != 12&ID != 13&ID != 24&ID != 25&ID != 36))

#------------------------------------------------------------------------------
#Is standard material results reproducible

stdcsv$Code <- paste(stdcsv$Code,"std1",sep="-")
stdcsv2$Code <- paste(stdcsv2$Code,"std2",sep="-")

stdcsv$run<- -1;stdcsv2$run <- -2; r1csv$run <- 1;r2csv$run <- 2;r3csv$run <- 3;r4csv$run <- 4;all <- rbind(stdcsv,stdcsv2,r1csv, r2csv,r3csv,r4csv)
refs<- rbind(stdcsv,stdcsv2, subset(r1csv,ID == 1|ID==13|ID==25),subset(r2csv,ID == 26|ID==37|ID==48),subset(r3csv,ID == 49|ID==73|ID==97),subset(r4csv,ID == 62|ID==74|ID==86))
#refs2<-rbind(stdcsv, subset(r1csv,ID == 1|ID==13|ID==25),subset(r2csv,ID == 26|ID==37|ID==48),subset(r3csv,ID == 49|ID==73|ID==97),subset(r4csv,ID == 62|ID==74|ID==86))
refs$run <- as.factor(refs$run);refs$Code <- as.factor(refs$Code)
sd(refs$TNC)

sum<- summaryBy(TNC~Code+run, data=refs, FUN=c(mean,sd)); sum$Code <- as.factor(sum$Code); sum$run<- as.factor(sum$run)
sum <- sum[order(sum$run),]

#lablist<-sum$Code
lablist<- c(rep("LH",10),rep("JED",10),rep("RO",10),rep("AV",22))
means <- sum$TNC.mean
stdev<- sum$TNC.sd
windows(30,20);palette(rainbow(6))
b<- barplot(sum$TNC.mean, axes=F, axisnames=F, ylim=c(0,25), col=sum$run, main="TNC (%)")
axis(2, at= seq(0,25,by=5))
axis(1, labels=lablist, at=b, las=2, cex.axis=0.8)
box()
segments(b, means - stdev, b, means + stdev, lwd=2)
segments(b - 0.1, means - stdev, b + 0.1, means - stdev, lwd=2)
segments(b - 0.1, means + stdev, b + 0.1, means + stdev, lwd=2)

#Sugars
sugarefs <- subset(refs, Class == "su")
sum<- summaryBy(corrglumgg~Code+run, data=sugarefs, FUN=c(mean,sd)); sum$Code <- as.factor(sum$Code); sum$run<- as.factor(sum$run)
sum <- sum[order(sum$run),]
sd(sugarefs$corrglumgg)
#lablist<-sum$Code
lablist<- c(rep("LH",10),rep("JED",10),rep("RO",10),rep("AV",22))
means <- sum$corrglumgg.mean
stdev<- sum$corrglumgg.sd
windows(30,20);palette(rainbow(6))
b<- barplot(means, axes=F, axisnames=F, ylim=c(0,130), col=sum$run, main="Sugars (mg/g)")
axis(2, at= seq(0,130,by=10))
axis(1, labels=lablist, at=b, las=2, cex.axis=0.8)
box()
segments(b, means - stdev, b, means + stdev, lwd=2)
segments(b - 0.1, means - stdev, b + 0.1, means - stdev, lwd=2)
segments(b - 0.1, means + stdev, b + 0.1, means + stdev, lwd=2)

#Starch
starchrefs <- subset(refs, Class == "st")
sum<- summaryBy(corrglumgg~Code+run, data=starchrefs, FUN=c(mean,sd)); sum$Code <- as.factor(sum$Code); sum$run<- as.factor(sum$run)
sum <- sum[order(sum$run),]
mean(starchrefs$corrglumgg)
#lablist<-sum$Code
lablist<- c(rep("LH",10),rep("JED",10),rep("RO",10),rep("AV",22))
means <- sum$corrglumgg.mean
stdev<- sum$corrglumgg.sd
windows(30,20);palette(rainbow(6))
b<- barplot(means, axes=F, axisnames=F, ylim=c(0,130), col=sum$run, main="Starch (mg/g)")
axis(2, at= seq(0,130,by=10))
axis(1, labels=lablist, at=b, las=2, cex.axis=0.8)
box()
segments(b, means - stdev, b, means + stdev, lwd=2)
segments(b - 0.1, means - stdev, b + 0.1, means - stdev, lwd=2)
segments(b - 0.1, means + stdev, b + 0.1, means + stdev, lwd=2)
#-------------------------
#- nested random effects model (JED 25-5-2015)
library(nlme)
refs$runfactor <- as.factor(refs$run)
refs$Code <- as.factor(refs$Code)

fm.su <- lme(corrglumgg~1,random=~1|runfactor/Code/rep, data=subset(refs,Class=="su"))

fm.st <- lme(corrglumgg~1,random=~1|runfactor/Code/rep, data=subset(refs,Class=="st"))

fm.su
fm.st

#Su: 
Formula: ~1 | runfactor
StdDev:    5.024716
Formula: ~1 | Code %in% runfactor
StdDev:    3.869472
Formula: ~1 | rep %in% Code %in% runfactor
StdDev:    7.703713 2.899509

#St
Formula: ~1 | runfactor
StdDev:    6.391506
Formula: ~1 | Code %in% runfactor
StdDev:     6.805545
Formula: ~1 | rep %in% Code %in% runfactor
StdDev:    5.105788 2.523027
#-------------------------

#Have alook at the other samples
samples<-rbind(r1csv[!(r1csv$ID %in% c(1,13,25)),],r2csv[!(r2csv$ID %in% c(26,37,48)),],r3csv[!(r3csv$ID %in% c(49,73,97)),],r4csv[!(r4csv$ID %in% c(62,74,86)),])

mass <- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/TNCMass.csv")
mass <- rename(mass,c("ID"="Code"));mass <- rename(mass,c("Order"="ID"));
mass$ch <- as.factor(mass$ch); mass$leaf..branch<-as.factor(mass$leaf..branch); mass$ID<-as.factor(mass$ID)
dat<- merge(mass, samples, by="ID")
dat1 <- dat[,-10];dat1<-rename(dat1, c("Code.x"="Code","leaf..branch"="leaf"))
write.csv(dat1,"W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/WTC3_TNC_20150603_L1.csv")

dat2<- read.csv("W:/WorkingData/GHS39/GLAHD/Share/Data/TNC/WTC3_TNC_20150603_L1.csv")
dat2$date<- as.Date(as.character(dat2$date, format="%yyyy/%mm/%dd"))
dat2$ch<-as.factor(dat2$ch)

summar<-summaryBy(TNC~date+ch, data=dat2, FUN=c(mean,sd))
summar<-summaryBy(corrglumgg~date+ch, data=subset(dat2, Class=="st"), FUN=c(mean,sd))

ggplot(summar, aes(x=date, y=corrglumgg.mean , group = ch)) +
  geom_line() +
  geom_point() +
  coord_cartesian(ylim = c(0, 250))

# Set overall shapes and line type
ggplot(summar, aes(x=date, y=TNC.mean, group = ch)) +
  geom_line(linetype="dashed",  # Dashed line
            size = 1.5) +       # Thicker line
  geom_point(shape = 0,         # Hollow squares
             size = 4)          # Large points

# Condition shapes and line type on variable cond
ggplot(summar, aes(x=date, y=TNC.mean, group = ch)) +
  geom_line(aes(linetype=ch), # Line type depends on cond
            size = 1.5) +       # Thicker line
  geom_point(aes(shape=ch),   # Shape depends on cond
             size = 4)          # Large points

su2<-subset(refs2, Class=="st")
sd(su2$corrglumgg)               
