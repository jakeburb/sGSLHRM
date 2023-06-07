#####################################################
################Factors Influencing Recruitment of Atlantic Herring##################

rm(list = ls())

library(reshape2)
library(tidyverse)
library(ggfortify)
library(MASS)
library(mgcv)
library(performance)
library(pscl)
library(gslea)
library(dplyr)
library(MuMIn)
library(gulf)
library(ggplot2)
library(scales)
library(gridExtra)

cat("\f")
clg()

current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))

fp <- getwd()

#We will load all data then conduct the modelling for spring and fall spwaning stocks seperately
#####################################SPRING###########################
#################
# Load SSB and Recruitment Data
#################
##SSB
mctmp <- read.table("mcoutSSBa.dat", sep="")

dim(mctmp)
dim(mctmp)[1]
ssb.5 <- apply(mctmp[26:dim(mctmp)[1],], 2, median)/1000
rm(mctmp)
##ssb.5 is median SSB estimates
#Need to append this to data frame below with recruitment
#first remove last 2 years of ssb data
ssbN <- head(ssb.5, -2)     
ssbN
ssbAb<-ssbN*1000

##Load in recruitment values (aged 2 lagged 2 years) based on 2022 Stock Assessment
springrecN<-read.csv("SpringR.csv")
plot(springrecN$year,springrecN$Rec, type="l")
springrecN[which.max(springrecN$Rec),]
springrecN[which.min(springrecN$Rec),]
#create dataframe with ssb and recruits
springrecr<-data.frame(springrecN$year,springrecN$Rec, springrecN$UL, springrecN$LL, ssbN)
#compute recruitment rate by year, times ssb to get it in kg
recrate<-springrecr$springrecN.Rec/(springrecr$ssbN*1000)
recrate
SRRUL<-springrecr$springrecN.UL/(springrecr$ssbN*1000)
SRRLL<-springrecr$springrecN.LL/(springrecr$ssbN*1000)
#create new data frame with recruitment rate and year
recrateT<-data.frame(springrecN$year,recrate,SRRUL,SRRLL, ssbAb)
colnames(recrateT)<-c("year","recrate", "RRUL", "RRLL", "SSB")
recrateT[which.max(recrateT$recrate),]
recrateT[which.min(recrateT$recrate),]
springrecN<-subset(springrecN, year>1977) 

########Load age structure data, proportion of age group compared to total age 4+
Sagestruct<-read.csv("SpringAgeProp.csv")
Sagestruct<-subset(Sagestruct, year>1977)

cor(springrecN$Rec,Sagestruct$X4)
cor(springrecN$Rec,Sagestruct$X5)
cor(springrecN$Rec,Sagestruct$X6)
cor(springrecN$Rec,Sagestruct$X7)
cor(springrecN$Rec,Sagestruct$X8)
cor(springrecN$Rec,Sagestruct$X9)
cor(springrecN$Rec,Sagestruct$X10)
cor(springrecN$Rec,Sagestruct$X11)

##Look at ages8 to 11+
Sage8to11<- rowSums(Sagestruct[,c(6,7,8,9)])
Sage8to11T<-data.frame(Sagestruct$year,Sage8to11)
colnames(Sage8to11T)<-c("year","Page8to11")
Sage8to11T<-subset(Sage8to11T, year>1977) 
excel(Sage8to11T)
mean(Sage8to11T$Page8to11)
Sage8to11T[which.max(Sage8to11T$Page8to11),]
Sage8to11T[which.min(Sage8to11T$Page8to11),]

cor(springrecN$Rec,Sage8to11T$Page8to11)
cor.test(springrecN$Rec,Sage8to11T$Page8to11)
plot(Sage8to11T$Page8to11, springrecN$Rec)
summary(lm(springrecN$Rec~Sage8to11T$Page8to11))

#Load weight at age data
SWatA<-read.csv("SpringWatA.csv")
SWatA[is.na(SWatA)] = 0
SWatA<-subset(SWatA, year>1977)

cor(springrecN$Rec,SWatA$X4)
cor(springrecN$Rec,SWatA$X5)
cor(springrecN$Rec,SWatA$X6)
cor(springrecN$Rec,SWatA$X7)
cor(springrecN$Rec,SWatA$X8)
cor(springrecN$Rec,SWatA$X9)
cor(springrecN$Rec,SWatA$X10)
cor(springrecN$Rec,SWatA$X11)

SWatASum<-rowSums(SWatA[,c(2,3,4,5,6,7,8,9)])
SWatASumT<-data.frame(Sagestruct$year,SWatASum)
colnames(SWatASumT)<-c("year","SWatAtot")
mean(SWatASumT$SWatAtot)
SWatASumT[which.max(SWatASumT$SWatAtot),]
SWatASumT[which.min(SWatASumT$SWatAtot),]


cor(springrecN$Rec, SWatASumT$SWatAtot)
cor.test(springrecN$Rec, SWatASumT$SWatAtot)
plot(SWatASumT$SWatAtot,log(springrecN$Rec))

summary(lm(log(springrecN$Rec)~SWatASumT$SWatAtot))


###Load abiotic and biotic variables available in gslea package
plank.var= vars.f(variable.type="planktonic")
formattable::formattable(plank.var)

physical.var= vars.f(variable.type="physical")
formattable::formattable(physical.var)

Variables<-EA.query.f(years=2001:2019, variables=c("amplitude","calanus.finmarchicus.annual","calanus.finmarchicus.early_summer",
                                                   "calanus.finmarchicus.fall", "calanus.hyperboreus.annual","calanus.hyperboreus.early_summer","calanus.hyperboreus.fall",
                                                   "chl0_100.annual", "chl0_100.early_summer", "chl0_100.fall", "chl0_100.late_summer","ci.civ.cfin.early_summer",
                                                   "ci.civ.cfin.fall", "civ.glac.fall","civ.hyp.early_summer","civ.hyp.fall","cold.annual", "cold.early_summer", "cold.fall","duration",
                                                   "dw2_t.annual","dw2_t.early_summer", "dw2_t.fall", "largecal.annual","largecal.early_summer", "largecal.fall", "magnitude",
                                                   "non.copepods.annual","non.copepods.early_summer","non.copepods.fall","pseudocalanus.annual", "pseudocalanus.early_summer",
                                                   "pseudocalanus.fall","smallcal.annual","smallcal.early_summer","smallcal.fall", "start","warm.annual", "warm.early_summer",
                                                   "warm.fall", "decrease.10", "decrease.12", "first.ice", "ice.duration", "ice.max", "last.ice", "sst", "sst.anomaly",
                                                   "sst.month10", "sst.month11", "sst.month5", "sst.month6", "sst.month7", "sst.month8", "sst.month9", "start.10",
                                                   "start.12","t.shallow"), EARs=c(3,5,6))

##Compute mean of each variable by year
meanvars = Variables %>% group_by(year, variable) %>% summarise(value = mean(value))
meanvars = reshape2::dcast(meanvars, year~variable)

#Merge variables from GLSEA package with recruitment data and other variables (age strcuture and weight at age)
vardat<- merge(meanvars, springrecN, by = "year")
vardatB<- merge(vardat, Sage8to11T, by = "year")
vardatC<- merge(vardatB, SWatASumT, by = "year")
vardatD<- merge(vardatC, recrateT, by = "year")
vardatD[is.na(vardatD)]<-0
##look at correlations between variables and Recruitment
cormatrix<-cor(vardatD) 

#Plot each variable that has high correlation or is hypothesized to be important (based on my conceptual diagram) with recruitment to briefly view the general relationship. 
plot(vardatD$ci.civ.cfin.early_summer,vardatD$Rec)
plot(vardatD$civ.hyp.early_summer,vardatD$Rec)
plot(vardatD$civ.glac.fall,vardatD$Rec)
plot(vardatD$cold.annual,vardatD$Rec)
plot(vardatD$dw2_t.annual,vardatD$Rec)
plot(vardatD$largecal.early_summer, vardatD$Rec)
plot(vardatD$sst.month6,vardatD$Rec)
plot(vardatD$sst.month8, vardatD$Rec)
plot(vardatD$start.10, vardatD$Rec)
plot(vardatD$warm.annual, vardatD$Rec)
##Create new dataframe with only variables that will be used in models
vardatES<-vardatD[ , c("year","ci.civ.cfin.early_summer", "civ.hyp.early_summer", "cold.annual","dw2_t.annual", "largecal.early_summer",
                          "sst.month6","sst.month8","start.10","warm.annual","SWatAtot","Rec","recrate", "Page8to11")]
vardatES
#Look breifly at correlation among variables 
cormatrixSN<-cor(vardatES)

plot(vardatES$SSB, vardatES$Rec)
##Compute mean and variance of recruitment
mean(vardatES$Rec)
var(vardatES$Rec)
#####################################
######Recruitment Modelling#########
###################################
##Recruitment as response
fullmodelSR<-glm.nb(Rec ~  ci.civ.cfin.early_summer + civ.hyp.early_summer + dw2_t.annual + 
                       largecal.early_summer +  sst.month6 + sst.month8 + start.10 + warm.annual + SWatAtot + Page8to11,    na.action = "na.fail", link = "log", data = vardatES)
summary(fullmodelSR)
modelcombs<-dredge(fullmodelSR)
100*with(summary(fullmodelSR), 1 - deviance/null.deviance)

##Best model
bestmodS<-glm.nb(Rec ~ sst.month8 + start.10 + warm.annual, link = "log", data = vardatES)
summary(bestmodS)
AICc(bestmodS)
100*with(summary(bestmodS), 1 - deviance/null.deviance)
check_model(bestmodS)
#Compute dispersion for best model
E4 <- resid(bestmodS, type = "pearson")
N  <- nrow(vardatES)
p  <- length(coef(bestmodS)) + 1  #the +1 is due to k
Dispersion <- sum(E4^2) / (N - p)
Dispersion

##First generate predictions#
p<- predict(bestmodS,  type = "link", se = T)
p
exp(p$fit)
##Compile data frame for plot##
spring_recruitsE<-data.frame(vardatES$year,vardatES$Rec, (exp(p$fit)), (exp(p$fit + 1.96 * p$se.fit)), (y= exp(p$fit - 1.96 * p$se.fit)))
colnames(spring_recruitsE)<-c("year","Rec","pred","ULpred","LLpred")
spring_recruitsE

####Now check full model but without biological variables in 2001-2019 as per reviewer comments
##Recruitment as response
fullmodelSRNB<-glm.nb(Rec ~   sst.month6 + sst.month8 + start.10 + SWatAtot + Page8to11,    na.action = "na.fail", link = "log", data = vardatES)
summary(fullmodelSRNB)
modelcombsNB<-dredge(fullmodelSRNB)
100*with(summary(fullmodelSRNB), 1 - deviance/null.deviance)
excel(modelcombsNB)

########Recruitment rate as response##############
fullmodelSRR<-glm(log(recrate) ~ ci.civ.cfin.early_summer + civ.hyp.early_summer + dw2_t.annual +
                    largecal.early_summer +  sst.month6 + sst.month8 + start.10 + warm.annual + SWatAtot + Page8to11,   na.action = "na.fail", data = vardatES)
summary(fullmodelSRR)
modelcombsRR<-dredge(fullmodelSRR)
100*with(summary(fullmodelSRR), 1 - deviance/null.deviance)
##Best model
bestmodRRS<-glm(log(recrate) ~ sst.month8 + dw2_t.annual, data = vardatES )
summary(bestmodRRS)
AICc(bestmodRRS)
100*with(summary(bestmodRRS), 1 - deviance/null.deviance)
check_model(bestmodRRS)


#Compute dispersion for best model
E4 <- resid(bestmodRRS, type = "pearson")
N  <- nrow(vardatES)
p  <- length(coef(bestmodRRS)) + 1  #the +1 is due to k
DispersionRR <- sum(E4^2) / (N - p)
DispersionRR



##First generate predictions#
pRR<- predict(bestmodRRS,  type = "link", se = T)
pRR
exp(pRR$fit)
##Compile data frame for plot##
spring_rrE<-data.frame(vardatES$year,vardatES$recrate, (exp(pRR$fit)), (exp(pRR$fit + 1.96 * pRR$se.fit)), (y= exp(pRR$fit - 1.96 * pRR$se.fit)))
colnames(spring_rrE)<-c("year","recrate","pred","ULpred","LLpred")
spring_rrE

####Now check full model but without biological variables in 2001-2019 as per reviewer comments
########Recruitment rate as response##############
fullmodelSRRNB<-glm(log(recrate) ~ sst.month6 + sst.month8 + start.10 + SWatAtot + Page8to11,   na.action = "na.fail", data = vardatES)
summary(fullmodelSRRNB)
modelcombsRRNB<-dredge(fullmodelSRRNB)
100*with(summary(fullmodelSRRNB), 1 - deviance/null.deviance)
excel(modelcombsRRNB)
#########################################################################
###########################FALL##########################################
#########################################################################
#################
# Load SSB and Recruits
#################
#SSB
mctmpF <- read.table("mcoutSSBaT.dat", sep="")

dim(mctmpF)
dim(mctmpF)[1]
ssb.5F <- apply(mctmpF[26:dim(mctmpF)[1],], 2, median)/1000
rm(mctmpF)
##ssb.5 is median SSB estimates
#Need to append this to data frame below with recruitment
#first remove last 2 years of ssb data
ssbNF <- head(ssb.5F, -2)     
ssbNF
ssbAbF<-ssbNF*1000 
##Load in recruitment values (aged 2 lagged 2 years) based on 2022 Stock Assessment
FallrecN<-read.csv("FallR.csv")
plot(FallrecN$year,FallrecN$Rec, type="l")
plot(ssbAbF, FallrecN$Rec)
FallrecN[which.max(FallrecN$Rec),]
FallrecN[which.min(FallrecN$Rec),]
#create dataframe with ssb and recruits
Fallrecr<-data.frame(FallrecN$year,FallrecN$Rec,FallrecN$UL,FallrecN$LL, ssbNF)
#compute recruitment rate by year, times ssb to get it in kg
recrateF<-Fallrecr$FallrecN.Rec/(Fallrecr$ssbN*1000)
recrateF
FRRUL<-Fallrecr$FallrecN.UL/(Fallrecr$ssbN*1000)
FRRLL<-Fallrecr$FallrecN.LL/(Fallrecr$ssbN*1000)
#create new data frame with recruitment rate and year
recrateFT<-data.frame(FallrecN$year,recrateF,FRRUL,FRRLL, ssbAbF)
colnames(recrateFT)<-c("year","recrate","RRUL", "RRLL", "ssb")
plot(recrateFT$year,recrateFT$recrate, type="l")
recrateFT[which.max(recrateFT$recrate),]
recrateFT[which.min(recrateFT$recrate),]

FallrecN<-subset(FallrecN, year>1977) 
#Load age structure data, proportion of age group compared to total age 4+
Fagestruct<-read.csv("FallAgeProp.csv")

cor(FallrecN$Rec,Fagestruct$X4)
cor(FallrecN$Rec,Fagestruct$X5)
cor(FallrecN$Rec,Fagestruct$X6)
cor(FallrecN$Rec,Fagestruct$X7)
cor(FallrecN$Rec,Fagestruct$X8)
cor(FallrecN$Rec,Fagestruct$X9)
cor(FallrecN$Rec,Fagestruct$X10)
cor(FallrecN$Rec,Fagestruct$X11)

Fage8to11<- rowSums(Fagestruct[,c(6,7,8,9)])
Fage8to11

Fage8to11T<-data.frame(Fagestruct$year,Fage8to11)
colnames(Fage8to11T)<-c("year","Page8to11")
Fage8to11T<-subset(Fage8to11T, year>1977)
excel(Fage8to11T)
mean(Fage8to11T$Page8to11)
Fage8to11T[which.max(Fage8to11T$Page8to11),]
Fage8to11T[which.min(Fage8to11T$Page8to11),]
Fage8to11T<-subset(Fage8to11T, year>1977) 

cor(FallrecN$Rec,Fage8to11T$Page8to11)
cor.test(FallrecN$Rec,Fage8to11T$Page8to11)
plot(Fage8to11T$Page8to11, FallrecN$Rec)
summary(lm(FallrecN$Rec~Fage8to11T$Page8to11))


##Load weight at age
FWatA<-read.csv("FallWatA.csv")
FWatA[is.na(FWatA)] = 0
FWatA<-subset(FWatA, year>1977)
  
cor(FallrecN$Rec,FWatA$X4)
cor(FallrecN$Rec,FWatA$X5)
cor(FallrecN$Rec,FWatA$X6)
cor(FallrecN$Rec,FWatA$X7)
cor(FallrecN$Rec,FWatA$X8)
cor(FallrecN$Rec,FWatA$X9)
cor(FallrecN$Rec,FWatA$X10)
cor(FallrecN$Rec,FWatA$X11)
Fagestruct<-subset(Fagestruct, year>1977)

FWatASum<-rowSums(FWatA[,c(2,3,4,5,6,7,8,9)])
FWatASumT<-data.frame(Fagestruct$year,FWatASum)
colnames(FWatASumT)<-c("year","FWatAtot")
FWatASumT<-subset(FWatASumT, year>1977)
mean(FWatASumT$FWatAtot)
FWatASumT[which.max(FWatASumT$FWatAtot),]
FWatASumT[which.min(FWatASumT$FWatAtot),]

cor(FallrecN$Rec, FWatASum)
plot(FWatASum,log(FallrecN$Rec))

summary(lm(log(FallrecN$Rec)~FWatASum))



###Now make dataframe with all variables pulled from GSLEA package and fall recruitment
vardatF<-merge(meanvars, FallrecN, by = "year")
vardatFB<- merge(vardatF, Fage8to11T, by = "year")
vardatFC<- merge(vardatFB, FWatASumT, by = "year")
vardatFD<- merge(vardatFC, recrateFT, by = "year")
vardatFD[is.na(vardatFD)]<-0

###Briefly look at correlation among variables
Fcormat<-cor(vardatFD)

vardatFE<-vardatFD[ , c("year","ci.civ.cfin.fall","civ.glac.fall", "civ.hyp.early_summer", "dw2_t.annual", "largecal.annual",
                       "non.copepods.annual","sst.month8","sst.month10","start.10","warm.annual","FWatAtot","Page8to11", "Rec","recrate")]
vardatFE

Fcordil<-cor(vardatFE)

#####################################
######Recruitment Modelling#########
###################################
###Recruitment models Fall###
fullmodelF<-glm.nb(Rec ~ ci.civ.cfin.fall + civ.glac.fall + dw2_t.annual + largecal.annual + 
                     sst.month8 + sst.month10 + start.10 + warm.annual + FWatAtot + Page8to11, data=vardatFE,na.action = "na.fail", link="log" )
summary(fullmodelF)
FmodelcombsR<-dredge(fullmodelF)
100*with(summary(fullmodelF), 1 - deviance/null.deviance)
excel(FmodelcombsR)

###Best model
FbestmodelR<-glm.nb(Rec ~ ci.civ.cfin.fall + dw2_t.annual + largecal.annual, link = "log", data = vardatFE)
summary(FbestmodelR)
AICc(FbestmodelR)
100*with(summary(FbestmodelR), 1 - deviance/null.deviance)
check_model(FbestmodelR)
#Compute dispersion for best model
E4 <- resid(FbestmodelR, type = "pearson")
N  <- nrow(vardatFS)
p  <- length(coef(FbestmodelR)) + 1  #the +1 is due to k
Dispersion <- sum(E4^2) / (N - p)
Dispersion

##First generate predictions#
pF<- predict(FbestmodelR,  type = "link", se = T)
pF
exp(pF$fit) 
##Compile data frame for ggplot plot##
fall_recruitsE<-data.frame(vardatFE$year,vardatFE$Rec, (exp(pF$fit)), (exp(pF$fit + 1.96 * pF$se.fit)), (y= exp(pF$fit - 1.96 * pF$se.fit)))
colnames(fall_recruitsE)<-c("year","Rec","pred","ULpred","LLpred")
fall_recruitsE

#2001-2019 recruits but with No Biological Variables
###Recruitment models Fall###
fullmodelFNB<-glm.nb(Rec ~ sst.month8 + sst.month10 + start.10 + FWatAtot + Page8to11, data=vardatFE,na.action = "na.fail", link="log" )
summary(fullmodelFNB)
FmodelcombsRNB<-dredge(fullmodelFNB)
100*with(summary(fullmodelFNB), 1 - deviance/null.deviance)
excel(FmodelcombsRNB)


####Now model with Recruitment Rate####################

FFullmodRR<-glm(log(recrate)~ci.civ.cfin.fall + civ.glac.fall + dw2_t.annual + largecal.annual + 
                  sst.month8 + sst.month10 + start.10 + warm.annual + FWatAtot + Page8to11, data=vardatFE, na.action = "na.fail")
summary(FFullmodRR)
FRRmodcombos<-dredge(FFullmodRR)
excel(FRRmodcombos)

###Best model
FbestmodRR<-glm(log(recrate)~ci.civ.cfin.fall + warm.annual, data=vardatFE)
summary(FbestmodRR)
AICc(FbestmodRR)
100*with(summary(FbestmodRR), 1 - deviance/null.deviance)
check_model(FbestmodRR)
#Compute dispersion for best model
E4 <- resid(FbestmodRR, type = "pearson")
N  <- nrow(vardatFE)
p  <- length(coef(FbestmodRR)) + 1  #the +1 is due to k
Dispersion <- sum(E4^2) / (N - p)
Dispersion


##First generate predictions#
pFRR<- predict(FbestmodRR,  type = "link", se = T)
pFRR
exp(pFRR$fit) 
##Compile data frame for plot##
fall_rrE<-data.frame(vardatFE$year,vardatFE$recrate, (exp(pFRR$fit)), (exp(pFRR$fit + 1.96 * pFRR$se.fit)), (y= exp(pFRR$fit - 1.96 * pFRR$se.fit)))
colnames(fall_rrE)<-c("year","recrate","pred","ULpred","LLpred")
fall_rrE

####Now check full model but without biological variables in 2001-2019 as per reviewer comments
########Recruitment rate as response##############
FFullmodRRNB<-glm(log(recrate)~ sst.month8 + sst.month10 + start.10 +  FWatAtot + Page8to11, data=vardatFE, na.action = "na.fail")
summary(FFullmodRRNB)
FRRmodcombosNB<-dredge(FFullmodRRNB)
excel(FRRmodcombosNB)


##############################################################################################
##########COMPILE RECRUITS AND RECRUITMENT RATE FROM ENTIRE TIMESERIES FOR FIGURE 1#####################
#############################################################################################
#Compile with environmental factors for consistency and to facilitate modelling over entire time series if eer needed (not done here)
##Start by loading relevant variables from across the time series from gslea package
##NOTE SST data doesnt start till 1982, therefore only take 1982-2019
VariablesT<-EA.query.f(years=1978:2019, variables=c( "decrease.10", "decrease.12", "first ice", "ice.duration", "ice.max", "last.ice", "sst", "sst.anomaly",
                                                   "sst.month10", "sst.month11", "sst.month5", "sst.month6", "sst.month7", "sst.month8", "sst.month9", "start.10",
                                                   "start.12","t.shallow"), EARs=c(3,5,6))

meanvarsT<- VariablesT %>% group_by(year, variable) %>% summarise(value = mean(value))
meanvarsT<- reshape2::dcast(meanvarsT, year~variable)
meanvarsT

#######################SPRING###############################
##Merge enviromental variables with recruits, recruitment rate, age strcuture and weight at age from above##

#Merge variables from GLSEA package with recruitment data and other variables (age strcuture and weight at age)
vardatST<- merge(meanvarsT, springrecN, by = "year")
vardatSTB<- merge(vardatST, Sage8to11T, by = "year")
vardatSTC<- merge(vardatSTB, SWatASumT, by = "year")
vardatSTD<- merge(vardatSTC, recrateT, by = "year")

cormatST<-cor(vardatSTD)
excel(cormatST)

########Recruits##########
##Compile data frame for ggplot plot##
spring_recruits<-data.frame(vardatSTD$year,vardatSTD$Rec,vardatSTD$UL,vardatSTD$LL)
colnames(spring_recruits)<-c("year","Rec","UL","LL")
spring_recruits

########Recruitment Rate###########
##Compile data frame for plot##
spring_rr<-data.frame(vardatSTD$year,vardatSTD$recrate,vardatSTD$RRUL,vardatSTD$RRLL)
colnames(spring_rr)<-c("year","recrate","RRUL","RRLL")
spring_rr


########################################################################
###############################FALL#####################################
#######################################################################
##Important/merge fall data
vardatFT<-merge(meanvarsT, FallrecN, by = "year")
vardatFBT<- merge(vardatFT, Fage8to11T, by = "year")
vardatFCT<- merge(vardatFBT, FWatASumT, by = "year")
vardatFDT<- merge(vardatFCT, recrateFT, by = "year")

cormatFT<-cor(vardatFDT)

##Compile data frame for plot##
fall_recruits<-data.frame(vardatFDT$year,vardatFDT$Rec,vardatFDT$UL, vardatFDT$LL)
colnames(fall_recruits)<-c("year","Rec","UL","LL")
fall_recruits

##Compile data frame for plot##
fall_rr<-data.frame(vardatFDT$year,vardatFDT$recrate, vardatFDT$RRUL, vardatFDT$RRLL)
colnames(fall_rr)<-c("year","recrate","RRUL", "RRLL")
fall_rr

#############Plotting Results################
#######FIGURE 2###################
years <- c(1978:2019)

output <- "Figure spring and fall recruit and recruitment rate_R1.png"

png(output, width = 3500, height = 3500, res = 300)

par(mfrow=c(2,2),mar=c(2,3.5,0.5,1),oma=c(4,4,0.5,1))

# spring recruit

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(1977, max(years)+1),
     ylim = c(0, 2000))
axis(1, at = seq(1975, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(1977, 2020, by = 5), label = F)
#axis(1, at = seq(1980, 2020, by = 2), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 2000, by = 50), tck = -0.008, label = F)
axis(2, at = seq(0, 2000, by = 250), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(spring_recruits$LL/1000,spring_recruits$UL[length(years):1]/1000), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,spring_recruits$Rec/1000,type='l',lwd=2)

plot_label <- "a)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

ylabel <- "Recruits (millions)"

mtext(ylabel, outer = F, side = 2, line = 4, cex = 2, font = 2)

# fall recruit

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(1977, max(years)+1),
     ylim = c(0, 3100))
axis(1, at = seq(1975, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(1975, 2020, by = 5), label = F)
#axis(1, at = seq(1980, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 3100, by = 50), tck = -0.008, label = F)
axis(2, at = seq(0, 3100, by = 500), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(fall_recruits$LL/1000,fall_recruits$UL[length(years):1]/1000), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,fall_recruits$Rec/1000,type='l',lwd=2)

plot_label <- "b)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

# spring rr

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(1977, max(years)+1),
     ylim = c(0, 35))
axis(1, at = seq(1975, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(1975, 2020, by = 5), label = F)
axis(1, at = seq(1975, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 35, by = 1), tck = -0.008, label = F)
axis(2, at = seq(0, 35, by = 5), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(spring_rr$RRLL,spring_rr$RRUL[length(years):1]), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,spring_rr$recrate,type='l',lwd=2)

plot_label <- "c)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

ylabel <- "Recruitment Rate (recruits/kg SSB)"
mtext(ylabel, outer = F, side = 2, line = 4, cex = 2, font = 2)

# fall recruitment rate

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(1977, max(years)+1),
     ylim = c(0, 25))
axis(1, at = seq(1975, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(1975, 2020, by = 5), label = F)
axis(1, at = seq(1975, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 25, by = 1), tck = -0.008, label = F)
axis(2, at = seq(0, 25, by = 2), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(fall_rr$RRLL,fall_rr$RRUL[length(years):1]), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,fall_rr$recrate,type='l',lwd=2)

plot_label <- "d)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

xlabel <-  "Year"
mtext(xlabel, outer = T, side = 1, line = 2, cex = 2, font = 2)

dev.off()



##########Figure 3#############
years <- c(2001:2019)

output <- "Figure spring and fall recruit and recruitment rate MODEL01-19_R1.png"

png(output, width = 3500, height = 3500, res = 300)

par(mfrow=c(2,2),mar=c(2,3.5,0.5,1),oma=c(4,4,0.5,1))

# spring recruit

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(2000, max(years)+1),
     ylim = c(0, 1050))
axis(1, at = seq(2000, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(2000, 2020, by = 5), label = F)
#axis(1, at = seq(1980, 2020, by = 2), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 2000, by = 50), tck = -0.008, label = F)
axis(2, at = seq(0, 2000, by = 250), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(spring_recruitsE$LLpred/1000,spring_recruitsE$ULpred[length(years):1]/1000), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,spring_recruitsE$pred/1000,type='l',lwd=2)
points(years,spring_recruitsE$Rec/1000,pch=19)

plot_label <- "a)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

ylabel <- "Recruits (millions)"

mtext(ylabel, outer = F, side = 2, line = 4, cex = 2, font = 2)

# fall recruit

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(2000, max(years)+1),
     ylim = c(0, 3000))
axis(1, at = seq(2000, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(2000, 2020, by = 5), label = F)
#axis(1, at = seq(1980, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 3100, by = 50), tck = -0.008, label = F)
axis(2, at = seq(0, 3100, by = 500), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(fall_recruitsE$LLpred/1000,fall_recruitsE$ULpred[length(years):1]/1000), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,fall_recruitsE$pred/1000,type='l',lwd=2)
points(years,fall_recruitsE$Rec/1000,pch=19)

plot_label <- "b)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

# spring rr

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(2000, max(years)+1),
     ylim = c(0, 22))
axis(1, at = seq(2000, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(2000, 2020, by = 5), label = F)
axis(1, at = seq(2000, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 35, by = 1), tck = -0.008, label = F)
axis(2, at = seq(0, 35, by = 5), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(spring_rrE$LLpred,spring_rrE$ULpred[length(years):1]), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,spring_rrE$pred,type='l',lwd=2)
points(years,spring_rrE$recrate,pch=19)

plot_label <- "c)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)
legend("topright", c("Observed", "Predicted"), lty = c(NA, 1), lwd = c(NA, 2), pch = c(19, NA), cex = 1.5)
box()

ylabel <- "Recruitment Rate (recruits/kg SSB)"
mtext(ylabel, outer = F, side = 2, line = 4, cex = 2, font = 2)

# fall recruitment rate

plot(NULL, xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", xlim = c(2000, max(years)+1),
     ylim = c(0, 10))
axis(1, at = seq(2000, max(years)+1, by = 1), tck = -0.008, label = F)
axis(1, at = seq(2000, 2020, by = 5), label = F)
axis(1, at = seq(2000, 2020, by = 5), tck = 0, cex.axis = 1.5)
axis(2, at = seq(0, 20, by = 1), tck = -0.008, label = F)
axis(2, at = seq(0, 20, by = 2), cex.axis = 1.5, las = 2)

polygon(c(years,years[length(years):1]),c(fall_rrE$LLpred,fall_rrE$ULpred[length(years):1]), col =  rgb(195, 195, 195, alpha = 125, max=255), border = NA)
lines(years,fall_rrE$pred,type='l',lwd=2)
points(years,fall_rrE$recrate,pch=19)

plot_label <- "d)"

legend("topleft", plot_label, bty = "n", cex = 1.5, text.font = 2)

box()

xlabel <-  "Year"
mtext(xlabel, outer = T, side = 1, line = 2, cex = 2, font = 2)

dev.off()