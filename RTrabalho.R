#Libarys

#install.packages("lme4", type = "source")

library(readr)

library(lme4)
library(lmerTest)
library(glmmTMB)

library(ggplot2)

library(fitdistrplus)
library(car)
library(emmeans)

#install.packages("Matrix")

#oo <- options(repos = "https://cran.r-project.org/")
#install.packages("Matrix")
#install.packages("lme4")
#options(oo)

#detach("package:lmerTest", unload = TRUE)
#detach("package:lme4", unload = TRUE)
#detach("package:Matrix", unload = TRUE)


#import the used excel tables for Fecundity assay

fecdata = read.table("Total-Fec-Mean.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",", na.strings = "NA")
fecdata

#summary of the data imported
summary(fecdata)
str(fecdata)

#change the order of the factors for better vizualization on statistical analysis
fecdata$Selection2<-factor(fecdata$Selection, levels=c("PT", "Cg-PT"))
fecdata

fecdata$EnvDensity2<-factor(fecdata$EnvDensity, levels=c("L", "H"))
fecdata



#model for our first assay
str(fecdata)
lm_fec = lmer(MeanOfEggs_F ~ Selection2 * EnvDensity2 + (1|AP) + (1|AP:Selection2) + (1|AP:EnvDensity2) + (1|AP:Selection2:EnvDensity2), data = fecdata)

summary(lm_fec)
Anova(lm_fec)

#models for mating duration (normal linear regression (used)) 

duration_matings = lmer(MATING_BEGINNING~Selection2 + (1|AP) + (1|Cage) + (1|AP:Selection2) , data = mattimedata)

summary(duration_matings)

Anova(duration_matings)

#main plots for our assay
#plot(fecdata$Selection, fecdata$MeanOfEggs_F, xlab = "Selection", ylab = "Nr of eggs per female", main = "Nr of eggs vs selection and Density")
#abline(coef(lm_fec), col = "green")

#boxplot(fecdata$MeanOfEggs_F ~ fecdata$Selection*fecdata$EnvDensity, xlab = "Selection vs Density", ylab = "Nr of eggs per female")


#Import our tables for Matings assay

mattimedata = read.table("Mating_times.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",")
mattimedata

#summary of the data
str(mattimedata)
summary(mattimedata)

matdata = read.table("Matings.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",")
matdata

#summary of the data
str(matdata)
summary(matdata)


#change the order o the factors for better vizualization of the statistical analysis
matdata$Selection2<-factor(matdata$Selection, levels=c("PT", "Cg-PT"))
matdata
mattimedata$Selection2<-factor(mattimedata$Selection, levels=c("PT", "Cg-PT"))
matdata


#histogramas

hist(mattimedata$MATING_DURATION) #normal
hist(matdata$Matings) #not normal
hist(matdata$Courts) #normal


hist(mattimedata$MATING_BEGINNING)

#distribuição de selection vs matings

boxplot(matdata$Courts ~ matdata$Selection) #ver com media


#distribuição de selection vs duração

boxplot(mattimedata$MATING_BEGINNING ~ mattimedata$Selection)

descdist(matdata$Matings, boot = 100, discrete = TRUE)
descdist(matdata$Courts, boot = 100, discrete = TRUE)


#models (matings)

#poisson
glm_matings = glmmTMB(Matings ~ Selection2 + (1|AP)+ (1|AP:Selection2), data = matdata, family = "poisson")

#negatic binominal
glm_matings_nb <- glmmTMB(Matings ~ Selection2 + (1|AP)+ (1|AP:Selection2), data = matdata,family=nbinom2, ziformula = ~1)


summary(glm_matings_nb)
summary(glm_matings)

anova(glm_matings, glm_matings_nb)

AIC(glm_matings,glm_matings_nb, k=2) # poisson is better

summary(glm_matings)

Anova(glm_matings)



#normal linear model Duration
lme_matings = lmer(MATING_DURATION ~ Selection2 + (1|AP)+ (1|AP:Selection), data = mattimedata)
summary(lme_matings)

anova(lme_matings)

###########normal linear model Beginning

lme_matings_b = lmer(MATING_BEGINNING ~ Selection2 + (1|AP)+ (1|AP:Selection), data = mattimedata)
summary(lme_matings_b)

anova(lme_matings_b)


#### modelos (Courts)

#poisson
glm_courts = glmmTMB(Courts ~ Selection + (1|AP)+ (1|AP:Selection), data = matdata, family = "poisson")

#negative binominal
glm_courts_nb <- glmmTMB(Courts ~ Selection + (1|AP)+ (1|AP:Selection), data = matdata, family=nbinom1, ziformula = ~1)

summary(glm_courts_nb)
summary(glm_courts)

anova(glm_courts, glm_courts_nb)

AIC(glm_courts,glm_courts_nb, k=2)

Anova(glm_courts)

#### Time

#distribuição de selection vs matings

boxplot(mattimedata$MATING_BEGINNING ~ mattimedata$Selection)

#distribuição

descdist(mattimedata$MATING_BEGINNING, boot = 100, discrete = TRUE)


#organizing factors
mattimedata$Selection2<-factor(mattimedata$Selection, levels=c("PT", "Cg-PT"))
mattimedata

#models for mating duration (normal linear regression (used)) 

duration_matings = lmer(MATING_BEGINNING~Selection2 + (1|AP) + (1|Cage) + (1|AP:Selection2) , data = mattimedata)

summary(duration_matings)

anova(duration_matings)




#histogramas, dispersão de dados 


#F4

ggplot(data = fecdata, aes(x = MeanOfEggsF_4, y = after_stat(density), fill = EnvDensity )) +
  geom_histogram() +
  facet_wrap(~ Selection*AP) +
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female in day 4") +
  ggtitle("Dispersal of data in Day 4") 
  

#F5

ggplot(data = fecdata, aes(x = MeanOfEggsF_5, y = after_stat(density), fill = EnvDensity )) +
  geom_histogram() +
  facet_wrap(~ Selection*AP) +
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female in day 5") +
  ggtitle("Dispersal of data in Day 5") 
  #geom_density(aes(group = Selection, color = Selection), linewidth = 1) 


#F6

ggplot(data = fecdata, aes(x = MeanOfEggsF_6, y = after_stat(density), fill = EnvDensity)) +
  geom_histogram() +
  facet_wrap(~ Selection*AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female in day 6") +
  ggtitle("Dispersal of data in Day 6") 
#  geom_density(color = "green", linewidth = 1)


#F9

  ggplot(data = fecdata, aes(x = MeanOfEggsF_9, y = after_stat(density), fill = EnvDensity)) +
  geom_histogram() +
  facet_wrap(~ Selection*AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female in day 9") +
  ggtitle("Dispersal of data in Day 9") 


#F_Total

ggplot(data = fecdata, aes(x = MeanOfEggs_F, y = after_stat(density), fill = EnvDensity)) +
  geom_histogram() +
  facet_wrap(~ Selection*AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female") +
  ggtitle("Dispersal of data of eggs per female") 


#boxplots para a fecundity por dias 

#F4

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggsF_4, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Data dispersion in day 4")

#F5

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggsF_5, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Data dispersion in day 5")

#F6

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggsF_6, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Data dispersion in day 6")

#F9

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggsF_9, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Data dispersion in day 9")


#boxplots para a fecundity media 

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggs_F, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection)+
  #facet_wrap(~ Rack)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Total of days")


ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggs_F, fill = Selection)) +
  geom_boxplot() +
  facet_wrap(~ Selection) +
  facet_wrap(~ AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Data dispersion")

#Matings

#histogramas, dispersão dos dados

#Mating Count 

ggplot(data = matdata, aes(x = Matings, y = after_stat(density), fill = Selection)) +
  geom_histogram() +
  facet_wrap(~ AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Matings in all populations") +
  ggtitle("Dispersal of data in Day 4")
#  geom_density(color = "green", linewidth = 1)


#Court Count 

ggplot(data = matdata, aes(x = Courts, y = after_stat(density), fill = Selection)) +
  geom_histogram() +
  facet_wrap(~ AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Density")+
  xlab("Number of Eggs per female in day 4") +
  ggtitle("Dispersal of data in Day 4")
#  geom_density(color = "green", linewidth = 1)


#boxplots para dados de mating duration


ggplot(mattimedata, aes(x = Selection, y = MATING_DURATION, fill = Selection)) +
  geom_boxplot() +
 # facet_wrap(~ AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Duraction of a mating (seconds)")+
  xlab("Selection Environment") +
  ggtitle("Total of days")



#####
#adictional test of corting data 

courtdata = read.table("Courtstime.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",")
courtdata
str(courtdata)
summary(courtdata)


ggplot(courtdata, aes(x = Selection, y = COURTSHIP_BEGINNING, fill = AP)) +
  geom_boxplot() +
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Duraction of a mating (seconds)")+
  xlab("Selection Environment") +
  ggtitle("Data dispersion")


#########################

#fecundity assay

ggplot(fecdata, aes(x = EnvDensity, y = MeanOfEggs_F, fill = Selection)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = Selection), position = position_dodge(width = 0.75)) +
  #facet_wrap(~ Rack)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs per female")+
  xlab("Environment") +
  ggtitle("Total of days")


#MATING_BEGINNING

ggplot(mattimedata, aes(x = Selection, y = MATING_BEGINNING, fill = Selection)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  # facet_wrap(~ AP)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Beginning of a mating (seconds)")+
  xlab("Selection Environment") +
  ggtitle("Total of days")


#how to citate and pick references form libraries
#packageVersion('fitdistrplus')

#citation("glmmTMB")

#citation()

#RStudio.Version()

#R.version

