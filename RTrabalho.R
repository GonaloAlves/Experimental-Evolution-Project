library(readr)

library(lme4)
library(lmerTest)
library(glmmTMB)

library(ggplot2)

library(fitdistrplus)
library(car)
library(emmeans)


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


#Import our tables for Matings assay

mattimedata = read.table("Mating_times.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",")
mattimedata

#summary of the data
str(mattimedata)
summary(mattimedata)


#Import our tables for Matings assay

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





#########################plots#############

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
  ggtitle("")


