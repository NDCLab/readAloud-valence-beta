# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2023-06-30

# INPUTS
# data/df: behavioral data, for each participant on each passage, with relevant participant information and trial-level stimulus information

# OUTPUTS
# TBD

# NOTES TO DO
# drop 150086 as only completed 12 of 20 passages and low accuracy

### SECTION 1: SETTING UP
library(dplyr)
library(lme4)
library(lmerTest)
library(interactions)

#visualization tools
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(colorspace)
library(colorblindr)

#set up date for output file naming
today <- Sys.Date()
today <- format(today, "%Y%m%d")

#set up directories for input/output data
data <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/readAloudBetaData_20230630.csv'
out_path <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/'

#read in data
df <- read.csv(data)

#organize participant demographic variables
df$sex <- as.factor(df$sex)
df$pronouns <- as.factor(df$pronouns)
df$ethnic <- as.factor(df$ethnic)
df$socclass <- as.factor(df$socclass)

#extract demo stats
summary(df$age)
sd(df$age)
summary(df$sex)/18
summary(df$sex)/18 / (nrow(df)/18)
summary(df$pronouns)/18
summary(df$pronouns)/18 / (nrow(df)/18)
summary(df$ethnic)/18
summary(df$ethnic)/18 / (nrow(df)/18)
summary(df$socclass)/18
summary(df$socclass)/18 / (nrow(df)/18)

#remove participants who were not engaged in the task
#TBD, ex. 150222

#remove participants whose challenge question accuracy was below 50% (chance = 25%)
dfTrim <- df
dfTrim <- dfTrim %>%
  group_by(id) %>%
  mutate(challengeAvgSub = mean(challengeACC)) %>%
  ungroup

dfTrim <- subset(dfTrim, challengeAvgSub>0.5)
length(unique(df$id)) - length(unique(dfTrim$id)) #number of participants removed

#calculate average accuracy
mean(dfTrim$challengeAvgSub)
sd(dfTrim$challengeAvgSub)


### SECTION 2: PASSAGE-LEVEL TRIMMING
passage_no_before_trimming <- nrow(dfTrim)

#remove passages with high omissions (participant did not complete reading)
##vegas 150013

passage_no_after_trim1 <- nrow(dfTrim)
passage_no_before_trimming - passage_no_after_trim1 #number of passages trimmed
(passage_no_before_trimming - passage_no_after_trim1) / passage_no_before_trimming #percentage of passages trimmed


### SECTION 3: ORGANIZE DATA FOR MODELING
errorDat <- dfTrim

#modify contrasts for categorical predictors
contrasts(errorDat$sex) <- contr.sum(2) #male: -1, female: +1

#center continuous predictors
errorDat$age_gmc <- errorDat$age - mean(errorDat$age)
errorDat$bfne_gmc <- errorDat$bfne - mean(errorDat$bfne)
errorDat$phq8_gmc <- errorDat$phq8 - mean(errorDat$phq8)
errorDat$scaaredTotal_gmc <- errorDat$scaaredTotal - mean(errorDat$scaaredTotal)
errorDat$scaaredGA_gmc <- errorDat$scaaredGA - mean(errorDat$scaaredGA)
errorDat$scaaredSoc_gmc <- errorDat$scaaredSoc - mean(errorDat$scaaredSoc)
errorDat$sps_gmc <- errorDat$sps - mean(errorDat$sps)
errorDat$lenSyll_gmc <- errorDat$lenSyll - mean(errorDat$lenSyll)
errorDat$lenWord_gmc <- errorDat$lenWord - mean(errorDat$lenWord)
errorDat$avgSyllPerWord_gmc <- errorDat$avgSyllPerWord - mean(errorDat$avgSyllPerWord)

#extract demo stats
errorDatStats <- subset(errorDat, !duplicated(errorDat$id))
summary(errorDatStats$age)
sd(errorDatStats$age)
summary(errorDatStats$sex)
summary(errorDatStats$sex) / length(unique(errorDatStats$id))
summary(errorDatStats$pronouns)
summary(errorDatStats$pronouns) / length(unique(errorDatStats$id))
summary(errorDatStats$ethnic)
summary(errorDatStats$ethnic) / length(unique(errorDatStats$id))
summary(errorDatStats$socclass)
summary(errorDatStats$socclass) / length(unique(errorDatStats$id))


### SECTION 4: MODEL RESULTS
#misprod_rate x bfne
model1 <- lmerTest::lmer(misprod_rate ~ bfne_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model1)

#misprod_rate x scaaredSoc
model2 <- lmerTest::lmer(misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model2)

#misprod_rate x sps
model3 <- lmerTest::lmer(misprod_rate ~ sps_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model3)

#hesitation_rate x bfne
model4 <- lmerTest::lmer(hesitation_rate ~ bfne_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model4)

#hesitation_rate x scaaredSoc
model5 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model5)

#hesitation_rate x sps
model6 <- lmerTest::lmer(hesitation_rate ~ sps_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model6)

#words_with_misprod_rate x bfne
model7 <- lmerTest::lmer(words_with_misprod_rate ~ bfne_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model7)

#words_with_misprod_rate x scaaredSoc
model8 <- lmerTest::lmer(words_with_misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model8)

#words_with_misprod_rate x sps
model9 <- lmerTest::lmer(words_with_misprod_rate ~ sps_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model9)

#words_with_hes_rate x bfne
model10 <- lmerTest::lmer(words_with_hes_rate ~ bfne_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model10)

#words_with_hes_rate x scaaredSoc
model11 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model11)

#words_with_hes_rate x sps
model12 <- lmerTest::lmer(words_with_hes_rate ~ sps_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(model12)