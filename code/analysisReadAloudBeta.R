# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2023-06-14

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
data <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/readAloudBetaData_20230614.csv'
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


### SECTION 2: INITIAL DATA TRIMMING
passage_no_before_trimming <- nrow(dfTrim)

#insert passage-level trimming here

passage_no_after_trim1 <- nrow(dfTrim)
passage_no_before_trimming - passage_no_after_trim1 #number of passages trimmed
(passage_no_before_trimming - passage_no_after_trim1) / passage_no_before_trimming #percentage of passages trimmed


### SECTION 3: TRANSITION DATA TO LONG FORMAT
errorDat <- data.frame(matrix(ncol=17, nrow=0))
colnames(errorDat) <- c("passage", "id",
                        "sex", "pronouns", "age", "ethnic", "socclass",
                        "bfne", "phq8", "scaaredTotal", "scaaredGA", "scaaredSoc", "sps",
                        "lenSyll", "lenWord", "avgSyllPerWord",
                        "errors")
for(j in 1:nrow(dfTrim)){
  passage <- dfTrim$passage[j]
  id <- dfTrim$id[j]
  sex <- as.character(dfTrim$sex[j])
  pronouns <- as.character(dfTrim$pronouns[j])
  age <- dfTrim$age[j]
  ethnic <- as.character(dfTrim$ethnic[j])
  socclass <- as.character(dfTrim$socclass[j])
  bfne <- dfTrim$bfne[j]
  phq8 <- dfTrim$phq8[j]
  scaaredTotal <- dfTrim$scaaredTotal[j]
  scaaredGA <- dfTrim$scaaredGA[j]
  scaaredSoc <- dfTrim$scaaredSoc[j]
  sps <- dfTrim$sps[j]
  lenSyll <- dfTrim$lenSyll[j]
  lenWord <- dfTrim$lenWord[j]
  avgSyllPerWord <- dfTrim$avgSyllPerWord[j]
  errors <- dfTrim$errors[j]
  errorDat[nrow(errorDat) + 1,] <-c(passage, id, sex, pronouns, age, ethnic, socclass, bfne, phq8, scaaredTotal, scaaredGA, scaaredSoc, sps, lenSyll, lenWord, avgSyllPerWord, errors)
}

#organize data types
errorDat$sex <- as.factor(errorDat$sex)
errorDat$pronouns <- as.factor(errorDat$pronouns)
errorDat$age <- as.numeric(errorDat$age)
errorDat$ethnic <- as.factor(errorDat$ethnic)
errorDat$socclass <- as.factor(errorDat$socclass)
errorDat$bfne <- as.numeric(errorDat$bfne)
errorDat$phq8 <- as.numeric(errorDat$phq8)
errorDat$scaaredTotal <- as.numeric(errorDat$scaaredTotal)
errorDat$scaaredGA <- as.numeric(errorDat$scaaredGA)
errorDat$scaaredSoc <- as.numeric(errorDat$scaaredSoc)
errorDat$sps <- as.numeric(errorDat$sps)
errorDat$lenSyll <- as.numeric(errorDat$lenSyll)
errorDat$lenWord <- as.numeric(errorDat$lenWord)
errorDat$avgSyllPerWord <- as.numeric(errorDat$avgSyllPerWord)
errorDat$errors <- as.numeric(errorDat$errors)

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
errorDat$errors_gmc <- errorDat$errors - mean(errorDat$errors)

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
modelErr <- lmerTest::lmer(errors ~ scaaredTotal_gmc * lenSyll_gmc + (1|id) + (1|passage), data=errorDat, REML=TRUE)
summary(modelErr)