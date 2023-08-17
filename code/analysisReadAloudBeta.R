# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2023-08-16

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

# ```
# Warning in install.packages :
# package ‘colorblindr’ is not available for this version of R
#
# A version of this package for your version of R might be available elsewhere,
# see the ideas at
# https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages
# ```

#set up date for output file naming
today <- Sys.Date()
today <- format(today, "%Y%m%d")

#set up directories for input/output data
# data <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/readAloudBetaData_20230630.csv'
# data <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData_20230810.csv'
# data <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData_20230815.csv'
data <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData_20230816.csv'
to_omit <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/input/passages-to-omit_20230810.csv'
# out_path <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/'
out_path <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/derivatives/'

#read in data
df <- read.csv(data, row.names = NULL) # output of prep script
passage_omissions_df <- read.csv(to_omit, row.names = NULL) # hand-crafted list of participant x passage entries to exclude, based on coder comments

#organize participant demographic variables
df$sex <- as.factor(df$sex)
df$pronouns <- as.factor(df$pronouns)
df$ethnic <- as.factor(df$ethnic)
df$socclass <- as.factor(df$socclass)


# FIXME: stats on df and dfTrim for reading speed are NA

#extract demo stats

# all these values are just in case they're useful - not needed per se for later
# steps of the logic in this script
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
dfTrim <- df

# removal based on coder notes of audible distractions, others in the room, etc.:
# 150015
# 150208
# 150245 had many passages that were entirely or near-entirely inaudible; the
# rest were dropped too under the assumption that the audible ones too would be
# too faint to identify errors in
dfTrim <- subset(dfTrim, !(id %in% c(150015, 150208, 150245)))

#remove participants whose challenge question accuracy was below 50% (chance = 25%)
dfTrim <- dfTrim %>%
  group_by(id) %>%
  mutate(challengeAvgSub = mean(challengeACC)) %>%
  ungroup

dfTrim <- subset(dfTrim, challengeAvgSub>0.5)
length(unique(df$id)) - length(unique(dfTrim$id)) #number of participants removed due to distraction or low accuracy

#calculate average accuracy
mean(dfTrim$challengeAvgSub)
sd(dfTrim$challengeAvgSub)

# calculate average speed
mean(dfTrim$timePerSyllable, na.rm=TRUE)
sd(dfTrim$timePerSyllable, na.rm=TRUE)
mean(dfTrim$timePerWord, na.rm=TRUE)
sd(dfTrim$timePerWord, na.rm=TRUE)


### SECTION 2: PASSAGE-LEVEL TRIMMING
passage_no_before_trimming <- nrow(dfTrim)

#remove passages with high omissions (participant did not complete reading) or other problems (someone else is in the room, etc.)
# e.g. vegas 150013
dfTrim <- anti_join(dfTrim,
                    passage_omissions_df,
                    by = join_by(id == participant, passage == passage))


passage_no_after_trim1 <- nrow(dfTrim)
passage_no_before_trimming - passage_no_after_trim1 #number of passages trimmed
(passage_no_before_trimming - passage_no_after_trim1) / passage_no_before_trimming #percentage of passages trimmed


# band-aid fix: remove passages without reading speed data so that we can run
# our analyses on them nonetheless
# todo

# these are the only four passages without reading time data...
# and incidentally? well, see their comments here...
c(150013, "vegas")      # N.B.: 161 omitted syllables of 318 total in passage
c(150022, "depression") # N.B.: 160 omitted syllables of 362 total in passage
c(150083, "caramel")    # N.B.: only one of four passages to have >= 5% of syllables omitted
c(150083, "cars")       # N.B.: only one of four passages to have >= 5% of syllables omitted


dfTrim <- filter(dfTrim, !is.na(timePerSyllable))
# itself, but without ones for which we have no reading data
# this ends up only dropping 0083, caramel - the other three already end up
# getting dropped based on other criteria



"TODO"

passage_no_after_trim2 <- nrow(dfTrim)
passage_no_after_trim1 - passage_no_after_trim2 #number of passages trimmed
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_after_trim1 #percentage of passages trimmed of last bunch
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_before_trimming #percentage of passages trimmed of whole





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

# LS additions 8/11/23
# errorDat$errors <- errorDat$errors - mean(errorDat$errors)
# errorDat$correction <- errorDat$corrections - mean(errorDat$corrections)
errorDat$error_rate <- errorDat$errors / errorDat$lenSyll
errorDat$correction_rate <- errorDat$corrections / errorDat$lenSyll

errorDat$timePerSyllable_gmc <- errorDat$timePerSyllable - mean(errorDat$timePerSyllable)
errorDat$timePerWord_gmc <- errorDat$timePerWord - mean(errorDat$timePerWord)

# todo center avg reading speed, probably


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

# Reading speed stats (ls additions 8/16/23)
summary(errorDatStats$timePerSyllable)
summary(errorDatStats$timePerWord)

summary(errorDatStats$timePerSyllable_gmc)
summary(errorDatStats$timePerWord_gmc)

### SECTION 4: MODEL RESULTS
#misprod_rate x bfne
model1 <- lmerTest::lmer(misprod_rate ~ bfne_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model1)

#misprod_rate x scaaredSoc
model2 <- lmerTest::lmer(misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model2)

#misprod_rate x sps
model3 <- lmerTest::lmer(misprod_rate ~ sps_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model3)

#hesitation_rate x bfne
model4 <- lmerTest::lmer(hesitation_rate ~ bfne_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model4)

#hesitation_rate x scaaredSoc
model5 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model5)

#hesitation_rate x sps
model6 <- lmerTest::lmer(hesitation_rate ~ sps_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model6)

#words_with_misprod_rate x bfne
model7 <- lmerTest::lmer(words_with_misprod_rate ~ bfne_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model7)

#words_with_misprod_rate x scaaredSoc
model8 <- lmerTest::lmer(words_with_misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model8)

#words_with_misprod_rate x sps
model9 <- lmerTest::lmer(words_with_misprod_rate ~ sps_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model9)

#words_with_hes_rate x bfne
model10 <- lmerTest::lmer(words_with_hes_rate ~ bfne_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(model10)

#words_with_hes_rate x scaaredSoc
model11 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(model11)

#words_with_hes_rate x sps
model12 <- lmerTest::lmer(words_with_hes_rate ~ sps_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(model12)


#### supplemental analyses
# see notes






# glmer(accuracy ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, family="binomial")
# "f_" : follow-up

# Accuracy/comprehension as explained by social anxiety: scaaredSoc
f_model1 <- glmer(challengeACC ~ scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model1)

# Accuracy/comprehension as explained by social anxiety: bfne
f_model2 <- glmer(challengeACC ~ bfne_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model2)


# Accuracy/comprehension as explained by social anxiety: sps
f_model3 <- glmer(challengeACC ~ sps_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model3)


# Accuracy/comprehension as explained by disfluencies: hesitations per syllable
f_model4 <- glmer(challengeACC ~ hesitation_rate + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model4)

# Accuracy/comprehension as explained by disfluencies: hesitations per word
f_model5 <- glmer(challengeACC ~ words_with_hes_rate + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model5)


# Accuracy/comprehension as explained by errors: misproductions per syllable
f_model6 <- glmer(challengeACC ~ misprod_rate + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model6)

# Accuracy/comprehension as explained by errors: misproductions per word
f_model7 <- glmer(challengeACC ~ words_with_misprod_rate + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model7)




# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with scaared
f_model8 <- glmer(challengeACC ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model8)

# Accuracy/comprehension as explained by disfluencies: hesitations per word with scaared
f_model9 <- glmer(challengeACC ~ words_with_hes_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model9)


# Accuracy/comprehension as explained by errors: misproductions per syllable with scaared
f_model10 <- glmer(challengeACC ~ misprod_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model10)

# Accuracy/comprehension as explained by errors: misproductions per word with scaared
f_model11 <- glmer(challengeACC ~ words_with_misprod_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model11)



# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with bfne
f_model12 <- glmer(challengeACC ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model12)

# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per word with bfne
f_model13 <- glmer(challengeACC ~ words_with_hes_rate * bfne_gmc + (1|id) + (1|passage),
                  data=errorDat, family = "binomial")
summary(f_model13)


# Accuracy/comprehension as explained by errors *and* SA: misproductions per syllable with bfne
f_model14 <- glmer(challengeACC ~ misprod_rate * bfne_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model14)

# Accuracy/comprehension as explained by errors *and* SA: misproductions per word with bfne
f_model15 <- glmer(challengeACC ~ words_with_misprod_rate * bfne_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model15)



# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with sps
f_model16 <- glmer(challengeACC ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model16)

# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per word with sps
f_model17 <- glmer(challengeACC ~ words_with_hes_rate * sps_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model17)


# Accuracy/comprehension as explained by errors *and* SA: misproductions per syllable with sps
f_model18 <- glmer(challengeACC ~ misprod_rate * sps_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model18)

# Accuracy/comprehension as explained by errors *and* SA: misproductions per word with sps
f_model19 <- glmer(challengeACC ~ words_with_misprod_rate * sps_gmc + (1|id) + (1|passage),
                   data=errorDat, family = "binomial")
summary(f_model19)



# Now, misproduction-hesitation relationships

# Errors as explained by disfluency: rate of misproduced syllables from rate of hesitated syllables
f_model20 <- lmerTest::lmer(misprod_rate ~ hesitation_rate + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model20) # ***

# Errors as explained by disfluency: rate of misproduced words from rate of hesitated words
f_model21 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model21) # ***


# Errors as explained by disfluency: rate of misproduced words from rate of hesitated syllables
f_model22 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model22) # ***



# Now, misproduction-hesitation interactions with social anxiety

# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and scaared
f_model23 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model23)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and scaared
f_model24 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model24)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and scaared
f_model25 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model25)


# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and bfne
f_model26 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model26)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and bfne
f_model27 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * bfne_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model27)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and bfne
f_model28 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model28)


# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and sps
f_model29 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model29)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and sps
f_model30 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * sps_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model30)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and sps
f_model31 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model31)



# Now: see if reading speed plays into it

# Does scaaredSoc predict reading speed?
# syllable level
rs_model_1 <- lmerTest::lmer(timePerSyllable_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_1)

# word level
rs_model_2 <- lmerTest::lmer(timePerWord_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_2)


rs_model_1_bfne <- lmerTest::lmer(timePerSyllable_gmc ~ bfne_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_1_bfne)

# word level
rs_model_2_bfne <- lmerTest::lmer(timePerWord_gmc ~ bfne_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_2_bfne)

rs_model_1_sps <- lmerTest::lmer(timePerSyllable_gmc ~ sps_gmc + (1|id) + (1|passage),
                                  data=errorDat, REML=TRUE)
summary(rs_model_1_sps)

# word level
rs_model_2_bfne <- lmerTest::lmer(timePerWord_gmc ~ sps_gmc + (1|id) + (1|passage),
                                  data=errorDat, REML=TRUE)
summary(rs_model_2_sps)



# Does scaaredSoc predict reading speed?
# syllable level
rs_model_3 <- lmerTest::lmer(timePerSyllable ~ scaaredSoc_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_3)

# word level
rs_model_4 <- lmerTest::lmer(timePerWord ~ scaaredSoc_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_4)






# What happens when we control for age?
#hesitation_rate x scaaredSoc
age_model1 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(age_model1)

#words_with_hes_rate x scaaredSoc
age_model2 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(age_model2)


# And now ->> check out work
# Does our hesitation ~ scaaredSoc finding hold with reading speed controlled for?
# syllable level
rs_model_3 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + timePerSyllable_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_3)

# word level
rs_model_4 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + timePerWord_gmc + (1|id) + (1|passage),
                             data=errorDat, REML=TRUE)
summary(rs_model_4)


# todo: check models 1 and 2 when hesitation rate is held still





# LS ideas:
# error rate ~ SA
# correction rate ~ error rate * SA

# comprehension ~ correction rate
# comprehension ~ SA * correction rate
# comprehension ~ error rate
# comprehension ~ SA * error rate
# comprehension ~ error rate * correction rate
# comprehension ~ error rate * correction rate * SA

# comprehension ~ hes rate * error rate
# comprehension ~ hes rate * correction rate

# comprehension ~ hes rate * error rate * SA
# comprehension ~ hes rate * correction rate * SA
# comprehension ~ hes rate * error rate * correction rate * SA

# error rate      ~ hes rate * SA
# correction rate ~ hes rate * SA, control for error rate

# corr_ : corrections rate
# err_  : errors rate


# error rate as explained by social anxiety
fls_model_err_scaared <- lmer(error_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                                   data=errorDat, REML=TRUE)
summary(fls_model_err_scaared)

# correction rate as explained by social anxiety
fls_model_corr_scaared <- lmer(correction_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
                              data=errorDat, REML=TRUE)
summary(fls_model_corr_scaared)

# correction rate as explained by social anxiety and errors
fls_model_corr_scaared_err <- lmer(correction_rate ~ scaaredSoc_gmc * error_rate + (1|id) + (1|passage),
                                   data=errorDat, REML=TRUE)
summary(fls_model_corr_scaared_err)

# comprehension ~ correction rate
fls_model_comp_corr <- glmer(challengeACC ~ correction_rate + (1|id) + (1|passage),
                             data=errorDat, family = "binomial")
summary(fls_model_comp_corr)

# comprehension ~ SA * correction rate
fls_model_comp_scaared_corr <- glmer(challengeACC ~ scaaredSoc_gmc * correction_rate + (1|id) + (1|passage),
                                     data=errorDat, family = "binomial")
summary(fls_model_comp_scaared_corr) # *

# comprehension ~ error rate
fls_model_comp_err <- glmer(challengeACC ~ error_rate + (1|id) + (1|passage),
                             data=errorDat, family = "binomial")
summary(fls_model_comp_err)

# comprehension ~ SA * error rate
fls_model_comp_scaared_err <- glmer(challengeACC ~ scaaredSoc_gmc * error_rate + (1|id) + (1|passage),
                                    data=errorDat, family = "binomial")
summary(fls_model_comp_scaared_err)

# comprehension ~ error rate * correction rate
fls_model_comp_err_corr <- glmer(challengeACC ~ error_rate * correction_rate + (1|id) + (1|passage),
                                 data=errorDat, family = "binomial")
summary(fls_model_comp_err_corr)

# comprehension ~ hes rate * error rate
fls_model_comp_hes_err <- glmer(challengeACC ~ hesitation_rate * error_rate + (1|id) + (1|passage),
                                data=errorDat, family = "binomial")
summary(fls_model_comp_hes_err)

# comprehension ~ hes rate * correction rate
fls_model_comp_hes_corr <- glmer(challengeACC ~ hesitation_rate * correction_rate + (1|id) + (1|passage),
                                 data=errorDat, family = "binomial")
summary(fls_model_comp_hes_corr) # ***


# Check the main one of interest, comprehension as explained by SA and
# corrections (see fls_model_comp_scaared_corr above), with all three tests
f_model_comp_scaared_corr <- glmer(challengeACC ~ scaaredSoc_gmc * correction_rate + (1|id) + (1|passage),
                                   data=errorDat, family = "binomial")
summary(f_model_comp_scaared_corr) # *

f_model_comp_bfne_corr <- glmer(challengeACC ~ bfne_gmc * correction_rate + (1|id) + (1|passage),
                                data=errorDat, family = "binomial")
summary(f_model_comp_bfne_corr) # N.S.

f_model_comp_sps_corr <- glmer(challengeACC ~ sps_gmc * correction_rate + (1|id) + (1|passage),
                               data=errorDat, family = "binomial")
summary(f_model_comp_sps_corr) # N.S., but p for interaction is at 0.0698


