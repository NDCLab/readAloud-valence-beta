# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2024-08-16

# INPUTS
# data/df: behavioral data, for each participant on each passage, with relevant participant information and trial-level stimulus information

# OUTPUTS
# models

# NOTES TO DO
# gmc

# Data dict

# errorDatMisprodHes:
#
#   our errorDat dataframe, just without the misprod-sequencing columns (which
#   we'll add in piecemeal by different names later)

# First, look at a given misproduction and check for nearby hesitations
#
# hes_position:
#
#   for long-form dataframes counting misproductions, this indicates whether the
#   relevant count is the number of hesitations before (0) or after (1) those
#   misproductions being counted in that row
#
#
# misprod_tally:
#
#   conversely, in long-form dataframes counting misproductions, this column
#   actually tracks how many misproductions there are in that reading
#   (participant x passage) that have a hesitation in the relevant relative
#   position


# justMisprodWithHesBefore:
#
#   this is the dataframe with every (participant x passage) reading, counting
#   the number of misproductions with a nearby preceding hesitation
#
#   i.e., for each reading, it counts the number of times (misprod_tally) that a
#   hesitation comes before a misproduction -- so for every entry, hes_position = 0
#
#
# justMisprodWithHesAfter
#
#   similarly, this is the dataframe with every (participant x passage) reading,
#   counting the number of misproductions with a nearby following hesitation
#
#   i.e., for each reading, it counts the number of times (misprod_tally) that a
#   hesitation comes after a misproduction -- so for every entry, hes_position = 1
#
#
# errorDatLongMisprodWithRelHes:
#
#   this is the long-form dataframe, with two rows per reading (participant x
#   passage): one for each position for a relative hesitation. i.e. this stacks
#   the two dataframes that respectively have (1) every passage, with a count of
#   misproductions for hes_position = 0, and (2) every passage, with a count of
#   misproductions for hes_position = 1



# Then, look at a given hesitation and check for nearby misproductions

# misprod_position:
#
#   for long-form dataframes counting hesitations, this indicates whether the
#   relevant count is the number of misproductions before (0) or after (1) those
#   hesitations being counted in that row
#
#
# hes_tally:
#
#   conversely, in long-form dataframes counting hesitations, this column
#   actually tracks how many hesitations there are in that reading (participant
#   x passage) that have a misproduction in the relevant relative position


# justHesWithMisprodBefore:
#
#   this is the dataframe with every (participant x passage) reading, counting
#   the number of hesitations with a nearby preceding misproduction
#
#   i.e., for each reading, it counts the number of times (hes_tally) that a
#   misproduction comes before a hesitation -- so for every entry,
#   misprod_position = 0
#
#
# justHesWithMisprodAfter
#
#   similarly, this is the dataframe with every (participant x passage) reading,
#   counting the number of hesitations with a nearby following misproduction
#
#   i.e., for each reading, it counts the number of times (hes_tally) that a
#   misproduction comes after a hesitation -- so for every entry,
#   misprod_position = 1
#
#
# errorDatLongHesWithRelMisprod:
#
#   this is the long-form dataframe, with two rows per reading (participant x
#   passage): one for each position for a relative misproduction, i.e. this
#   stacks the two dataframes that respectively have (1) every passage, with a
#   count of hesitations for misprod_position = 0, and (2) every passage, with a
#   count of hesitations for misprod_position = 1


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
data <- '/home/luc/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData_20230825.csv'
to_omit <- '/home/luc/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/input/passages-to-omit_20230810.csv'

#read in data
df <- read.csv(data, row.names = NULL) # output of prep script
passage_omissions_df <- read.csv(to_omit, row.names = NULL) # hand-crafted list of participant x passage entries to exclude, based on coder comments

#organize participant demographic variables
df$sex <- as.factor(df$sex)
df$pronouns <- as.factor(df$pronouns)
df$ethnic <- as.factor(df$ethnic)
df$socclass <- as.factor(df$socclass)

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
mean(dfTrim$challengeAvgSub) # bad?: weights participants by the number of passages they read
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

passage_no_after_trim2 <- nrow(dfTrim)
passage_no_after_trim1 - passage_no_after_trim2 #number of passages trimmed
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_after_trim1 #percentage of passages trimmed of last bunch
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_before_trimming #percentage of passages trimmed of whole





### SECTION 3: ORGANIZE DATA FOR MODELING
errorDat <- dfTrim


#modify contrasts for categorical predictors
contrasts(errorDat$sex) <- contr.sum(2) #male: -1, female: +1
# now verify:
contrasts(errorDat$sex)
# see below for how I've handled challengeACC- a special case as it is used both
# as a predictor and as an outcome


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

errorDat$timePerSyllable_gmc <- errorDat$timePerSyllable - mean(errorDat$timePerSyllable)
errorDat$timePerWord_gmc <- errorDat$timePerWord - mean(errorDat$timePerWord)

# likewise, for error rates

# First we'll encapsulate centering so we don't write each column name out
# manually three times. This way a typo in a column name can't screw it up. (It
# throws an error instead, rather than allowing in the incorrect data.)
gmc <- function(vec) { vec - mean(vec) }

add_gmc <- function(df, col) { # makes a new column: centered version of variable
  vec <- pull(df, {{col}}) # what's the data we're centering?

  df %>%
    mutate("{{col}}_gmc" := gmc(vec)) # ex. challengeAvgSub->challengeAvgSub_gmc
}

errorDat <-
  errorDat %>%
  add_gmc(words_with_misprod_rate) %>%
  add_gmc(words_with_hes_rate) %>%
  add_gmc(challengeAvgSub)


# sanity check
for (col_index in which(stringr::str_detect(colnames(errorDat), ".*_gmc"))) {
  colname <- names(errorDat[col_index])
  col <- errorDat[[col_index]]
  print(colname)

  avg <- mean(col)
  print(paste('  mean:', avg,
              '  rounded mean:', round(avg, digits = 10)))
}


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


# Questionnaires stats (ls additions 9/6/24)
print("SCAARED Social")
summary(errorDatStats$scaaredSoc)
sd(errorDatStats$scaaredSoc)

print("BFNE")
summary(errorDatStats$bfne)
sd(errorDatStats$bfne)

print("SPS")
summary(errorDatStats$sps)
sd(errorDatStats$sps)


### SECTION 3.1: Correct data for contrasts and numerics according to whether
#                they are predictors or outcomes in the models to follow
errorDatPredictorsOutcomes <- errorDat # separate them

# we will use the following columns as predictors:
# scaaredSoc_gmc
# interaction: scaaredSoc_gmc and words_with_hes_rate_gmc
# interaction: scaaredSoc_gmc and words_with_hes_rate
# interaction: scaaredSoc_gmc and words_with_misprod_rate

# the following as outcomes:
# challengeACC
# words_with_misprod_rate_gmc

# and the following as both:
# words_with_hes_rate
# words_with_hes_rate_gmc
# words_with_misprod_rate

#  Note also ~ scaaredSoc_gmc + age_gmc : I am not sure what we call/categorize
#  age here as

# Among the above, the ONLY binary column is challengeACC. So we just make two
# versions. So

# make predictor version of challengeACC:
# s/b -1 +1 factors

# make outcome version of challengeACC:
# s/b numeric


# verify current status:
class(errorDatPredictorsOutcomes$challengeACC) # -> "integer"

# the following lines are to repair unexpected behavior in the model below
errorDatPredictorsOutcomes$challengeACC_predictor <- as.numeric(errorDat$challengeACC)
errorDatPredictorsOutcomes$challengeACC_outcome <- as.numeric(errorDat$challengeACC)

# to be safe:
errorDatPredictorsOutcomes$challengeACC_predictor <-
  replace(
    errorDatPredictorsOutcomes$challengeACC_predictor,
    which(errorDatPredictorsOutcomes$challengeACC_predictor == 0),
    -1
  )

# now make it a factor
errorDatPredictorsOutcomes$challengeACC_predictor <- as.factor(errorDatPredictorsOutcomes$challengeACC_predictor)

# current status
contrasts(errorDatPredictorsOutcomes$challengeACC_predictor) # -> -1 (incorrect): 0, 1 (correct): 1

# set it to -1 rather than 0
contrasts(errorDatPredictorsOutcomes$challengeACC_predictor) <- rev(contr.sum(2)) # fix
contrasts(errorDatPredictorsOutcomes$challengeACC_predictor) # -> -1 (incorrect): -1, 1 (correct): 1

# prevent accidental use:
errorDatPredictorsOutcomes$challengeACC <- NULL


### SECTION 3.2: preparing for misprod-hes sequential analyses

# ignore the misprod-hes columns for now
errorDatMisprodHes <- select(errorDat, !contains("_syllables"))

# First: look at a given misproduction and check for nearby hesitations
justMisprodWithHesBefore <- cbind(errorDatMisprodHes,
                                  hes_position = 0, # "before",
                                  misprod_tally = errorDat$misprod_with_hes_in_previous_syllables)

justMisprodWithHesAfter <- cbind(errorDatMisprodHes,
                                 hes_position = 1, # "after",
                                 misprod_tally = errorDat$misprod_with_hes_in_next_syllables)


# stack the ones before and the ones after as rows of a single df (my attempt at long form)
errorDatLongMisprodWithRelHes <- rbind(justMisprodWithHesBefore, justMisprodWithHesAfter)

# track the binary relative position as a factor
errorDatLongMisprodWithRelHes$hes_position <- as.factor(errorDatLongMisprodWithRelHes$hes_position)

# Then: look at a given hesitation and check for nearby misproductions
justHesWithMisprodBefore <- cbind(errorDatMisprodHes,
                                  misprod_position = 0, # "before",
                                  hes_tally = errorDat$hes_with_misprod_in_previous_syllables)

justHesWithMisprodAfter <- cbind(errorDatMisprodHes,
                                 misprod_position = 1, # "after",
                                 hes_tally = errorDat$hes_with_misprod_in_next_syllables)

# stack the ones before and the ones after as rows of a single df (my attempt at long form)
errorDatLongHesWithRelMisprod <- rbind(justHesWithMisprodBefore, justHesWithMisprodAfter)

# track the binary relative position as a factor
errorDatLongHesWithRelMisprod$misprod_position <- as.factor(errorDatLongHesWithRelMisprod$misprod_position)


### SECTION 4: MODEL RESULTS
# for every model involving comprehension accuracy, rather than errorDat we use
# errorDatPredictorsOutcomes, which differentiates how comprehension accuracy is
# represented as a predictor versus as an outcome

#misprod_rate x bfne
# model1 <- lmerTest::lmer(misprod_rate ~ bfne_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model1)

#misprod_rate x scaaredSoc
# model2 <- lmerTest::lmer(misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model2)

#misprod_rate x sps
# model3 <- lmerTest::lmer(misprod_rate ~ sps_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model3)

#hesitation_rate x bfne
# model4 <- lmerTest::lmer(hesitation_rate ~ bfne_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model4)

#hesitation_rate x scaaredSoc
# model5 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model5)

#hesitation_rate x sps
# model6 <- lmerTest::lmer(hesitation_rate ~ sps_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model6)

#words_with_misprod_rate x bfne
# model7 <- lmerTest::lmer(words_with_misprod_rate ~ bfne_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model7)

#words_with_misprod_rate x scaaredSoc
# model8 <- lmerTest::lmer(words_with_misprod_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model8)

# fix: gmc
model8_center <- lmerTest::lmer(words_with_misprod_rate_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model8_center)

#words_with_misprod_rate x sps
# model9 <- lmerTest::lmer(words_with_misprod_rate ~ sps_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model9)

#words_with_hes_rate x bfne
# model10 <- lmerTest::lmer(words_with_hes_rate ~ bfne_gmc + (1|id) + (1|passage),
#                           data=errorDat, REML=TRUE)
# summary(model10)

# ! words_with_hes_rate x scaaredSoc
# model11 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                           data=errorDat, REML=TRUE)
# summary(model11)

# fix: gmc
model11_center <- lmerTest::lmer(words_with_hes_rate_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(model11_center)

#words_with_hes_rate x sps
# model12 <- lmerTest::lmer(words_with_hes_rate ~ sps_gmc + (1|id) + (1|passage),
#                           data=errorDat, REML=TRUE)
# summary(model12)


#### supplemental analyses

# glmer(accuracy ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, family="binomial")
# "f_" : follow-up

# Accuracy/comprehension as explained by social anxiety: scaaredSoc

# outcome is binary 0/1 numeric
f_model1 <- glmer(challengeACC_outcome ~ scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model1)

# fix: gmc (same column but with -1 (incorrect) and +1 (correct))
# errorDat$challengeACC <- replace(errorDat$challengeACC, which(errorDat$challengeACC == 0), -1)
# confirm:
# unique(errorDatPredictorsOutcomes$challengeACC_outcome) # -> (incorrect:) 0, (correct:) 1
# f_model1_center <- glmer(challengeACC ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                          data=errorDatPredictorsOutcomes, family = "binomial")
# summary(f_model1_center)


# Accuracy/comprehension as explained by social anxiety: bfne
# f_model2 <- glmer(challengeACC ~ bfne_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model2)


# Accuracy/comprehension as explained by social anxiety: sps
# f_model3 <- glmer(challengeACC ~ sps_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model3)


# Accuracy/comprehension as explained by disfluencies: hesitations per syllable
# f_model4 <- glmer(challengeACC ~ hesitation_rate + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model4)

# Accuracy/comprehension as explained by disfluencies: hesitations per word
# f_model5 <- glmer(challengeACC ~ words_with_hes_rate + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model5)

# fix: gmc

# outcome is binary 0/1 numeric
f_model5_center <- glmer(challengeACC_outcome ~ words_with_hes_rate_gmc + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model5_center)



# Accuracy/comprehension as explained by errors: misproductions per syllable
# f_model6 <- glmer(challengeACC ~ misprod_rate + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model6)

# Accuracy/comprehension as explained by errors: misproductions per word
f_model7 <- glmer(challengeACC_outcome ~ words_with_misprod_rate + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model7)




# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with scaared
# f_model8 <- glmer(challengeACC ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model8)

# Accuracy/comprehension as explained by disfluencies: hesitations per word with scaared
# f_model9 <- glmer(challengeACC_outcome ~ words_with_hes_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDatPredictorsOutcomes, family = "binomial")
# summary(f_model9)

# fix: gmc
f_model9_center <- glmer(challengeACC_outcome ~ words_with_hes_rate_gmc * scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model9_center)


# Accuracy/comprehension as explained by errors: misproductions per syllable with scaared
# f_model10 <- glmer(challengeACC ~ misprod_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model10)

# Accuracy/comprehension as explained by errors: misproductions per word with scaared
f_model11 <- glmer(challengeACC_outcome ~ words_with_misprod_rate * scaaredSoc_gmc + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model11)



# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with bfne
# f_model12 <- glmer(challengeACC ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model12)

# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per word with bfne
# f_model13 <- glmer(challengeACC ~ words_with_hes_rate * bfne_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model13)


# Accuracy/comprehension as explained by errors *and* SA: misproductions per syllable with bfne
# f_model14 <- glmer(challengeACC ~ misprod_rate * bfne_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model14)

# Accuracy/comprehension as explained by errors *and* SA: misproductions per word with bfne
# f_model15 <- glmer(challengeACC ~ words_with_misprod_rate * bfne_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model15)



# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with sps
# f_model16 <- glmer(challengeACC ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model16)

# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per word with sps
# f_model17 <- glmer(challengeACC ~ words_with_hes_rate * sps_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model17)


# Accuracy/comprehension as explained by errors *and* SA: misproductions per syllable with sps
# f_model18 <- glmer(challengeACC ~ misprod_rate * sps_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model18)

# Accuracy/comprehension as explained by errors *and* SA: misproductions per word with sps
# f_model19 <- glmer(challengeACC ~ words_with_misprod_rate * sps_gmc + (1|id) + (1|passage),
#                    data=errorDat, family = "binomial")
# summary(f_model19)



# Now, misproduction-hesitation relationships

# Errors as explained by disfluency: rate of misproduced syllables from rate of hesitated syllables
# f_model20 <- lmerTest::lmer(misprod_rate ~ hesitation_rate + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model20) # ***

# Errors as explained by disfluency: rate of misproduced words from rate of hesitated words
# f_model21 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model21) # ***

# fix: gmc
f_model21_center <- lmerTest::lmer(words_with_misprod_rate_gmc ~ words_with_hes_rate_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model21_center) # ***

# Errors as explained by disfluency: rate of misproduced words from rate of hesitated syllables
# f_model22 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model22) # ***



# Now, misproduction-hesitation interactions with social anxiety

# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and scaared
# f_model23 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model23)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and scaared
# f_model24 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model24)

# fix: gmc
f_model24_center <- lmerTest::lmer(words_with_misprod_rate_gmc ~ words_with_hes_rate_gmc * scaaredSoc_gmc + (1|id) + (1|passage),
                                   data=errorDat, REML=TRUE)
summary(f_model24_center)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and scaared
# f_model25 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model25)


# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and bfne
# f_model26 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model26)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and bfne
# f_model27 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * bfne_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model27)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and bfne
# f_model28 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * bfne_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model28)


# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and sps
# f_model29 <- lmerTest::lmer(misprod_rate ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model29)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and sps
# f_model30 <- lmerTest::lmer(words_with_misprod_rate ~ words_with_hes_rate * sps_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model30)


# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and sps
# f_model31 <- lmerTest::lmer(words_with_misprod_rate ~ hesitation_rate * sps_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model31)



# Now: see if reading speed plays into it

# Does scaaredSoc predict reading speed?
# syllable level
# rs_model_1 <- lmerTest::lmer(timePerSyllable_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_1)

# word level
# rs_model_2 <- lmerTest::lmer(timePerWord_gmc ~ scaaredSoc_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_2)


# rs_model_1_bfne <- lmerTest::lmer(timePerSyllable_gmc ~ bfne_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_1_bfne)

# word level
# rs_model_2_bfne <- lmerTest::lmer(timePerWord_gmc ~ bfne_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_2_bfne)

# rs_model_1_sps <- lmerTest::lmer(timePerSyllable_gmc ~ sps_gmc + (1|id) + (1|passage),
#                                   data=errorDat, REML=TRUE)
# summary(rs_model_1_sps)

# word level
# rs_model_2_bfne <- lmerTest::lmer(timePerWord_gmc ~ sps_gmc + (1|id) + (1|passage),
#                                   data=errorDat, REML=TRUE)
# summary(rs_model_2_sps)


# What happens when we control for age?
#hesitation_rate x scaaredSoc
# age_model1 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(age_model1)

#words_with_hes_rate x scaaredSoc
age_model2 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(age_model2)


# And now ->> check our work
# Does our hesitation ~ scaaredSoc finding hold with reading speed controlled for?
# syllable level
# rs_model_3 <- lmerTest::lmer(hesitation_rate ~ scaaredSoc_gmc + timePerSyllable_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_3)

# word level
# rs_model_4 <- lmerTest::lmer(words_with_hes_rate ~ scaaredSoc_gmc + timePerWord_gmc + (1|id) + (1|passage),
#                              data=errorDat, REML=TRUE)
# summary(rs_model_4)


# misprod-hes ordering

# Is the number of hesitations adjacent to misproductions in a particular
# reading predicted by the

# Does the position of misproductions relative to hesitations


# we have a number of occurrences of a misproduction in a particular position
# relative to a passage's hesitations. does knowing the position (before/after)
# predict the number of these sequences we have?

# does misproduction location relative to a hesitation predict how many
# instances we get in a particular reading?

# hes_with_rel_misprod_model_1 <- lmerTest::lmer(hes_tally ~ misprod_position + (1|id) + (1|passage),
#                                                data=errorDatLongHesWithRelMisprod, REML=TRUE)
# summary(hes_with_rel_misprod_model_1)
#
# misprod_with_rel_hes_model_1 <- lmerTest::lmer(misprod_tally ~ hes_position + (1|id) + (1|passage),
#                                                data=errorDatLongMisprodWithRelHes, REML=TRUE)
# summary(misprod_with_rel_hes_model_1)

## does it interact with SA?
# hes_with_rel_misprod_model_3 <- lmerTest::lmer(hes_tally ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|passage),
#                                                data=errorDatLongHesWithRelMisprod, REML=TRUE)
# # summary(hes_with_rel_misprod_model_3)
#
# misprod_with_rel_hes_model_4 <- lmerTest::lmer(misprod_tally ~ hes_position * scaaredSoc_gmc + (1|id) + (1|passage),
#                                                data=errorDatLongMisprodWithRelHes, REML=TRUE)
# summary(misprod_with_rel_hes_model_4)



# Word frequency analysis
# Does a passage's average word frequency predict participants' hesitation rate or misproduction rate?
# wordfreq_model_1 <- lmerTest::lmer(hesitation_rate ~ avgWordFreq + (1|id) + (1|passage),
#                                    data=errorDat, REML=TRUE)
# summary(wordfreq_model_1)
# wordfreq_model_2 <- lmerTest::lmer(misprod_rate ~ avgWordFreq + (1|id) + (1|passage),
#                                    data=errorDat, REML=TRUE)
# summary(wordfreq_model_2)

# Do social anxiety and frequency interact to predict hesitation rate or misproduction rate?
# wordfreq_model_3 <- lmerTest::lmer(hesitation_rate ~ avgWordFreq * scaaredSoc_gmc + (1|id) + (1|passage),
#                                    data=errorDat, REML=TRUE)
# summary(wordfreq_model_3)
#
# wordfreq_model_4 <- lmerTest::lmer(misprod_rate ~ avgWordFreq * scaaredSoc_gmc + (1|id) + (1|passage),
#                                    data=errorDat, REML=TRUE)
# summary(wordfreq_model_4)
