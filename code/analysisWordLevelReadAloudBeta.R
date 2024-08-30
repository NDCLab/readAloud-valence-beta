# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2024-08-17

# INPUTS
# data/df: behavioral data, for each participant on each passage, with relevant participant information and trial-level stimulus information

# OUTPUTS
# models

# NOTES TO DO
# center

# Data dict

# errorDatMisprodHes:
#
#   our errorDat dataframe, just without the misprod-sequencing columns (which
#   we'll add in piecemeal by different names later)

# First, look at a given misproduction and check for nearby hesitations
#
# hes_position:
#
#   for long-form dataframes looking for misproductions, this indicates whether
#   the relevant binary value is tracking whether a hesitation comes before (0)
#   or after (1) the misproduction being counted in that row
#
#
# misprod_in_adjacent_window:
#
#   conversely, in long-form dataframes looking at misproductions, this column
#   actually tracks whether or not there is a misproduction in the relevant
#   relative position (to a hesitation), within the specified window


# justMisprodWithHesBefore:
#
#   this is the dataframe with every (participant x passage x word) read,
#   checking for each word whether there is a misproduction with a nearby
#   preceding hesitation
#
#   i.e., for each reading, it only looks at rows where a hesitation comes
#   before a misproduction -- so for every entry, hes_position = 0
#
#
# justMisprodWithHesAfter
#
#   similarly, this is the dataframe with every (participant x passage x word)
#   read, checking each time whether there is a misproduction with a nearby
#   following hesitation
#
#   i.e., for each reading, it only looks at rows where a hesitation comes after
#   a misproduction -- so for every entry, hes_position = 1
#
#
# errorDatLongMisprodWithRelHes:
#
#   this is the long-form dataframe, with two rows per word read (participant x
#   passage x word): one for each position for a relative hesitation, i.e. this
#   stacks the two dataframes that respectively have (1) every word, where
#   misprod_in_adjacent_window corresponds to hes_position = 0, and (2) every
#   word, where misprod_in_adjacent_window corresponds to hes_position = 1



# Then, look at a given hesitation and check for nearby misproductions

# misprod_position:
#
#   for long-form dataframes looking for hesitations, this indicates whether the
#   relevant binary value is tracking whether a misproduction comes before (0)
#   or after (1) the hesitation being counted in that row
#
#
# hes_in_adjacent_window:
#
#   conversely, in long-form dataframes looking at hesitations, this column
#   actually tracks whether or not there is a hesitation in the relevant
#   relative position (to a misproduction), within the specified window


# justHesWithMisprodBefore:
#
#   this is the dataframe with every (participant x passage x word) read,
#   checking for each word whether there is a hesitation with a nearby
#   preceding misproduction
#
#   i.e., for each reading, it only looks at rows where a misproduction comes
#   before a hesitation -- so for every entry, misprod_position = 0
#
#
# justHesWithMisprodAfter
#
#   similarly, this is the dataframe with every (participant x passage x word)
#   read, checking for each word whether there is a hesitation with a nearby
#   following misproduction
#
#   i.e., for each reading, it only looks at rows where a misproduction comes
#   after a hesitation -- so for every entry, misprod_position = 1
#
#
# errorDatLongHesWithRelMisprod:
#
#   this is the long-form dataframe, with two rows per word read (participant x
#   passage x word): one for each position for a relative misproduction, i.e.
#   this stacks the two dataframes that respectively have (1) every word, where
#   hes_in_adjacent_window corresponds to misprod_position = 0, and (2) every
#   word, where hes_in_adjacent_window corresponds to misprod_position = 1


### SECTION 1: SETTING UP
library(dplyr)
library(purrr)
library(lme4)
library(lmerTest)
library(interactions)

#visualization tools
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(colorspace)
library(effects)
# library(colorblindr)

#set up date for output file naming
today <- Sys.Date()
today <- format(today, "%Y%m%d")

#set up directories for input/output data
data <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData-wordLevel_20240130.csv'
to_omit <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/input/passages-to-omit_20230810.csv'
# out_path <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/'
out_path <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/'

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
# from prep script:
demo_cols <- c("age", "sex", "pronouns", "ethnic", "socclass")
sample_size <- select(df, "id") %>% unique %>% nrow


summary_unique <- function(df, key, column, f = summary) {
  # Consider every value of `column` but without duplicates within any unique
  # identifiers in `key`. Summarize (or `f`) what remains.
  # unique(select(df, column, key))[[column]] %>% f
  df %>%
    select({{column}}, {{key}}) %>%
    unique %>%
    pull({{column}}) %>%
    f
}

demo_cols %>% # as totals, for each of our demographic columns
  map(\(col) summary_unique(df, id, col)) # `map` not `for`: return, not print

demo_cols %>% # as a percent (but that doesn't mean anything for min/max/median)
  discard(~. == "age") %>% # (...so drop age, which `summary` formats that way)
  map(\(col) summary_unique(df, id, col) / sample_size)

summary_unique(df, id, "age", f = sd) # also do stdev for age


# now remove participants who were not engaged in the task

# removal based on coder notes of audible distractions, others in the room, etc.:
# 150015 and 150208
# 150245 had many passages that were entirely or near-entirely inaudible; the
# rest were dropped too under the assumption that the audible ones too would be
# too faint to identify errors in
dfTrim <- subset(df, !(id %in% c(150015, 150208, 150245)))

#remove participants whose challenge question accuracy was below 50% (chance = 25%)
dfTrim <- dfTrim %>%
  group_by(id) %>%
  mutate(challengeAvgSub = mean(challengeACC)) %>%
  ungroup %>%
  subset(challengeAvgSub>0.5)

length(unique(df$id)) - length(unique(dfTrim$id)) #number of participants removed due to distraction or low accuracy

#calculate average accuracy
summary_unique(dfTrim, "id", "challengeAvgSub")
summary_unique(dfTrim, "id", "challengeAvgSub", f = sd)


### SECTION 2: PASSAGE-LEVEL TRIMMING
passage_no_before_trimming <- summary_unique(dfTrim, "id", "passage", f = length)

#remove passages with high omissions (participant did not complete reading) or other problems (someone else is in the room, etc.)
# e.g. vegas 150013
dfTrim <- anti_join(dfTrim,
                    passage_omissions_df,
                    by = join_by(id == participant, passage == passage))


passage_no_after_trim1 <- summary_unique(dfTrim, "id", "passage", f = length)
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

passage_no_after_trim2 <- summary_unique(dfTrim, "id", "passage", f = length)
passage_no_after_trim1 - passage_no_after_trim2 #number of passages trimmed
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_after_trim1 #percentage of passages trimmed of last bunch
(passage_no_after_trim1 - passage_no_after_trim2) / passage_no_before_trimming #percentage of passages trimmed of whole

# report state of demographic data after all trimming is complete
demo_cols %>% # as totals, for each of our demographic columns
  map(\(col) summary_unique(dfTrim, "id", col))

analyzed_sample_size <- select(dfTrim, "id") %>% unique %>% nrow

demo_cols %>% # as a percent (but that doesn't mean anything for min/max/median)
  discard(~. == "age") %>% # (...so drop age, which `summary` formats that way)
  map(\(col) summary_unique(dfTrim, "id", col) / analyzed_sample_size)

summary_unique(dfTrim, "id", "age", f = sd) # also do stdev for age



### SECTION 3: ORGANIZE DATA FOR MODELING
errorDat <- dfTrim

#modify contrasts for categorical predictors
contrasts(errorDat$sex) <- contr.sum(2) #male: -1, female: +1

# changes here!! TODO

binary_to_plus_minus_contrast <- function(vec) {
  vec2 <- vec %>% replace(which(vec == 0), -1) %>% as.factor()

  contrasts(vec2) <- rev(contr.sum(2))

  return(vec2)
}

#errorDat <-
#  mutate(errorDat, across(challengeACC, where(is.logical), binary_to_plus_minus_factor))
errorDat <-
  errorDat %>%
  mutate(across(where(is.logical), binary_to_plus_minus_contrast)) %>%
  mutate(challengeACC = binary_to_plus_minus_contrast(challengeACC))


DEBUG <- FALSE # show examples proving this works
if (DEBUG) {
  sandboxDat <- errorDat

  sandboxDat %>% select(where(is.logical)) %>% colnames # see what's there
  sandboxDat %>% select(where(is.factor)) %>% colnames # see what's there

  sandboxDat <-
    sandboxDat %>%
    mutate(across(where(is.logical), binary_to_plus_minus_contrast)) %>%
    mutate(challengeACC = binary_to_plus_minus_contrast(challengeACC))


  sandboxDat %>% select(where(is.logical)) %>% colnames # see what's there
  sandboxDat %>% select(where(is.factor)) %>% colnames # see what's there
  # as.factor(sandboxDat$misprod) %>% contrasts()

  # double check:
  sample_words <- slice_sample(sandboxDat, n = 20) #you can run this until a mix
  sample_words %>% pull(misprod)
  sample_words %>% pull(misprod) %>% as.factor

  identical(as.factor(sample_words$misprod), sample_words$misprod) # TRUE
  class(as.factor(sample_words$misprod))
  class(sample_words$misprod)
  contrasts(as.factor(sample_words$misprod))
  contrasts(sample_words$misprod)

  # confirm contrasts are as you expect
  sandboxDat %>% select(where(is.factor)) %>% map(contrasts)

  ### this gives the following

  #
  # $pronouns
  # other she/her they/them undisclosed
  # he/him          0       0         0           0
  # other           1       0         0           0
  # she/her         0       1         0           0
  # they/them       0       0         1           0
  # undisclosed     0       0         0           1
  #
  # $ethnic
  # AA LX M UND W
  # A    0  0 0   0 0
  # AA   1  0 0   0 0
  # LX   0  1 0   0 0
  # M    0  0 1   0 0
  # UND  0  0 0   1 0
  # W    0  0 0   0 1
  #
  # $socclass
  # middle poor working
  # affluent      0    0       0
  # middle        1    0       0
  # poor          0    1       0
  # working       0    0       1
  #
  ###

  # # alternate version
  # sandboxDat <- errorDat
  # #sandboxDat$challengeACC <- zero_to_neg_one(sandboxDat$challengeACC)
  #
  # sandboxDat$challengeACC <- binary_to_plus_minus_factor(sandboxDat$challengeACC)
  # # confirm:
  # class(sandboxDat$challengeACC) # -> "factor"
  #
  # # confirm:
  # contrasts(sandboxDat$challengeACC) # -> -1 (i.e. no error): 0, 1 (i.e. error): 1
  #
  # # convert all relevant columns
  # testDat <-
  #   mutate(sandboxDat, across(where(is.logical), binary_to_plus_minus_factor))
  # testDat %>% select(where(is.logical)) # nothing
  #
  # # confirm contrasts are as you expect
  # testDat %>% select(where(is.factor)) %>% map(contrasts)


}


#center continuous predictors
errorDat$age_gmc <- errorDat$age - mean(errorDat$age)
errorDat$bfne_gmc <- errorDat$bfne - mean(errorDat$bfne)
errorDat$phq8_gmc <- errorDat$phq8 - mean(errorDat$phq8)
errorDat$scaaredTotal_gmc <- errorDat$scaaredTotal - mean(errorDat$scaaredTotal)
errorDat$scaaredGA_gmc <- errorDat$scaaredGA - mean(errorDat$scaaredGA)
errorDat$scaaredSoc_gmc <- errorDat$scaaredSoc - mean(errorDat$scaaredSoc)
errorDat$sps_gmc <- errorDat$sps - mean(errorDat$sps)

# First we'll encapsulate centering so we don't write each column name out
# manually three times. This way a typo in a column name can't screw it up. (It
# throws an error instead, rather than allowing in the incorrect data.)
gmc <- function(vec) { vec - mean(vec) }

add_gmc <- function(df, col) { # makes a new column: centered version of variable
  vec <- pull(df, {{col}}) # what's the data we're centering?

  df %>%
    mutate("{{col}}_gmc" := gmc(vec)) # ex. challengeAvgSub->challengeAvgSub_gmc
}

sandboxDat %>% add_gmc(log10frequency) %>% slice_sample(n = 50) %>%  as.data.frame()


errorDat <-
  errorDat %>%
  add_gmc(log10frequency)


# sanity check
for (col_index in which(stringr::str_detect(colnames(errorDat), ".*_gmc"))) {
  colname <- names(errorDat[col_index])
  col <- errorDat[[col_index]]
  print(colname)

  avg <- mean(col)
  print(paste('  mean:', avg,
              '  rounded mean:', round(avg, digits = 10)))
}

### SECTION 3.5: preparing for misprod-hes sequential analyses

# ignore the misprod-hes columns for now
errorDatMisprodHes <- select(errorDat, !contains("any_"))

# First: look at a given misproduction and check for nearby hesitations
justMisprodWithHesBefore <- cbind(errorDatMisprodHes,
                                  hes_position = 0, # "before",
                                  misprod_in_adjacent_window = errorDat$misprod_with_any_prior_hesitation)

justMisprodWithHesAfter <- cbind(errorDatMisprodHes,
                                 hes_position = 1, # "after",
                                 misprod_in_adjacent_window = errorDat$misprod_with_any_upcoming_hesitation)


# stack the ones before and the ones after as rows of a single df (my attempt at long form)
errorDatLongMisprodWithRelHes <- rbind(justMisprodWithHesBefore, justMisprodWithHesAfter)

# track the binary relative position as a factor
errorDatLongMisprodWithRelHes$hes_position <- as.factor(errorDatLongMisprodWithRelHes$hes_position)

# Then: look at a given hesitation and check for nearby misproductions
justHesWithMisprodBefore <- cbind(errorDatMisprodHes,
                                  misprod_position = 0, # "before",
                                  hes_in_adjacent_window = errorDat$hesitation_with_any_prior_misprod)

justHesWithMisprodAfter <- cbind(errorDatMisprodHes,
                                 misprod_position = 1, # "after",
                                 hes_in_adjacent_window = errorDat$hesitation_with_any_upcoming_misprod)

# stack the ones before and the ones after as rows of a single df (my attempt at long form)
errorDatLongHesWithRelMisprod <- rbind(justHesWithMisprodBefore, justHesWithMisprodAfter)

# track the binary relative position as a factor
errorDatLongHesWithRelMisprod$misprod_position <- as.factor(errorDatLongHesWithRelMisprod$misprod_position)



### SECTION 4: MODEL RESULTS
#misprod x bfne
# model1 <- lmerTest::lmer(misprod ~ bfne_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model1)

#misprod x scaaredSoc
model2 <- lmerTest::lmer(misprod ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model2)

#misprod x scaaredSoc control for word
model2.5 <- lmerTest::lmer(misprod ~ scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                         data=errorDat, REML=TRUE)
summary(model2.5)
# tldr 2/28/24 SA indv. do not misproduce more/less


#misprod x sps
# model3 <- lmerTest::lmer(misprod ~ sps_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model3)

#hesitation x bfne
# model4 <- lmerTest::lmer(hesitation ~ bfne_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model4)

#hesitation x scaaredSoc
model5 <- lmerTest::lmer(hesitation ~ scaaredSoc_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model5)

# hesitation x scaaredSoc, control for word
model5.5 <- lmerTest::lmer(hesitation ~ scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                         data=errorDat, REML=TRUE)
summary(model5.5)

# results are similar
# tldr 2/28/24: SA indv.s DO hesitate more
# todo plot?, maybe?


#hesitation x sps
# model6 <- lmerTest::lmer(hesitation ~ sps_gmc + (1|id) + (1|passage),
#                          data=errorDat, REML=TRUE)
# summary(model6)


#### supplemental analyses
# see notes

# Now, misproduction-hesitation relationships

# Errors as explained by disfluency: rate of misproduced syllables from rate of hesitated syllables
f_model20 <- lmerTest::lmer(misprod ~ hesitation + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model20) # ***
# should we (1|word) here?
# no... that's not rate, is it? 3/27/24
# that's whether a word being hesitated on predicts it (that same word) being misproduced
# right?

# Errors as explained by disfluency: rate of misproduced words from rate of hesitated words
# f_model21 <- lmerTest::lmer(words_with_misprod ~ words_with_hes + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model21) # ***
#
#
# # Errors as explained by disfluency: rate of misproduced words from rate of hesitated syllables
# f_model22 <- lmerTest::lmer(words_with_misprod ~ hesitation + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model22) # ***
# commented out 2/28/24 bc word level isn't applicable here; _all_ of this is at
# word level
# NB my * comments here (this section of models at least) are out of date


# Now, misproduction-hesitation interactions with social anxiety

# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and scaared
f_model23 <- lmerTest::lmer(misprod ~ hesitation * scaaredSoc_gmc + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model23)

# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated words and scaared
# f_model24 <- lmerTest::lmer(words_with_misprod ~ words_with_hes * scaaredSoc_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model24)
# commented out 2/28/24 bc word level isn't applicable here; _all_ of this is at
# word level



# Errors as explained by disfluency and SA: rate of misproduced words from rate of hesitated syllables and scaared
# f_model25 <- lmerTest::lmer(words_with_misprod ~ hesitation * scaaredSoc_gmc + (1|id) + (1|passage),
#                             data=errorDat, REML=TRUE)
# summary(f_model25)
# commented out 2/28/24 bc word level isn't applicable here; _all_ of this is at
# word level



# What happens when we control for age?
#hesitation x scaaredSoc
age_model1 <- lmerTest::lmer(hesitation ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(age_model1)

#words_with_hes x scaaredSoc
# age_model2 <- lmerTest::lmer(words_with_hes ~ scaaredSoc_gmc + age_gmc + (1|id) + (1|passage),
#                           data=errorDat, REML=TRUE)
# summary(age_model2)
# commented out 2/28/24 bc word level isn't applicable here; _all_ of this is at
# word level


# misprod-hes ordering

# we have a number of occurrences of a misproduction in a particular position
# relative to a passage's hesitations. does knowing the position (before/after)
# predict the number of these sequences we have?

# does misproduction location relative to a hesitation predict how many
# instances we get in a particular reading?

hes_with_rel_misprod_model_1 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|passage),
                                               data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_1) # n.s., 0.271

misprod_with_rel_hes_model_1 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|passage),
                                               data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_1) # n.s., 0.108

## does it interact with SA?
hes_with_rel_misprod_model_3 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|passage),
                                               data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_3) # n.s.

misprod_with_rel_hes_model_4 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|passage),
                                               data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_4) # n.s.

# what if we control for word?
hes_with_rel_misprod_model_1.5 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|passage) + (1|word),
                                               data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_1.5) # n.s., sameish

misprod_with_rel_hes_model_1.5 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|passage) + (1|word),
                                               data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_1.5) # ., 0.0974

## does it interact with SA?
hes_with_rel_misprod_model_3.5 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                               data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_3.5) # n.s.

misprod_with_rel_hes_model_4.5 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                               data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_4.5) # n.s.

# and if we ignore passage?
hes_with_rel_misprod_model_1.6 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|word),
                                                 data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_1.6) # n.s., sameish

misprod_with_rel_hes_model_1.6 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|word),
                                                 data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_1.6) # made no difference, as you might expect

## does it interact with SA?
hes_with_rel_misprod_model_3.6 <- lmerTest::lmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|word),
                                                 data=errorDatLongHesWithRelMisprod, REML=TRUE)
summary(hes_with_rel_misprod_model_3.6) # ""

misprod_with_rel_hes_model_4.6 <- lmerTest::lmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|word),
                                                 data=errorDatLongMisprodWithRelHes, REML=TRUE)
summary(misprod_with_rel_hes_model_4.6) # ""



# Word frequency analysis with words absent from corpus dropped
# Does a word's frequency predict hesitation on that word?
errorDatAttestedFreqs <- filter(errorDat, log10frequency > 0)
wordfreq_model_1 <- lmerTest::lmer(hesitation ~ log10frequency + (1|id) + (1|passage) + (1|word),
                                   data=errorDatAttestedFreqs, REML=TRUE)
summary(wordfreq_model_1)
wordfreq_model_2 <- lmerTest::lmer(misprod ~ log10frequency + (1|id) + (1|passage) + (1|word),
                                   data=errorDatAttestedFreqs, REML=TRUE)
summary(wordfreq_model_2)


# Do social anxiety and frequency interact to predict hesitation rate or misproduction rate?
wordfreq_model_3 <- lmerTest::lmer(hesitation ~ log10frequency * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                   data=errorDatAttestedFreqs, REML=TRUE)
summary(wordfreq_model_3) # looks good!

wordfreq_model_4 <- lmerTest::lmer(misprod ~ log10frequency * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                   data=errorDatAttestedFreqs, REML=TRUE)
summary(wordfreq_model_4) # tldr no?



summary(errorDatAttestedFreqs$log10frequency)


# Word frequency analysis with words absent from corpus set to corpus median
# subtlexus_minimum = 0.301 # or: # subtlexus %>% select(Lg10WF) %>% min
subtlexus_median = 1 # or: # median(subtlexus$Lg10WF)
errorDat$log10frequency_with_absents_as_median <- case_match(
  errorDat$log10frequency,
  0 ~ subtlexus_median,
  .default = errorDat$log10frequency)

compare_freq <- data.frame(cbind(old = errorDat$log10frequency,
                                 new = errorDat$log10frequency_with_absents_as_median))

filter(compare_freq, old != new) %>% # confirm it worked as expected
  filter(old != 0 | new != subtlexus_median) %>%
  nrow == 0 # TRUE


# Does a word's frequency predict hesitation on that word?
wordfreq_model_with_absents_as_median_1 <- lmerTest::lmer(hesitation ~ log10frequency_with_absents_as_median + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
summary(wordfreq_model_with_absents_as_median_1)

# "" misprod on that word?
wordfreq_model_with_absents_as_median_2 <- lmerTest::lmer(misprod ~ log10frequency_with_absents_as_median + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
summary(wordfreq_model_with_absents_as_median_2)


# Do social anxiety and frequency interact to predict hesitation rate or misproduction rate?
wordfreq_model_with_absents_as_median_3 <- lmerTest::lmer(hesitation ~ log10frequency_with_absents_as_median * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
summary(wordfreq_model_with_absents_as_median_3)

wordfreq_model_with_absents_as_median_4 <- lmerTest::lmer(misprod ~ log10frequency_with_absents_as_median * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                     data=errorDat, REML=TRUE)
summary(wordfreq_model_with_absents_as_median_4)


# Now use the median!

# Do frequency and the presence of a hesitation interact to predict the presence of a misproduction?
wordfreq_model_5 <- lmerTest::lmer(misprod ~ log10frequency_with_absents_as_median * hesitation + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
summary(wordfreq_model_5)

# Do frequency, SA, and the presence of a hesitation all interact to predict the presence of a misproduction?
wordfreq_model_6 <- lmerTest::lmer(misprod ~ log10frequency_with_absents_as_median * hesitation * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
summary(wordfreq_model_6)


# todo check comprehension
word_level_comprehension_model_0 <-
  lmerTest::lmer(challengeACC ~ scaaredSoc + (1|id) + (1|passage) + (1|word),
                 data=errorDat, REML = TRUE)
summary(word_level_comprehension_model_0)

word_level_comprehension_model_1 <-
  lmerTest::lmer(challengeACC ~ misprod + (1|id) + (1|passage) + (1|word),
                 data=errorDat, REML = TRUE)
summary(word_level_comprehension_model_1)

word_level_comprehension_model_2 <-
  lmerTest::lmer(challengeACC ~ hesitation + (1|id) + (1|passage) + (1|word),
                 data=errorDat, REML = TRUE)
summary(word_level_comprehension_model_2)

word_level_comprehension_model_3 <-
  lmerTest::lmer(challengeACC ~ misprod * scaaredSoc + (1|id) + (1|passage) + (1|word),
                 data=errorDat, REML = TRUE)
summary(word_level_comprehension_model_3)

word_level_comprehension_model_4 <-
  lmerTest::lmer(challengeACC ~ hesitation * scaaredSoc + (1|id) + (1|passage) + (1|word),
                 data=errorDat, REML = TRUE)
summary(word_level_comprehension_model_4)

# effects
# eff <- effect("log10frequency_with_absents_as_median", wordfreq_model_with_absents_as_median_1)
# plot(eff, se = TRUE, rug = FALSE, xlab = "log10frequency_with_absents_as_median", ylab = "hesitation", col.points = "red", col.lines = "blue", lty = 1)
# eff_noabsents <- effect("log10frequency", wordfreq_model_1)
# plot(eff_noabsents, se = TRUE, rug = FALSE, xlab = "log10frequency", ylab = "hesitation", col.points = "red", col.lines = "blue", lty = 1)

plot_lmer <- function(model, predictor, outcome, xlab = predictor, ...) {
  # NB `outcome` will not catch your mistake; it's just a label
  eff <- effect(predictor, model)
  #dots <- substitute(alist(...)) # eval(substitute(alist(...)))
  #print(dots)
  #print(names(dots))
  plot(eff, se = TRUE, rug = FALSE, xlab = xlab, ylab = outcome,
       col.points = "red", col.lines = "blue", lty = 1,
       ...)
}

# as in
plot_lmer(wordfreq_model_1, "log10frequency", "hesitation")
plot_lmer(wordfreq_model_with_absents_as_median_1, "log10frequency_with_absents_as_median", "hesitation")
plot_lmer(wordfreq_model_with_absents_as_median_2, "log10frequency_with_absents_as_median", "misprod")
plot_lmer(wordfreq_model_2, "log10frequency", "misprod")


# plot_lmer(wordfreq_model_with_absents_as_median_1,
#           "log10frequency_with_absents_as_median",
#           "Hesitation rate\n(% of all words)",
#           main = "Hesitation rate and word familiarity",
#           xlab = expression(
#             atop("log"['10']*" word frequency",
#                  "(lower = rarer)")))

plot_lmer(model2.5,
          "scaaredSoc_gmc",
          "Misproduction rate\n(% of all words)",
          main = "Speech errors and social anxiety",
          xlab = expression("SCAARED-Social Score"))

plot_lmer(model5.5,
          "scaaredSoc_gmc",
          "Hesitation rate\n(% of all words)",
          main = "Hesitation and social anxiety",
          xlab = expression("SCAARED-Social Score"))


# hesitation ~ wf x SA
interact_plot(model = wordfreq_model_3,
              pred = log10frequency, modx = scaaredSoc_gmc, interval = TRUE)

interact_plot(model = wordfreq_model_with_absents_as_median_3,
              pred = log10frequency_with_absents_as_median,
              modx = scaaredSoc_gmc, interval = TRUE)

interact_plot(model = wordfreq_model_with_absents_as_median_3,
              pred = log10frequency_with_absents_as_median,
              modx = scaaredSoc_gmc,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")), #'testerx',
              y.label = 'Hesitation rate\n(% of all words)',
              legend.main = "SCAARED-Social score")

# misprod ~ wf x SA
interact_plot(model = wordfreq_model_4,
              pred = log10frequency, modx = scaaredSoc_gmc, interval = TRUE)

interact_plot(model = wordfreq_model_with_absents_as_median_4,
              pred = log10frequency_with_absents_as_median, modx = scaaredSoc_gmc, interval = TRUE)

interact_plot(model = wordfreq_model_with_absents_as_median_4,
              pred = log10frequency_with_absents_as_median,
              modx = scaaredSoc_gmc,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")), #'testerx',
              y.label = 'Misproduction rate\n(% of all words)',
              legend.main = "SCAARED-Social score")

# misprod ~ hes x wf
interact_plot(model = wordfreq_model_5,
              pred = log10frequency_with_absents_as_median, modx = hesitation, interval = TRUE)

interact_plot(model = wordfreq_model_6,
              pred = log10frequency_with_absents_as_median,
              modx = hesitation, mod2 = scaaredSoc_gmc,
              interval = TRUE)

# summary(errorDat$log10frequency_with_absents_as_median)

