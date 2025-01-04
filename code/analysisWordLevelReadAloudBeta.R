# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2024-12-23

# INPUTS
# data/df: behavioral data, for each participant on each passage, with relevant participant information and trial-level stimulus information

# OUTPUTS
# models, plots

# NOTES TO DO
# plot models

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
analyzed_sample_size

demo_cols %>% # as a percent (but that doesn't mean anything for min/max/median)
  discard(~. == "age") %>% # (...so drop age, which `summary` formats that way)
  map(\(col) summary_unique(dfTrim, "id", col) / analyzed_sample_size)

summary_unique(dfTrim, "id", "age", f = sd) # also do stdev for age



### SECTION 3: ORGANIZE DATA FOR MODELING
DEBUG <- FALSE # show examples proving this works
errorDat <- dfTrim


### SECTION 3.1: Correct data for contrasts and numerics according to whether
#                they are predictors or outcomes in the models to follow

#modify contrasts for categorical predictors
contrasts(errorDat$sex) <- contr.sum(2) #male: -1, female: +1

errorDatPredictorsOutcomes <- errorDat # separate them

# given limited testing, it appears to function as intended
differentiate_predictor_and_outcome <- function(df, colname, keep_col = FALSE) {
  # Return a copy of df with two new columns.
  # colname should be binary data
  # {{colname}}_predictor is colname's data as -1 and 1 -- a factor
  # {{colname}}_outcome   is colname's data as  0 and 1 -- a numeric

  # Original column is removed by default to prevent accidental use

  vec <- pull(df, colname)

  # Fail if not applicable to the data passed
  stopifnot(class(vec) %in% c("integer", "logical", "numeric"))

  numeric_data <- as.numeric(vec)
  # to be safe, we'll set even the predictor's numeric data to +1 and -1 for now
  factor_data <- replace(as.numeric(vec), which(vec == 0), -1)

  factor_data <- as.factor(factor_data) # actually make it a factor
  contrasts(factor_data) <- rev(contr.sum(2)) # and set -1, +1 contrasts

  if (DEBUG) {
    print(vec); print(numeric_data); print(factor_data)
    print(contrasts(factor_data))
  } #show it


  df <- mutate(df, "{colname}_predictor" := factor_data,
                   "{colname}_outcome"   := numeric_data)

  if(!keep_col) df <- select(df, -colname)
  return(df)
}

split_many_predictors_and_outcomes <- function(df, colnames) {
  acc <- df
  for(col in colnames) {
    # print({{col}})
    # print(class({{col}}))
    acc <- differentiate_predictor_and_outcome(acc, col)
  }
  return(acc)
}


# binary_to_plus_minus_contrast <- function(vec) {
#   vec2 <- vec %>% replace(which(vec == 0), -1) %>% as.factor()
#
#   contrasts(vec2) <- rev(contr.sum(2))
#
#   return(vec2)
# }

relevant_columns <-
  errorDat %>% select(where(is.logical), challengeACC) %>% colnames


# e.g.
errorDatPredictorsOutcomes <-
  errorDatPredictorsOutcomes %>%
  split_many_predictors_and_outcomes(relevant_columns)


#errorDat <-
#  mutate(errorDat, across(challengeACC, where(is.logical), binary_to_plus_minus_factor))
# errorDat <-
#   errorDat %>%
#   mutate(across(where(is.logical), binary_to_plus_minus_contrast)) %>%
#   mutate(challengeACC = binary_to_plus_minus_contrast(challengeACC))
#

if (DEBUG) {
  sandboxDat <- errorDatPredictorsOutcomes

  sandboxDat %>% select(where(is.logical)) %>% colnames # see what's there
  sandboxDat %>% select(where(is.factor)) %>% colnames # see what's there
  sandboxDat %>% select(where(is.numeric)) %>% colnames # see what's there

  # double check:
  sample_words <- slice_sample(sandboxDat, n = 20) #you can run this until a mix
  sample_words %>% pull(misprod) # gone!
  sample_words %>% pull(misprod_predictor)
  sample_words %>% pull(misprod_outcome)

  class(sample_words$misprod_predictor)
  class(sample_words$misprod_outcome)
  contrasts(sample_words$misprod_predictor)
  contrasts(sample_words$misprod_outcome) # -> fails

  # confirm contrasts are as you expect
  sandboxDat %>% select(where(is.factor)) %>% map(contrasts)

  ### this gives the following
  #
  # $sex
  # [,1]
  # female    1
  # male     -1
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
  # $misprod_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $ins_dup_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $omit_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $word_stress_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $filled_pause_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $hesitation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $elongation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $corrected_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $any_upcoming_hesitation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $misprod_with_any_upcoming_hesitation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $any_prior_misprod_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $hesitation_with_any_prior_misprod_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $any_upcoming_misprod_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $hesitation_with_any_upcoming_misprod_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $any_prior_hesitation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $misprod_with_any_prior_hesitation_predictor
  # [,1]
  # -1   -1
  # 1     1
  #
  # $challengeACC_predictor
  # [,1]
  # -1   -1
  # 1     1
  ###
}

### SECTION 3.2: mean-center continuous predictors

# now prevent ourselves from using data that isn't explicitly set up to be
# predictor or outcome: use the version that omitted the ambiguous columns
errorDatBackup <- errorDat
errorDat <- errorDatPredictorsOutcomes



# z-score continuous predictors
errorDat$age_z <- scale(errorDat$age, center = TRUE, scale = TRUE)
errorDat$bfne_z <- scale(errorDat$bfne, center = TRUE, scale = TRUE)
errorDat$phq8_z <- scale(errorDat$phq8, center = TRUE, scale = TRUE)
errorDat$scaaredTotal_z <- scale(errorDat$scaaredTotal, center = TRUE, scale = TRUE)
errorDat$scaaredGA_z <- scale(errorDat$scaaredGA, center = TRUE, scale = TRUE)
errorDat$scaaredSoc_z <- scale(errorDat$scaaredSoc, center = TRUE, scale = TRUE)
errorDat$sps_z <- scale(errorDat$sps, center = TRUE, scale = TRUE)

# First we'll encapsulate centering so we don't write each column name out
# manually three times. This way a typo in a column name can't screw it up. (It
# throws an error instead, rather than allowing in the incorrect data.)
gmc <- function(vec) { vec - mean(vec) }

add_gmc <- function(df, col) { # makes a new column: centered version of variable
  vec <- pull(df, {{col}}) # what's the data we're centering?

  df %>%
    mutate("{{col}}_gmc" := gmc(vec)) # ex. challengeAvgSub->challengeAvgSub_gmc
}

if (DEBUG) {
  sandboxDat %>%
  add_gmc(log10frequency) %>%
  slice_sample(n = 50) %>%
  as.data.frame()
}

# errorDat <-
#   errorDat %>%
#   add_gmc(log10frequency)

errorDat$log10frequency_z <- scale(errorDat$log10frequency, center = TRUE, scale = TRUE)[,1]

# sanity check
if (DEBUG) {
  for (col_index in which(stringr::str_detect(colnames(errorDat), ".*_gmc"))) {
    colname <- names(errorDat[col_index])
    col <- errorDat[[col_index]]
    print(colname)

    avg <- mean(col)
    print(paste('  mean:', avg,
                '  rounded mean:', round(avg, digits = 10)))
  }
}

### SECTION 3.3: preparing for misprod-hes sequential analyses

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
# again, now prevent ourselves from using data that isn't explicitly set up to
# be predictor or outcome: use the version that omitted the ambiguous columns

# set up plotting in same area so mistakes are less likely

plot_glmer <- function(model, predictor, outcome, xlab = predictor, ...) {
  # NB `outcome` will not catch your mistake; it's just a label
  eff <- effect(predictor, model)
  plot(eff, se = TRUE, rug = FALSE, xlab = xlab, ylab = outcome,
       col.points = "red", col.lines = "blue", lty = 1,
       ...)
}

# template:
# plot_glmer(model,
#            predictor = '',
#            outcome = '',
#            xlab = '',
#            main = '')


# as in
if (DEBUG) {
  plot_glmer(wordfreq_model_1, "log10frequency", "hesitation")
  plot_glmer(wordfreq_model_with_absents_as_median_1, "log10frequency_with_absents_as_median", "hesitation")
  plot_glmer(wordfreq_model_with_absents_as_median_2, "log10frequency_with_absents_as_median", "misprod")
  plot_glmer(wordfreq_model_2, "log10frequency", "misprod")
}

# these versions fail - as intended: object 'misprod' not found
if (DEBUG) {
  #misprod x scaaredSoc
  old_model2 <- lmerTest::lmer(misprod ~ scaaredSoc_gmc + (1|id) + (1|passage),
                           data=errorDat, REML=TRUE)
  summary(old_model2)

  #misprod x scaaredSoc control for word
  old_model2.5 <- lmerTest::lmer(misprod ~ scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                 data=errorDat, REML=TRUE)
}

if(DEBUG) { # compare to how it was/would've been before splitting outcome and predictor data
  wrong_model2.5 <- lmerTest::lmer(misprod_predictor ~ scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, REML=TRUE)
  # Error in mkRespMod(fr, REML = REMLpass) : response must be numeric
  summary(wrong_model2.5)
}



# n.s.
model2.5_z_scored_logistic <- glmer(misprod_outcome ~ scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                           data=errorDat, family = "binomial")
summary(model2.5_z_scored_logistic)

plot_glmer(model2.5_z_scored_logistic,
           predictor = 'scaaredSoc_z',
           outcome = 'Probability of misproduction \n(word-level)',
           xlab = 'SCAARED-Social Score\n(z-scored)',
           main = 'Social Anxiety Severity and Item-Level Misproductions')


# hesitation x scaaredSoc_z, control for word, binary outcome
# scaaredSoc_z  0.21491    0.07715   2.786  0.00534 **
model5.5_z_scored_logistic <- glmer(hesitation_outcome ~ scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                                    data=errorDat, family = "binomial")
summary(model5.5_z_scored_logistic)


# corrected
plot_glmer(model5.5_z_scored_logistic,
           predictor = 'scaaredSoc_z',
           outcome = 'Probability of hesitation\n(word level)',
           xlab = 'SCAARED-Social Score\n(z-scored)',
           main = 'Social Anxiety Severity and Item-Level Hesitations')




#### supplemental analyses
# see notes

# Now, misproduction-hesitation relationships

# Errors as explained by disfluency: rate of misproduced syllables from rate of hesitated syllables
# hesitation_predictor1  0.39511    0.03292   12.00   <2e-16 ***
f_model20_logistic <- glmer(misprod_outcome ~ hesitation_predictor + (1|id) + (1|passage),
                            data=errorDat, family = "binomial")
summary(f_model20_logistic)

effect(hesitation_predictor, f_model20_logistic) %>%
  plot(se = TRUE, rug = FALSE, xlab = "Hesitation (presence/absence)", ylab = "Probability of misproduction\n(word-level)",
       col.points = "red", col.lines = "blue", lty = 1)



plot_glmer(f_model20_logistic,
           predictor = 'hesitation_predictor',
           outcome = 'Probability of misproduction\n(word-level)',
           xlab = 'Hesitation',
           main = 'Item-Level Misproductions and Hesitations')


# should we (1|word) here?
# Errors as explained by disfluency: rate of misproduced syllables from rate of hesitated syllables, control for word
# hesitation_predictor1  0.09129    0.03670   2.488   0.0129 *
f_model20.5_logistic <- glmer(misprod_outcome ~ hesitation_predictor + (1|id) + (1|passage) + (1|word),
                            data=errorDat, family = "binomial")
summary(f_model20.5_logistic)

plot_glmer(f_model20.5_logistic,
           predictor = 'hesitation_predictor',
           outcome = 'Probability of misproduction\n(word-level)',
           xlab = 'Hesitation',
           main = 'Item-Level Misproductions and Hesitations')


# Now, misproduction-hesitation interactions with social anxiety

# Errors as explained by disfluency and SA: rate of misproduced syllables from rate of hesitated syllables and scaared
# hesitation_predictor1               0.396320   0.033051  11.991   <2e-16 ***
# scaaredSoc_z                        0.006507   0.078505   0.083    0.934
# hesitation_predictor1:scaaredSoc_z -0.013847   0.031155  -0.444    0.657
f_model23_z_scored_logistic <- glmer(misprod_outcome ~ hesitation_predictor * scaaredSoc_z + (1|id) + (1|passage),
                            data=errorDat, family = "binomial")
summary(f_model23_z_scored_logistic)


interact_plot(model = f_model23_z_scored_logistic,
              pred = scaaredSoc_z,
              modx = hesitation_predictor,
              interval = TRUE,
              x.label = "SCAARED-Social score\n(z-scored)",
              y.label = expression('Probability of misproduction (word-level)'),
              legend.main = "Presence/absence of hesitation (word-level)",
              main.title = "Social Anxiety Severity and Item-Level Hesitations and Misproductions") +
   theme(plot.title = element_text(hjust = 0.5))

interact_plot(model = f_model23_z_scored_logistic,
              pred = hesitation_predictor,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = "Presence/absence of hesitation (word-level)",
              y.label = expression('Probability of misproduction (word-level)'),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Social Anxiety Severity and Item-Level Hesitations and Misproductions") +
  theme(plot.title = element_text(hjust = 0.5))






# Errors as explained by disfluency and SA: misproduction from hesitation and scaared_z, control for word
# hesitation_predictor1               0.09261    0.03686   2.512    0.012 *
# scaaredSoc_z                        0.01159    0.08326   0.139    0.889
# hesitation_predictor1:scaaredSoc_z -0.01278    0.03332  -0.384    0.701
f_model23.5_z_scored_logistic <- glmer(misprod_outcome ~ hesitation_predictor * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                            data=errorDat, family = "binomial")
summary(f_model23.5_z_scored_logistic)

interact_plot(model = f_model23.5_z_scored_logistic,
              pred = scaaredSoc_z,
              modx = hesitation_predictor,
              interval = TRUE,
              x.label = "SCAARED-Social score\n(z-scored)",
              y.label = expression('Probability of misproduction (word-level)'),
              legend.main = "Presence/absence of hesitation (word-level)",
              main.title = "Social Anxiety Severity and Item-Level Hesitations and Misproductions") +
  theme(plot.title = element_text(hjust = 0.5))

interact_plot(model = f_model23.5_z_scored_logistic,
              pred = hesitation_predictor,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = "Presence/absence of hesitation (word-level)",
              y.label = expression('Probability of misproduction (word-level)'),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Social Anxiety Severity and Item-Level Hesitations and Misproductions") +
  theme(plot.title = element_text(hjust = 0.5))



# What happens when we control for age?
# hesitation x scaaredSoc_z, cf. model5.5
# scaaredSoc_z  0.19651    0.07384   2.661  0.00779 **
# age_z         0.16195    0.07340   2.206  0.02737 *
age_model1_z_scored_logistic <- glmer(hesitation_outcome ~ scaaredSoc_z + age_z + (1|id) + (1|passage),
                         data=errorDat, family = "binomial")
summary(age_model1_z_scored_logistic)

plot_glmer(age_model1_z_scored_logistic,
           predictor = 'scaaredSoc_z',
           outcome = 'Probability of hesitation\n(word level)',
           xlab = 'SCAARED-Social score\n(z-scored)',
           main = 'Social Anxiety Severity and Item-Level Hesitations\n(accounting for age)')


# What happens when we control for age?
# hesitation x scaaredSoc_z, control for word, cf. model5.5
# scaaredSoc_z  0.20149    0.07361   2.737   0.0062 **
# age_z         0.16721    0.07311   2.287   0.0222 *
age_model1.5_z_scored_logistic <- glmer(hesitation_outcome ~ scaaredSoc_z + age_z + (1|id) + (1|passage) + (1|word),
                                      data=errorDat, family = "binomial")
summary(age_model1.5_z_scored_logistic)

plot_glmer(age_model1.5_z_scored_logistic,
           predictor = 'scaaredSoc_z',
           outcome = 'Probability of hesitation\n(word level)',
           xlab = 'SCAARED-Social score\n(z-scored)',
           main = 'Social Anxiety Severity and Item-Level Hesitations\n(accounting for age)')


# misprod-hes ordering

# we have a number of occurrences of a misproduction in a particular position
# relative to a passage's hesitations. does knowing the position (before/after)
# predict the number of these sequences we have?

# does misproduction location relative to a hesitation predict how many
# instances we get in a particular reading?

# TODO unfixed per earlier predictor/outcome differentiation
if (FALSE) {
  hes_with_rel_misprod_model_1 <- glmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|passage),
                                                 data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_1) # n.s., 0.271

  misprod_with_rel_hes_model_1 <- glmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|passage),
                                                 data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_1) # n.s., 0.108

  ## does it interact with SA?
  hes_with_rel_misprod_model_3 <- glmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|passage),
                                                 data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_3) # n.s.

  misprod_with_rel_hes_model_4 <- glmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|passage),
                                                 data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_4) # n.s.

  # what if we control for word?
  hes_with_rel_misprod_model_1.5 <- glmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|passage) + (1|word),
                                                 data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_1.5) # n.s., sameish

  misprod_with_rel_hes_model_1.5 <- glmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|passage) + (1|word),
                                                 data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_1.5) # ., 0.0974

  ## does it interact with SA?
  hes_with_rel_misprod_model_3.5 <- glmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                                 data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_3.5) # n.s.

  misprod_with_rel_hes_model_4.5 <- glmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
                                                 data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_4.5) # n.s.

  # and if we ignore passage?
  hes_with_rel_misprod_model_1.6 <- glmer(hes_in_adjacent_window ~ misprod_position + (1|id) + (1|word),
                                                   data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_1.6) # n.s., sameish

  misprod_with_rel_hes_model_1.6 <- glmer(misprod_in_adjacent_window ~ hes_position + (1|id) + (1|word),
                                                   data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_1.6) # made no difference, as you might expect

  ## does it interact with SA?
  hes_with_rel_misprod_model_3.6 <- glmer(hes_in_adjacent_window ~ misprod_position * scaaredSoc_gmc + (1|id) + (1|word),
                                                   data=errorDatLongHesWithRelMisprod, family = "binomial")
  summary(hes_with_rel_misprod_model_3.6) # ""

  misprod_with_rel_hes_model_4.6 <- glmer(misprod_in_adjacent_window ~ hes_position * scaaredSoc_gmc + (1|id) + (1|word),
                                                   data=errorDatLongMisprodWithRelHes, family = "binomial")
  summary(misprod_with_rel_hes_model_4.6) # ""
}


# # Word frequency analysis with words absent from corpus dropped
# # Does a word's frequency predict hesitation on that word?
# errorDatAttestedFreqs <- filter(errorDat, log10frequency > 0)
# wordfreq_model_1 <- glmer(hesitation_outcome ~ log10frequency_gmc + (1|id) + (1|passage) + (1|word),
#                                    data=errorDatAttestedFreqs, family = "binomial")
# summary(wordfreq_model_1)
# wordfreq_model_2 <- glmer(misprod_outcome ~ log10frequency_gmc + (1|id) + (1|passage) + (1|word),
#                                    data=errorDatAttestedFreqs, family = "binomial")
# summary(wordfreq_model_2)
#
#
# # Do social anxiety and frequency interact to predict hesitation rate or misproduction rate?
# wordfreq_model_3 <- glmer(hesitation_outcome ~ log10frequency_gmc * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
#                                    data=errorDatAttestedFreqs, family = "binomial")
# summary(wordfreq_model_3) # looks good!
#
# wordfreq_model_4 <- glmer(misprod_outcome ~ log10frequency_gmc * scaaredSoc_gmc + (1|id) + (1|passage) + (1|word),
#                                    data=errorDatAttestedFreqs, family = "binomial")
# summary(wordfreq_model_4) # tldr no?
#
# summary(errorDatAttestedFreqs$log10frequency)


# Word frequency analysis with words absent from corpus set to corpus median
subtlexus_median = 1 # or: # median(subtlexus$Lg10WF)
errorDat$log10frequency_with_absents_as_median <- case_match(
  errorDat$log10frequency,
  0 ~ subtlexus_median,
  .default = errorDat$log10frequency)

if (DEBUG) {
  compare_freq <- data.frame(cbind(old = errorDat$log10frequency,
                                   new = errorDat$log10frequency_with_absents_as_median))

  filter(compare_freq, old != new) %>% # confirm it worked as expected
    filter(old != 0 | new != subtlexus_median) %>%
    nrow == 0 # TRUE
}

# Mean center... then make sure we don't accidentally use the wrong one
# errorDat <- errorDat %>%
#   add_gmc(log10frequency_with_absents_as_median) %>%
#   select(-log10frequency_with_absents_as_median)
errorDat$log10frequency_with_absents_as_median_z <- scale(errorDat$log10frequency_with_absents_as_median, center = TRUE, scale = TRUE)[,1]
errorDat <- errorDat %>% select(-log10frequency_with_absents_as_median)

# Does a word's frequency predict hesitation on that word?
# log10frequency_with_absents_as_median_z -0.44253    0.04150  -10.66   <2e-16 ***
wordfreq_model_with_absents_as_median_1_z_scored_logistic <- glmer(hesitation_outcome ~ log10frequency_with_absents_as_median_z + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, family = "binomial")
summary(wordfreq_model_with_absents_as_median_1_z_scored_logistic)

plot_glmer(wordfreq_model_with_absents_as_median_1_z_scored_logistic,
           predictor = 'log10frequency_with_absents_as_median_z',
           outcome = 'Probability of hesitation\n(word-level)',
           xlab = expression(
             atop("log"['10']*" word frequency",
                  "(lower = rarer)")),
           main = 'Item-level Word Frequency and Hesitations')


# Does a word's frequency predict misprod on that word?
# log10frequency_with_absents_as_median_z -0.68143    0.05109  -13.34   <2e-16 ***
wordfreq_model_with_absents_as_median_2_z_scored_logistic <- glmer(misprod_outcome ~ log10frequency_with_absents_as_median_z + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, family = "binomial")
summary(wordfreq_model_with_absents_as_median_2_z_scored_logistic)

plot_glmer(wordfreq_model_with_absents_as_median_2_z_scored_logistic,
           predictor = 'log10frequency_with_absents_as_median_z',
           outcome = 'Probability of misproduction\n(word-level)',
           xlab = expression(
             atop("log"['10']*" word frequency",
                  "(lower = rarer)")),
           main = 'Item-level Word Frequency and Misproductions')


# Do social anxiety and frequency interact to the presence of a hesitation?
# log10frequency_with_absents_as_median_z              -0.444116   0.041562 -10.686  < 2e-16 ***
# scaaredSoc_z                                          0.217698   0.077599   2.805  0.00502 **
# log10frequency_with_absents_as_median_z:scaaredSoc_z  0.009726   0.013909   0.699  0.48441
wordfreq_model_with_absents_as_median_3_z_scored_logistic <- glmer(hesitation_outcome ~ log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                                   data=errorDat, family = "binomial")
summary(wordfreq_model_with_absents_as_median_3_z_scored_logistic)

# hesitation ~ wf x SA
interact_plot(model = wordfreq_model_with_absents_as_median_3_z_scored_logistic,
              pred = log10frequency_with_absents_as_median_z,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")),
              y.label =  expression(
                atop("Probability of hesitation",
                     "(word-level)")),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Item-Level Word Frequency, Social Anxiety Severity, and Item-Level Hesitations") +
  theme(plot.title = element_text(hjust = 0.5))


# Do social anxiety and frequency interact to the presence of a misproduction?
# log10frequency_with_absents_as_median_z              -0.681946   0.051084 -13.349   <2e-16 ***
# scaaredSoc_z                                          0.007296   0.078917   0.092   0.9263
# log10frequency_with_absents_as_median_z:scaaredSoc_z -0.032234   0.016058  -2.007   0.0447 *
wordfreq_model_with_absents_as_median_4_z_scored_logistic <- glmer(misprod_outcome ~ log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                                     data=errorDat, family = "binomial")
summary(wordfreq_model_with_absents_as_median_4_z_scored_logistic)

interact_plot(model = wordfreq_model_with_absents_as_median_4_z_scored_logistic,
              pred = log10frequency_with_absents_as_median_z,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")),
              y.label =  expression(
                atop("Probability of misproduction",
                     "(word-level)")),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Item-Level Word Frequency, Social Anxiety Severity, and Item-Level Misproductions") +
  theme(plot.title = element_text(hjust = 0.5))



# Do frequency and the presence of a hesitation interact to predict the presence of a misproduction?
# log10frequency_with_absents_as_median_z                       -0.73044    0.06304 -11.587   <2e-16 ***
# hesitation_predictor1                                          0.02546    0.05489   0.464    0.643
# log10frequency_with_absents_as_median_z:hesitation_predictor1 -0.05939    0.04113  -1.444    0.149
wordfreq_model_with_absents_as_median_5_z_scored_logistic <- glmer(misprod_outcome ~ log10frequency_with_absents_as_median_z * hesitation_predictor + (1|id) + (1|passage) + (1|word),
                                                                   data=errorDat, family="binomial")
summary(wordfreq_model_with_absents_as_median_5_z_scored_logistic)

interact_plot(model = wordfreq_model_with_absents_as_median_5_z_scored_logistic,
              pred = log10frequency_with_absents_as_median_z,
              modx = hesitation_predictor,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")), #'testerx',
              y.label = expression('Probability of misproduction on a given word'),
              legend.main = "Hesitation presence/absence",#\n(z-scored)",
              modx.values = factor(c(-1, 1)), # implicit, but specifying s.t. labels are guaranteed to align with the right value
              modx.labels = c("no hesitation", "hesitation"),
              main.title = "Item-Level Word Frequency, Hesitations, and Misproductions") +
  theme(plot.title = element_text(hjust = 0.5))



# Do frequency, SA, and the presence of a hesitation all interact to predict the presence of a misproduction?
#

wordfreq_model_with_absents_as_median_6_z_scored_logistic <- glmer(misprod_outcome ~ log10frequency_with_absents_as_median_z * hesitation_predictor * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                                                                   data=errorDat, family="binomial")
# does not converge
summary(wordfreq_model_with_absents_as_median_6_z_scored_logistic)


# interact_plot(model = wordfreq_model_with_absents_as_median_6_z_scored,
#               pred = log10frequency_with_absents_as_median_z,
#               modx = hesitation_predictor,
#               mod2 = scaaredSoc_z,
#               interval = TRUE,
#               x.label = expression(
#                 atop("log"['10']*" word frequency",
#                      "(lower = rarer)")), #'testerx',
#               y.label = expression('Probability of misproduction (word-level)'),
#               legend.main = "Hesitation presence/absence",#\n(z-scored)",
#               modx.values = factor(c(-1, 1)), # implicit, but specifying s.t. labels are guaranteed to align with the right value
#               modx.labels = c("no hesitation", "hesitation"),
#               mod2.values = "mean-plus-minus",
#               mod2.labels = c("Mean SCAARED-Social score - 1 SD", "Mean SCAARED-Social score", "Mean SCAARED-Social score + 1 SD"),
#               main.title = "Social Anxiety Severity and Item-Level Word Frequency, Hesitations, and Misproductions")


# try again, method 1: raise # of iterations
wordfreq_model_with_absents_as_median_6_z_scored_logistic_bobyqa <- # does work
  update(wordfreq_model_with_absents_as_median_6_z_scored_logistic,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(wordfreq_model_with_absents_as_median_6_z_scored_logistic_bobyqa)
interact_plot(model = wordfreq_model_with_absents_as_median_6_z_scored_logistic_bobyqa,
              pred = log10frequency_with_absents_as_median_z,
              modx = hesitation_predictor,
              mod2 = scaaredSoc_z,
              interval = TRUE,
              x.label = expression(
                atop("log"['10']*" word frequency",
                     "(lower = rarer)")), #'testerx',
              y.label = expression('Probability of misproduction (word-level)'),
              legend.main = "Hesitation presence/absence",#\n(z-scored)",
              modx.values = factor(c(-1, 1)), # implicit, but specifying s.t. labels are guaranteed to align with the right value
              modx.labels = c("no hesitation", "hesitation"),
              mod2.values = "mean-plus-minus",
              mod2.labels = c("Mean SCAARED-Social score - 1 SD", "Mean SCAARED-Social score", "Mean SCAARED-Social score + 1 SD"),
              main.title = "Social Anxiety Severity and Item-Level Word Frequency, Hesitations, and Misproductions")


# try again, method 2: remove random effect
wordfreq_model_with_absents_as_median_6_z_scored_logistic_no_psg <-
  glmer(misprod_outcome ~ log10frequency_with_absents_as_median_z * hesitation_predictor * scaaredSoc_z + (1|id) + (1|word),
        data=errorDat, family="binomial")
# does not converge
# summary(wordfreq_model_with_absents_as_median_6_z_scored_logistic_no_psg)

# try again: raise # of iterations
# log10frequency_with_absents_as_median_z                                    -0.73250    0.06288 -11.649   <2e-16 ***
# hesitation_predictor1                                                       0.02969    0.05514   0.538    0.590
# scaaredSoc_z                                                               -0.00530    0.09156  -0.058    0.954
# log10frequency_with_absents_as_median_z:hesitation_predictor1              -0.05852    0.04132  -1.416    0.157
# log10frequency_with_absents_as_median_z:scaaredSoc_z                       -0.02604    0.03864  -0.674    0.500
# hesitation_predictor1:scaaredSoc_z                                         -0.01228    0.05203  -0.236    0.813
# log10frequency_with_absents_as_median_z:hesitation_predictor1:scaaredSoc_z  0.00722    0.03866   0.187    0.852
wordfreq_model_with_absents_as_median_6_z_scored_logistic_no_psg_bobyqa <-
  update(wordfreq_model_with_absents_as_median_6_z_scored_logistic_no_psg,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(wordfreq_model_with_absents_as_median_6_z_scored_logistic_no_psg_bobyqa)


# Comprehension models
# "boundary (singular) fit: see help('isSingular')" for all of the below models

# Social anxiety to predict comprehension accuracy
# n.s.
# scaaredSoc_z -0.05235    0.18485  -0.283    0.777
word_level_comprehension_model_0_logistic <-
  glmer(challengeACC_outcome ~ scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                 data=errorDat, family = "binomial")
summary(word_level_comprehension_model_0_logistic)

plot_glmer(
  word_level_comprehension_model_0_logistic,
  predictor = 'scaaredSoc_z',
  outcome = 'Probability of accurate comprehension\n(word-level)',
  xlab = 'SCAARED-Social score\n(z-scored)',
  main = 'Social Anxiety Severity and Item-Level Comprehension Accuracy'
)


# Misproduction to predict comprehension accuracy
# n.s.
# misprod_predictor1 0.009751   0.023084   0.422    0.673
word_level_comprehension_model_1_logistic <-
  glmer(challengeACC_outcome ~ misprod_predictor + (1|id) + (1|passage) + (1|word),
                 data=errorDat, family = "binomial")
summary(word_level_comprehension_model_1_logistic)

plot_glmer(
  word_level_comprehension_model_1_logistic,
  predictor = 'misprod_predictor',
  outcome = 'Probability of accurate comprehension\n(word-level)',
  xlab = 'Presence/absence of misproduction (word-level)',
  main = 'Item-Level Misproductions and Comprehension Accuracy'
)

# Hesitation to predict comprehension accuracy
# n.s.
# hesitation_predictor1  0.03030    0.02167   1.398    0.162
word_level_comprehension_model_2_logistic <-
  glmer(challengeACC_outcome ~ hesitation_predictor + (1|id) + (1|passage) + (1|word),
                 data=errorDat, family = "binomial")
summary(word_level_comprehension_model_2_logistic)

plot_glmer(
  word_level_comprehension_model_2_logistic,
  predictor = 'hesitation_predictor',
  outcome = 'Probability of accurate comprehension\n(word-level)',
  xlab = 'Presence/absence of hesitation (word-level)',
  main = 'Item-Level Hesitations and Comprehension Accuracy'
)

# Interaction between misproduction and social anxiety to predict comprehension accuracy
# n.s.
# misprod_predictor1               0.00975    0.02309   0.422    0.673
# scaaredSoc_z                    -0.04265    0.19825  -0.215    0.830
# misprod_predictor1:scaaredSoc_z  0.01022    0.02119   0.482    0.630
word_level_comprehension_model_3_logistic <-
  glmer(challengeACC_outcome ~ misprod_predictor * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                 data=errorDat, family = "binomial")
summary(word_level_comprehension_model_3_logistic)

interact_plot(model = word_level_comprehension_model_3_logistic,
              pred = misprod_predictor,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = 'Presence/absence of misproduction (word-level)',
              y.label = 'Probability of accurate comprehension\n(word-level)',
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Item-Level Misproductions, Comprehension Accuracy, and Social Anxiety Severity") +
  theme(plot.title = element_text(hjust = 0.5))

# alternative version
interact_plot(model = word_level_comprehension_model_3_logistic,
              pred = scaaredSoc_z,
              modx = misprod_predictor,
              interval = TRUE,
              x.label = 'SCAARED-Social score\n(z-scored)',
              y.label = 'Probability of accurate comprehension\n(word-level)',
              legend.main = "Presence/absence of misproduction (word-level)",
              main.title = "Item-Level Misproductions, Comprehension Accuracy, and Social Anxiety Severity") +
  theme(plot.title = element_text(hjust = 0.5))


# Interaction between hesitation and social anxiety to predict comprehension accuracy
# n.s.
# hesitation_predictor1               0.031587   0.022160   1.425    0.154
# scaaredSoc_z                       -0.058026   0.193953  -0.299    0.765
# hesitation_predictor1:scaaredSoc_z -0.005708   0.020126  -0.284    0.777
word_level_comprehension_model_4_logistic <-
  glmer(challengeACC_outcome ~ hesitation_predictor * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
                 data=errorDat, family = "binomial")
summary(word_level_comprehension_model_4_logistic)

interact_plot(model = word_level_comprehension_model_4_logistic,
              pred = hesitation_predictor,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = 'Presence/absence of hesitation (word-level)',
              y.label = 'Probability of accurate comprehension\n(word-level)',
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Item-Level Hesitations, Comprehension Accuracy, and Social Anxiety Severity") +
  theme(plot.title = element_text(hjust = 0.5))

# alternative version
interact_plot(model = word_level_comprehension_model_4_logistic,
              pred = scaaredSoc_z,
              modx = hesitation_predictor,
              interval = TRUE,
              x.label = 'SCAARED-Social score\n(z-scored)',
              y.label = 'Probability of accurate comprehension\n(word-level)',
              legend.main = "Presence/absence of hesitation (word-level)",
              main.title = "Item-Level Hesitations, Comprehension Accuracy, and Social Anxiety Severity") +
  theme(plot.title = element_text(hjust = 0.5))
