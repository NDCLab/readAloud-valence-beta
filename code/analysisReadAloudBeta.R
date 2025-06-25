# readAloud-valence-beta Reading Task Analyses
# Authors: Luc Sahar, Jessica M. Alexander
# Last Updated: 2025-05-17

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
library(interactions)
library(colorspace)
library(insight); library(sjPlot) # tables
library(effects)
library(xml2) # for saving tables to disk
library(htmlTable) # for descriptive table
# library(colorblindr)
library(MetBrewer)
library(RColorBrewer)

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
data <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/readAloudBetaData_20230825.csv'
to_omit <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/input/passages-to-omit_20230810.csv'
results_path <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/results/'

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
DEBUG <- FALSE
errorDat <- dfTrim


#modify contrasts for categorical predictors
contrasts(errorDat$sex) <- rev(contr.sum(2)) # female: -1, male: +1
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

errorDat$age_z <- scale(errorDat$age, center = TRUE, scale = TRUE)
errorDat$bfne_z <- scale(errorDat$bfne, center = TRUE, scale = TRUE)
errorDat$phq8_z <- scale(errorDat$phq8, center = TRUE, scale = TRUE)
errorDat$scaaredTotal_z <- scale(errorDat$scaaredTotal, center = TRUE, scale = TRUE)
errorDat$scaaredGA_z <- scale(errorDat$scaaredGA, center = TRUE, scale = TRUE)
errorDat$scaaredSoc_z <- scale(errorDat$scaaredSoc, center = TRUE, scale = TRUE)
errorDat$sps_z <- scale(errorDat$sps, center = TRUE, scale = TRUE)
errorDat$lenSyll_z <- scale(errorDat$lenSyll, center = TRUE, scale = TRUE)
errorDat$lenWord_z <- scale(errorDat$lenWord, center = TRUE, scale = TRUE)
errorDat$avgSyllPerWord_z <- scale(errorDat$avgSyllPerWord, center = TRUE, scale = TRUE)

errorDat$timePerSyllable_z <- scale(errorDat$timePerSyllable, center = TRUE, scale = TRUE)
errorDat$timePerWord_z <- scale(errorDat$timePerWord, center = TRUE, scale = TRUE)

errorDat$words_with_misprod_rate_z <- scale(errorDat$words_with_misprod_rate, center = TRUE, scale = TRUE)
errorDat$words_with_hes_rate_z <- scale(errorDat$words_with_hes_rate, center = TRUE, scale = TRUE)
errorDat$challengeAvgSub_z <- scale(errorDat$challengeAvgSub, center = TRUE, scale = TRUE)



# to check if gmc and scale return the same outputs
if (DEBUG) {
  res <-
    errorDat %>%
    mutate(base_R_scale = scale(words_with_misprod_rate,
                                center = TRUE,
                                scale = FALSE)) %>%
    select(words_with_misprod_rate_gmc, base_R_scale)

  print(res) # scope it out

  print(
    all.equal(res$words_with_misprod_rate_gmc, # 'manual' way we'd done it before
              res$base_R_scale[,1]))           # base R built-in way
  # -> TRUE

  rm(res)

  # sanity check: did we get means of zero
  for (col_index in which(stringr::str_detect(colnames(errorDat), ".*_gmc"))) {
    colname <- names(errorDat[col_index])
    col <- errorDat[[col_index]]
    print(colname)

    avg <- mean(col)
    print(paste('  mean:', avg,
                '  rounded mean:', round(avg, digits = 10)))
  }
}

# make sure we don't accidentally use the wrong version
errorDat <- select(errorDat, -ends_with('_gmc'))

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


### SECTION 4: MODEL RESULTS AND PLOTS
# for every model involving comprehension accuracy, rather than errorDat we use
# errorDatPredictorsOutcomes, which differentiates how comprehension accuracy is
# represented as a predictor versus as an outcome

# generic plotting wrapper:
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

# visual details
rwe_palette <- brewer.pal(4, "Purples")
rwe_palette <- colorRampPalette(rwe_palette)(17)
rwe_palette <- rwe_palette[4:17]

# helpers

#helper functions
digit_display <- function(number){
  if(abs(number)<0.001){
    x <- sprintf("%.4f", number)
  }else{
    x <- sprintf("%.3f", number)
  }
  return (x) #this is a string
}

tinyps <- function(pval){
  if(pval < 0.001){
    x = "< 0.001"
  }else {
    tmp <- round(pval,3)
    x = paste("= ", as.character(tmp), sep="")
  }
  return (x) #this is a string
}


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

# fix: z-score
model8_z_scored <- lmerTest::lmer(words_with_misprod_rate_z ~ scaaredSoc_z + (1|id) + (1|passage),
                         data=errorDat, REML=TRUE)
summary(model8_z_scored)

plot_lmer(model8_z_scored,
          predictor = 'scaaredSoc_z',
          outcome = 'Rate of misproductions per word\n(z-scored)',
          xlab = 'SCAARED-Social Score\n(z-scored)',
          main = 'Social Anxiety Severity and Rate of Misproduction')


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

# fix: z-score

model11_z_scored <- lmerTest::lmer(words_with_hes_rate_z ~ scaaredSoc_z + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(model11_z_scored)

plot_lmer(model11_z_scored,
          predictor = 'scaaredSoc_z',
          outcome = 'Rate of hesitations per word\n(z-scored)',
          xlab = 'SCAARED-Social Score\n(z-scored)',
          main = 'Social Anxiety Severity and Rate of Hesitation')

# alt
plot_model(model11_z_scored,
           type = "pred",
           terms = "scaaredSoc_z",
           colors = "mediumpurple2"
         # colors = rwe_palette
) + theme(
  plot.title = element_text(size = 18),
  text = element_text(size = 16),
) + # geom_line(size = 2) +
 theme_bw() +
# scale_color_manual(values = c("plum3")) +
# scale_color_manual(values = rwe_palette)
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_line(linewidth = 0.6, linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.6, linetype = 'dashed', color = '#bbbbbb60'),
        axis.ticks.x = element_blank()
        #  plot.background = element_rect(color = '#ffffff')
        # rect = element_blank() #element_rect(fill = 'white')
  ) +
  labs(title = 'Social Anxiety Severity and Rate of Hesitation') +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('SCAARED-Social Score\n(z-scored)') +
  ylab('Rate of hesitations per word\n(z-scored)') +
 # scale_x_continuous(breaks = -1:5) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1, 1.5, 2)) +
  theme(plot.title = element_text(hjust = 0.05))


# Jess' version
plot_fig_2 <- function() {
  coefsmodel11z <- summary(model11_z_scored)$coef
  cis <- confint(model11_z_scored)
  b0 <- coefsmodel11z[1]
  b1 <- coefsmodel11z[2]
  se <- coefsmodel11z[4]

  #bootstrap ci ribbon
  iterations = 1000
  a <- tibble(i=rep(1:iterations,))
  a <- mutate(a, intercept=NA, beta=NA)
  for(i in 1:nrow(a)){
    rows <- sample(1:nrow(errorDat), nrow(errorDat), replace=TRUE)
    df <- errorDat[rows, c('id', 'passage', 'scaaredSoc_z', 'words_with_hes_rate_z')]
    mdl <- lme4::lmer(words_with_hes_rate_z ~ scaaredSoc_z + (1|id) + (1|passage),
                      data=df, REML=TRUE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    a[i,2] <- lme4::fixef(mdl)[1]
    a[i,3] <- lme4::fixef(mdl)[2]
  }


  #create df for annotation
  label_text <- data.frame(
    label = c(paste("\u03b2 = ", digit_display(b1),
                    "\nSE = ", digit_display(se),
                    "\nCI = [", digit_display(cis[5,1]), " - ", digit_display(cis[5,2]), "]",
                    "\np ", tinyps(coefsmodel11z[10]), sep="")),
    scaaredSoc_z = c(-1.1),
    #words_with_hes_rate_z = c(4.5)) #location for plot with all datapoints
    words_with_hes_rate_z = c(0.75)) #location for plot with limited y-axis

  #plot
  p <- ggplot(errorDat, aes(x=scaaredSoc_z, y=words_with_hes_rate_z)) +
    geom_jitter(aes(color=factor(scaaredSoc)), alpha=0.5, width=0.05, show.legend=FALSE) +
    scale_color_manual(values=rwe_palette)

  for(i in 1:nrow(a)){ #add bootstrapped lines to show confidence interval
    p <- p + geom_abline(intercept=as.numeric(a[i,2]), slope=as.numeric(a[i,3]), color=rwe_palette[3], alpha=0.1)
  }

  p <- p + geom_abline(intercept=b0, slope=b1, color=rwe_palette[14], linewidth=1) +
    guides(color=FALSE, shape=FALSE) +
    geom_label(data=label_text, aes(x=scaaredSoc_z, y=words_with_hes_rate_z, label=label), size=3) +
    ylim(-0.9, 0.9) + #remove this line for plot with all datapoints
    theme_bw() +
    theme(plot.title = element_text(size=18, hjust=0.05, face='bold'),
          text = element_text(size=16),
          panel.border = element_blank(),
          panel.grid = element_line(linewidth=0.6, linetype='dashed'),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(linewidth=0.6, linetype='dashed', color='#bbbbbb60'),
          axis.ticks.x = element_blank()) +
    labs(title="Social Anxiety Symptoms × Hesitation Rate",
         x="SCAARED-Social Score\n(z-scored)",
         y="Rate of Hesitations\n(per word, z-scored)")
  return(p)
}

#save file (adjust width/height as needed)
ggsave(file.path(outpath, "fig2.jpg"), plot=plot_fig_2(), width=8, height=5, units="in")




# same, controlling for sex:
sex_model11_z_scored <- lmerTest::lmer(words_with_hes_rate_z ~ scaaredSoc_z + sex + (1|id) + (1|passage),
                                   data=errorDat, REML=TRUE)
summary(sex_model11_z_scored)

# same, controlling for age:
age_model11_z_scored <- lmerTest::lmer(words_with_hes_rate_z ~ scaaredSoc_z + age_z + (1|id) + (1|passage),
                                       data=errorDat, REML=TRUE)
summary(age_model11_z_scored)

#words_with_hes_rate x sps
# model12 <- lmerTest::lmer(words_with_hes_rate ~ sps_gmc + (1|id) + (1|passage),
#                           data=errorDat, REML=TRUE)
# summary(model12)


#### supplemental analyses

# glmer(accuracy ~ scaaredSoc_gmc + (1|id) + (1|passage), data=errorDat, family="binomial")
# "f_" : follow-up

# Accuracy/comprehension as explained by social anxiety: scaaredSoc

# outcome is binary 0/1 numeric
f_model1_z_scored <- glmer(challengeACC_outcome ~ scaaredSoc_z + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model1_z_scored)

plot_lmer(f_model1_z_scored,
          predictor = 'scaaredSoc_z',
          outcome = 'Average comprehension accuracy',
          xlab = 'SCAARED-Social Score\n(z-scored)',
          main = 'Social Anxiety Severity and Comprehension Accuracy')

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
f_model5_z_scored <- glmer(challengeACC_outcome ~ words_with_hes_rate_z + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model5_z_scored)

plot_lmer(f_model5_z_scored,
          predictor = 'words_with_hes_rate_z',
          outcome = 'Average comprehension accuracy',
          xlab = 'Rate of hesitations per word\n(z-scored)',
          main = 'Rate of Hesitation and Comprehension Accuracy')



# Accuracy/comprehension as explained by errors: misproductions per syllable
# f_model6 <- glmer(challengeACC ~ misprod_rate + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model6)

# Accuracy/comprehension as explained by errors: misproductions per word
f_model7_z_scored <- glmer(challengeACC_outcome ~ words_with_misprod_rate_z + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model7_z_scored)




# Accuracy/comprehension as explained by disfluencies *and* SA: hesitations per syllable with scaared
# f_model8 <- glmer(challengeACC ~ hesitation_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model8)

# Accuracy/comprehension as explained by disfluencies: hesitations per word with scaared
# f_model9 <- glmer(challengeACC_outcome ~ words_with_hes_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDatPredictorsOutcomes, family = "binomial")
# summary(f_model9)

# fix: z score
f_model9_z_scored <- glmer(challengeACC_outcome ~ words_with_hes_rate_z * scaaredSoc_z + (1|id) + (1|passage),
                         data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model9_z_scored)


interact_plot(model = f_model9_z_scored,
              pred = words_with_hes_rate_z,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = expression(
                atop("Rate of hesitations per word", "(z-scored)")),
              y.label = #expression(
                #atop(
                "Average comprehension accuracy",# "(z-scored)")),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Rate of Hesitation, Social Anxiety Severity, and Comprehension Accuracy")


# Accuracy/comprehension as explained by errors: misproductions per syllable with scaared
# f_model10 <- glmer(challengeACC ~ misprod_rate * scaaredSoc_gmc + (1|id) + (1|passage),
#                   data=errorDat, family = "binomial")
# summary(f_model10)

# Accuracy/comprehension as explained by errors: misproductions per word with scaared
f_model11_z_scored <- glmer(challengeACC_outcome ~ words_with_misprod_rate_z * scaaredSoc_z + (1|id) + (1|passage),
                  data=errorDatPredictorsOutcomes, family = "binomial")
summary(f_model11_z_scored)



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

# fix: z score
f_model21_z_scored <- lmerTest::lmer(words_with_misprod_rate_z ~ words_with_hes_rate_z + (1|id) + (1|passage),
                            data=errorDat, REML=TRUE)
summary(f_model21_z_scored) # ***



plot_lmer(f_model21_z_scored,
          predictor = 'words_with_hes_rate_z',
          outcome = 'Rate of misproductions per word\n(z-scored)',
          xlab = 'Rate of hesitations per word\n(z-scored)',
          main = 'Rate of Hesitation and Rate of Misproduction')



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

# fix: z-score
f_model24_z_scored <- lmerTest::lmer(words_with_misprod_rate_z ~ words_with_hes_rate_z * scaaredSoc_z + (1|id) + (1|passage),
                                   data=errorDat, REML=TRUE)
summary(f_model24_z_scored)


interact_plot(model = f_model24_z_scored,
              pred = words_with_hes_rate_z,
              modx = scaaredSoc_z,
              interval = TRUE,
              x.label = expression(
                atop("Rate of hesitations per word", "(z-scored)")),
              y.label = expression(
                atop('Rate of misproductions per word', '(z-scored)')),
              legend.main = "SCAARED-Social score\n(z-scored)",
              main.title = "Rate of Hesitation, Social Anxiety Severity, and Rate of Misproduction")

# alt: don't plot interaction; up font

plot_model(f_model24_z_scored,
           type = "pred",
           terms = "words_with_hes_rate_z",
           colors = "#4b9bc7"
) + theme(
  plot.title = element_text(size = 18),
  text = element_text(size = 16),
) + # geom_line(size = 2) +
  theme_bw() +
  scale_color_manual(values = c("blue")) +
  theme(plot.title = element_text(size = 18),
        text = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid = element_line(linewidth = 0.6, linetype = 'dashed'),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(linewidth = 0.6, linetype = 'dashed', color = '#bbbbbb60'),
        axis.ticks.x = element_blank()
      #  plot.background = element_rect(color = '#ffffff')
        # rect = element_blank() #element_rect(fill = 'white')
        ) +
  labs(title = 'Rate of Hesitation and Rate of Misproduction') +
  theme(plot.title = element_text(face = 'bold')) +
  xlab('Rate of hesitations per word\n(z-scored)') +
  ylab('Rate of misproductions per word\n(z-scored)') +
  scale_x_continuous(breaks = -1:5) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) +
  theme(plot.title = element_text(hjust = 0.05))

# Jess' version, wip
plot_fig_3 <- function() {
  # determine degrees of purple needed for this variable
  rwe_palette_custom <- brewer.pal(4, "Purples")
  number_of_values <-
    pull(errorDat, words_with_hes_rate_z) %>%
    unique %>%
    length

  rwe_palette_custom <- colorRampPalette(rwe_palette_custom)(number_of_values+3)
  rwe_palette_custom <- rwe_palette_custom[4:(number_of_values+3)]

  coefsmodel11z <- summary(f_model24_z_scored)$coef
  cis <- confint(f_model24_z_scored)
  b0 <- coefsmodel11z[1]
  b1 <- coefsmodel11z[2]
  se <- coefsmodel11z[4]

  #bootstrap ci ribbon
  iterations = 1000
  a <- tibble(i=rep(1:iterations,))
  a <- mutate(a, intercept=NA, beta=NA)
  for(i in 1:nrow(a)){
    rows <- sample(1:nrow(errorDat), nrow(errorDat), replace=TRUE)
    df <- errorDat[rows, c('id', 'passage', 'words_with_hes_rate_z', 'words_with_misprod_rate_z')]
    mdl <- lme4::lmer(words_with_misprod_rate_z ~ words_with_hes_rate_z + (1|id) + (1|passage),
                      data=df, REML=TRUE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
    a[i,2] <- lme4::fixef(mdl)[1]
    a[i,3] <- lme4::fixef(mdl)[2]
  }


  #create df for annotation
  label_text <- data.frame(
    label = c(paste("\u03b2 = ", digit_display(b1),
                    "\nSE = ", digit_display(se),
                    "\nCI = [", digit_display(cis[5,1]), " - ", digit_display(cis[5,2]), "]",
                    "\np ", tinyps(coefsmodel11z[10]), sep="")),
    words_with_hes_rate_z = c(-1.1),
    words_with_misprod_rate_z = c(0.75)) #location for plot with limited y-axis

  #plot
  p <- ggplot(errorDat, aes(x=words_with_hes_rate_z, y=words_with_misprod_rate_z)) +
    geom_jitter(aes(color=factor(words_with_hes_rate_z)), alpha=0.5, width=0.05, show.legend=FALSE) +
    scale_color_manual(values=rwe_palette_custom)

  for(i in 1:nrow(a)){ #add bootstrapped lines to show confidence interval
    p <- p + geom_abline(intercept=as.numeric(a[i,2]), slope=as.numeric(a[i,3]), color=rwe_palette_custom[3], alpha=0.1)
  }

  p <- p + geom_abline(intercept=b0, slope=b1, color=rwe_palette_custom[number_of_values], linewidth=1) +
    guides(color=FALSE, shape=FALSE) +
    geom_label(data=label_text, aes(x=words_with_hes_rate_z, y=words_with_misprod_rate_z, label=label), size=3) +
    ylim(-0.9, 0.9) + #remove this line for plot with all datapoints
    theme_bw() +
    theme(plot.title = element_text(size=18, hjust=0.05, face='bold'),
          text = element_text(size=16),
          panel.border = element_blank(),
          panel.grid = element_line(linewidth=0.6, linetype='dashed'),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(linewidth=0.6, linetype='dashed', color='#bbbbbb60'),
          axis.ticks.x = element_blank()) +
    labs(title="Hesitation Rate × Misproduction Rate",
         x="Rate of Hesitations\n(per word, z-scored)",
         y="Rate of Misproductions\n(per word, z-scored)")
  return(p)
}
ggsave(file.path(outpath, "fig3.jpg"), plot=plot_fig_3(), width=8, height=5, units="in")



# same, controlling for sex:
sex_f_model24_z_scored <-  lmerTest::lmer(words_with_misprod_rate_z ~ words_with_hes_rate_z * scaaredSoc_z + sex + (1|id) + (1|passage),
                                          data=errorDat, REML=TRUE)
summary(sex_f_model24_z_scored)

# same, controlling for age:
age_f_model24_z_scored <-  lmerTest::lmer(words_with_misprod_rate_z ~ words_with_hes_rate_z * scaaredSoc_z + age_z + (1|id) + (1|passage),
                                          data=errorDat, REML=TRUE)
summary(age_f_model24_z_scored)



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
age_model2_z_scored <- lmerTest::lmer(words_with_hes_rate_z ~ scaaredSoc_z + age_z + (1|id) + (1|passage),
                          data=errorDat, REML=TRUE)
summary(age_model2_z_scored)


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

# Generate descriptive stats as table
errorDat %>%
  select(words_with_hes_rate, words_with_misprod_rate, challengeAvgSub, scaaredSoc, bfne, sps) %>%
  mutate(
    across(words_with_hes_rate:challengeAvgSub,
           ~ 100 * .)) %>%
  reframe(
    across(words_with_hes_rate:sps,
           \(x) { c(mean(x), sd(x)) } )) %>%
  mutate(
    across(words_with_hes_rate:challengeAvgSub,
           ~ paste0(round(., digits = 2), '%'))) %>%
  mutate(
    across(scaaredSoc:sps, ~ round(., 2))) %>%
  as.data.frame %>%
 `rownames<-`(c("Mean", "SD")) %>%
 `colnames<-`(c("Hesitation rate", "Misproduction rate", "Comprehension accuracy", "SCAARED Social", "BFNE", "SPS")) %>%
  htmlTableWidget()

# Generate tables for core analyses
SAVE_TABLES_TO_DISK <- FALSE # assume don't want to overwrite; change if desired

write_table_html_to_disk <- function(model_table, namespace, table_name = deparse(substitute(model_table))) {
  # reads variable name as default output name unless specified otherwise
  table_basename <- paste0(table_name, '.html')
  table_dirname <- paste(sep = '/', results_path, namespace)
  table_out_path <- paste(sep = '/', table_dirname, table_basename)
  table_html <- read_html(model_table$page.complete)

  if (SAVE_TABLES_TO_DISK) {
    fs::dir_create(table_dirname)
    write_html(table_html, table_out_path)
  }
}
tables_core_passage_level <- 'tables'
tables_non_sig_passage_level <- 'tables/non-sig'

table_model11_z_scored <-
  tab_model(model11_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c("Intercept", "SCAARED social (z-scored)"),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_model11_z_scored

write_table_html_to_disk(
  table_model11_z_scored,
  tables_core_passage_level
)

table_f_model24_z_scored <-
  tab_model(f_model24_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c("Intercept",
                            "Hesitation rate",
                            "SCAARED social",
                            "Hesitation rate × SCAARED social"),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_f_model24_z_scored

write_table_html_to_disk(
  table_f_model24_z_scored,
  tables_core_passage_level
)

# non significant core analyses:
table_model8_z_scored <-
  tab_model(model8_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c("Intercept", "SCAARED social (z-scored)"),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_model8_z_scored

write_table_html_to_disk(
  table_model8_z_scored,
  tables_non_sig_passage_level
)

table_f_model1_z_scored <-
  tab_model(f_model1_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
            transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "z", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c("Intercept", "SCAARED social (z-scored)"),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_f_model1_z_scored

write_table_html_to_disk(
  table_f_model1_z_scored,
  tables_non_sig_passage_level
)


# Generate tables for control analyses
ctrl_tables_sex_passage_level <- 'tables/control-analyses/sex'
ctrl_tables_age_passage_level <- 'tables/control-analyses/age'
ctrl_tables_final <- 'tables/control-alt'

table_age_model11_z_scored <-
  tab_model(age_model11_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c("Intercept", "SCAARED social (z-scored)", "Age"),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_age_model11_z_scored

write_table_html_to_disk(
  table_age_model11_z_scored,
  ctrl_tables_age_passage_level
)

table_sex_model11_z_scored <-
  tab_model(sex_model11_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
         #  transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            pred.labels = c(
              "Intercept",
              "SCAARED social (z-scored)",
              "Male sex"
            ),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_sex_model11_z_scored

write_table_html_to_disk(
  table_sex_model11_z_scored,
  ctrl_tables_sex_passage_level
)

table_age_f_model24_z_scored <-
  tab_model(age_f_model24_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            show.r2 = FALSE,
            show.icc = FALSE,
            show.ngroups = FALSE,
            show.dev = FALSE,
            show.re.var = FALSE,
            show.obs = FALSE,
            pred.labels = c(
              "Intercept",
              "Hesitation rate (z-scored)",
              "SCAARED social (z-scored)",
              "Age",
              "Hesitation rate (z-scored) × SCAARED social (z-scored)"
            ),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_age_f_model24_z_scored

write_table_html_to_disk(
  table_age_f_model24_z_scored,
  ctrl_tables_final
)


table_sex_f_model24_z_scored <-
  tab_model(sex_f_model24_z_scored,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
          # transform = exp, # make sure it standardizes betas
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "t", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            show.r2 = FALSE,
            show.icc = FALSE,
            show.ngroups = FALSE,
            show.dev = FALSE,
            show.re.var = FALSE,
            show.obs = FALSE,
            pred.labels = c(
              "Intercept",
              "Hesitation rate",
              "SCAARED social",
              "Male sex",
              "Hesitation rate × SCAARED social"
            ),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); table_sex_f_model24_z_scored

write_table_html_to_disk(
  table_sex_f_model24_z_scored,
  ctrl_tables_final
)


# save.image("~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/RWE-passage-level-with-zscoring-may-13-25.RData")
