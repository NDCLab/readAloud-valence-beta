# assumes code from analysisWordLevelReadAloudBeta.R has been loaded into session

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
library(insight); library(sjPlot) # tables
library(effects)
library(xml2) # for saving tables to disk
# library(colorblindr)


FIRSTRUN=FALSE
if (FIRSTRUN) { revisedDat <- errorDat }

# first, we make a model with one observation per word
#

# Does the PRESENCE (not position) of an adjacent hesitation predict a misproduction?
revisedDat <- revisedDat %>%
  mutate(
    any_adjacent_hesitation =
      as.integer(any_prior_hesitation_outcome | any_upcoming_hesitation_outcome))

revisedDat <- differentiate_predictor_and_outcome(revisedDat, "any_adjacent_hesitation")

# remove first and last 5 words of each passage
revisedDatNoNAs <- revisedDat %>% filter(!is.na(any_adjacent_hesitation_outcome))

hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median <-
  glmer(misprod_outcome ~ any_adjacent_hesitation_predictor * log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|passage) + (1|word),
        data=revisedDatNoNAs, family = "binomial")

# because it does not converge:
hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_bobyqa <-
  update(hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e11)))

# still fails to converge, so borrowing from the other models, we remove remove
# the passage random intercept and try again
hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_no_psg <-
  glmer(misprod_outcome ~ any_adjacent_hesitation_predictor * log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|word),
        data=revisedDatNoNAs, family = "binomial")

# because this too does not converge:
hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa <-
  update(hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_no_psg,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e9)))
