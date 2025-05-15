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

# now plot
interact_plot(model = hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa,
              pred = log10frequency_with_absents_as_median_z,
              modx = any_adjacent_hesitation_predictor,
              interval = TRUE,
             #fixme x.label = "SCAARED-Social score\n(z-scored)",
              y.label =  expression(
                atop("Probability of misproduction",
                     "(word-level)")),
              # legend.main = expression(
              #   atop("Hesitation position",
              #        "(before or after misproduction in question)")),
              # modx.values = factor(c(-1, 1)), # implicit, but specifying s.t. labels are guaranteed to align with the right value
              # modx.labels = c("preceding misproduction", "following hes_position_predictor"),
              # main.title = 'Item-Level Misproduction, Hesitation Position, and Social Anxiety'
) # + theme(plot.title = element_text(hjust = 0.5))


# second model
# two observations per word, once forward and once back
# new predictor variable is NA when there is no adjacent hesitation
# the predictor is 1 when there IS an adjacent hesitation and it is AFTER the syllable in question
# the predictor is -1 when there IS an adjacent hesitation and it is BEFORE the syllable in question
# having each word as two observations accounts for the possibility that a given word might have a hesitation on both sides

# so I think we double all rows
# in theory, each row corresponds to one instance of predictor = 1 and predictor = -1
# but for half of them, if any_prior_hesitation is false, we set the predictor to NA
# and for the other half, if any_upcoming_hesitation, we set the predictor to NA

# similar to what we did before, we can create two new variant dataframes and join them
# one has hes_position set to NA
# then in that dataframe, for any row where any_prior_hesitation is TRUE, we change hes_position to -1

# in the second dataframe, we set hes_position to NA again
# then for any row in that dataframe, for any row where any_upcoming_hesitation is TRUE, we change hes_position to 1

# then we stack the two sets of rows on top of each other

revisedDatLookBackwards <- revisedDatNoNAs
revisedDatLookBackwards$hesitation_position_predictor <- if_else(
  as.logical(revisedDatLookBackwards$any_prior_hesitation_outcome),
  true = -1, # meaning before the syllable in question
  false = NA)

revisedDatLookBackwards$hesitation_position_predictor


revisedDatLookForwards <- revisedDatNoNAs
revisedDatLookForwards$hesitation_position_predictor <- if_else(
  as.logical(revisedDatLookForwards$any_upcoming_hesitation_outcome),
  true = 1, # meaning after the syllable in question
  false = NA)
revisedDatLookForwards$hesitation_position_predictor

revisedDatWithPosition <- rbind(revisedDatLookBackwards, revisedDatLookForwards)


