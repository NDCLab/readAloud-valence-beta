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
              x.label = expression(
                atop("log"['10']*" word frequency",
                             "(lower = more rare)")),
              y.label =  expression(
                atop("Probability of misproduction",
                     "(word-level)")),
              legend.main = # expression(
              #  atop(
               "Adjacent hesitation (within 5 words)",
              #        "(before or after misproduction in question)")),
              modx.values = factor(c(-1, 1)), # implicit, but specifying s.t. labels are guaranteed to align with the right value
              modx.labels = c("no adjacent hesitation", "adjacent hesitation"),
              main.title = 'Item-Level Misproduction, Adjacent Hesitation, and Word Frequency'
) + theme(plot.title = element_text(hjust = 0.5))

# stats table
hesitation_adjacent_misproduction_table_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa <-
  tab_model(hesitation_adjacent_misproduction_model_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa,
            show.est = TRUE, # estimates
            string.est = "β", # ...which it will transform by default, (i.e. standardize, by way of `exp`). In the manuscript we clarify that standardized is the default when not explicitly specified otherwise. This is easy to verify, e.g. exp(-0.44) gives 0.6440364, which is exactly what the table will print by default. (where -0.44 is the "raw" estimate shown in the output of summary(model_name))
            show.se = TRUE,
            string.se = "SE",
            std.response = TRUE, # trying with this: don't restandardize what we already did
            show.stat = TRUE, # test statistic
            string.stat = "z", # determined per model output
            show.df = TRUE, # degrees of freedom, Kenward-Rogers approximation
            p.style = "numeric_stars",
            show.r2 = FALSE,
            show.icc = FALSE,
            show.ngroups = FALSE,
            show.dev = FALSE,
            show.re.var = FALSE,
            show.obs = FALSE,
            pred.labels = c("Intercept",
                            "Presence of adjacent hesitation",
                            "Word frequency",
                            "SCAARED social (z-scored)",
                            "Presence of adjacent hesitation × Word frequency",
                            "Presence of adjacent hesitation × SCAARED social (z-scored)",
                            "Word frequency × SCAARED social (z-scored)",
                            "Presence of adjacent hesitation × Word frequency × SCAARED social (z-scored)"
            ),
            dv.labels = "",
            col.order = c("est", "se", "df.error", "ci", "ci.inner", "ci.outer",
                          "stat", "p", "est", "response.level")); hesitation_adjacent_misproduction_table_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa

write_table_html_to_disk(
  hesitation_adjacent_misproduction_table_1_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa,
  tables_core_pre_post
)


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
revisedDatLookBackwards$hesitation_position <- if_else(
  as.logical(revisedDatLookBackwards$any_prior_hesitation_outcome),
  true = -1, # meaning before the syllable in question
  false = NA)

revisedDatLookBackwards$hesitation_position


revisedDatLookForwards <- revisedDatNoNAs
revisedDatLookForwards$hesitation_position <- if_else(
  as.logical(revisedDatLookForwards$any_upcoming_hesitation_outcome),
  true = 1, # meaning after the syllable in question
  false = NA)
revisedDatLookForwards$hesitation_position

revisedDatWithPosition <- rbind(revisedDatLookBackwards, revisedDatLookForwards)

revisedDatWithPosition <-
  differentiate_predictor_and_outcome(revisedDatWithPosition, "hesitation_position")
revisedDatWithPosition$hesitation_position_outcome <- NULL # because unused

# set up model 2
hesitation_adjacent_misproduction_model_2_logistic_wordfreq_with_absents_as_median_no_psg <- # starting sans passage because prior model with fewer predictors still did not initially converge
  glmer(misprod_outcome ~ any_adjacent_hesitation_predictor * hesitation_position_predictor * log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|word),
        data=revisedDatWithPosition, family = "binomial")
# fails, try again without the earlier predictor, because it is already encompassed by the new positional predictor

hesitation_adjacent_misproduction_model_3_logistic_wordfreq_with_absents_as_median_no_psg <- # starting sans passage because prior model with fewer predictors still did not initially converge
  glmer(misprod_outcome ~ hesitation_position_predictor * log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|word),
        data=revisedDatWithPosition, family = "binomial")

# raise # of iterations
hesitation_adjacent_misproduction_model_3_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa <-
  update(hesitation_adjacent_misproduction_model_3_logistic_wordfreq_with_absents_as_median_no_psg,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e9)))

# now plot
interact_plot(model = hesitation_adjacent_misproduction_model_3_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa,
              pred = log10frequency_with_absents_as_median_z,
              modx = hesitation_position_predictor,
              interval = TRUE,
              #fixme x.label = "SCAARED-Social score\n(z-scored)",
              y.label =  expression(
                atop("Probability of misproduction",
                     "(word-level)")),
)

# third model attempt
# two observations per word, once forward and once back
# hesitationLocation, which is identical to hes_position in original batch of models
# (always -1 for looking backwards, always +1 for looking forwards)
# I want to name it direction_searched_for_potential_hesitations, because it is defined even when there is no adjacent hesitation

# when direction_searched_for_potential_hesitations is +1,
# adjacent_hesitation_present_in_direction_looked is equal to any_upcoming_hesitation

# when direction_searched_for_potential_hesitations is -1,
# adjacent_hesitation_present_in_direction_looked is equal to any_prior_hesitation


# similar to what we did before, we can create two new variant dataframes and join them
# then we stack the two sets of rows on top of each other

revisedDatLookBackwardsPositional <- revisedDatNoNAs
revisedDatLookForwardsPositional <- revisedDatNoNAs

revisedDatLookBackwardsPositional$direction_searched_for_potential_hesitation <- -1 # before word in question
revisedDatLookForwardsPositional$direction_searched_for_potential_hesitation <- 1 # after word in question

# for lookback  instances of each original observation (word), adjacent_hesitation_present_in_direction_looked is just any_prior_hesitation
# for lookahead instances of each original observation (word), adjacent_hesitation_present_in_direction_looked is just any_upcoming_hesitation

revisedDatLookBackwardsPositional$adjacent_hesitation_present_in_direction_looked <- revisedDatLookBackwardsPositional$any_prior_hesitation_predictor
revisedDatLookForwardsPositional$adjacent_hesitation_present_in_direction_looked  <- revisedDatLookForwardsPositional$any_upcoming_hesitation_predictor

# combine
revisedDatWithPositionAlt <- rbind(revisedDatLookBackwardsPositional, revisedDatLookForwardsPositional)


# want to explore?
revisedDatWithPositionAlt %>%
  slice_sample(n = 200) %>%
  select(direction_searched_for_potential_hesitation,
         adjacent_hesitation_present_in_direction_looked,
         any_upcoming_hesitation_predictor,
         any_prior_hesitation_predictor,
         misprod_outcome,
         any_adjacent_hesitation_predictor) %>%
  View()

# then we would just have to differentiate predictors and outcomes, but
# unnecessary because all the new variables created will only ever be predictors

# set up model 3
hesitation_adjacent_misproduction_model_4_logistic_wordfreq_with_absents_as_median_no_psg <- # starting sans passage because prior model with fewer predictors still did not initially converge
  glmer(misprod_outcome ~ adjacent_hesitation_present_in_direction_looked * direction_searched_for_potential_hesitation * log10frequency_with_absents_as_median_z * scaaredSoc_z + (1|id) + (1|word),
        data=revisedDatWithPositionAlt, family = "binomial") # does not converge

# raise # of iterations
hesitation_adjacent_misproduction_model_4_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa <-
  update(hesitation_adjacent_misproduction_model_4_logistic_wordfreq_with_absents_as_median_no_psg,
         control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1e9)))

# now plot at least initially
interact_plot(model = hesitation_adjacent_misproduction_model_4_logistic_wordfreq_with_absents_as_median_no_psg_bobyqa,
              pred = adjacent_hesitation_present_in_direction_looked,
              modx = direction_searched_for_potential_hesitation,
              interval = TRUE,
              #fixme x.label = "SCAARED-Social score\n(z-scored)",
              y.label =  expression(
                atop("Probability of misproduction",
                     "(word-level)")),
)
