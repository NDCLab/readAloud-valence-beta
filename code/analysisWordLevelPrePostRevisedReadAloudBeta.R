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

revisedDat <- revisedDat %>%
  mutate(
    any_adjacent_hesitation =
      as.integer(any_prior_hesitation_outcome | any_upcoming_hesitation_outcome))

revisedDat <- differentiate_predictor_and_outcome(revisedDat, "any_adjacent_hesitation")

# first, we make a model with one observation per word
#
