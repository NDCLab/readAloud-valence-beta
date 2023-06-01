# Reading in CSVs of preprocessed error data and participant information,
# writing a new CSV including the preprocessed data as well as comprehension,
# demographic, reading, and passage information
#
# Luc Sahar and Jessica M. Alexander -- NDCLab, Florida International University
# last updated 5/31/23

# NB passages "sun" and "broccoli" as coded contain errors. Namely, broccoli had
# "iodized _table_ counteracts" instead of the intended "table salt", and sun
# showed participants "_empower_ individuals" whereas it was coded as "enable"

library(readxl) # read_xlsx
library(stringr) # str_extract
library(dplyr) # most things
library(purrr) # map, map_df; generally good to have
library(lubridate) # now
library(readr) # write_csv

ext_default = 'csv'
tz_default = "America/New_York"
date_format_default = "%Y%m%d_%I%M%P"


build_output_filename <- function(label, ext = ext_default, timezone = tz_default, date_format = date_format_default) {
  # `label` may include the destination directory, if different from the working directory when the script is run
  current_datetime <- now(timezone) %>% format(date_format)
  paste(label, '_', current_datetime, '.', ext, sep = "")
}

collapse_by_participant <- function(filename_in, filename_out) {
  by_participant <- read_csv(filename_in) %>%
    unique %>% # dedup
    group_by(id) %>% summarize(across(misprod:total_uncorrected_errors, sum)) # summarize by participant, across all passages
    # TODO change the columns selected, once more have been added to the output of the preproc script

  write_csv(by_participant, filename_out)
  return(filename_out)
}


# base = "~/Documents/ndclab/analysis-sandbox/github-structure-mirror/readAloud-valence-dataset/derivatives/preprocessed"
base = "/home/data/NDClab/datasets/readAloud-valence-dataset/derivatives/preprocessed"
preprocessed_summary_filename = "TODO"
collapsed_filename = build_output_filename(label = paste(base, "disfluencies_subject", sep='/'))

collapse_by_participant(preprocessed_summary_file, collapsed_filename)