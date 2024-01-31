# readAloud-valence-beta Analysis Preparation
# Authors: Luc Sahar and Jessica M. Alexander
# NDCLab, Florida International University
# Last updated: 2024-01-30

# INPUTS
# data/df: behavioral (error-related) data, for each participant for each word
# accDat:  comprehension question accuracy (0/1) for each participant for each
#          passage
#            as of 1/30/24 this is alphabetical by *passage* now (more easily
#            automated, less stateful, less hard coding)
# readDat: stimuli characteristics (by passage half) - 2024-01-30 - now NA?
# redcap:  participant data, incl. demographics and responses + scored factors
#          for questionnaires:
#   bfne (brief fear of negative evaluation): bfne_b_scrdTotal (fear of negative
#         evaluation total)
#   phq8 (patient health questionnaire): phq8_scrdTotal (depression scale total)
#   scaared, total (screen for adult anxiety disorders): scaared_b_scrdTotal
#        (total anxiety)
#   scaared, social (screen for adult anxiety disorders): scaared_b_scrdSoc
#        (social phobias)
#   scaared, general (screen for adult anxiety disorders): scaared_b_scrdGA
#        (general anxiety)
#   sps (social phobia scale): sias6sps6_b_scrdSPS

# OUTPUTS
# dfTrim: for each passage, for each participant, details on:
#   participant behavior: reading errors made, comprehension question accuracy
#   passage characteristics: length (syllable and word), average syllables per
#                            word
#   participant data: demographics, language history, mood and mood disorder
#                     scores


### SECTION 1: SETTING UP
library(readxl) # read_xlsx
library(stringr) # str_extract
library(dplyr) # most things
library(purrr) # map, map_df; generally good to have
library(lubridate) # now
library(readr) # write_csv

#set up defaults for output file naming
today <- Sys.Date() %>% format("%Y%m%d")

#set up directories for input/output
# main_dataset <- '/Users/jalexand/github/readAloud-valence-dataset/'
# main_analyses <- '/Users/jalexand/github/readAloud-valence-beta/'
# out_path <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/'
main_dataset <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-dataset/'
main_analyses <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/'
out_path <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/derivatives/'


#load input files
# data <- paste(main_dataset, 'derivatives/preprocessed/disfluencies_subject-x-passage_20230616_1229pm.csv', sep="", collapse=NULL)
data <- "~/Documents/ndclab/rwe-analysis-sandbox/rwe-dataset/derivatives/disfluencies_subject-x-passage-x-word_20240126_0253pm.csv"
accDat_path <- paste(main_dataset,'derivatives/preprocessed/readAloud_passage-level_summary_20220812.csv', sep="", collapse=NULL)
readDat_path <- paste(main_dataset, 'derivatives/analysisStimuli_readDat_20230614.csv', sep="", collapse=NULL)
redcap_path <- paste(main_dataset,'derivatives/preprocessed/202201v0readAloudval_SCRD_2022-06-20_1019.csv', sep="", collapse=NULL)
agedat_path <- paste(main_dataset,'derivatives/preprocessed/202201v0readAloudval_SCRD_2022-06-20_1019_ageonly.csv', sep="", collapse=NULL)
speedDat_path <- paste(main_dataset, "derivatives/preprocessed/valence-timing/timingpitch_subject-by-passage_2022-09-09.csv", sep="", collapse=NULL)
freqDat_path <- paste(main_analyses, "derivatives/prepWordFreq_readDat20230825.csv", sep="")
scaffolds_path <- paste(main_dataset, 'code/scaffolds.xlsx', sep="", collapse=NULL)
# c(data, accDat_path, readDat_path, redcap_path, agedat_path, speedDat_path, scaffolds_path) %>% fs::as_fs_path() %>% fs::is_file()
# ✅: all TRUE

all_passages <- excel_sheets(scaffolds_path)
df <- read.csv(data, na.strings='NA')
redcap <- read.csv(redcap_path, na.strings='NA') #participant questionnaire responses
agedat <- read.csv(agedat_path, na.strings='NA') #participant age information
readDat <- read.csv(readDat_path, na.strings='N') #passage-level characteristics from analysisStimuli.R
accDat <- read.csv(accDat_path, na.strings='NA', check.names=FALSE) #passage level accuracy for each subject
accDat$passage <- all_passages  #rename passages with short-name
speedDat <- read.csv(speedDat_path, na.strings='NA')
freqDat <- read.csv(freqDat_path, na.strings = 'NA')

#organize data types
# df[,3:30] <- sapply(df[,3:30],as.numeric)
# not applicable for now, we'll see

#add missing passages for 150086 so that nrow is divisible by 20
# passages_read <- df$passage[which(df$id=="150086")]
# tempdf <- data.frame(matrix(nrow=0, ncol=ncol(df)))
# colnames(tempdf) <- colnames(df)
# for(passage in 1:length(all_passages)){
#   if(all_passages[passage] %in% passages_read){next}else{
#     tempdf[nrow(tempdf) + 1,] <- c("150086", all_passages[passage], rep(NA, 30))
#   }
# }
# df <- rbind(df, tempdf)

### SECTION 2: BUILD DEMOGRAPHIC DATA DF
demoDat <- redcap[,c(1,5)]
# demoDat_imperative <- redcap[,c(1,5)]; demoDat_dplyr <- demoDat_imperative
# as we refactor, our test case is: all.equal(demoDat_imperative, demoDat_dplyr)
# this has been confirmed to work for with every new column added

#biological sex: replace numerical values with text description
demoDat$sex <- case_match(redcap$demo_b_sex_s1_r1_e1,
                          1 ~ 'male', 2 ~ 'female', 3 ~ 'intersex',
                          4 ~ 'other', 5 ~ 'unknown', .default = 'undisclosed')

#preferred pronouns: replace numerical values with text description
# for(b in 1:nrow(redcap)){ ... }
# 8 lines, 466 chars; boilerplatey

# try rewriting as case_when
demoDat$pron <- case_match (redcap$demo_b_pronouns_s1_r1_e1,
                            1 ~ "she/her", 2 ~ "he/him", 3 ~ "they/them",
                            5 ~ "other", .default = "undisclosed")
# `.default` catches both NA and everything else


#ethnicity affiliation: map to text description
demoDat$ethnic <- case_when( # first, check for self-identifying 2+ ethnicities
  redcap %>% select(matches("demo_b_ethnic_*")) %>% rowSums >= 2 ~ 'M', # multi
  redcap$demo_b_ethnic_s1_r1_e1___1 == 1 ~ 'AI', #american indian/alaska native
  redcap$demo_b_ethnic_s1_r1_e1___2 == 1 ~ 'A', #asian
  redcap$demo_b_ethnic_s1_r1_e1___3 == 1 ~ 'AA', #african american
  redcap$demo_b_ethnic_s1_r1_e1___4 == 1 ~ 'LX', #hispanic/latinx
  redcap$demo_b_ethnic_s1_r1_e1___5 == 1 ~ 'ME', #middle eastern
  redcap$demo_b_ethnic_s1_r1_e1___6 == 1 ~ 'PI', #pacific islander
  redcap$demo_b_ethnic_s1_r1_e1___7 == 1 ~ 'W', #white
  redcap$demo_b_ethnic_s1_r1_e1___8 == 1 ~ 'O', #other
  .default = 'UND' #undisclosed
)

#social class affiliation: replace numerical values with text description
demoDat$socclass <- case_match(redcap$demo_b_socclass_s1_r1_e1,
                                     1 ~ "poor", 2 ~ "working", 3 ~ "middle",
                                     4 ~ "affluent", .default = "undisclosed")

#communication disorders diagnoses: sum across childhood, adolescence, and adulthood
# nb. there was a typo in the old version: no adult diagnoses were being checked
# because the column name did not exist and `sum` with the df$col syntax did not
# catch that
demoDat$commdis <- select(redcap, matches("demo_b_comdis.*e1")) %>% rowSums

#language history: transfer directly
demoDat$eng <- redcap$demo_b_eng_s1_r1_e1 #participant monolingualism
demoDat$langhis <- redcap$demo_b_langhis_s1_r1_e1 #participant language history
demoDat$ageen <- redcap$demo_b_ageen_s1_r1_e1 #participant age of English acquisition
demoDat$profen <- redcap$demo_b_profen_s1_r1_e1 #participant English proficiency

#mood and mood disorders: transfer directly
demoDat$bfne <- redcap$bfne_b_scrdTotal
demoDat$phq8 <- redcap$phq8_scrdTotal #phq8 depression scale
demoDat$scaaredTotal <- redcap$scaared_b_scrdTotal #scaared total anxiety
demoDat$scaaredGA <- redcap$scaared_b_scrdGA #scaared general anxiety
demoDat$scaaredSoc <- redcap$scaared_b_scrdSoc #scaared social phobias
demoDat$sps <- redcap$sias6sps6_b_scrdSPS #sps social phobia scale

#age: pull from separate file
demoDat <- left_join(demoDat, # can't just assign: matching matters given new df
                     select(agedat, record_id, age = info_age_s1_r1_e1))


### SECTION 3: SET UP DERIVED FIELDS FOR SPEED ANALYSES
speedDat$readingTime <- speedDat$readEnd - speedDat$readStart
df <- left_join(df, speedDat, by = c("id", "passage")) # now reading timestamps and duration are looped into df


### SECTION 4: BUILD TRIAL-LEVEL DF (ADD DEMODAT, READDAT, ACCDAT, AND FREQDAT to DF)
# this takes over an hour and a half to run - todo: reimplement for efficiency
# (worthwhile given it's 200,000 rows to loop over)
for(i in 1:nrow(df)){
  subject <- df$id[i] #extract subject number for matching
  passage <- df$passage[i] #extract passage name for matching

  #production errors of interest
  # misprod = raw misproduction errors
  # hesitation = raw hesitations
  # words_with_misprod = distinct words with misproduction errors
  # words_with_hes = distinct words with pre-word or word-internal hesitation
  # misprod_rate = rate of raw misproduction errors
  # hesitation_rate = rate of raw hesitations
  # words_with_misprod_rate = rate of word-level misproduction errors
  # words_with_hes_rate = rate of word-level hesitations

  #extract passage characteristics from readDat
  df$lenSyll[i] <- sum(readDat$lengthSyll[which(readDat$passage==passage)]) #length of passage in syllables
  df$lenWord[i] <- sum(readDat$lengthWord[which(readDat$passage==passage)]) #length of passage in words
  df$avgSyllPerWord[i] <- df$lenSyll[i]/df$lenWord[i]

  #extract passage characteristics from freqDat
  df$avgWordFreq[i] <- freqDat$avgFreq[which(freqDat$passage==passage)]

  #extract participant accuracy from accDat
  df$challengeACC[i] <- accDat[match(passage, accDat$passage), as.character(subject)] #passage-specific challenge question accuracy for subject

  #extract participant demographics from demoDat
  df$sex[i] <- demoDat$sex[match(df$id[i], demoDat$record_id)]   #participant biological sex
  df$pronouns[i] <- demoDat$pron[match(df$id[i], demoDat$record_id)]   #participant preferred pronouns
  df$age[i] <- demoDat$age[match(df$id[i], demoDat$record_id)] #participant age
  df$ethnic[i] <- demoDat$ethnic[match(df$id[i], demoDat$record_id)] #participant ethnic group affiliation
  df$socclass[i] <- demoDat$socclass[match(df$id[i], demoDat$record_id)]   #participant social class identification
  df$eng[i] <- demoDat$eng[match(df$id[i], demoDat$record_id)] #participant multilingualism (0=EN only, 1=EN+another)
  df$langhis[i] <- demoDat$langhis[match(df$id[i], demoDat$record_id)] #participant language learning history (1=EN first, 2=other first, 3=EN+other same, 4=something else)
  df$ageen[i] <- demoDat$ageen[match(df$id[i], demoDat$record_id)] #participant age of English acquisition (if not L1)
  df$profen[i] <- demoDat$profen[match(df$id[i], demoDat$record_id)] #participant English proficiency (1=Native, 2=Advanced, 3=Intermediate, 4=Elementary, 5=Not proficient)
  df$commdis[i] <- demoDat$commdis[match(df$id[i], demoDat$record_id)] #participant communication disorder history (0=none, 1+=diagnoses to review)
  df$bfne[i] <- demoDat$bfne[match(df$id[i], demoDat$record_id)] #participant fear of negative evaluation
  df$phq8[i] <- demoDat$phq8[match(df$id[i], demoDat$record_id)] #participant depression
  df$scaaredTotal[i] <- demoDat$scaaredTotal[match(df$id[i], demoDat$record_id)] #participant overall anxiety
  df$scaaredGA[i] <- demoDat$scaaredGA[match(df$id[i], demoDat$record_id)] #participant general anxiety
  df$scaaredSoc[i] <- demoDat$scaaredSoc[match(df$id[i], demoDat$record_id)] #participant social phobias (scaared)
  df$sps[i] <- demoDat$sps[match(df$id[i], demoDat$record_id)] #participant social phobias (sias6sps6)
}
# succeeded given changes above it; to be revised

#organize participant demographic variables
df$sex <- as.factor(df$sex)
df$pronouns <- as.factor(df$pronouns)
df$ethnic <- as.factor(df$ethnic)
df$socclass <- as.factor(df$socclass)

# compute speed
# nb these are passage level data but are tracked by word
df$timePerSyllable <- df$readingTime / df$lenSyll
df$timePerWord     <- df$readingTime / df$lenWord


### SECTION 5: CROSS-CHECK ALL PARTICIPANTS MET INCLUSION CRITERIA
# note: given the time required to annotated errors, only participants who met
#       inclusion criteria were annotated
sum(df$eng==1 & df$langhis %in% c(2,4) & df$ageen>6) #confirm all subjects monolingual English OR natively bilingual OR learned English before age 6
sum(df$commdis>0) #confirm no subject diagnosed with any communication disorder
filter(df, profen > 3) %>% select(id) %>% unique %>% nrow #one remaining subject (150060) rates own English proficiency as not "elementary" or "not proficient", but reads fluidly and achieved 80% accuracy on challenge questions, so not excluded

#extract age and sex stats
# all these values are just in case they're useful - not needed per se for later
# steps of the logic in this script
summary_unique <- function(df, key, column, f = summary) {
  unique(select(df, column, key))[[column]] %>% f
}

c("age", "sex", "pronouns", "ethnic", "socclass") %>%
  map(\(col) summary_unique(df, "id", col)) # `map` not `for`: return, not print

summary_unique(df, "id", "age", f = sd) # also do stdev for age


### SECTION 6: TRIM PASSAGES DUE TO EXPERIMENTER ERROR
dfTrim <- subset(df, !(df$passage=='broccoli')) #remove broccoli passage due to typo in the last sentence as presented on-screen to participants
dfTrim <- subset(dfTrim, !(dfTrim$passage=='sun')) #remove sun passage due to error in coding Excel


### SECTION 7: OUTPUT DATAFRAME
out_filename <- paste(out_path, "readAloudBetaData-wordLevel_", today, ".csv", sep="", collapse=NULL)
write.csv(dfTrim, out_filename, row.names = FALSE)
out_filename

# collapse_by_participant <- function(filename_in, filename_out) {
#   by_participant <- read_csv(filename_in) %>%
#     unique %>% # dedup
#     group_by(id) %>% summarize(across(misprod:total_uncorrected_errors, sum)) # summarize by participant, across all passages
#     # TODO change the columns selected, once more have been added to the output of the preproc script
#
#   write_csv(by_participant, filename_out)
#   return(filename_out)
# }
#
#
# # base = "~/Documents/ndclab/analysis-sandbox/github-structure-mirror/readAloud-valence-dataset/derivatives/preprocessed"
# base = "/home/data/NDClab/datasets/readAloud-valence-dataset/derivatives/preprocessed"
# preprocessed_summary_filename = "TODO"
# collapsed_filename = "TODO"
#
# collapse_by_participant(preprocessed_summary_file, collapsed_filename)
