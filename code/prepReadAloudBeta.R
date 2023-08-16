# readAloud-valence-beta Analysis Preparation
# Authors: Luc Sahar and Jessica M. Alexander -- NDCLab, Florida International University
# Last updated: 2023-08-15

# INPUTS
# data/df: behavioral (error-related) data, for each participant on each passage
# accDat: comprehension question accuracy (0/1) for each participant for each passage
# readDat: stimuli characteristics (by passage half)
# redcap: participant data, incl. demographics and responses + scored factors for questionnaires:
  # bfne (brief fear of negative evaluation): bfne_b_scrdTotal (fear of negative evaluation total)
  # phq8 (patient health questionnaire): phq8_scrdTotal (depression scale total)
  # scaared, total (screen for adult anxiety disorders): scaared_b_scrdTotal (total anxiety)
  # scaared, social (screen for adult anxiety disorders): scaared_b_scrdSoc (social phobias)
  # scaared, general (screen for adult anxiety disorders): scaared_b_scrdGA (general anxiety)
  # sps (social phobia scale): sias6sps6_b_scrdSPS

# OUTPUTS
# dfTrim: for each passage, for each participant, details on:
  # participant behavior: reading errors made, comprehension question accuracy
  # passage characteristics: length (syllable and word), average syllables per word
  # participant data: demographics, language history, mood and mood disorder scores


### SECTION 1: SETTING UP
library(readxl) # read_xlsx
library(stringr) # str_extract
library(dplyr) # most things
library(purrr) # map, map_df; generally good to have
library(lubridate) # now
library(readr) # write_csv

#set up defaults for output file naming
# ext_default = 'csv'
# tz_default = "America/New_York"
# date_format_default = "%Y%m%d_%I%M%P"
#
# build_output_filename <- function(label, ext = ext_default, timezone = tz_default, date_format = date_format_default) {
#   # `label` may include the destination directory, if different from the working directory when the script is run
#   current_datetime <- now(timezone) %>% format(date_format)
#   paste(label, '_', current_datetime, '.', ext, sep = "")
# }
today <- Sys.Date()
today <- format(today, "%Y%m%d")

#set up directories for input/output
# main_dataset <- '/Users/jalexand/github/readAloud-valence-dataset/'
# main_analyses <- '/Users/jalexand/github/readAloud-valence-beta/'
# out_path <- '/Users/jalexand/github/readAloud-valence-beta/derivatives/'
main_dataset <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-dataset/'
main_analyses <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/'
out_path <- '/home/luc/Documents/ndclab/analysis-sandbox/rwe-analysis/derivatives/'


#load input files
# data <- paste(main_dataset, 'derivatives/preprocessed/disfluencies_subject-x-passage_20230616_1229pm.csv', sep="", collapse=NULL)
data <- "/home/luc/Documents/ndclab/analysis-sandbox/output-csvs/disfluencies_subject-x-passage_20230616_1229pm.csv"
accDat_path <- paste(main_dataset,'derivatives/preprocessed/readAloud_passage-level_summary_20220812.csv', sep="", collapse=NULL)
readDat_path <- paste(main_dataset, 'derivatives/analysisStimuli_readDat_20230614.csv', sep="", collapse=NULL)
redcap_path <- paste(main_dataset,'derivatives/preprocessed/202201v0readAloudval_SCRD_2022-06-20_1019.csv', sep="", collapse=NULL)
agedat_path <- paste(main_dataset,'derivatives/preprocessed/202201v0readAloudval_SCRD_2022-06-20_1019_ageonly.csv', sep="", collapse=NULL)
speedDat_path <- paste(main_dataset, "derivatives/preprocessed/valence-timing/timingpitch_subject-by-passage_2022-09-09.csv", sep="", collapse=NULL)

# c(data, accDat_path, readDat_path, redcap_path, agedat_path, speedDat_path) %>% fs::as_fs_path() %>% fs::is_file()
# ✅: all TRUE


df <- read.csv(data, na.strings='NA')
redcap <- read.csv(redcap_path, na.strings='NA') #participant questionnaire responses
agedat <- read.csv(agedat_path, na.strings='NA') #participant age information
readDat <- read.csv(readDat_path, na.strings='N') #passage-level characteristics from analysisStimuli.R
accDat <- read.csv(accDat_path, na.strings='NA', check.names=FALSE) #passage level accuracy for each subject
accDat$passage <- c("dams", "flying", "bats", "broccoli", "realty", "bees", "dogshow", "dolphins", "icefishing",
                    "cars", "vegas", "sun", "caramel", "congo", "antarctica", "depression", "skunkowl", "grizzly",
                    "mantis", "dentist")        #rename passages with short-name
speedDat <- read.csv(speedDat_path, na.strings='NA')

#organize data types
df[,3:30] <- sapply(df[,3:30],as.numeric)

#add missing passages for 150086 so that nrow is divisible by 20
passages_read <- df$passage[which(df$id=="150086")]
all_passages <- unique(df$passage)
tempdf <- data.frame(matrix(nrow=0, ncol=ncol(df)))
colnames(tempdf) <- colnames(df)
for(passage in 1:length(all_passages)){
  if(all_passages[passage] %in% passages_read){next}else{
    tempdf[nrow(tempdf) + 1,] <- c("150086", all_passages[passage], rep(NA, 28))
  }
}
df <- rbind(df, tempdf)

### SECTION 2: BUILD DEMOGRAPHIC DATA DF
demoDat <- redcap[,c(1,5)]
#biological sex: replace numerical values with text description
for(a in 1:nrow(redcap)){
  if(is.na(redcap$demo_b_sex_s1_r1_e1[a])){demoDat$sex[a] <- 'undisclosed'}
  else if(redcap$demo_b_sex_s1_r1_e1[a]==1){demoDat$sex[a] <- 'male'}
  else if(redcap$demo_b_sex_s1_r1_e1[a]==2){demoDat$sex[a] <- 'female'}
  else if(redcap$demo_b_sex_s1_r1_e1[a]==3){demoDat$sex[a] <- 'intersex'}
  else if(redcap$demo_b_sex_s1_r1_e1[a]==4){demoDat$sex[a] <- 'other'}
  else if(redcap$demo_b_sex_s1_r1_e1[a]==5){demoDat$sex[a] <- 'unknown'}
  else{demoDat$sex[a] <- 'undisclosed'}
}

#preferred pronouns: replace numerical values with text description
for(b in 1:nrow(redcap)){
  if(is.na(redcap$demo_b_pronouns_s1_r1_e1[b])){demoDat$pron[b] <- 'undisclosed'}
  else if(redcap$demo_b_pronouns_s1_r1_e1[b]==1){demoDat$pron[b] <- 'she/her'}
  else if(redcap$demo_b_pronouns_s1_r1_e1[b]==2){demoDat$pron[b] <- 'he/him'}
  else if(redcap$demo_b_pronouns_s1_r1_e1[b]==3){demoDat$pron[b] <- 'they/them'}
  else if(redcap$demo_b_pronouns_s1_r1_e1[b]==5){demoDat$pron[b] <- 'other'}
  else{demoDat$pron[b] <- 'undisclosed'}
}

#ethnicity affiliation: map to text description
for(c in 1:nrow(redcap)){
  if(redcap$demo_b_ethnic_s1_r1_e1___1[c]==1){demoDat$ethnic[c] <- 'AI'} #american indian/alaska native
  else if(redcap$demo_b_ethnic_s1_r1_e1___2[c]==1){demoDat$ethnic[c] <- 'A'} #asian
  else if(redcap$demo_b_ethnic_s1_r1_e1___3[c]==1){demoDat$ethnic[c] <- 'AA'} #african american
  else if(redcap$demo_b_ethnic_s1_r1_e1___4[c]==1){demoDat$ethnic[c] <- 'LX'} #hispanic/latinx
  else if(redcap$demo_b_ethnic_s1_r1_e1___5[c]==1){demoDat$ethnic[c] <- 'ME'} #middle eastern
  else if(redcap$demo_b_ethnic_s1_r1_e1___6[c]==1){demoDat$ethnic[c] <- 'PI'} #pacific islander
  else if(redcap$demo_b_ethnic_s1_r1_e1___7[c]==1){demoDat$ethnic[c] <- 'W'} #white
  else if(redcap$demo_b_ethnic_s1_r1_e1___8[c]==1){demoDat$ethnic[c] <- 'O'} #other
  else{demoDat$ethnic[c] <- 'UND'} #undisclosed
}

#social class affiliation: replace numerical values with text description
for(d in 1:nrow(redcap)){
  if(is.na(redcap$demo_b_socclass_s1_r1_e1[d])){demoDat$socclass[d] <- 'undisclosed'}
  else if(redcap$demo_b_socclass_s1_r1_e1[d]==1){demoDat$socclass[d] <- 'poor'}
  else if(redcap$demo_b_socclass_s1_r1_e1[d]==2){demoDat$socclass[d] <- 'working'}
  else if(redcap$demo_b_socclass_s1_r1_e1[d]==3){demoDat$socclass[d] <- 'middle'}
  else if(redcap$demo_b_socclass_s1_r1_e1[d]==4){demoDat$socclass[d] <- 'affluent'}
  else{demoDat$socclass[d] <- 'undisclosed'}
}

#communication disorders diagnoses: sum across childhood, adolescence, and adulthood
for(e in 1:nrow(redcap)){
  demoDat$commdis[e] <- sum(redcap$demo_b_comdiskid_s1_r1_e1[e],
                            redcap$demo_b_comdisteen_s1_r1_e1[e],
                            redcap$demo_b_comdisad_s1_r1_e[e])
}

#language history: transfer directly
for(f in 1:nrow(redcap)){
  demoDat$eng[f] <- redcap$demo_b_eng_s1_r1_e1[match(demoDat$record_id[f], redcap$record_id)] #participant monolingualism
  demoDat$langhis[f] <- redcap$demo_b_langhis_s1_r1_e1[match(demoDat$record_id[f], redcap$record_id)] #participant language history
  demoDat$ageen[f] <- redcap$demo_b_ageen_s1_r1_e1[match(demoDat$record_id[f], redcap$record_id)] #participant age of English acquisition
  demoDat$profen[f] <- redcap$demo_b_profen_s1_r1_e1[match(demoDat$record_id[f], redcap$record_id)] #participant English proficiency
}

#mood and mood disorders: transfer directly
for(g in 1:nrow(redcap)){
  demoDat$bfne[g] <- redcap$bfne_b_scrdTotal[match(demoDat$record_id[g], redcap$record_id)] #bfne total score
  demoDat$phq8[g] <- redcap$phq8_scrdTotal[match(demoDat$record_id[g], redcap$record_id)] #phq8 depression scale
  demoDat$scaaredTotal[g] <- redcap$scaared_b_scrdTotal[match(demoDat$record_id[g], redcap$record_id)] #scaared total anxiety
  demoDat$scaaredGA[g] <- redcap$scaared_b_scrdGA[match(demoDat$record_id[g], redcap$record_id)] #scaared general anxiety
  demoDat$scaaredSoc[g] <- redcap$scaared_b_scrdSoc[match(demoDat$record_id[g], redcap$record_id)] #scaared social phobias
  demoDat$sps[g] <- redcap$sias6sps6_b_scrdSPS[match(demoDat$record_id[g], redcap$record_id)] #sps social phobia scale
}

#age: pull from separate file
for(h in 1:nrow(demoDat)){
  demoDat$age[h] <- agedat$info_age_s1_r1_e1[match(demoDat$record_id[h], agedat$record_id)]
}


### SECTION 3: SET UP DERIVED FIELDS FOR SPEED ANALYSES
speedDat$readingTime <- speedDat$readEnd - speedDat$readStart
speedDat$id <- as.character(speedDat$id) # so we can join and it doesn't complain about type comparison
df <- left_join(df, speedDat, by = c("id", "passage")) # now reading timestamps and duration are looped into df


### SECTION 4: BUILD TRIAL-LEVEL DF (ADD DEMODAT, READDAT, and ACCDAT to DF)
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

#organize participant demographic variables
df$sex <- as.factor(df$sex)
df$pronouns <- as.factor(df$pronouns)
df$ethnic <- as.factor(df$ethnic)
df$socclass <- as.factor(df$socclass)

# compute speed
df$timePerSyllable <- df$readingTime / df$lenSyll
df$timePerWord     <- df$readingTime / df$lenWord


### SECTION 5: CROSS-CHECK ALL PARTICIPANTS MET INCLUSION CRITERIA
#note: given the time required to annotated errors, only participants who met inclusion criteria were annotated
#sum(df$eng==1 & df$langhis %in% c(2,4) & df$ageen>6) #confirm all subjects monolingual English OR natively bilingual OR learned English before age 6
#sum(df$commdis>0) #confirm no subject diagnosed with any communication disorder
#sum(df$profen>3, na.rm=TRUE)/20 #one remaining subject (150060) rates own English proficiency as not "elementary" or "not proficient", but reads fluidly and achieved 80% accuracy on challenge questions, so not excluded

#extract age and sex stats

# all these values are just in case they're useful - not needed per se for later
# steps of the logic in this script

summary(df$age) #age range and mean
sd(df$age) #age standard deviation
summary(df$sex)/20 #number of participants by sex
summary(df$sex)/20 / (nrow(df)/20) #percentage of participants by sex
summary(df$pronouns)/20 #number of participants by preferred pronoun
summary(df$pronouns)/20 / (nrow(df)/20) #percentage of participants by preferred pronoun
summary(df$ethnic)/20 #number of participants by ethnic affiliation
summary(df$ethnic)/20 / (nrow(df)/20) #percentage of participants by ethnic affiliation
summary(df$socclass)/20 #number of participants by social class affiliation
summary(df$socclass)/20 / (nrow(df)/20) #percentage of participants by social class affiliation


### SECTION 6: TRIM PASSAGES DUE TO EXPERIMENTER ERROR
dfTrim <- subset(df, !(df$passage=='broccoli')) #remove broccoli passage due to typo in the last sentence as presented on-screen to participants
dfTrim <- subset(dfTrim, !(dfTrim$passage=='sun')) #remove sun passage due to error in coding Excel


### SECTION 7: OUTPUT DATAFRAME
write.csv(dfTrim, paste(out_path, "readAloudBetaData_", today, ".csv", sep="", collapse=NULL), row.names = FALSE)

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
# collapsed_filename = build_output_filename(label = paste(base, "disfluencies_subject", sep='/'))
#
# collapse_by_participant(preprocessed_summary_file, collapsed_filename)
