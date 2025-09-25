# readAloud-valence-beta Analysis Preparation for word frequency
# Authors: Jessica M. Alexander and Luc Sahar -- NDCLab, Florida International University
# Last updated: 2023-08-25

# Input:
# - SUBTLEXus corpus
# - Scaffolds xlsx, used only as the simplest way to enumerate the passages
# - Stimulus characteristics for our passages

# Output:
# - a CSV, with one row per passage, showing the average word frequency (in
#   SUBTLEXus) for the words in that passage

### SECTION 1: SETTING UP
library(readxl) #load excel files
library(textstem) #lemmatize_words function

#set up date for output file naming
today <- Sys.Date()
today <- format(today, "%Y%m%d")

#set up directories for input/output
# main_dataset <- '/Users/jalexand/github/readAloud-valence-dataset/'
main_dataset <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-dataset/'
main_analysis <- '~/Documents/ndclab/rwe-analysis-sandbox/rwe-analysis/'
# out_path_readDat <- '/Users/jalexand/github/readAloud-valence-dataset/derivatives/'
out_path_readDat <- paste(main_analysis, 'derivatives/', sep = '')

#load input files
scaffolds <- paste(main_dataset, 'code/scaffolds.xlsx', sep="", collapse=NULL)
readAloudStimChar <- paste(main_dataset, 'materials/readAloud-ldt/stimuli/readAloud/readAloud-stimuli_characteristics.xlsx', sep="", collapse=NULL)
SUBList <- paste(main_dataset, 'materials/readAloud-ldt/stimuli/resources/SUBTLEXus74286wordstextversion.txt', sep="", collapse=NULL) #downloaded from https://www.ugent.be/pp/experimentele-psychologie/en/research/documents/subtlexus on 06/13/2022

# maybe include passDat? it has something about word data, and seems to have flesch too?
subtlexus <- read.table(SUBList, header=TRUE)
subtlexus$Word <- tolower(subtlexus$Word) #make all entries in SUBTLEXUS lower-case
passages <-excel_sheets(scaffolds) # antarctica ... vegas

# Jess' suggestions
# L72-124 from analysisStimuli.R

### SECTION 2: SET UP PASSAGE LIST AND PREPARE LEMMAS FOR ACCESSING FREQUENCY INFO
#create manual mapping of words to SUBTLEXUS corpus when lemma doesn't automatically match
manualLemma <- data.frame(matrix(ncol=2, nrow=0)) #initialize a table to hold words without a frequency match
colnames(manualLemma) <- c("stimWord",
                           "lemma")
for(i in 1:length(passages)){
  passage <- passages[i]
  passageDat <- read_xlsx(readAloudStimChar, sheet=passage, skip=1, na="#") #read in passage data
  passWords <- passageDat[,1:2] #pull word list
  for(a in 1:length(passWords$stimWord)){        #correct apostrophes (curly to straight)
    string <- passWords$stimWord[a]
    passWords$stimWord[a] <- gsub("’", "'",string)
  }
  passWords$stimWord <- tolower(passWords$stimWord) #shift word list to lowercase to match SUBTLEXUS
  passWords$lemma <- lemmatize_words(passWords$stimWord) #lemmatize word list
  passWords$freq <- rep(0, nrow(passWords)) #initialize new column to hold word frequency data
  for(f in 1:nrow(passWords)){        #add log word frequency from SUBTLEXUS corpus for each word lemma in the passage
    passWords$freq[f] <- subtlexus$Lg10WF[match(passWords$lemma[f], subtlexus$Word)]
  }
  noFreqTable <- subset(passWords, is.na(passWords$freq))[,2:3] #extract words that did not get a frequency match
  manualLemma <- rbind(manualLemma, noFreqTable) #bind to running table of words without a frequency match
}
manualLemma <- manualLemma[!duplicated(manualLemma$stimWord),] #remove duplicate rows
manualLemma <- manualLemma[order(manualLemma$stimWord),] #alphabetize
manualLemma$newLemma <- c(rep(0,nrow(manualLemma))) #initialize column to hold manual lemmas extracted by researcher
for(lemma in 1:nrow(manualLemma)){       #manually adjust possessives by dropping "'s" (except for "it's")
  string <- manualLemma$lemma[lemma]
  if(substr(string, nchar(string)-1, nchar(string))=="'s" & string!="it's"){
    manualLemma$newLemma[lemma] <- substr(string, 0, nchar(string)-2)
  }
}
#manually adjust plurals
manualLemma[match("brittles", manualLemma$stimWord),3] <- "brittle"
#manually adjust adjectives
manualLemma[match("club-like", manualLemma$stimWord),3] <- "club"
manualLemma[match("in-flight", manualLemma$stimWord),3] <- "flight"
manualLemma[match("mid-", manualLemma$stimWord),3] <- "middle"
#no manual adjustment of the following categories:
#compound words with no obvious "primary" lemma: ccc, don't, it's, long-term
#proper nouns: delano, nissan
#words simply not in the SUBTLEXUS database: ecotourism, hydropower, jetsetter, megabat, microbats, photoreceptor, plumicorn, powertrain, spearer, trinocular
#ordinal numbers: nineteeth, second, twentieth
manualLemma <- subset(manualLemma, manualLemma$newLemma!=0) #drop words without a manual mapping
manualLemma$freq <- rep(0, nrow(manualLemma))       #add log word frequency from SUBTLEXUS corpus
for(f in 1:nrow(manualLemma)){
  manualLemma$freq[f] <- subtlexus$Lg10WF[match(manualLemma$newLemma[f], subtlexus$Word)]
}

# Jess' suggestions
# L170-189 from analysisStimuli.R

### SECTION 3: BUILD AND OUTPUT READDAT MATRIX
readDat <- data.frame(matrix(ncol=2, nrow=0))
colnames(readDat) <- c("passage", "avgFreq")

#calculate characteristics per passage half
for(j in 1:length(passages)){
  passage <- passages[j]
  passageDat <- read_xlsx(readAloudStimChar, sheet=passage, skip=1, na="#") #read in passage word list

  #extract passage word list
  passWords <- passageDat[,1:2] #pull word list
  for(a in 1:length(passWords$stimWord)){        #correct apostrophes (curly to straight)
    string <- passWords$stimWord[a]
    passWords$stimWord[a] <- gsub("’", "'",string)
  }
  passWords$stimWord <- tolower(passWords$stimWord) #shift word list to lowercase to match SUBTLEXUS
  passWords$lemma <- lemmatize_words(passWords$stimWord) #lemmatize word list

  #add frequency data
  passWords$freq <- rep(0, nrow(passWords)) #initialize column to hold frequency data
  for(f in 1:nrow(passWords)){
    if(!is.na(subtlexus$Lg10WF[match(passWords$stimWord[f], subtlexus$Word)])){           #automatic matching of full word to SUBTLEXUS corpus
      passWords$freq[f] <- subtlexus$Lg10WF[match(passWords$stimWord[f], subtlexus$Word)]
    } else if(!is.na(subtlexus$Lg10WF[match(passWords$lemma[f], subtlexus$Word)])){       #automatic matching of lemma to SUBTLEXUS corpus
      passWords$freq[f] <- subtlexus$Lg10WF[match(passWords$lemma[f], subtlexus$Word)]
    } else if (passWords$lemma[f] %in% manualLemma$lemma){                                #check manual table if no auto-match
      passWords$freq[f] <- manualLemma$freq[match(passWords$lemma[f], manualLemma$lemma)]
    } else {passWords$freq[f] <- median(subtlexus$Lg10WF)}                                #impute median value for non-matches
  }


  #calculate average frequency
  freqAvg <- mean(passWords$freq) #calculate mean frequency for whole passage

  #create vectors for each passage and add to readDat
  vectorWholePassage <- c(passage, freqAvg) # instead of vectorPre and vectorPost, as in analysisStimuli.R
  readDat[nrow(readDat) + 1,] <- c(vectorWholePassage)
}

#organize data types
readDat$avgFreq <- as.numeric(readDat$avgFreq)

#output readDat
out_path_readDat_csv <- paste(out_path_readDat, 'prepWordFreq_readDat', today, '.csv', sep="", collapse=NULL)
write.csv(readDat, out_path_readDat_csv)
out_path_readDat_csv # print out to terminal so it's easy to open!
