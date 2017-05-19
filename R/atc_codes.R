# atc_codes.R
# read in drugbank1.tsv and sort out atc codes.
# This is Daniel Himmelstein's file, it contains 7759 entries with 8 variables. You will need to change the filepath
# to where you have saved it on your computer.Its two years old so not a recent image of drugbank but so far the best we can do.
# https://github.com/dhimmel/drugbank/tree/gh-pages/data

library(dplyr)
library(tidyr)

DB <- file.path('C://R-files//drugbank//','drugbank1.tsv') %>% read.delim(na.strings='',sep='\t',header=T,comment.char="#",stringsAsFactors = FALSE)

# we can only use drugs with ATC codes
drugsATC <- DB %>%
  filter(!is.na(atc_codes))

no_drugs <- drugsATC %>% distinct(drugbank_id)
head(no_drugs)
paste("There are", nrow(no_drugs), "DRUGS with valid ATC codes in our database.")
dim(drugsATC)

# ENSURE YOU USE drugsATC structure and NOT DB structure!!!!!
# Create 1st, 2nd, 3rd & 4th level ATC codes
atc1 <- rep("", times = nrow(drugsATC)) # empty string vector
atc2 <- rep("", times = nrow(drugsATC)) # empty string vector
atc3 <- rep("", times = nrow(drugsATC)) # empty string vector
atc4 <- rep("", times = nrow(drugsATC)) # empty string vector
for(i in 1:nrow(drugsATC)){
  atc1[i] <- substr(drugsATC[i,5], start=1, stop=1)
  atc2[i] <- substr(drugsATC[i,5], start=1, stop=3)
  atc3[i] <- substr(drugsATC[i,5], start=1, stop=4)
  atc4[i] <- substr(drugsATC[i,5], start=1, stop=5)
}
drugsATC <- cbind(drugsATC,atc1,atc2,atc3,atc4)  # add to DB structure

no_atc1 = drugsATC %>% distinct(atc1)
no_atc2 = drugsATC %>% distinct(atc2)
no_atc3 = drugsATC %>% distinct(atc3)
no_atc4 = drugsATC %>% distinct(atc4)
paste("There are", nrow(no_atc1), "unique L1 ATC codes in our database.")
paste("There are", nrow(no_atc2), "unique L2 ATC codes in our database.")
paste("There are", nrow(no_atc3), "unique L3 ATC codes in our database.")
paste("There are", nrow(no_atc4), "unique L4 ATC codes in our database.")


