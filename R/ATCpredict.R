# ATCpredict.R
#
# You need to go to the top menu of RStudio to let R know where you store your data.
# Session -> Set Working Directory -> Choose Directory

library(dplyr)
library(tidyr)
library(ChemmineR)
library(drugbankR)
library(RSQLite)
library(XML)


setwd("C:/R-files/ATC")
ATCdata <- read.table('ATCfile1.csv', header=TRUE, sep=',')
names(ATCdata)
length(unique(ATCdata$GenericName))
length(unique(ATCdata$ATCCode1))
head(ATCdata)

# downloaded from drugbank, no ATC codes though! :(
setwd("C:/R-files/ATC")
ATClinks <- read.table('drug_links.csv', header=TRUE, sep=',')

# This sdf file contains all the drug structures for each drug in database
# 23 structures are invalid.
DBsdf <- read.SDFset("C:\\R-files\\ATC\\DBstructures.sdf")
valid <- validSDF(DBsdf); 
DBsdf <- DBsdf[valid]
length(DBsdf)

blockmatrixDB <- datablock2ma(datablocklist=datablock(DBsdf)) # Converts data block to matrix 
numchar <- splitNumChar(blockmatrix=blockmatrixDB) # Splits to numeric and character matrix 
data1 <- numchar[2]
data1 <- as.data.frame(data1)
names(data1) <- gsub("numMA.", "", names(data1)) # gets rid of "numMA." which is prefixed to column names
dim(data1)


# This is daniel himmelstein's file
DB <- read.table('C:\\R-files\\drugbank\\drugbank.tsv', header=TRUE, sep='\t')

## I downloaded the drugbank.xml file and then
## convert drugbank database (xml file) into dataframe:
## below code takes 20 minutes on my machine so only do it the once!! So commented out.

## drugbank_df <- dbxml2df(xmlfile="drugbank.xml")





