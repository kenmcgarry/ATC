# ATCpredict1.R
# Ken McGarry, Feb 2017
# You need to go to the top menu of RStudio to let R know where you store your data.
# Session -> Set Working Directory -> Choose Directory

library(dplyr)
library(tidyr)
library(ChemmineR)
library(caret)
library(nnet)
library(neuralnet)
library(KEGGREST)

# speed up with multithreading cores - ONLY use if you have MICROSOFT R Open installed  and NOT the usual R!!!!
# https://mran.microsoft.com/open/
#mycores <- getMKLthreads()
#setMKLthreads(mycores)

# This is Daniel Himmelstein's file, it contains 7759 entries with 8 variables. You will need to change the filepath
# to where you have saved it on your computer.Its two years old so not a recent image of drugbank but so far the best we can do.
# https://github.com/dhimmel/drugbank/tree/gh-pages/data

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


# This sdf file contains the drug structures, 1,885 drugs in database
# 23 structures are invalid, these are removed.The SDF stands for "Structure-Data File" format
#DBsdf <- read.SDFset("C:\\R-files\\ATC\\DBstructures.sdf") # from Daniel Himmelstein
DBsdf <- read.SDFset("C:\\R-files\\ATC\\all_structures.sdf") # from Drugbank 1/2/2017

valid <- validSDF(DBsdf); 
DBsdf <- DBsdf[valid]
length(DBsdf)

blockmatrixDB <- datablock2ma(datablocklist=datablock(DBsdf)) # Converts data block to matrix 
## Assign drugbank IDs from datablock to cid slot

cid(DBsdf) <- as.character(blockmatrixDB[,"DRUGBANK_ID"])

numchar <- splitNumChar(blockmatrix=blockmatrixDB) # Splits to numeric and character matrix 
data1 <- numchar[2] # numchar[1] contains numbers, numchar[2] contains character version
data1 <- as.data.frame(data1,stringsasfactors=FALSE)
names(data1) <- gsub("charMA.", "", names(data1)) # gets rid of "charMA." which is prefixed to column names
dim(data1)

names(data1)[names(data1)=="DRUGBANK_ID"] <- "drugbank_id"  # must rename for dplyr inner_join to recognise it

# Now need to match up SDF drugnames (data1) with drugbank ID's for class label dataframe (drugsATC)
ID <- drugsATC %>% distinct(drugbank_id) # get list unique Id's
chemstructure <- filter(data1, drugbank_id %in% ID[1:nrow(ID),]) # see if they exist in data1

# Merge the chemical data in data1 with those drugs in DB (for atc names)
drugsATC <- inner_join(drugsATC, chemstructure, by = "drugbank_id")

# ==== Create the chemical FINGERPRINTS =====
apset <- sdf2ap(DBsdf) # fist create atom pairs
#fpset <- desc2fp(apset, descnames=1024, type="FPset")  
fpset <- desc2fp(apset, descnames=256, type="FPset") # convert atom-pairs into fingerprints
fpdrugs <- fpset[drugsATC$drugbank_id,]  # only keep those drug ATC and chem structures we actually have

fpma <- as.matrix(fpdrugs)  # convert fingerprints to matrix

#rownames(fpma) <- drugsATC$drugbank_id  # keep drugbank names
drugbank_id <- drugsATC$drugbank_id 
atc1 <- drugsATC$atc1
atc2 <- drugsATC$atc2
atc3 <- drugsATC$atc3
atc4 <- drugsATC$atc4

atcLabels1 <- class.ind(atc1) # great little function for creating binary 1 of n class labels!!
fpdataL1 <- data.frame(drugbank_id, fpma, atcLabels1)

# Start building classifiers using fpdata, we have four levels of ATC codes matched up with data
# for 1,334 drugs with chemical data (fingerprints) and ATC codes. ATC codes are the class labels.
# library(randomForest)

n <- nrow(fpdataL1)
index <- sample(1:n, size = round(0.80*n), replace=FALSE)

#L1net <- nnet(fpma[index,],atcLabels1[index,],size=12,rang=0.1,decay=0.1,maxit=2000,MaxNWts=9000000)
#predicted <- (round(as.matrix(predict(L1net, fpma[-index,]))))

actual <- (as.matrix(atcLabels1[-index,]))
#rownames(actual) <- sort(no_atc1$atc1) # ensure both have same drugbank id's
#colnames(predicted) <- sort(no_atc1$atc1) # replaces numbers with the ATC codes
colnames(actual) <- sort(no_atc1$atc1) # replaces numbers with the ATC codes

# NEURAL NETWORK TRAINING - MLP (MULTI-LAYER PERCEPTRONS)
# appending the labels to the training data for the neuralnet package requirements
output <- atcLabels1[index,]
testout <- atcLabels1[-index,]
colnames(testout)<-colnames(atcLabels1)
colnames(output)<-colnames(atcLabels1)
output.names<-colnames(atcLabels1)
# Must rename the input variables as the long names cause neuralnet problems i.e.
# Error in terms.formula(formula) : invalid model formula in ExtractVars
j<-ncol(fpma)
colnames(fpma) <- paste("v", 1:j, sep = "")
input.names<-colnames(fpma)

trainingData <- as.data.frame(fpma[index,])
trainingData <- as.data.frame(cbind(output,trainingData))
testData <- as.data.frame(fpma[-index,]) # Note: 

formula1 <- as.formula(paste(paste(output.names,collapse=" + ")," ~",paste(input.names,collapse=" + ")))
L1net <- neuralnet(formula1,data=trainingData,hidden=c(25,25),lifesign="full",linear.output=F,algorithm='rprop+',learningrate=0.1,rep=1)

prednn <- neuralnet::compute(L1net,testData)
cm <- test.cl(actual,round(prednn$net.result),colnames(atcLabels1))

# Define some basic statistics that will be needed to compute the evaluation metrics.See website
# http://blog.revolutionanalytics.com/2016/03/com_class_eval_metrics_r.html
showClass("FPset") 
ni <- sum(cm) # number of instances
nc <- nrow(cm) # number of classes
diag <- diag(cm) # number of correctly classified instances per class 
rowsums <- apply(cm, 1, sum) # number of instances per class
colsums <- apply(cm, 2, sum) # number of predictions per class
p <- rowsums / ni # distribution of instances over the actual classes
q <- colsums / ni # distribution of instances over the predicted classes
classaccuracy <- round((diag/rowsums)*100,digits=2) # accuracy for each class %
s <- matrix(0, nrow = 2, ncol = 2)
avgAccuracy <- sum(diag(s)) / sum(s)
accuracy <- sum(diag) / ni 

accuracy
classaccuracy

precision <- diag / colsums 
recall <- diag / rowsums 
f1 <- 2 * precision * recall / (precision + recall) #f1 is the harmonic mean (or a weighted average) of precision and recall.

data.frame(precision, recall, f1) 

test.cl <- function(Actual, NNpred, thenames) {
  colnames(Actual) <- thenames
  actual <- max.col(Actual)
  predicted <- max.col(NNpred)
  x<-table(actual, predicted)
  colnames(x) <- thenames
  rownames(x) <-thenames
  return(x)
}


oneVsAll <- lapply(1 : nc,
                  function(i){
                    v = c(cm[i,i],
                          rowsums[i] - cm[i,i],
                          colsums[i] - cm[i,i],
                          n-rowsums[i] - colsums[i] + cm[i,i]);
                    return(matrix(v, nrow = 2, byrow = T))})

s = matrix(0, nrow = 2, ncol = 2)
for(i in 1 : nc){s = s + oneVsAll[[i]]}

avgAccuracy <- sum(diag(s)) / sum(s)

# save the chemcical and drug data ATC code for later.
save(drugsATC, file = "C:\\R-files\\ATC\\ATC.RData")



