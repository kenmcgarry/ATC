# som_drugs_atc.R
# Using the kohonen self-organising feature map for prediction of drug
# anatomical, theraputic & chemical codes (ATC). The intention is to
# better understand the pharmaco-theraputic nature of drugs for:
# 1. candidates drugs for repositioning; 2. predicting side-effects 3. protein-targets
# 27/2/2017
#     source("http://bioconductor.org/biocLite.R")
#     biocLite("tidyr")

library(caret)
#library(kernlab)
library(dplyr)
library(tidyr)
library(kohonen)
library(RColorBrewer)

# This code assumes that the drugsATC data structure has already been created using
# ATCpredict1.R, so load it in.

load("C:\\R-files\\ATC\\ATC.RData")

drugDATA <- drugsATC[,c(10,22,23,25,29,30,31,33:34,37:40)]  # take only the columns we need 
rownames(drugDATA)<- drugsATC[,2]   # keep the drug names for reference thru rownames.

# unfortunately drugDATA is all of type factors, use this little function to convert to integers!
drugDATA[] <- lapply(drugDATA[,-1], function(x)
  as.integer(levels(x))[x])

drugDATA[,1] <- drugsATC[,10]  # get the factors for classlabels back (ATC1)
drugDATA <- na.omit(drugDATA) # remove rows with NA's, if present as we cannot train SOM

ATClabels <- drugDATA[,1]
drugDATA <- scale(drugDATA[,2:ncol(drugDATA)])
drugDATA <- cbind.data.frame(ATClabels,drugDATA)  # get the factors for ATC labels back. 

#------------- ATC --------------------
smp_size <- floor(0.90 * nrow(drugDATA)) ## 90% of the sample size for training the model
set.seed(1234)       ## set the seed to make the partition reproducible
rownames(drugDATA)<-NULL

train_ind <- sample(seq_len(nrow(drugDATA)), size = smp_size)
Xtraining <- scale(drugDATA[train_ind,2:13 ])  # we have "smp_size" amount of training data
Xtest     <-  scale(drugDATA[-train_ind,2:13 ] , center = attr(Xtraining, "scaled:center"),scale = attr(Xtraining, "scaled:scale")) 
trainingdata <- list(measurements=Xtraining,ATC=ATClabels[train_ind])
testdata     <- list(measurements = Xtest,                ATC=ATClabels[-train_ind])
testdata     <- list(measurements = Xtest)

mygrid = somgrid(15, 15, "hexagonal")
som.atc <- supersom(trainingdata, grid = mygrid)
som.prediction <- predict(som.atc, newdata = testdata)

table(ATClabels[-train_ind], som.prediction$predictions$ATC)
summary(som.atc)
#------------------------------------------


plot(som.atc, type="changes")
plot(som.atc, type="count")
plot(som.atc, type="dist.neighbours")
plot(som.atc, type="codes", main = c("Codes X", "Codes Y"))
plot(som.atc, type="counts")

## palette suggested by Leo Lopes
coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}

plot(som.atc, type="quality", palette.name = coolBlueHotRed)
#plot(som.atc, type="mapping", labels = Xtraining[,1], col = Xtraining[,1]+1, main = "ATC codes mapping plot")

plot(som.atc, type="mapping",labels = (ATClabels), main = "ATC codes mapping plot")

## NEEDS TO BE DEBUGGED AGAIN FROM THIS POINT ONWARDS.
## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(som.atc$codes$measurements)), 10)
# plot these results:
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')
plot(som.atc, type="mapping", bgcol = pretty_palette[som_cluster], main = "Drug ATC codes") 
add.cluster.boundaries(som.atc, som_cluster)

Umat <- plot(som.atc, type="dist.neighbours", main = "SOM neighbour distances")
## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(som.atc, "codes")), 5)
add.cluster.boundaries(som.atc, som.hc)




# SUPERVISED SOM MAPPINGS
som.xyf <- xyf(data = as.matrix(Xtraining[,-1]),Y = factor(Xtraining[,1]),grid = somgrid(15, 15, "hexagonal"))
par(mfrow = c(1, 2))
plot(som.xyf, type = "counts", main = "Supervised SOM: counts")
plot(som.xyf, type = "quality",main = "Supervised SOM: mapping quality")
plot(som.xyf, type="codes", main = c("Codes X", "Codes Y"))

som.spred <- predict(som.xyf,   newdata = as.matrix(test[,-1]), trainX = as.matrix(train[,-1]), trainY = factor(train[,1]))
table(test[,1], som.prediction$prediction)







