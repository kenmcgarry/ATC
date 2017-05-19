library(dplyr)
library(tidyr)
library(kohonen)

#---------------- WINE ----------------
data(wines)
training <- sample(nrow(wines), 120)
Xtraining <- scale(wines[training, ])
Xtest     <- scale(wines[-training, ], center = attr(Xtraining, "scaled:center"),scale = attr(Xtraining, "scaled:scale"))
trainingdata <- list(measurements = Xtraining, vintages = vintages[training])
testdata     <- list(measurements = Xtest,     vintages = vintages[-training])
mygrid = somgrid(5, 5, "hexagonal")
som.wines <- som(Xtraining, grid = mygrid)
som.prediction <- predict(som.wines, newdata = Xtest)
table(vintages[-training], som.prediction$predictions[["vintages"]])

#------------- ATC --------------------
smp_size <- floor(0.90 * nrow(drugDATA)) ## 90% of the sample size for training the model
set.seed(123)       ## set the seed to make the partition reproducible
rownames(drugDATA)<-NULL

train_ind <- sample(seq_len(nrow(drugDATA)), size = smp_size)
Xtraining <- scale(drugDATA[train_ind,2:13 ])  # we have "smp_size" amount of training data
Xtest     <-  scale(drugDATA[-train_ind,2:13 ] , center = attr(train, "scaled:center"),scale = attr(train, "scaled:scale")) 
trainingdata <- list(measurements=Xtraining,ATC=ATClabels[train_ind])
testdata     <- list(measurements = Xtest,ATC=ATClabels[-train_ind])

mygrid = somgrid(15, 15, "hexagonal")
som.atc <- supersom(Xtraining, grid = mygrid)
som.prediction <- predict(som.atc, newdata = testdata)

table(ATClabels[-trainingdata], som.prediction$predictions[["ATClabels"]])
table(testdata$ATC, som.prediction$predictions$ATC)
summary(som.atc)
#------------------------------------------



