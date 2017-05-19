
 x<-confusionMatrix(iris$Species, sample(iris$Species))

newPrior <- c(.05, .8, .15)
names(newPrior) <- levels(iris$Species)

confusionMatrix(iris$Species, sample(iris$Species))

expected <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 0))
predicted <- factor(c(1, 0, 0, 1, 0, 0, 1, 1, 1, 0))
results <- confusionMatrix(data=predicted, reference=expected)
print(results)

# try this.
colnames(predicted) <- sort(no_atc1$atc1) # replaces numbers with the ATC codes
colnames(actual) <- sort(no_atc1$atc1) # replaces numbers with the ATC codes
colnames(prednn$net.result) <- sort(no_atc1$atc1)
prednn$net.result
actual

input.matrix <- data.matrix(input)
input.matrix.normalized <- normalize(input.matrix)

colnames(input.matrix.normalized) = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")
rownames(input.matrix.normalized) = colnames(input.matrix.normalized)

confusion <- as.data.frame(as.table(input.matrix.normalized))
