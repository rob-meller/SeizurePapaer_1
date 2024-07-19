#### Prediction Modeling:  Seizure project
## Edited by RMELLER
## REFS https://topepo.github.io/caret/index.html
## REFS #https://www.machinelearningplus.com/machine-learning/caret-package/

# Intital steps 
setwd("/path/to/project_folder")

getwd()
# Data are in a folder called "DATA"
# Create a folder called "PLOTS" for results

library(readr)
library(corrplot)
library(caret)
library(glmnet)
library(stringr)
library(plyr)
library(dplyr)
library(skimr)
library(caretEnsemble)
library(reshape2)
library(pROC)
library(caretSDM)
library(rminer)
library(MLeval)

##Load Datasheet (TMM gene expression values calculated in Partek Genomics Studio)
df_0 <- read.table("DATA/SEIZURE_GENE_TMM.orig.txt", sep = "\t", header = TRUE); typeofdata="Gene tmm"

#Remove S0994 unknown etiology
df1 <- df_0[!df_0$Sample=="S994",]
df1[1:10,1:20];dim(df1)# 89 samples

## remove discharge
df2 = df1[!df1$Time=="Disc",]
df2[1:10,1:20];dim(df2)# 62 samples

# Create test and training data frames 
## SEIZ samples
df3 = df2[(df2$"EEG.change"=="YES"),]
df3[1:10,1:20];dim(df3)# 54 samples

## PNES Samples
df4 = df2[(df2$"EEG.change"=="NO"),]
df4[1:8,1:20];dim(df4)# 8 samples

# To address error message, convert df to data frame object 
train_df <- as.data.frame(df3)
test_df <- data.frame(df4[,])
train_df <- train_df[,c(5,19:8025)]
train_df[1:10,1:20]; dim(train_df)

##Pull in data sheet of sig genes..
genes=read.csv("DATA/60 SEIZ4-6hvsAll_genelist.csv")
genes[1:20,1:5]; dim(genes) # 60x26 # Genes are in colum Gene.Symbol

## Now subset gene expression values based on list of sig genes.
train_df_0  <- (train_df[,names(train_df) %in% genes$"Gene.Symbol"])
head(train_df_0);dim(train_df_0)
# now add back Time info
dfz=cbind(train_df$Time, train_df_0[,])
# rename col1
colnames(dfz)[1] <- "Time"
head(dfz); dim(dfz)
## rename back to train_df
train_df = dfz

## Data needs to be a binomial/ classifier model!!! recode data to  Y/ N
train_df$Time <- recode(train_df$Time, 'Base' = 'N', '4-6h'='Y', 
	.default=NA_character_)
head(train_df); dim(train_df)

# better recode the test data
test_df$Time <- recode(test_df$Time, 'Base' = 'N', '4-6h'='Y', 
	.default=NA_character_)
head(test_df); dim(test_df)

# quick assessment..
skim(train_df)

#############################################################################
#                        PREDICTION MODELING						    #
#############################################################################

##Preprocessing
preProcess_range_model <- preProcess(train_df, method=c('scale', 'center'))
train_df <- predict(preProcess_range_model, newdata = train_df)

# process the test data too
test_df <- predict(preProcess_range_model, newdata = test_df)

# Review expression of DEGs
featurePlot(x = train_df[, 2:60], 
            y = as.factor(train_df$Time), 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

featurePlot(x = train_df[, 2:60], 
            y = as.factor(train_df$Time), 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))


## Let's try a series of models...

# Define the training control
fitControl <- trainControl(
    method = 'cv',                   # k-fold cross validation
    number = 10,                      # number of folds
    savePredictions = 'final',       # saves predictions for optimal tuning parameter
    classProbs = TRUE,                  # should class probabilities be returned
    summaryFunction=twoClassSummary  # results summary function
) 


# Train the model using earth
set.seed(1234); model_mars = train(Time~ ., 
		data=train_df, method='earth', tuneLength=20, trControl = fitControl)
varimp_mars <- varImp(model_mars)
res_mars <- evalm(model_mars)

# Train the model using glm
set.seed(1234); model_glm = train(Time~ ., 
		data=train_df, method='glm', tuneLength=20, trControl = fitControl)
varimp_glm <- varImp(model_glm)
res_glm <- evalm(model_glm)

# Train the model using rf
set.seed(1234); model_rf = train(Time~ ., 
		data=train_df, method='rf', tuneLength=20, trControl = fitControl, importance=TRUE)
varimp_rf <- varImp(model_rf)
res_rf <- evalm(model_rf)

# Train the model using svm
set.seed(1234); model_svmRadial = train(Time~ ., 
		data=train_df, method='svmRadial', tuneLength=10, trControl = fitControl, importance=TRUE)
svm.imp <-Importance(model_svmRadial, data=train)
res_SVM <- evalm(model_svmRadial)

# Train the model using naive Bayes
set.seed(1234); model_NB = train(Time~ ., 
		data=train_df, method='naive_bayes', tuneLength=10, trControl = fitControl)
varimp_NB <- varImp(model_NB)
res_NB <- evalm(model_NB)

# Train the model using Gradient Boost
set.seed(1234); model_GBM = train(Time~ ., 
		data=train_df, method='gbm', tuneLength=10, trControl = fitControl)
varimp_GBM <- summary.gbm(model_GBM)
res_GBM <- evalm(model_GBM)

# Train the model using decision trees
set.seed(1234); model_DT = train(Time~ ., 
		data=train_df, method='rpart', tuneLength=10, trControl = fitControl)
varimp_DT <- varImp(model_DT)
res_DT <- evalm(model_DT)

##  Then compare model performances using resample()
models_compare <- resamples(list(GLM=model_glm, RF=model_rf,
					NaiveBayes=model_NB,MARS=model_mars, 
					GBM=model_GBM,SVM=model_svmRadial, DT=model_DT))

# Summary of the models performances
sink("1.SEIZURE_Training set.txt")
cat("Summary of Training using ROC metrics")
print(summary(models_compare))
sink()

# Plot ROCs
test4 <- evalm(list(GLM=model_glm, RF=model_rf,
			NaiveBayes=model_NB,MARS=model_mars, 
			GBM=model_GBM,SVM=model_svmRadial, DT=model_DT),
			gnames=c('glm', 'rf', 'NB', 'mars', 'GBM', 'svm', 'DT'),rlinethick=0.8,fsize=8,
plots='r')

tiff(file="PLOTS/1.Seizure_train__ROC.tiff", unit= "in", res = 300, width = 6, height = 6)
test4$roc
dev.off()

# Draw box plots to compare models
tiff(file="PLOTS/1.Seizure_train_Box_comp).tiff", unit= "in", res = 300, width = 6, height = 6)
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)
dev.off()

#############################################################################
#                        TESTING		         				    #
#############################################################################
# Predict on test.  Data print to summary sheet
sink("2.SEIZURE model testing.txt")
cat("Summary of Testing predictions using confusion matrix/n
cat("Reference testing population = 4 baseline PNES, 4 post PNES pateints./n")
cat("Reference training population = 27 baseline, 27 post SEIZURE pateints./n")

# model_mars
cat("MARS-model on test\n")
predicted_test <- predict(model_mars, newdata=test_df)
cat("Predicted values\n")
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
cat("MARS-model on train\n")
#train
predicted_train=predict(model_mars, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), 
		data = predicted_train, mode='everything', positive='Y'))

#model_glm
cat("glm model on test\n")
predicted_test <- predict(model_glm, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y')) 
#train
cat("glm model on train\n")
predicted_train=predict(model_glm, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_rf
cat("rf-model on test\n")
predicted_test <- predict(model_rf, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y')) 
#train
cat("rf model on train\n")
predicted_train=predict(model_rf, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_NB
cat("Naive-Bayes model on test\n")
predicted_test <- predict(model_NB, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))

#train
cat("Naive-Bayes model on train\n")
predicted_train=predict(model_NB, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))			

#model_GBM
cat("GBM model on test\n")
predicted_test <- predict(model_GBM, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("GBM model on train\n")
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_svmRadial
cat("svmRadial model on test\n")
predicted_test <- predict(model_svmRadial, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("svmRadial model on train\n")
predicted_train=predict(model_svmRadial, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_DT
cat("Decision Tree model on test\n")
predicted_test <- predict(model_DT, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("Decision Tree model on train\n")
predicted_train=predict(model_DT, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))
sink()

## Multiple using caretEnsembl.
# recall control parameters from above
trainControl <- trainControl(method='repeatedcv', 
                             number=10, 
                             repeats=3,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

algorithmList <- c('rf', 'gbm', 'earth','svmRadial','naive_bayes', 'rpart')

# Then run
set.seed(1234); models <- caretList(Time~., 
		data=train_df, 
		trControl=trainControl,
		methodList=algorithmList) 
results <- resamples(models)
summary(results)

# plot results as a correlation
xyplot(resamples(models))
data <- modelCor(resamples(models))
data
 
# creating correlation matrix
library (corrplot)
corrplot(data, method = 'color', order = 'alphabet')
# prefer this one!

# Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

# Combine the predictions of multiple models to form a final prediction
# Create the trainControl
stackControl <- trainControl(method="repeatedcv", 
                             number=20, 
                             repeats=3,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

## Combining models into single glm
set.seed(1234);stack.glm <- caretStack(models, method="glm", metric="ROC", trControl=stackControl)
print(stack.glm)

#test_df <- predict(preProcess_range_model, newdata = test_df)
# Predict on testData
predicted_test <- predict(stack.glm, newdata=test_df)
confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, mode='everything', positive='Y') 

predicted_train=predict(stack.glm, newdata=train_df)
confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, mode='everything', positive='Y')

###################################################################################
###################################################################################
###Now create a model on the PNES data

# Create test and training data frame 
## PNES samples and controls
df5 = df2[!(df2$"EEG.change"=="YES" & df2$Time=="4-6h"),]
df5[1:10,1:20];dim(df5)# 35 samples

df6 = df2[(df2$"EEG.change"=="YES"),]
df6[1:8,1:20];dim(df6)# 54 samples

#Paritition (split) and create index matrix of selected values 
#index <- createDataPartition(df$NIHSS_Delta, p=.7, list=FALSE, times=1)

# To address error message, convert df to data frame object 
train_df <- as.data.frame(df5)
test_df <- data.frame(df6[,])
train_df <- train_df[,c(5,20:12572)]
train_df[1:10,1:20]; dim(train_df)

##Pull in data sheet of sig genes..
genes=read.csv("DATA/60 PNES4-6hvsAll_genelist.csv")
genes[1:20,1:5]; dim(genes) # 60x26 # Genes are in colum Gene.Symbol


## Now subset gene expression values based on inclusion in list.
train_df_0  <- (train_df[,names(train_df) %in% genes$"Gene.Symbol"])
head(train_df_0);dim(train_df_0)
# now add back Time info
dfz=cbind(train_df$Time, train_df_0[,])
# rename col1
colnames(dfz)[1] <- "Time"
head(dfz); dim(dfz)
## rename back to train_df
train_df = dfz

## Data needs to be a binomial/ classifier model!!! recode data to A and B or Y/ N
train_df$Time <- recode(train_df$Time, 'Base' = 'N', '4-6h'='Y', 
	.default=NA_character_)
head(train_df); dim(train_df)

# better recode the test data
test_df$Time <- recode(test_df$Time, 'Base' = 'N', '4-6h'='Y', 
	.default=NA_character_)
head(test_df); dim(test_df)

# quick assessment..
skim(train_df)


#############################################################################
#                        PREDICTION MODELING				          #
#############################################################################
##Preprocessing
preProcess_range_model <- preProcess(train_df, method=c('scale', 'center'))
train_df <- predict(preProcess_range_model, newdata = train_df)

# process the test data too
test_df <- predict(preProcess_range_model, newdata = test_df)


#preProcess_range_model <- preProcess(train_df, method='pca')
#train_df <- predict(preProcess_range_model, newdata = train_df)
dim(train_df)

featurePlot(x = train_df[, 2:51], 
            y = as.factor(train_df$Time), 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

featurePlot(x = train_df[, 2:51], 
            y = as.factor(train_df$Time), 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))


## Lets try all of these together...
set.seed(1234)

# Define the training control
fitControl <- trainControl(
    method = 'cv',                   # k-fold cross validation
    number = 10,                      # number of folds
    savePredictions = 'final',       # saves predictions for optimal tuning parameter
    classProbs = TRUE,                  # should class probabilities be returned
    summaryFunction=twoClassSummary  # results summary function
) 

# Train the model using earth
set.seed(1234); model_mars = train(Time~ ., 
		data=train_df, method='earth', tuneLength=20, trControl = fitControl)
# Train the model using glm
set.seed(1234); model_glm = train(Time~ ., 
		data=train_df, method='glm', tuneLength=20, trControl = fitControl)
# Train the model using rf
set.seed(1234); model_rf = train(Time~ ., 
		data=train_df, method='rf', tuneLength=20, trControl = fitControl)
# Train the model using svm
set.seed(1234); model_svmRadial = train(Time~ ., 
		data=train_df, method='svmRadial', tuneLength=10, trControl = fitControl)
# Train the model using naive Bayes
set.seed(1234); model_NB = train(Time~ ., 
		data=train_df, method='naive_bayes', tuneLength=10, trControl = fitControl)
# Train the model using Gradient Boost
set.seed(1234); model_GBM = train(Time~ ., 
		data=train_df, method='gbm', tuneLength=10, trControl = fitControl)

# Train the model using decision trees
set.seed(1234); model_DT = train(Time~ ., 
		data=train_df, method='rpart', tuneLength=10, trControl = fitControl)

##  Then compare
# Compare model performances using resample()
models_compare <- resamples(list(GLM=model_glm, RF=model_rf,
					NaiveBayes=model_NB,MARS=model_mars, 
					GBM=model_GBM,SVM=model_svmRadial, DT=model_DT))
# Summary of the models performances

sink("3.PNES_Training set.txt")
cat("Summary of Training using ROC metrics")
print(summary(models_compare))
sink()

# Plot ROCs
test4 <- evalm(list(GLM=model_glm, RF=model_rf,
			NaiveBayes=model_NB,MARS=model_mars, 
			GBM=model_GBM,SVM=model_svmRadial, DT=model_DT),
			gnames=c('glm', 'rf', 'NB', 'mars', 'GBM', 'svm', 'DT'),rlinethick=0.8,fsize=8,
plots='r')

tiff(file="PLOTS/3.PNES_train__ROC.tiff", unit= "in", res = 300, width = 6, height = 6)
test4$roc
dev.off()

# Draw box plots to compare models
tiff(file="PLOTS/4.PNES_train_Box_comp).tiff", unit= "in", res = 300, width = 6, height = 6)
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)
dev.off()

#############################################################################
#                        TESTING		         				    #
#############################################################################
# Predict on test.  Data print to summary sheet
sink("4.PNES model testing.txt")
cat("Summary of Testing predictions using confusion matrix/n
cat("Reference testing population = 27 baseline, 27 post SEIZURE pateints ./n")
cat("Reference training population = 4 baseline PNES, 4 post PNES patients./n")

# model_mars
cat("MARS-model on test\n")
predicted_test <- predict(model_mars, newdata=test_df)
cat("Predicted values\n")
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
cat("MARS-model on train\n")
#train
predicted_train=predict(model_mars, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), 
		data = predicted_train, mode='everything', positive='Y'))

#model_glm
cat("glm model on test\n")
predicted_test <- predict(model_glm, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y')) 
#train
cat("glm model on train\n")
predicted_train=predict(model_glm, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_rf
cat("rf-model on test\n")
predicted_test <- predict(model_rf, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y')) 
#train
cat("rf model on train\n")
predicted_train=predict(model_rf, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_NB
cat("Naive-Bayes model on test\n")
predicted_test <- predict(model_NB, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))

#train
cat("Naive-Bayes model on train\n")
predicted_train=predict(model_NB, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))			

#model_GBM
cat("GBM model on test\n")
predicted_test <- predict(model_GBM, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("GBM model on train\n")
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_svmRadial
cat("svmRadial model on test\n")
predicted_test <- predict(model_svmRadial, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("svmRadial model on train\n")
predicted_train=predict(model_svmRadial, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))

#model_DT
cat("Decision Tree model on test\n")
predicted_test <- predict(model_DT, newdata=test_df)
print(predicted_test)
print(confusionMatrix(reference = as.factor(test_df$Time), data = predicted_test, 
		mode='everything', positive='Y'))
#train
cat("Decision Tree model on train\n")
predicted_train=predict(model_DT, newdata=train_df)
print(predicted_train)
print(confusionMatrix(reference = as.factor(train_df$Time), data = predicted_train, 
		mode='everything', positive='Y'))
sink()


####################################################################################################################################






