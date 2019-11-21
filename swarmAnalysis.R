setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
swarm <- read.csv(file = "sAreaTable.csv")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
swarm
tiff(filename = "histSwarm.tiff")
hist(swarm[, 2])
dev.off()

shapiroRes <- list()
for(i in seq_along(levels(swarm[, 1]))){
        vec <- swarm[which(swarm[, 1] == levels(swarm[, 1])[i]), 2]
        if(length(vec) >= 3){
                shapiroRes[[i]] <- shapiro.test(vec) 
        }else{
                shapiroRes[[i]] <- "Not enough samples"
        }
}
names(shapiroRes) <- levels(swarm[, 1])
shapiroRes <- shapiroRes[!unlist(lapply(shapiroRes, function(f) (f == "Not enough samples")[1]))]

shapPVals <- sapply(shapiroRes, function(f) f$p.value)

sum(shapPVals > 0.05)
sum(shapPVals <= 0.05)

# We assume normal distribution


swarmMeans <- c()
for(i in seq_along(levels(swarm[, 1]))){
        vec <- swarm[which(swarm[, 1] == levels(swarm[, 1])[i]), 2]
        swarmMeans <- c(swarmMeans, mean(vec))
}
names(swarmMeans) <- levels(swarm[, 1])

####################################################################################################################
#   Logistic regression
####################################################################################################################
if (!require(caret)) install.packages('caret')
library(caret)

# Binarize swarming according to histogram: if greater than -2 they are swarmers.

binarizeSwarm <- function(x, threshold = -2){
        swarmBin <- x
        swarmBin[swarmMeans > threshold] <- "Swarmer"
        swarmBin[swarmMeans < threshold] <- "nonSwarmer"
        swarmBin <- as.factor(swarmBin)
        return(swarmBin)
}

swarmBin <- binarizeSwarm(swarmMeans, threshold = -2)
names(swarmBin) <- gsub(".*_", names(swarmBin), replacement = "")
strainNames <- unique(gsub("\\_.*|(PA14).*", rownames(ccmn_norm_mets_good_old), rep = "\\1"))
# W70332 is mispelled
strainNames[strainNames == "W70322"] <- "W70332"
# Remove strains that we don't have in metabolomic data
swarmBin <- swarmBin[names(swarmBin) %in% strainNames]

swarmBin <- swarmBin[match(strainNames, names(swarmBin))]
swarmBin <- swarmBin[!is.na(swarmBin)]

notInSwarmData <- strainNames[!strainNames %in% names(swarmBin)]

rownames(ccmn_norm_mets_good_old)[grep("W70322", rownames(ccmn_norm_mets_good_old))] <- paste("W70332", 
                                                                                              grep("W70322", rownames(ccmn_norm_mets_good_old)), 
                                                                                              sep = "_")

ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[-grep(notInSwarmData, rownames(ccmn_norm_mets_good_old)), ]

ccmn_norm_mets_good_old <- cbind.data.frame(ccmn_norm_mets_good_old, 
                                            swarmBin[match(gsub("\\_.*|(PA14).*", 
                                                                rownames(ccmn_norm_mets_good_old), 
                                                                rep = "\\1"), 
                                                           names(swarmBin))])

colnames(ccmn_norm_mets_good_old)[ncol(ccmn_norm_mets_good_old)] <- "swarmData"
colnames(ccmn_norm_mets_good_old) <- make.names(colnames(ccmn_norm_mets_good_old))

set.seed(107)
inTrain <- createDataPartition(
        y = ccmn_norm_mets_good_old$swarmData,
        p = 0.5,
        list = F
)
training <- ccmn_norm_mets_good_old[inTrain,]
testing <- ccmn_norm_mets_good_old[-inTrain,]
nrow(training)
nrow(testing)

ctrl <- trainControl(method = "repeatedcv", 
                     repeats = 3,
                     classProbs = T,
                     summaryFunction = twoClassSummary)

bstLogFit <- train(
        swarmData~.,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)


bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$swarmData)#-->Good accuracy, but it doesn't beat No Info Rate

#test importance of variables

importance <- varImp(bstLogFit, scale=F)
print(importance)
plot(importance)
#repeat with 3 more important variables
bstLogFit <- train(
        swarmData~X4.Aminobutanoate.3 +
                X2.Acetyl.aminoadipate +
                X2.5.Dioxopentanoate.2 +
                X2.Amino.6.oxopimelate +
                X5.Oxo.proline.1,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$swarmData)#-->Good accuracy, but it doesn't beat No Info Rate


#Apply recursive feature elimination to evaluate the importance of variables
if(!require(randomForest)) install.packages("randomForest")
library(randomForest)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)

results <- rfe(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) -1)], 
               ccmn_norm_mets_good_old[, ncol(ccmn_norm_mets_good_old)], 
               sizes = c(1:(ncol(ccmn_norm_mets_good_old) -1)), 
               rfeControl = control)

print(results)

# Accuracy higher than 90% with 16 metabolites. Those are:

# Acetylhomoserine.2 
# Acetylhomoserine
# X2.Aminoadipate.1
# Tartronate.semialdehyde
# Formyl.methionine            
# Dimethylmalate.1
# Erithro.3.Methylmalate.2
# d.ala.d.ala.2
# X2.aminoadipate.6.semialdehyde
# N.acetyl.2.4.diaminobutanoate.1
# Anthranilate.1
# X3.Methyl.2.oxopentanoate.1
# Histidine
# X3.Amino.isobutanoate.1
# Citrulline                   
# glutamate 

predictors(results)

plot(results, type=c("g", "o"))

bstLogFit <- train(
        swarmData~Acetylhomoserine.2 +
                Acetylhomoserine +
                X2.Aminoadipate.1 +
                Tartronate.semialdehyde +
                Formyl.methionine +
                Dimethylmalate.1 +
                Dimethylmalate.1 +
                Erithro.3.Methylmalate.2 +
                d.ala.d.ala.2 +
                X2.aminoadipate.6.semialdehyde +
                N.acetyl.2.4.diaminobutanoate.1 +
                Anthranilate.1 +
                X3.Methyl.2.oxopentanoate.1 +
                Histidine +
                X3.Amino.isobutanoate.1 +
                Citrulline +
                glutamate,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$swarmData)
