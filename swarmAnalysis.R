setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
swarm <- read.csv(file = "sAreaTable.csv")
#Load functions
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/diffMetAnal_functions.R")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
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
getSwarmMeans <- function(swarmDat){
        swarmMeans <- c()
        for(i in seq_along(levels(swarmDat[, 1]))){
                vec <- swarmDat[which(swarmDat[, 1] == levels(swarmDat[, 1])[i]), 2]
                swarmMeans <- c(swarmMeans, mean(vec))
        }
        names(swarmMeans) <- levels(swarmDat[, 1])
        return(swarmMeans)
}

swarmMeans <- getSwarmMeans(swarm)

####################################################################################################################
####################################################################################################################
#   Logistic regression
####################################################################################################################
####################################################################################################################

####################################################################################################################
# With metabolomic data 
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

swarmBin <- binarizeSwarm(swarmMeans, threshold = -1.7)
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

# Remove ambiguous mets

#ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[, !is.na(dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)])]
#colnames(ccmn_norm_mets_good_old) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
#ccmn_norm_mets_good_old <- rmAmbig(ccmn_norm_mets_good_old)


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

bstLogFitPre <- train(
        swarmData~.,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)


bstClassesPre <- predict(bstLogFitPre, newdata = testing)
head(bstClassesPre)

confusionMatrix(data = bstClassesPre, testing$swarmData)#-->Good accuracy, but it doesn't beat No Info Rate

#test importance of variables

importance <- varImp(bstLogFitPre, scale=F)
print(importance)
plot(importance)
formImp <- as.formula(paste("swarmData", paste(rownames(importance$importance)[1:5], collapse = " + "), sep = " ~ "))
#repeat with 3 more important variables
bstLogFitVarImp <- train(
        formImp,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClassesVarImp <- predict(bstLogFitVarImp, newdata = testing)
head(bstClassesVarImp)

confusionMatrix(data = bstClassesVarImp, testing$swarmData)#-->Good accuracy, but it doesn't beat No Info Rate


#Apply recursive feature elimination to evaluate the importance of variables
if(!require(randomForest)) install.packages("randomForest")
library(randomForest)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)

results <- rfe(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) -1)], 
               ccmn_norm_mets_good_old[, ncol(ccmn_norm_mets_good_old)], 
               sizes = c(1:(ncol(ccmn_norm_mets_good_old) -1)), 
               rfeControl = control)

print(results)


predictors(results)

plot(results, type=c("g", "o"))

rfeResult <- dictionary$Consensus[match(predictors(results), make.names(dictionary$`Old Data Names`))]
rfeResult <- rfeResult[!is.na(rfeResult)]
names(rfeResult) <- rfeResult
rfeResult <- rmAmbig(rfeResult)
formRFE <- as.formula(paste("swarmData", paste(make.names(rfeResult), collapse = " + "), sep = " ~ "))
rfeResultKEGGIDs <- dictionary$`KEGG IDs`[match(rfeResult, dictionary$Consensus)]

bstLogFitRFE <- train(
        formRFE,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClassesRFE <- predict(bstLogFitRFE, newdata = testing)
head(bstClassesRFE)

confusionMatrix(data = bstClassesRFE, testing$swarmData) # --> Accuracy of 0.8 when using only unambiguous predictors. No beating NIR

# Do pathway enrichment 
if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)

# With ORA
allMets <- dictionary$`KEGG IDs`[!grepl("_?", dictionary$Consensus, fixed = T)][!is.na(dictionary$`KEGG IDs`[!grepl("_?", dictionary$Consensus, fixed = T)])]

ORA_rfeResults <- doORA(diffMetObjkt = rfeResultKEGGIDs,
                        allMetsObjkt = allMets,
                        org = "pae")

# With FELLA
graph <- buildGraphFromKEGGREST(
        organism = "pae",
        filter.path = c("01100", "01200", "01210", "01212", "01230")
)

buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = NULL,
        internalDir = TRUE,
        matrices = "none",
        normality = "diffusion",
        niter = 100)


alias2entrez <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
entrez2ec <- KEGGREST::keggLink("enzyme", "pae")
entrez2path <- KEGGREST::keggLink("pathway", "pae")
fella.data <- loadKEGGdata(
        #databaseDir = "created__2019-03-25;meta__pae_Release_89.0_03_23_Mar_19",
        internalDir = T,
        loadMatrix = "none"
)

fella.data

id.cpd <- getCom(fella.data, level = 5, format = "id") %>% names
id.rx <- getCom(fella.data, level = 4, format = "id") %>% names
id.ec <- getCom(fella.data, level = 3, format = "id") %>% names

analysis.rfeRes<- enrich(
        compounds = rfeResultKEGGIDs,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.rfeRes %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.rfeRes)

g_rfeSwarm <- generateResultsGraph(
        object = analysis.rfeRes,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_rfeSwarm

tiff("FELLA_clusts_rfeRes_swarm.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_rfeSwarm
        #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_rfeSwarm <- generateResultsTable(
        object = analysis.rfeRes,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)


# Same pathways than the ones obtained in the ORA with differential metabolites

# Lets make dendogram with swarm in label

if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

ccmnNormMeds <- getStrainMedian(ccmn_norm_mets_good_old[, -ncol(ccmn_norm_mets_good_old)])
ccmnNormMeds <- ccmnNormMeds[-1, ]
rownames(ccmnNormMeds)[1] <- "PA14"

dend <- ccmnNormMeds %>% dist(method = "euclidean") %>%
        hclust(method = "ward.D") 

df_ccmnNormMeds <- data.frame(Strain = rownames(ccmnNormMeds),
                              SwarmData = as.character(swarmBin))

dend_hcData <- dendro_data(dend, type="rectangle")

dend_hcData$labels <- merge(x = dend_hcData$labels,
                            y = df_ccmnNormMeds,  by.x = "label", by.y = "Strain")

swarmCols <- c("firebrick", "blue3")

ggplot() +
        geom_segment(data=segment(dend_hcData), aes(x=x, y=y, xend=xend, yend=yend)) +
        geom_text(data = label(dend_hcData), aes(x=x, y=y, label=label, colour = SwarmData, hjust=0), size=5) +
        geom_point(data = label(dend_hcData), aes(x=x, y=y), size=0.1, shape = 21) +
        coord_flip() +
        scale_y_reverse(expand=c(0.2, 0)) +
        scale_colour_brewer(palette = "Dark2") + 
        scale_color_manual(values = swarmCols)


ggsave(filename = "HCA_with_swarm.tiff", plot = last_plot(), device = "tiff", 
       path = NULL,
       scale = 1, width = 18, height = 40, units = "cm",
       dpi = 300, limitsize = F)

# Whith a threshold of -1.7 in log swarming area all the strains in cluster 2 are nonswarmers,
# while the 66.66% of the swarmer strains come from blood. 

####################################################################################################################
# With genomic data 
####################################################################################################################
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_enz_tab_filt.RData")


jaccGroup <- function(genePresAbsObjkt, threshold = 0.05){
        if(!require(proxy)) install.packages('proxy')
        library(proxy)
        geneTab <- sapply(genePresAbsObjkt, as.integer)
        rownames(geneTab) <- rownames(genePresAbsObjkt)
        colnames(geneTab) <- make.names(colnames(geneTab), unique = T)
        jaccDist <- proxy::dist(geneTab, by_rows = F, method = "Jaccard", diag = T)
        jaccHCA <- hclust(jaccDist, method = "complete")
        jaccGroups <- cutree(jaccHCA, h = threshold)
        groups <- c()
        for(i in 1:max(jaccGroups)){
                groups <- c(groups, names(sort(jaccGroups)[which(sort(jaccGroups) == i)])[1])
        }
        geneTabGrouped <- geneTab[, groups]
        return(geneTabGrouped)
}

geneEnzTabFiltGrouped <- jaccGroup(gene_enz_tab_filt)

rownames(geneEnzTabFiltGrouped)[rownames(geneEnzTabFiltGrouped) == "W70322"] <- "W70332"

geneEnzTabFiltGrouped <- geneEnzTabFiltGrouped[rownames(geneEnzTabFiltGrouped) %in% names(swarmBin), ]

geneEnzTabFiltGrouped <- cbind.data.frame(as.data.frame(geneEnzTabFiltGrouped), 
                                          swarmBin[names(swarmBin) %in% rownames(geneEnzTabFiltGrouped)])
colnames(geneEnzTabFiltGrouped)[ncol(geneEnzTabFiltGrouped)] <- "swarmData"

set.seed(108)
inTrainEnz <- createDataPartition(
        y = geneEnzTabFiltGrouped$swarmData,
        p = 0.5,
        list = F
)
trainingEnz <- geneEnzTabFiltGrouped[inTrainEnz,]
testingEnz <- geneEnzTabFiltGrouped[-inTrainEnz,]


bstLogFitPreEnz <- train(
        swarmData~.,
        data = trainingEnz,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClassesPreEnz <- predict(bstLogFitPreEnz, newdata = testingEnz)
head(bstClassesPreEnz)

confusionMatrix(data = bstClassesPreEnz, testingEnz$swarmData)

importanceEnz <- varImp(bstLogFitPreEnz, scale=F)
print(importanceEnz)
plot(importanceEnz)
formImpEnz <- as.formula(paste("swarmData", paste(rownames(importanceEnz$importance)[1:5], collapse = " + "), sep = " ~ "))
#repeat with 5 more important variables
bstLogFitVarImpEnz <- train(
        formImpEnz,
        data = trainingEnz,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClassesVarImpEnz <- predict(bstLogFitVarImpEnz, newdata = testingEnz)
head(bstClassesVarImpEnz)

confusionMatrix(data = bstClassesVarImpEnz, testingEnz$swarmData)


# RFE with genomic data (only accessory enzymatic genome)
resultsEnz <- rfe(geneEnzTabFiltGrouped[, 1:(ncol(geneEnzTabFiltGrouped) -1)], 
                  geneEnzTabFiltGrouped[, ncol(geneEnzTabFiltGrouped)], 
                  sizes = c(1:(ncol(geneEnzTabFiltGrouped) -1)), 
                  rfeControl = control)

print(resultsEnz)

plot(resultsEnz, type=c("g", "o"))

predictors(resultsEnz)

formRFEEnz <- as.formula(paste("swarmData", paste(predictors(resultsEnz), collapse = " + "), sep = " ~ "))

bstLogFitRFEEnz <- train(
        formRFEEnz,
        data = trainingEnz,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClassesRFEEnz <- predict(bstLogFitRFEEnz, newdata = testingEnz)
head(bstClassesRFEEnz)

confusionMatrix(data = bstClassesRFEEnz, testingEnz$swarmData) 


# The enzymatic genes obtained are symilar to the ones obtained when classifying the strains in the 
# two major clusters.
