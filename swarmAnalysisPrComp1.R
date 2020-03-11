setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
load("pcomp1AreaPercCircularity.RData")
swarm <- read.csv(file = "sAreaTable.csv")
swarmMeans <- pcomp1Means
#Load functions
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/swarmFunctions.R")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")

####################################################################################################################
####################################################################################################################
####################################################################################################################
# With metabolomic data 
####################################################################################################################
####################################################################################################################
####################################################################################################################

####################################################################################################################
####################################################################################################################
#   Logistic regression
####################################################################################################################
####################################################################################################################
if (!require(caret)) install.packages('caret')
library(caret)

# W70332 is mispelled
rownames(ccmn_norm_mets_good_old) <- gsub("W70322", replacement = "W70332", rownames(ccmn_norm_mets_good_old))
names(swarmMeans) <- gsub("W70322", replacement = "W70332", names(swarmMeans))

# Binarize swarming according to histogram: if greater than -2 they are swarmers.
swarmBin <- binarizeSwarm(swarmMeans, threshold = -0.25)
names(swarmBin) <- gsub(".*_", names(swarmBin), replacement = "")
strainNames <- unique(gsub("\\_.*|(PA14).*", rownames(ccmn_norm_mets_good_old), rep = "\\1"))


# Remove strains that we don't have in metabolomic data
swarmBin <- swarmBin[names(swarmBin) %in% strainNames]

swarmBin <- swarmBin[match(strainNames, names(swarmBin))]
swarmBin <- swarmBin[!is.na(swarmBin)]

#notInSwarmData <- strainNames[!strainNames %in% names(swarmBin)]


#ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[-grep(notInSwarmData, rownames(ccmn_norm_mets_good_old)), ]


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

save(rfeResultKEGGIDs, file = "rfeResultKEGGIDs.RData")

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

write.csv(ORA_rfeResults, file = "ORA_rfeResults.csv")
save(ORA_rfeResults, file = "ora_rfeResults.RData")

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

write.csv(tab_rfeSwarm, file = "tab_rfeSwarm.csv")
save(tab_rfeSwarm, file = "tab_rfeSwarm.RData")
# Same pathways than the ones obtained in the ORA with differential metabolites

# Lets make dendogram with swarm in label

if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

if(!require(ggdendro)) install.packages("ggdendro")
library(ggdendro)

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

ccmnNormMeds <- getStrainMedian(ccmn_norm_mets_good_old[, -ncol(ccmn_norm_mets_good_old)])

swarmForDend <- ccmn_norm_mets_good_old$swarmData
names(swarmForDend) <- gsub("_.*", rownames(ccmn_norm_mets_good_old), replacement = "")

swarmForDend <- swarmForDend[match(unique(names(swarmForDend)), names(swarmForDend))]

dend <- ccmnNormMeds %>% dist(method = "euclidean") %>%
        hclust(method = "ward.D") 

df_ccmnNormMeds <- data.frame(Strain = rownames(ccmnNormMeds),
                              SwarmData = as.character(swarmForDend))

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
####################################################################################################################
#   Multivariate analysis
####################################################################################################################
####################################################################################################################

if(!require(ropls)) BiocManager::install("ropls")
library(ropls)

swarmMeansRearr <- swarmMeans
names(swarmMeansRearr) <- gsub(".*_", names(swarmMeans), replacement = "")
swarmMeansRearr <- swarmMeansRearr[names(swarmMeansRearr) %in% strainNames]
swarmMeansRearr <- swarmMeansRearr[match(strainNames, names(swarmMeansRearr))]
swarmMeansRearr <- swarmMeansRearr[!is.na(swarmMeansRearr)]

save(swarmMeansRearr, file = "swarmMeansRearr.RData")

ccmn_norm_mets_good_old <- cbind.data.frame(ccmn_norm_mets_good_old, 
                                            swarmMeansRearr[match(gsub("\\_.*|(PA14).*", 
                                                                       rownames(ccmn_norm_mets_good_old), 
                                                                       rep = "\\1"), 
                                                                  names(swarmMeansRearr))])

colnames(ccmn_norm_mets_good_old)[ncol(ccmn_norm_mets_good_old)] <- "swarmQuant"

ccmnNormPCA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)])
swarmData <- ccmn_norm_mets_good_old$swarmData
tiff(filename = "pcaMetDataSwarmQual.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(ccmnNormPCA,
     typeVc = "x-score",
     parAsColFcVn = swarmData)
dev.off()

tiff(filename = "pcaMetDataSwarmQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(ccmnNormPCA,
     typeVc = "x-score",
     parAsColFcVn = ccmn_norm_mets_good_old$swarmQuant)
dev.off()

ccmnNormPLSQuant <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)], ccmn_norm_mets_good_old$swarmQuant)
ccmnNormPLSQual <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)], ccmn_norm_mets_good_old$swarmData)

tiff(filename = "ccmnNormOPLSDAQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
ccmnNormOPLSDAQuant <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)],
                            ccmn_norm_mets_good_old$swarmQuant, 
                            predI = 1, 
                            orthoI = NA#,
                            #subset = as.vector(inTrain)
)
dev.off()

jpeg(filename = "ccmnNormOPLSDAQuant.jpeg", res = 300, height = 3000, width = 3000, units = "px")
ccmnNormOPLSDAQuant <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)],
                            ccmn_norm_mets_good_old$swarmQuant, 
                            predI = 1, 
                            orthoI = NA#,
                            #subset = as.vector(inTrain)
)
dev.off()

tiff(filename = "ccmnNormOPLSDAQual.tiff", res = 300, height = 3000, width = 3000, units = "px")
ccmnNormOPLSDAQual <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)],
                           ccmn_norm_mets_good_old$swarmData, 
                           predI = 1, 
                           orthoI = NA#,
                           #subset = as.vector(inTrain),
)
dev.off()

table(ccmn_norm_mets_good_old$swarmData[as.vector(inTrain)], fitted(ccmnNormOPLSDAQual)[as.vector(inTrain)])

table(ccmn_norm_mets_good_old$swarmData[as.vector(inTrain)],
      predict(ccmnNormOPLSDAQual, ccmn_norm_mets_good_old[as.vector(inTrain), 1:(ncol(ccmn_norm_mets_good_old) - 2)]))

table(ccmn_norm_mets_good_old$swarmData[-as.vector(inTrain)],
      predict(ccmnNormOPLSDAQual, ccmn_norm_mets_good_old[-as.vector(inTrain), 1:(ncol(ccmn_norm_mets_good_old) - 2)]))
# Accuracy of 0.95 with the testing dataset

# Build a matrix of the loading values for each component
OPLSDAQuantLoads <- cbind(getLoadingMN(ccmnNormOPLSDAQuant)[, 1],
                          getLoadingMN(ccmnNormOPLSDAQuant, orthoL = T))
colnames(OPLSDAQuantLoads)[1] <- "p1"

OPLSDAQuantLoads <- getLoadingMN(ccmnNormOPLSDAQuant)

OPLSDAQualLoads <- cbind(getLoadingMN(ccmnNormOPLSDAQual)[, 1],
                         getLoadingMN(ccmnNormOPLSDAQual, orthoL = T))
colnames(OPLSDAQualLoads)[1] <- "p1"

OPLSDAQualLoads <- getLoadingMN(ccmnNormOPLSDAQual)

# Make barplot of loadings for OPLSDA quantitative. Previously we're gonna remove the ambiguous compounds.

OPLSDAQuantLoadsRmAmbig <- OPLSDAQuantLoads
rownames(OPLSDAQuantLoadsRmAmbig) <- dictionary$Consensus[match(rownames(OPLSDAQuantLoads), make.names(dictionary$`Old Data Names`))]
OPLSDAQuantLoadsRmAmbig <- OPLSDAQuantLoadsRmAmbig[!is.na(rownames(OPLSDAQuantLoadsRmAmbig)), ]
OPLSDAQuantLoadsRmAmbig <- OPLSDAQuantLoadsRmAmbig[-grep("_?", names(OPLSDAQuantLoadsRmAmbig), fixed = T)]
OPLSDAQuantLoadsRmAmbig <- cbind.data.frame(names(OPLSDAQuantLoadsRmAmbig), OPLSDAQuantLoadsRmAmbig)
colnames(OPLSDAQuantLoadsRmAmbig) <- c("metabolite", "loading")

ggplot(data = OPLSDAQuantLoadsRmAmbig,
       aes(x = reorder(metabolite, loading), y = loading)) +
  labs(x = "compound", y = "OPLS-DA Loading") +
  geom_bar(stat = "identity") +
  coord_flip()
ggsave(filename = "OPLSDAQuantLoadsBarplot.tiff")
ggsave(filename = "OPLSDAQuantLoadsBarplot.jpeg")
# Obtain the 3variables with most extreme loading values (positive and negative) of each component
extremValsOPLSDAQual <- getExtremVals(OPLSDAQualLoads, n = 14)
extremValsOPLSDAQuant <- getExtremVals(OPLSDAQuantLoads, n = 14)

OPLSDAQualResult <- dictionary$Consensus[match(extremValsOPLSDAQual$uniqueExtremeVars, make.names(dictionary$`Old Data Names`))]
OPLSDAQualResult <- OPLSDAQualResult[!is.na(OPLSDAQualResult)]
names(OPLSDAQualResult) <- OPLSDAQualResult
OPLSDAQualResult <- rmAmbig(OPLSDAQualResult)
OPLSDAQualResultKEGGIDs <- dictionary$`KEGG IDs`[match(OPLSDAQualResult, dictionary$Consensus)]
OPLSDAQualResultKEGGIDs <- OPLSDAQualResultKEGGIDs[!is.na(OPLSDAQualResultKEGGIDs)]
save(OPLSDAQualResultKEGGIDs, file = "OPLSDAQualResultKEGGIDs.RData")

OPLSDAQuantResult <- dictionary$Consensus[match(extremValsOPLSDAQuant$uniqueExtremeVars, make.names(dictionary$`Old Data Names`))]
OPLSDAQuantResult <- OPLSDAQuantResult[!is.na(OPLSDAQuantResult)]
names(OPLSDAQuantResult) <- OPLSDAQuantResult
OPLSDAQuantResult <- rmAmbig(OPLSDAQuantResult)
OPLSDAQuantResultKEGGIDs <- dictionary$`KEGG IDs`[match(OPLSDAQuantResult, dictionary$Consensus)]
OPLSDAQuantResultKEGGIDs <- OPLSDAQuantResultKEGGIDs[!is.na(OPLSDAQuantResultKEGGIDs)]
save(OPLSDAQuantResultKEGGIDs, file = "OPLSDAQuantResultKEGGIDs.RData")

OPLSDAQualResultTab <- cbind(dictionary$Consensus[match(OPLSDAQualResultKEGGIDs, dictionary$`KEGG IDs`)], OPLSDAQualResultKEGGIDs)
OPLSDAQuantResultTab <- cbind(dictionary$Consensus[match(OPLSDAQuantResultKEGGIDs, dictionary$`KEGG IDs`)], OPLSDAQuantResultKEGGIDs)
write.csv(OPLSDAQualResultTab, file = "OPLSDAQualResultTab.csv")
write.csv(OPLSDAQuantResultTab, file = "OPLSDAQuantResultTab.csv")

ORA_OPLSDAQual <- doORA(OPLSDAQualResultKEGGIDs, allMets, org = "pae")
ORA_OPLSDAQuant <- doORA(OPLSDAQuantResultKEGGIDs, allMets, org = "pae")

write.csv(ORA_OPLSDAQual, file = "ORA_OPLSDAQual.csv")
write.csv(ORA_OPLSDAQuant, file = "ORA_OPLSDAQuant.csv")
save(ORA_OPLSDAQual, file = "ORA_OPLSDAQual.RData")
save(ORA_OPLSDAQuant, file = "ORA_OPLSDAQuant.RData")

analysis.OPLSDAQual <- enrich(
        compounds = OPLSDAQualResultKEGGIDs,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.OPLSDAQual %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.OPLSDAQual)

g_OPLSDAQual <- generateResultsGraph(
        object = analysis.OPLSDAQual,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_OPLSDAQual

tiff("FELLA_OPLSDAQual.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_OPLSDAQual
  #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_OPLSDAQual <- generateResultsTable(
        object = analysis.OPLSDAQual,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_OPLSDAQual, file = "tab_OPLSDAQual.csv")
save(tab_OPLSDAQual, file = "tab_OPLSDAQual.RData")

analysis.OPLSDAQuant <- enrich(
        compounds = OPLSDAQuantResultKEGGIDs,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.OPLSDAQuant %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.OPLSDAQuant)

g_OPLSDAQuant <- generateResultsGraph(
        object = analysis.OPLSDAQuant,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_OPLSDAQuant

tiff("FELLA_OPLSDAQuant.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_OPLSDAQuant
  #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_OPLSDAQuant <- generateResultsTable(
        object = analysis.OPLSDAQuant,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

save(tab_OPLSDAQuant, file = "tab_OPLSDAQuant.RData")

write.csv(tab_OPLSDAQuant, file = "tab_OPLSDAQuant.csv")

grep("*", dictionary$`KEGG IDs`[!is.na(dictionary$`KEGG IDs`)], fixed = T)

dictionary$`KEGG IDs`[!is.na(dictionary$`KEGG IDs`)][-grep("*", dictionary$`KEGG IDs`[!is.na(dictionary$`KEGG IDs`)], fixed = T)]

# See overlaps
if(!require(Vennerable)) install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library(Vennerable)

overlapPathsORA <- list("ORA RFE" = ORA_rfeResults[, 1], 
                        "ORA OPLS-DA quantitative" = ORA_OPLSDAQuant[, 1], 
                        "ORA OPLS-DA qualitative" = ORA_OPLSDAQuant[, 1])

vennPathsORA <- Venn(Sets = overlapPathsORA)
tiff(filename = "overlapPathsORA.tiff", height = 1400, width = 1800, res = 300)
plot(vennPathsORA, doWeights = T, type = "circles")
dev.off()

overlapORA <- Reduce(intersect, overlapPathsORA)
save(overlapORA, file = "overlapORA.RData")

tab_rfeSwarm[grep("pathway", tab_rfeSwarm[, 2]), 3]

overlapPathsFELLA <- list("FELLA RFE" = tab_rfeSwarm[grep("pathway", tab_rfeSwarm[, 2]), 3], 
                          "FELLA OPLS-DA quantitative" = tab_OPLSDAQuant[grep("pathway", tab_OPLSDAQuant[, 2]), 3], 
                          "FELLA OPLS-DA qualitative" = tab_OPLSDAQual[grep("pathway", tab_OPLSDAQual[, 2]), 3])

vennPathsFELLA <- Venn(Sets = overlapPathsFELLA)
tiff(filename = "overlapPathsFELLA.tiff", height = 1400, width = 1800, res = 300)
plot(vennPathsFELLA, doWeights = T, type = "circles")
dev.off()

overlapFELLA <- Reduce(intersect, overlapPathsFELLA)
save(overlapFELLA, file = "overlapFELLA.RData")

tab_rfeSwarm[grep("pathway", tab_rfeSwarm[, 2]),]
tab_OPLSDAQuant[grep("pathway", tab_OPLSDAQuant[, 2]),]
tab_OPLSDAQual[grep("pathway", tab_OPLSDAQual[, 2]),]


pae00630Genes <- keggLink(source = c("map00630", "map00630"), target = "enzyme")
pae00630Genes <- gsub("ec:", pae00630Genes, replacement = "")

tab_rfeSwarm[tab_rfeSwarm$KEGG.id %in% pae00630Genes, ]

# both ORA and FELLA seem to output same pathways when using random metabolites. Let's do some random iterations to 
# see what pathways get enriched.

sampleMetsRFE <- list()
for(i in 1:1000){
      sampleMetsRFE[[i]] <- sample(allMets, length(rfeResultKEGGIDs))  
}

sampleMetsOPLSDA <- list()
for(i in 1:1000){
        sampleMetsOPLSDA[[i]] <- sample(allMets, length(OPLSDAQuantResultKEGGIDs))  
}

# Modify doORA function to speed it up (we need to do 1000 random iterations)

paePaths <- keggList("pathway", "pae")
totPathsInMets <- unique(unlist(sapply(allMets, keggLink, target = "pathway")))
totPathsInMets <- totPathsInMets[gsub("map", replacement = "pae", totPathsInMets) %in% names(paePaths)]
compsPerPathAllPaths <- sapply(totPathsInMets, keggLink, target = "compound")

randORARFE <- lapply(sampleMetsRFE, doORAMod, alpha = 2)
randORAOPLSDA <- lapply(sampleMetsOPLSDA, doORAMod, alpha = 2)
ORARFEReal <- doORA(rfeResultKEGGIDs, alpha = 2, allMetsObjkt = allMets, org = "pae")
ORAOPLSDAReal <- doORA(OPLSDAQuantResult, alpha = 2, allMetsObjkt = allMets, org = "pae")

getRandEnrichment <- function(randTest, realTest){
        pathsRandPermPval <- c()
        for(i in 1:nrow(realTest)){
                pathRandPvals <- sapply(randTest, function(x) x$p.values[i])
                pathPVal <- sum(pathRandPvals <= realTest$p.values[i])/length(randTest)
                pathsRandPermPval <- c(pathsRandPermPval, pathPVal)
        }
        names(pathsRandPermPval) <- realTest$Pathways
        return(pathsRandPermPval)
}

pathsRandPermPvalORARFE <- getRandEnrichment(randORARFE, ORARFEReal)
pathsRandPermPvalORAOPLSDA <- getRandEnrichment(randORAOPLSDA, ORAOPLSDAReal)

# with random values the probability of getting the same enriched pathways than with ours is very small.
# all pvalues much bigger than 0.05


####################################################################################################################
####################################################################################################################
####################################################################################################################
# With genomic data 
####################################################################################################################
####################################################################################################################
####################################################################################################################
library(gplots)

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_enz_tab_filt.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_tab_filt.RData")

rownames(gene_enz_tab_filt) <- gsub("PA14.jbx", "PA14", rownames(gene_enz_tab_filt))
rownames(gene_tab_filt) <- gsub("PA14.jbx", "PA14", rownames(gene_tab_filt))

# Group strains with a Jaccard similarity of 0.95 or above
geneEnzTabFiltGrouped <- jaccGroup(gene_enz_tab_filt, threshold = 0.05)$GroupedMat
geneGroups <- jaccGroup(gene_enz_tab_filt)$Groups
rownames(geneEnzTabFiltGrouped)[rownames(geneEnzTabFiltGrouped) == "W70322"] <- "W70332"
geneEnzTabFiltGrouped <- geneEnzTabFiltGrouped[rownames(geneEnzTabFiltGrouped) %in% names(swarmBin), ]

geneTabFiltGrouped <- jaccGroup(gene_tab_filt, threshold = 0.05)$GroupedMat
geneAllGroups <- jaccGroup(gene_tab_filt)$Groups
rownames(geneTabFiltGrouped)[rownames(geneTabFiltGrouped) == "W70322"] <- "W70332"
geneTabFiltGrouped <- geneTabFiltGrouped[rownames(geneTabFiltGrouped) %in% names(swarmBin), ]

geneEnzTabFiltGrouped <- cbind.data.frame(as.data.frame(geneEnzTabFiltGrouped), 
                                          swarmBin[names(swarmBin) %in% rownames(geneEnzTabFiltGrouped)])
colnames(geneEnzTabFiltGrouped)[ncol(geneEnzTabFiltGrouped)] <- "swarmData"

geneTabFiltGrouped <- cbind.data.frame(as.data.frame(geneTabFiltGrouped), 
                                       swarmBin[names(swarmBin) %in% rownames(geneTabFiltGrouped)])
colnames(geneTabFiltGrouped)[ncol(geneTabFiltGrouped)] <- "swarmData"

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

# Ungroup first 10 more important genes (according to RFE) grouped in jaccard
predGenes <- c()
for(i in 1:10){
        predGenes <- c(predGenes, names(which(geneGroups == 
                                                geneGroups[which(names(geneGroups) == 
                                                                   predictors(resultsEnz)[i])])))
}

rownames(gene_enz_tab_filt)[grep("W70322", rownames(gene_enz_tab_filt))] <- "W70332"

# Do a submatrix of GeneEnzTab with the unpacked 10 most important group according to RFE
geneEnzTabSign <- gene_enz_tab_filt[match(names(sort(swarmMeansRearr)), rownames(gene_enz_tab_filt)), 
                                    which(make.names(colnames(gene_enz_tab_filt)) %in% predGenes)]

geneEnzTabSign <- geneEnzTabSign[rownames(geneEnzTabSign) != "NA", ]

# Do heatmap of the submatrix

cols <- topo.colors(nrow(gene_enz_tab_filt) + 2)
cols <- cols[-c(1, 2)]#, grep("W70332", rownames(gene_enz_tab_filt)))]


tiff("geneEnzTabSign.tiff", width = 10000, height = 5000, units = "px", pointsize = 100)
heatmap.2(t(as.matrix(geneEnzTabSign)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(rownames(geneEnzTabSign), rownames(gene_enz_tab_filt))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of genes related to swarming (RFE)", margins = c(8, 60), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          cexRow = 1.75,
          cexCol = 1.75,
          
          #sepcolor = "black",
          #colsep=1:ncol(t(gene_tab)),
          #rowsep=1:nrow(t(gene_tab)),
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
dev.off()



geneEnzTabFiltGrouped <- cbind.data.frame(geneEnzTabFiltGrouped, swarmMeansRearr[match(rownames(geneEnzTabFiltGrouped), names(swarmMeansRearr))])
colnames(geneEnzTabFiltGrouped)[ncol(geneEnzTabFiltGrouped)] <- "swarmQuant"

geneTabFiltGrouped <- cbind.data.frame(geneTabFiltGrouped, swarmMeansRearr[match(rownames(geneTabFiltGrouped), names(swarmMeansRearr))])
colnames(geneTabFiltGrouped)[ncol(geneTabFiltGrouped)] <- "swarmQuant"

geneEnzTabPCA <- opls(geneEnzTabFiltGrouped[, 1:(ncol(geneEnzTabFiltGrouped) - 2)])
geneTabPCA <- opls(geneTabFiltGrouped[, 1:(ncol(geneTabFiltGrouped) - 2)])

tiff(filename = "pcaEnzDataSwarmQual.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(geneEnzTabPCA,
     typeVc = "x-score",
     parAsColFcVn = geneEnzTabFiltGrouped$swarmData)
dev.off()

tiff(filename = "pcaEnzDataSwarmQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(geneEnzTabPCA,
     typeVc = "x-score",
     parAsColFcVn = geneEnzTabFiltGrouped$swarmQuant)
dev.off()

jpeg(filename = "pcaEnzDataSwarmQuant.jpeg", res = 300, height = 3000, width = 3000, units = "px")
plot(geneEnzTabPCA,
     typeVc = "x-score",
     parAsColFcVn = geneEnzTabFiltGrouped$swarmQuant)
dev.off()

tiff(filename = "pcaGenesSwarmQual.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(geneTabPCA,
     typeVc = "x-score",
     parAsColFcVn = geneTabFiltGrouped$swarmData)
dev.off()

tiff(filename = "pcaGenesSwarmQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(geneTabPCA,
     typeVc = "x-score",
     parAsColFcVn = geneTabFiltGrouped$swarmQuant)
dev.off()

tiff(filename = "geneEnzTabOPLSDAQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
geneEnzTabOPLSDAQuant <- opls(geneEnzTabFiltGrouped[, 1:(ncol(geneEnzTabFiltGrouped) - 2)],
                              geneEnzTabFiltGrouped$swarmQuant, 
                              predI = 1, 
                              orthoI = NA#,
                            #subset = as.vector(inTrain)
)
dev.off()

jpeg(filename = "geneEnzTabOPLSDAQuant.jpeg", res = 300, height = 3000, width = 3000, units = "px")
geneEnzTabOPLSDAQuant <- opls(geneEnzTabFiltGrouped[, 1:(ncol(geneEnzTabFiltGrouped) - 2)],
                              geneEnzTabFiltGrouped$swarmQuant, 
                              predI = 1, 
                              orthoI = NA#,
                              #subset = as.vector(inTrain)
)
dev.off()

tiff(filename = "geneTabOPLSDAQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
geneTabOPLSDAQuant <- opls(geneTabFiltGrouped[, 1:(ncol(geneTabFiltGrouped) - 2)],
                           geneTabFiltGrouped$swarmQuant, 
                           predI = 1, 
                           orthoI = NA#,
                           #subset = as.vector(inTrain)
)
dev.off()

geneEnzTabOPLSDAQuantLoads <- cbind(getLoadingMN(geneEnzTabOPLSDAQuant)[, 1]#,
                                    #getLoadingMN(geneEnzTabOPLSDAQuant, orthoL = T)
                                    )
colnames(geneEnzTabOPLSDAQuantLoads)[1] <- "p1"

geneEnzTabExtremLoads <- getExtremVals(geneEnzTabOPLSDAQuantLoads, n = 5)


geneOPLSDALoads4Barplot <- cbind.data.frame(rownames(geneEnzTabOPLSDAQuantLoads), geneEnzTabOPLSDAQuantLoads[, 1])
colnames(geneOPLSDALoads4Barplot) <- c("gene", "loading")

ggplot(data = geneOPLSDALoads4Barplot,
       aes(x = reorder(gene, loading), y = loading)) +
  labs(x = "gene", y = "OPLS-DA Loading") +
  geom_bar(stat = "identity") +
  coord_flip()
ggsave(filename = "enzOPLSDAQuantLoadsBarplot.tiff", width = 500, height = 1000, units = "mm")


predGenesOPLSDA <- c()
for(i in seq_along(geneEnzTabExtremLoads$uniqueExtremeVars)){
  predGenesOPLSDA <- c(predGenesOPLSDA, names(which(geneGroups == 
                                                      geneGroups[which(names(geneGroups) == 
                                                                         geneEnzTabExtremLoads$uniqueExtremeVars[i])])))
}
predGenesOPLSDA

# Do a submatrix of geneEnzTab of the 10 groups with most extreme loading values of p1 from OPLS-DA

geneEnzTabSignOPLSDA <- gene_enz_tab_filt[match(names(sort(swarmMeansRearr)), rownames(gene_enz_tab_filt)), 
                                          which(make.names(colnames(gene_enz_tab_filt)) %in% predGenesOPLSDA)]

geneEnzTabSignOPLSDA <- geneEnzTabSignOPLSDA[rownames(geneEnzTabSignOPLSDA) != "NA", ]

save(geneEnzTabSignOPLSDA, file = "geneEnzTabSignOPLSDA.RData")

# Plot it

tiff("geneEnzTabSignOPLSDA.tiff", width = 10000, height = 20000, units = "px", pointsize = 100)
heatmap.2(t(as.matrix(geneEnzTabSignOPLSDA)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(rownames(geneEnzTabSignOPLSDA), rownames(gene_enz_tab_filt))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of genes related to swarming (OPLS-DA)", margins = c(8, 60), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          cexRow = 1.75,
          cexCol = 1.75,
          
          #sepcolor = "black",
          #colsep=1:ncol(t(gene_tab)),
          #rowsep=1:nrow(t(gene_tab)),
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
dev.off()

jpeg("geneEnzTabSignOPLSDA.jpeg", width = 2200, height = 3000, units = "px", pointsize = 100)
heatmap.2(t(as.matrix(geneEnzTabSignOPLSDA)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(rownames(geneEnzTabSignOPLSDA), rownames(gene_enz_tab_filt))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of genes related to swarming (OPLS-DA)", margins = c(5, 6), 
          keysize = 0.5,
          sepwidth = c(0.1, 0.05),
          cexRow = 0.5,
          cexCol = 0.5,
          
          #sepcolor = "black",
          #colsep=1:ncol(t(gene_tab)),
          #rowsep=1:nrow(t(gene_tab)),
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })
dev.off()

geneEnzTabSignOPLSDAHM <- heatmap.2(t(as.matrix(geneEnzTabSignOPLSDA)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
                                    density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
                                    col = c("blue", "red"), ColSideColors = cols[match(rownames(geneEnzTabSignOPLSDA), rownames(gene_enz_tab_filt))], notecol = NULL, trace = "none", xlab = "Strains", 
                                    ylab = "Genes", main = "Presence/absence of genes related to swarming (OPLS-DA)", margins = c(8, 30), 
                                    keysize = 1,
                                    sepwidth = c(0.1, 0.05),
                                    cexRow = 1.75,
                                    cexCol = 1.75,
                                    
                                    #sepcolor = "black",
                                    #colsep=1:ncol(t(gene_tab)),
                                    #rowsep=1:nrow(t(gene_tab)),
                                    key.xtickfun=function() {
                                      cex <- par("cex")*par("cex.axis")
                                      side <- 1
                                      line <- 0
                                      col <- par("col.axis")
                                      font <- par("font.axis")
                                      mtext("low", side=side, at=0, adj=0,
                                            line=line, cex=cex, col=col, font=font)
                                      mtext("high", side=side, at=1, adj=1,
                                            line=line, cex=cex, col=col, font=font)
                                      return(list(labels=FALSE, tick=FALSE))
                                    })

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")

# Let's see if any of these enzymatic genes appears in the fella subgraph

sum(tab_OPLSDAQuant$KEGG.id %in% OPLSDASignEnzsECnums)
tab_OPLSDAQuant[which(tab_OPLSDAQuant$KEGG.id %in% OPLSDASignEnzsECnums), ]

# EC 2.2.1.6 appears. The most relevant pathway where it appears, according to the compounds related with swarming, 
# is Valine Biosynthesis
dictEnzymes[which(dictEnzymes$ECnums == "2.2.1.6"), ]

valBiosComps <- gsub("cpd:", replacement = "", keggLink(target = "compound", source = "map00290"))
OPLSDAQuantResultTab[OPLSDAQuantResultTab[, 2] %in% valBiosComps, ]

OPLSDASignEnzsAnnot <- dictEnzymes$Annotation[match(colnames(geneEnzTabSignOPLSDA), dictEnzymes$Gene)]

signEnzymes <- cbind.data.frame(OPLSDASignEnzsAnnot, 
                                colnames(geneEnzTabSignOPLSDA), 
                                OPLSDASignEnzsECnums)

colnames(signEnzymes) <- c("Annotation", "Gene", "EC_num")


# Let's see if the genes that are in F34365 are contiguous in the genome
sort(dictEnzymes$F34365[dictEnzymes$Gene %in% signEnzymes$Gene])

# Do keggLink for each gene to see pathways were each gene appears. 
# Can we do an ORA? 
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

oraEnzs <- doORA(signEnzymes$EC_num, dictEnzymes$ECnums, org = "pae", alpha = 0.05, target = "enzyme")

tyrMetComps <- gsub("cpd:", replacement = "", keggLink("map00350", target = "compound"))

OPLSDAQuantResultTab[, 2][OPLSDAQuantResultTab[, 2] %in% tyrMetComps]

signEnzymes$EC_num[signEnzymes$EC_num %in% gsub("ec:", replacement = "", keggLink("map00350", target = "enzyme"))]
signEnzymes[signEnzymes$EC_num %in% gsub("ec:", replacement = "", keggLink("map00350", target = "enzyme")), ]

# Formyl methionine is the metabolite that is most correlated with swarming. It only appears in Cysteine and methionine 
# metabolism. Let's see if we have any gene in this pathway correlated with swarming.
cysMethMet <- gsub("ec:", replacement = "", keggLink("enzyme", "map00270"))
sum(signEnzymes$EC_num %in% cysMethMet)
dictEnzymes$Gene[which(dictEnzymes$ECnums %in% cysMethMet)]
signEnzymes[which(signEnzymes$EC_num %in% cysMethMet), ]
# Now without doing the gene grouping 

geneEnzTabFiltWSwarm <- cbind.data.frame(gene_enz_tab_filt, swarmMeansRearr[match(rownames(gene_enz_tab_filt), names(swarmMeansRearr))])
colnames(geneEnzTabFiltWSwarm)[ncol(geneEnzTabFiltWSwarm)] <- "swarmQuant"
geneEnzTabFiltWSwarm <- geneEnzTabFiltWSwarm[!is.na(geneEnzTabFiltWSwarm$swarmQuant), ]
dim(geneEnzTabFiltWSwarm)

tiff(filename = "geneEnzTabWOGroupOPLSDAQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
geneEnzTabWOGroupOPLSDAQuant <- opls(geneEnzTabFiltWSwarm[, 1:(ncol(geneEnzTabFiltWSwarm) - 1)],
                                     geneEnzTabFiltWSwarm$swarmQuant, 
                                     predI = 1, 
                                     orthoI = NA
)
dev.off()

geneEnzTabWOGroupOPLSDAQuantLoads <- getLoadingMN(geneEnzTabWOGroupOPLSDAQuant)


geneEnzTabWOGroupExtremLoads <- getExtremVals(geneEnzTabWOGroupOPLSDAQuantLoads, n = 5)


geneOPLSDALoads4BarplotWOGroup <- cbind.data.frame(rownames(geneEnzTabWOGroupOPLSDAQuantLoads), geneEnzTabWOGroupOPLSDAQuantLoads[, 1])
colnames(geneOPLSDALoads4BarplotWOGroup) <- c("gene", "loading")

ggplot(data = geneOPLSDALoads4BarplotWOGroup,
       aes(x = reorder(gene, loading), y = loading)) +
  labs(x = "gene", y = "OPLS-DA Loading") +
  geom_bar(stat = "identity") +
  coord_flip()
ggsave(filename = "enzWOGroupOPLSDAQuantLoadsBarplot.tiff", width = 500, height = 1000, units = "mm")

# Doesn't work, too many variables, turn the loadings very low. 

# Do OPLS-DA mixing metabolomic and presence/absence data
rownames(gene_enz_tab_filt)
rownames(ccmn_norm_mets_good_old)
mixedDat <- cbind.data.frame(ccmn_norm_mets_good_old, 
                             gene_enz_tab_filt[match(gsub("\\_.*|(PA14).*", 
                                                          rownames(ccmn_norm_mets_good_old), 
                                                          rep = "\\1"), 
                                                     rownames(gene_enz_tab_filt)), ])
mixedDat <- mixedDat[!is.na(mixedDat[, ncol(mixedDat)]), ]
mixedDat <- cbind.data.frame(mixedDat, swarmMeansRearr[match(gsub("\\_.*|(PA14).*", 
                                                                  rep = "\\1",
                                                                  rownames(mixedDat)), names(swarmMeansRearr))])

colnames(mixedDat)[ncol(mixedDat)] <- "swarmQuant"

mixedDat <- mixedDat[, -87]
mixedDat <- mixedDat[, -87]

mixedDatPCA  <- opls(mixedDat[, 1:(ncol(mixedDat) - 1)])
swarmData <- mixedDat$swarmQuant
tiff(filename = "pcamixedSwarmQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
plot(mixedDatPCA,
     typeVc = "x-score",
     parAsColFcVn = swarmData)
dev.off()

tiff(filename = "mixedOPLSDAQuant.tiff", res = 300, height = 3000, width = 3000, units = "px")
mixedOPLSDAQuant <- opls(mixedDat[, 1:(ncol(mixedDat) - 1)],
                         mixedDat$swarmQuant, 
                         predI = 1, 
                         orthoI = NA)
dev.off()

mixedLoadings <- getLoadingMN(mixedOPLSDAQuant)


# Valine, Leucine and Isoleucine Biosynthesis seems to be affected. Gene ilvB3 (codes for Acetolactate isozyme, EC 2.2.1.6)
# Is only in F34365, F22031 (swarmers) and T6313 (nonswarmer). Let's see if using only the metabolite abundances of
# the compounds that appear in this pathway the strains cluser according to how much they swarm. 
ccmn_norm_mets_good_old

ValLeuIleBioCpds <- gsub("cpd:", replacement = "", keggLink("compound", "map00290"))

metMat <- ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 2)]
dictionary$`KEGG IDs`[match(colnames(metMat), make.names(dictionary$`Old Data Names`))]

colnames(metMat) <- dictionary$`KEGG IDs`[match(colnames(metMat), make.names(dictionary$`Old Data Names`))]

valMetMat <- metMat[, colnames(metMat) %in% ValLeuIleBioCpds]


dictionary$Consensus[match(colnames(valMetMat), dictionary$`KEGG IDs`)]

groups <- unique(gsub("_.*", replacement = "", rownames(valMetMat)))
cols <- topo.colors(length(groups))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(valMetMat))))
for(ii in 1:nrow(valMetMat)) {
        selected <- which(groups == gsub("\\_.*", "", rownames(valMetMat))[ii])
        Cols[ii] <- cols[selected]
}
valMetMat <- getStrainMedian(valMetMat)
Cols <- unique(Cols)
tiff("heatmapValBiosynthMets.tiff", width = 3000, height = 1000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(valMetMat[, c("C00183", "C02504", "C02631", "C06032", "C02226")])), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          scale = "row",
          #Colv = colThick,
          #Rowv = rowThick,
          #colCol = colCols,
          cellnote = round(as.matrix(t(ccmnNormMets)), 2),
          notecex = 0.7,
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })

dev.off()

jpeg("heatmapValBiosynthMets.jpeg", width = 3000, height = 1000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(valMetMat[, c("C00183", "C02504", "C02631", "C06032", "C02226")])), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          scale = "none",
          #Colv = colThick,
          #Rowv = rowThick,
          #colCol = colCols,
          cellnote = round(as.matrix(t(valMetMat[, c("C00183", "C02504", "C02631", "C06032", "C02226")])), 2),
          notecex = 0.7,
          key.xtickfun=function() {
            cex <- par("cex")*par("cex.axis")
            side <- 1
            line <- 0
            col <- par("col.axis")
            font <- par("font.axis")
            mtext("low", side=side, at=0, adj=0,
                  line=line, cex=cex, col=col, font=font)
            mtext("high", side=side, at=1, adj=1,
                  line=line, cex=cex, col=col, font=font)
            return(list(labels=FALSE, tick=FALSE))
          })

dev.off()

pcaValBiosynthMets <- prcomp(valMetMat)

tiff("pcaSwarmMets.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaValBiosynthMets, 
             col.ind = c(ccmn_norm_mets_good_old$swarmQuant[1],
                         unique(ccmn_norm_mets_good_old$swarmQuant)),
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Log Swarming Area")
dev.off()

