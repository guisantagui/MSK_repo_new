setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/decTreesLogReg")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes gene_tab.csv (data of presence/absence of genes) generated in gene_presAbs.R script and trains     #################################
# decision trees and logistic regression models in order to select the features that are more useful for classifying   #################################
# the strains in the 4 major clusters that were observed during the differential metabolite analysis                   #################################
# (diffMet_analysis.R).                                                                                                #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

if (!require(tree)) install.packages('tree')
library(tree)


#Load data

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_tab.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_tab_filt.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/gene_enz_tab_filt.RData")

#filter gene_tab to get only genes that doesn't appear in all strains (remove core), and to remove PA14, because it is outside
#all clusters

genes_strains_filt <- gene_enz_tab_filt
genes_strains_filt <- sapply(genes_strains_filt, as.integer)
rownames(genes_strains_filt) <- rownames(gene_tab_filt)
colnames(genes_strains_filt) <- make.names(colnames(genes_strains_filt), unique = T)

if(!require(proxy)) install.packages('proxy')
library(proxy)

jacc_mat_genes_dist <- proxy::dist(genes_strains_filt, by_rows = F, method = "Jaccard", diag = T)

hca_jaccMat_genes <- hclust(jacc_mat_genes_dist, method = "complete")

tiff("hca_jaccGenes.tiff", height = 15000, width = 15000)
plot(hca_jaccMat_genes)
dev.off()

#Stablish a threshold of 0.95 similarity and group genes that are as similar as it. Those groups will be identified by the first member 
#the group
hca_jaccMat_group_genes <- cutree(hca_jaccMat_genes, h = 0.05)

gene_groups <- c()
for(i in 1:max(hca_jaccMat_group_genes)){
        gene_groups <- c(gene_groups, names(sort(hca_jaccMat_group_genes)[which(sort(hca_jaccMat_group_genes) == i)])[1])
}

genes_strains_filt_grouped <- genes_strains_filt[, gene_groups]

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/metClusts_oldGood.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/metClusts_new.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/metClusts_hyb.RData")

geneTabFilt_grouped_oldG <- cbind.data.frame(genes_strains_filt_grouped, Clusters_4 = metClusts_oldGood[-1])
geneTabFilt_grouped_oldG <- cbind.data.frame(geneTabFilt_grouped_oldG, 
                                             Clusters_2 = as.factor(gsub("2.*", 
                                                                         "Clust.2", 
                                                                         gsub("1.*", 
                                                                              "Clust.1", 
                                                                              geneTabFilt_grouped_oldG[, "Clusters_4"]))))

save(geneTabFilt_grouped_oldG, file = "geneTabFilt_grouped_oldG.RData")
write.csv(geneTabFilt_grouped_oldG, file = "geneTabFilt_grouped_oldG.csv")

geneTabFilt_grouped_new <- cbind.data.frame(genes_strains_filt_grouped, Clusters_4 = metClusts_new[-1])
geneTabFilt_grouped_new <- cbind.data.frame(geneTabFilt_grouped_new, 
                                            Clusters_2 = as.factor(gsub("2.*", 
                                                                        "Clust.2", 
                                                                        gsub("1.*", 
                                                                             "Clust.1", 
                                                                             geneTabFilt_grouped_new[, "Clusters_4"]))))

save(geneTabFilt_grouped_new, file = "geneTabFilt_grouped_new.RData")
write.csv(geneTabFilt_grouped_new, file = "geneTabFilt_grouped_new.csv")

geneTabFilt_grouped_hyb <- cbind.data.frame(genes_strains_filt_grouped, Clusters_4 = metClusts_hyb[-1])
geneTabFilt_grouped_hyb <- cbind.data.frame(geneTabFilt_grouped_hyb, 
                                            Clusters_2 = as.factor(gsub("2.*", 
                                                                        "Clust.2", 
                                                                        gsub("1.*", 
                                                                             "Clust.1", 
                                                                             geneTabFilt_grouped_hyb[, "Clusters_4"]))))

save(geneTabFilt_grouped_hyb, file = "geneTabFilt_grouped_hyb.RData")
write.csv(geneTabFilt_grouped_hyb, file = "geneTabFilt_grouped_hyb.csv")

########################################################################################################################################################
#
# Decision trees
#
########################################################################################################################################################

####
# Old Data
########################################################################################################################################################

tree.genes_old = tree(Clusters_2 ~ . - Clusters_4, data = geneTabFilt_grouped_oldG)
sum <- summary(tree.genes_old)
plot(tree.genes_old)
text(tree.genes_old, pretty = 0)
cv.genes_old = cv.tree(tree.genes_old, FUN = prune.misclass)
plot(cv.genes_old)

names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[1])]))
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[2])]))
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[3])]))

tree.genes2_old = tree(Clusters_4 ~ . - Clusters_2, data = geneTabFilt_grouped_oldG)
sum <- summary(tree.genes2_old)
plot(tree.genes2_old)
text(tree.genes2_old, pretty = 0)

names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[1])]))
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[2])]))
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          sum$used[3])]))

set.seed(101)
train_old=sample(1:nrow(geneTabFilt_grouped_oldG), 15)

tree.genes_old = tree(Clusters_2~.-Clusters_4, geneTabFilt_grouped_oldG, subset=train_old)
plot(tree.genes_old)
text(tree.genes_old, pretty=0)

tree.pred_old = predict(tree.genes_old, geneTabFilt_grouped_oldG[-train_old,], type="class")

with(geneTabFilt_grouped_oldG[-train_old,], table(tree.pred_old, Clusters_2))
5/11

####
# New Data
########################################################################################################################################################

tree.genes_new = tree(Clusters_2 ~ . - Clusters_4, data = geneTabFilt_grouped_new)
summary(tree.genes_new)
plot(tree.genes_new)
text(tree.genes_new, pretty = 0)
cv.genes_new = cv.tree(tree.genes_new, FUN = prune.misclass)
plot(cv.genes_new)

tree.genes2_new = tree(Clusters_4 ~ . - Clusters_2, data = geneTabFilt_grouped_new)
summary(tree.genes2_new)
plot(tree.genes2_new)
text(tree.genes2_new, pretty = 0)

set.seed(101)
train_new=sample(1:nrow(geneTabFilt_grouped_new), 15)

tree.genes_new = tree(Clusters_2~.-Clusters_4, geneTabFilt_grouped_new, subset=train_new)
plot(tree.genes_new)
text(tree.genes_new, pretty=0)

tree.pred_new = predict(tree.genes_new, geneTabFilt_grouped_new[-train_new,], type="class")

with(geneTabFilt_grouped_new[-train_new,], table(tree.pred_new, Clusters_2))

####
# Hybrid Data
########################################################################################################################################################

tree.genes_hyb = tree(Clusters_2 ~ . - Clusters_4, data = geneTabFilt_grouped_hyb)
summary(tree.genes_hyb)
plot(tree.genes_hyb)
text(tree.genes_hyb, pretty = 0)
cv.genes_hyb = cv.tree(tree.genes_hyb, FUN = prune.misclass)
plot(cv.genes_hyb)

tree.genes2_hyb = tree(Clusters_4 ~ . - Clusters_2, data = geneTabFilt_grouped_hyb)
summary(tree.genes2_hyb)
plot(tree.genes2_hyb)
text(tree.genes2_hyb, pretty = 0)

set.seed(101)
train_hyb=sample(1:nrow(geneTabFilt_grouped_hyb), 15)

tree.genes_hyb = tree(Clusters_2~.-Clusters_4, geneTabFilt_grouped_hyb, subset=train_hyb)
plot(tree.genes_hyb)
text(tree.genes_hyb, pretty=0)

tree.pred_hyb = predict(tree.genes_hyb, geneTabFilt_grouped_hyb[-train_hyb,], type="class")

with(geneTabFilt_grouped_hyb[-train_hyb,], table(tree.pred_hyb, Clusters_2))


########################################################################################################################################################
#
# Logistic Regression
#
########################################################################################################################################################

glm.genes.fit <- glm(Clusters_2 ~. -Clusters_4, data = geneTabFilt_grouped_oldG, family = binomial)
summary(glm.genes.fit)

glm.genes.probs <- predict(glm.genes.fit,type = "response")
glm.genes.probs[1:5]
glm.genes.pred <- ifelse(glm.genes.probs > 0.5, "Clust.2", "Clust.1")
table(glm.genes.pred, geneTabFilt_grouped_oldG$Clusters_2)
mean(glm.genes.pred == geneTabFilt_grouped_oldG$Clusters_2)

#Assigns correctly each strain to its cluster 100% of times when using all data
#Use a training set to see if there is overfitting
train <- rep(F, dim(geneTabFilt_grouped_oldG)[1])
train[sample(1:nrow(geneTabFilt_grouped_oldG), 15)] <- T

glm.genes.fit <- glm(Clusters_2 ~.-Clusters_4, 
                     data = geneTabFilt_grouped_oldG, 
                     family = binomial, 
                     subset = train)

glm.genes.probs <- predict(glm.genes.fit, 
                           newdata = geneTabFilt_grouped_oldG[!train,], 
                           type = "response")

glm.genes.pred <- ifelse(glm.genes.probs > 0.5, "Clust.2", "Clust.1")
table(glm.genes.pred, geneTabFilt_grouped_oldG$Clusters_2[!train])
mean(glm.genes.pred == geneTabFilt_grouped_oldG$Clusters_2[!train])
#There is overfitting, so we select a smaller number of variables

glm.genes.fit <- glm(Clusters_2 ~ hypothetical.protein.APECO1_2271 + Retron.type.RNA.directed.DNA.polymerase..EC.2.7.7.49. + 
                             Type.I.restriction.modification.system..restriction.subunit.R..EC.3.1.21.3., 
                     data = geneTabFilt_grouped_oldG, 
                     family = binomial, 
                     subset = train) 
## choose hypothetical.protein.APECO1_2271, Retron.type.RNA.directed.DNA.polymerase..EC.2.7.7.49. and 
## Type.I.restriction.modification.system..restriction.subunit.R..EC.3.1.21.3. because according to decision 
## trees those where the ones that better classify the strains
## in major clusters 1 and 2.

glm.genes.probs <- predict(glm.genes.fit, 
                           newdata = geneTabFilt_grouped_oldG[!train,], 
                           type = "response")
glm.genes.pred <- ifelse(glm.genes.probs > 0.5, "Clust.2", "Clust.1")
table(glm.genes.pred, geneTabFilt_grouped_oldG$Clusters_2[!train])
mean(glm.genes.pred == geneTabFilt_grouped_oldG$Clusters_2[!train]) ##-->> With hypothetical.protein.APECO1_2271,
## Retron.type.RNA.directed.DNA.polymerase..EC.2.7.7.49. and
## Type.I.restriction.modification.system..restriction.subunit.R..EC.3.1.21.3. we can predict to which of the 2 major 
## clusters a strain belongs with a high probability.

#See which genes are in those 3 groups of genes
#hypothetical.protein.APECO1_2271
hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == "hypothetical.protein.APECO1_2271")]
which(hca_jaccMat_group_genes == 219) #-->Only this gene

#Retron.type.RNA.directed.DNA.polymerase..EC.2.7.7.49.
hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == "Retron.type.RNA.directed.DNA.polymerase..EC.2.7.7.49.")]
which(hca_jaccMat_group_genes == 370) #-->Only this gene

#Type.I.restriction.modification.system..restriction.subunit.R..EC.3.1.21.3.
hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == "Type.I.restriction.modification.system..restriction.subunit.R..EC.3.1.21.3.")]
which(hca_jaccMat_group_genes == 188) #-->Only this gene


########################################################################################################################################################
#
# Now we repeat logistic regression step, but using caret R package, which contains many tools for training logistic regression models
#
########################################################################################################################################################

if (!require(caret)) install.packages('caret')
library(caret)

set.seed(107)
inTrain <- createDataPartition(
        y = geneTabFilt_grouped_oldG$Clusters_2,
        p = 0.5,
        list = F
)
training <- geneTabFilt_grouped_oldG[inTrain,]
testing <- geneTabFilt_grouped_oldG[-inTrain,]
nrow(training)
nrow(testing)

ctrl <- trainControl(method = "repeatedcv", 
                     repeats = 3,
                     classProbs = T,
                     summaryFunction = twoClassSummary)

bstLogFit <- train(
        Clusters_2~.-Clusters_4,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)


bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$Clusters_2)#-->Good accuracy, but it doesn't beat No Info Rate

#test importance of variables

importance <- varImp(bstLogFit, scale=F)
print(importance)
plot(importance)
#repeat with 3 more important variables
bstLogFit <- train(
        Clusters_2~RNA.polymerase.ECF.type.sigma.factor +
                Oxaloacetate.decarboxylase..EC.4.1.1.3. +
                Trans.aconitate.2.methyltransferase..EC.2.1.1.144. -Clusters_4,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$Clusters_2)#-->Good accuracy, but it doesn't beat No Info Rate

##Look if inside of the goups are other genes with biological meaning
for(i in 1:nrow(importance$importance)){
        print(names(which(hca_jaccMat_group_genes == 
                            hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                                  rownames(importance[[1]])[i])])))
}

# Similar to Agmatine deiminase gene is present in strains F34365, H5708, M1608, M37351, W45909, which belong to clusters
# 2.1, 1.2, 2.1, 1.2 & 1.1 respectively.

#Apply recursive feature elimination to evaluate the importance of variables
if(!require(randomForest)) install.packages("randomForest")
library(randomForest)
control <- rfeControl(functions=rfFuncs, method="cv", number=10)

results <- rfe(geneTabFilt_grouped_oldG[, 1:(ncol(geneTabFilt_grouped_oldG) -2)], 
               geneTabFilt_grouped_oldG[, ncol(geneTabFilt_grouped_oldG)], 
               sizes = c(1:(ncol(geneTabFilt_grouped_oldG) -2)), 
               rfeControl = control)

print(results)

predictors(results)

plot(results, type=c("g", "o"))

#Repeat with top5 genes

bstLogFit <- train(
        Clusters_2~Serine..glyoxylate.aminotransferase..EC.2.6.1.45. +
                NAD.P..dependent.glyceraldehyde.3.phosphate.dehydrogenase.archaeal..EC.1.2.1.59. +
                Thioredoxin.2..EC.1.8.1.8. +
                hypothetical.protein.APECO1_2271 +
                Sulfur.carrier.protein.ThiS.adenylyltransferase..EC.2.7.7.73. +
                -Clusters_4,
        data = training,
        method = "glm",
        preProc = c("center", "scale"),
        trControl = ctrl,
        metric = "ROC"
)

bstClasses <- predict(bstLogFit, newdata = testing)
head(bstClasses)

confusionMatrix(data = bstClasses, testing$Clusters_2)#-->Good accuracy, but it doesn't beat No Info Rate
which(colnames(genes_strains_filt_grouped) == "Similar.to.Agmatine.deiminase")

##Look if inside of the goups are other genes with biological meaning
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[1])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[2])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[3])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[4])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[5])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[6])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[7])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[8])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[9])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[10])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[11])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[12])]))#-->No more genes
names(which(hca_jaccMat_group_genes == 
                    hca_jaccMat_group_genes[which(names(hca_jaccMat_group_genes) == 
                                                          predictors(results)[13])]))#-->No more genes

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/presAbsC1_2Genes_old_filtEnz.RData")

mixedSign <- rbind(presAbsC1_2Genes_old_filtEnz, 
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Thioredoxin", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Arsenical pump-driving", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Dihydropteroate", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Aminoglycoside", 
                                                                                         colnames(gene_enz_tab_filt))[1]],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Aminoglycoside", 
                                                                                         colnames(gene_enz_tab_filt))[2]],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Homospermidine", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Mercuric ion reductase", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Sulfur carrier", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("restriction subunit R", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("methyltransferase subunit M", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Methyl-directed", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("Phosphate acetyltransferase (EC 2.3.1.8)", 
                                                                                         colnames(gene_enz_tab_filt))],
                   gene_enz_tab_filt[match(colnames(presAbsC1_2Genes_old_filtEnz), 
                                           rownames(gene_enz_tab_filt)), grep("1,3,6,8", 
                                                                                         colnames(gene_enz_tab_filt))])
rownames(mixedSign)[5:nrow(mixedSign)] <- c(colnames(gene_enz_tab_filt)[grep("Thioredoxin", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Arsenical pump-driving", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Dihydropteroate", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Aminoglycoside", colnames(gene_enz_tab_filt))[1]],
                                            colnames(gene_enz_tab_filt)[grep("Aminoglycoside", colnames(gene_enz_tab_filt))[2]],
                                            colnames(gene_enz_tab_filt)[grep("Homospermidine", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Mercuric ion reductase", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Sulfur carrier", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("restriction subunit R", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("methyltransferase subunit M", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Methyl-directed", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("Phosphate acetyltransferase (EC 2.3.1.8)", colnames(gene_enz_tab_filt))],
                                            colnames(gene_enz_tab_filt)[grep("1,3,6,8", colnames(gene_enz_tab_filt))])

mixedSign_old <- mixedSign        
save(mixedSign_old, file = "mixedSign_old.RData")        


load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/metClusts_oldGood.RData")
clustOrdC1_2_old <- names(sort(metClusts_oldGood[2:length(metClusts_oldGood)]))

if(!require(gplots)) install.packages('gplots')
library(gplots)

cols <- topo.colors(nrow(gene_enz_tab_filt))
cols <- cols[-1]
cols[1] <- '#000000'

tiff("minedDiffGenes_old_withRFE.tiff", width = 7000, height = 5000, units = "px", pointsize = 100)
heatmap.2(as.matrix(mixedSign_old), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(clustOrdC1_2_old, rownames(gene_tab))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of diff genes", margins = c(6, 50), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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

