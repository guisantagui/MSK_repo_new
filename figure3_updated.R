setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3")

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

if(!require(ropls)) BiocManager::install("ropls")
library(ropls)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(ggpubr)) install.packages("ggpubr")
library(ggpubr)

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/pcomp1AreaPercCircularity.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")

ccmn_norm_mets_good_old <- ccmnNormMets

rownames(ccmn_norm_mets_good_old) <- gsub("W70322", replacement = "W70332", rownames(ccmn_norm_mets_good_old))

strains <- unique(gsub("\\_.*|(PA14).*", 
                       rownames(ccmn_norm_mets_good_old), 
                       rep = "\\1"))

swarmMeans <- pcomp1Means*-1

ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[-25,]

names(swarmMeans) <- gsub("W70322", replacement = "W70332", names(swarmMeans))

swarmMeansFilt <- swarmMeans
swarmMeansFilt <- swarmMeansFilt[match(gsub("\\_.*|(PA14).*", 
                                            rownames(ccmn_norm_mets_good_old), 
                                            rep = "\\1"), names(swarmMeans))]


ccmn_norm_mets_good_old$swarmData <- swarmMeansFilt


# Add binarized swarming data
swarmBin <- rep(NA, nrow(ccmn_norm_mets_good_old))
swarmBin[which(swarmMeansFilt > 0.19)] <- 1
swarmBin[is.na(swarmBin)] <- 0
swarmBin <- as.factor(swarmBin)
ccmn_norm_mets_good_old$swarmDatBin <- swarmBin

# Add rhamnolipid production data

glycRhamn3Cats <- c(2, 2, 2, 2, 2, 0, 2, 0, 2, 1, 0, 1, 1, 2, 0, 2, 2, 2, 1, 1, 2, 0, 2, 0, 1, 2, 2, 2)
glycRhamn2Cats <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1)

#glycRhamn2Cats <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1)# let's try what happens if mid producers are considered non producers

names(glycRhamn3Cats) <- c(strains[1], strains)
names(glycRhamn2Cats) <- c(strains[1], strains)

rhamnSwarmMat <- data.frame(strain = names(glycRhamn2Cats),
                            rhamn2cats = glycRhamn2Cats,
                            rhamn3cats = glycRhamn3Cats,
                            swarmScore = swarmMeans[match(names(glycRhamn2Cats),
                                                          names(swarmMeans))])[-1, ]

write.csv(rhamnSwarmMat, file = "rhamnSwarmMat.csv")
save(rhamnSwarmMat, file = "rhamnSwarmMat.RData")


# 2nd version of file for Jinyuan, with strains not in metabolimic data
glycRhamn2Cats_2 <- glycRhamn2Cats
glycRhamn3Cats_2 <- glycRhamn3Cats

names(glycRhamn3Cats_2)[2] <- "F23197"
names(glycRhamn2Cats_2)[2] <- "F23197"

LBRhamn2Cats_2 <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1)
LBRhamn3Cats_2 <- c(2, 2, 2, 1, 2, 0, 2, 0, 2, 1, 0, 0, 0, 2, 0, 2, 1, 1, 1, 1, 0, 0, 2, 0, 0, 1, 2, 1)

glycRhamn3Cats_2[13] <- 0
glycRhamn2Cats_2[13] <- 0
names(glycRhamn3Cats_2)[13] <- "F5677"
names(glycRhamn2Cats_2)[13] <- "F5677"

names(LBRhamn2Cats_2) <- names(glycRhamn2Cats_2)
names(LBRhamn3Cats_2) <- names(glycRhamn3Cats_2)

glycRhamn3Cats_2 <- c(glycRhamn3Cats_2, 2, 2, 0)
glycRhamn2Cats_2 <- c(glycRhamn2Cats_2, 1, 1, 0)
names(glycRhamn3Cats_2)[(length(glycRhamn3Cats_2)-2):length(glycRhamn3Cats_2)] <- c("PAO1", "PA7", "M55212")
names(glycRhamn2Cats_2)[(length(glycRhamn2Cats_2)-2):length(glycRhamn2Cats_2)] <- c("PAO1", "PA7", "M55212")

LBRhamn3Cats_2 <- c(LBRhamn3Cats_2, 2, 2, 0)
LBRhamn2Cats_2 <- c(LBRhamn2Cats_2, 1, 1, 0)
names(LBRhamn3Cats_2)[(length(LBRhamn3Cats_2)-2):length(LBRhamn3Cats_2)] <- c("PAO1", "PA7", "M55212")
names(LBRhamn2Cats_2)[(length(LBRhamn2Cats_2)-2):length(LBRhamn2Cats_2)] <- c("PAO1", "PA7", "M55212")

glycRhamn2Cats_2 <- glycRhamn2Cats_2[order(names(glycRhamn2Cats_2))]
glycRhamn3Cats_2 <- glycRhamn3Cats_2[order(names(glycRhamn3Cats_2))]

LBRhamn2Cats_2 <- LBRhamn2Cats_2[order(names(LBRhamn2Cats_2))]
LBRhamn3Cats_2 <- LBRhamn3Cats_2[order(names(LBRhamn3Cats_2))]

rhamnSwarmMat_2 <- data.frame(strain = names(glycRhamn2Cats_2),
                              glyc_rhamn2cats = glycRhamn2Cats_2,
                              glyc_rhamn3cats = glycRhamn3Cats_2,
                              LBRhamn2Cats = LBRhamn2Cats_2,
                              LBRhamn3Cats = LBRhamn3Cats_2,
                              swarmScore = swarmMeans[match(names(glycRhamn2Cats_2),
                                                          names(swarmMeans))])

save(rhamnSwarmMat_2, file = "rhamnSwarmMat_2.RData")
write.csv(rhamnSwarmMat_2, file = "rhamnSwarmMat_2.csv")

rhamnSwarmMat_3 <- rbind.data.frame(rhamnSwarmMat_2,
                                    data.frame(strain = "M6075",
                                               glyc_rhamn2cats = 1,
                                               glyc_rhamn3cats = 1,
                                               LBRhamn2Cats = 1,
                                               LBRhamn3Cats = 1,
                                               swarmScore = -0.001025926))

rownames(rhamnSwarmMat_3)[nrow(rhamnSwarmMat_3)] <- "M6075"

rhamnSwarmMat_3 <- rhamnSwarmMat_3[order(as.character(rhamnSwarmMat_3$strain)), ]

write.csv(rhamnSwarmMat_3, file = "rhamnSwarmMat_3.csv")

glycRhamn3Cats <- glycRhamn3Cats[match(gsub("\\_.*|(PA14).*", 
                                            rownames(ccmn_norm_mets_good_old), 
                                            rep = "\\1"), 
                                       names(glycRhamn3Cats))]

glycRhamn2Cats <- glycRhamn2Cats[match(gsub("\\_.*|(PA14).*", 
                                            rownames(ccmn_norm_mets_good_old), 
                                            rep = "\\1"), 
                                       names(glycRhamn2Cats))]

ccmn_norm_mets_good_old$rhamn2Cats <- as.factor(glycRhamn2Cats)
ccmn_norm_mets_good_old$rhamn3Cats <- glycRhamn3Cats

ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[-grep("_?", colnames(ccmn_norm_mets_good_old), fixed = T)]

# Do PCA with swarming color code (figure S3)

PCAMetSwarms <- prcomp(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)])

pdf("PCAMetSwarms.pdf", height = 12, width = 12)
fviz_pca_ind(PCAMetSwarms, 
             col.ind = ccmn_norm_mets_good_old$swarmData,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Swarming Score")
dev.off()

##############################################################################################################################################
# OPLS-DA Swarming
##############################################################################################################################################

# OPLS-DA

# Qualitative and with partition

OPLSDA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 4)],
               ccmn_norm_mets_good_old$swarmDatBin, 
               predI = 1, 
               orthoI = NA,
               subset = "odd"
)

trainVi <- getSubsetVi(OPLSDA)
table(ccmn_norm_mets_good_old$swarmDatBin[trainVi], fitted(OPLSDA))

confusion <- table(ccmn_norm_mets_good_old$swarmDatBin[-trainVi],
                   predict(OPLSDA, ccmn_norm_mets_good_old[-trainVi, 1:(ncol(ccmn_norm_mets_good_old)-4)]))

(confusion[1, 1] + confusion[2, 2])/sum(confusion)*100 # --> 95% of accuracy

# Quantitative without partition
pdf(file = "OPLSDA_swarm.pdf", width = 10, height = 10)
OPLSDA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 4)],
               ccmn_norm_mets_good_old$swarmData, 
               predI = 1, 
               orthoI = 3,
               permI = 2000
)
dev.off()

# Obtain barplot

OPLSDALoads <- getLoadingMN(OPLSDA)


plot(getScoreMN(OPLSDA), ccmn_norm_mets_good_old$swarmData)

scoreYDF <- data.frame("scores" = getScoreMN(OPLSDA),
                       "swarmScore" = ccmn_norm_mets_good_old$swarmData)

ggscatter(scoreYDF, x = "p1", y = "swarmScore", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10, label.y = 3)) +
        stat_regline_equation(label.x = 3)

scorePredYDF <- data.frame("scores" = getScoreMN(OPLSDA),
                           "swarmScore" = predict(OPLSDA, ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)]))

ggscatter(scorePredYDF, x = "p1", y = "swarmScore", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10, label.y = 3)) +
        stat_regline_equation(label.x = 3)

plot(getScoreMN(OPLSDA), predict(OPLSDA, ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)]))

plot(scale(as.matrix(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)]) %*% getLoadingMN(OPLSDA)[, 1],
           center = T,
           scale = F)[, 1], 
     ccmn_norm_mets_good_old$swarmData)

# Remove ambiguous mets
if(length(grep("_?", rownames(OPLSDALoads), fixed = T)) == 0){
        namesOPLSDALoads <- rownames(OPLSDALoads)
        OPLSDALoads <- as.vector(OPLSDALoads)
        names(OPLSDALoads) <- namesOPLSDALoads
        
}else{
        OPLSDALoads <- OPLSDALoads[-grep("_?", rownames(OPLSDALoads), fixed = T), ]
}


defNames <- dictionary$definitiveNames[match(names(OPLSDALoads), dictionary$Consensus)]

names(OPLSDALoads)[!is.na(defNames)] <- defNames[!is.na(defNames)]

OPLSDALoads <- cbind.data.frame(names(OPLSDALoads), OPLSDALoads)
colnames(OPLSDALoads) <- c("metabolite", "loading")

ggplot(data = OPLSDALoads,
       aes(x = reorder(metabolite, loading), y = loading)) +
        labs(x = "Metabolites", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity") +
        coord_flip()
ggsave(filename = "OPLSDALoadsBarplot.pdf", limitsize = F, height = 15, width = 7)

# Let's do a fancier version of the OPLS-DA plot

OPLSDAScores <- cbind(getScoreMN(OPLSDA), getScoreMN(OPLSDA, orthoL = T)[, 1])

colnames(OPLSDAScores)[2] <- "o1"

OPLSDAScores <- as.data.frame(OPLSDAScores)

OPLSDAScores$swarmDat <- swarmMeansFilt


#low = "#003CFF", 
#mid = "#66FF00",
#high = "#FF0000",
#space = "Lab")

# get means of scores for plotting purposes

OPLSDAScoresMeans <- c()
for(i in seq_along(strains)){
        strMean <- apply(OPLSDAScores[gsub("\\_.*|(PA14).*", 
                                           rownames(OPLSDAScores), 
                                           rep = "\\1") %in% strains[i], ], 2, mean)
        OPLSDAScoresMeans <- rbind(OPLSDAScoresMeans, strMean)
}
rownames(OPLSDAScoresMeans) <- strains
OPLSDAScoresMeans <- as.data.frame(OPLSDAScoresMeans)


OPLSDAScorePlot<-ggplot(OPLSDAScores, aes(x=p1, y=o1, color=swarmDat)) + geom_point()


mid<-(max(OPLSDAScores$swarmDat) - min(OPLSDAScores$swarmDat))/2 + min(OPLSDAScores$swarmDat)
OPLSDAScorePlot+scale_color_gradient2(midpoint=mid, 
                          low="#003CFF",
                          mid="#66FF00",
                          high="#FF0000",
                          space ="Lab",
                          limits = c(min(OPLSDAScores$swarmDat),
                                     max(OPLSDAScores$swarmDat))) +
        geom_text(label=rownames(OPLSDAScores), nudge_y = 0.3) + 
        labs(col="Swarming Score") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")

ggsave("OPLSDAFancy_swarm.pdf", height = 10, width = 10)

OPLSDAScoreMeansPlot<-ggplot(OPLSDAScoresMeans, aes(x=p1, y=o1, color=swarmDat)) + 
        geom_point() +
        scale_color_gradient2(midpoint=mid, 
                              low="#003CFF",
                              mid="#66FF00",
                              high="#FF0000",
                              space ="Lab",
                              limits = c(min(OPLSDAScores$swarmDat),
                                         max(OPLSDAScores$swarmDat))) +
        geom_text(label=rownames(OPLSDAScoresMeans), nudge_y = 0.3) + 
        labs(col="Swarming Score") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("t1 (7%)") +
        ylab("to1")

OPLSDAScoreMeansPlot

ggsave("OPLSDAFancy_swarm_means.pdf", height = 10, width = 10)



# Spearman correlation test to find metabolites significatively related to swarming (univariate)
spearCor <- p.adjust(apply(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)],
                           2,
                           function(x) cor.test(x=x, 
                                                y = ccmn_norm_mets_good_old$swarmData,
                                                method = "spearman")$p.value),
                     method = "BH")

spearCorSign <- spearCor[spearCor < 0.05]

# get VIP values

VIP <- getVipVn(OPLSDA)

VIPsign <- VIP[VIP > 1]

statMetrics <- data.frame(metabolites = names(spearCor), spearCor = spearCor, VIP = VIP)

vipSign2Add <- statMetrics[statMetrics$spearCor < 0.07 & statMetrics$VIP > 1, ]

# Add to the significative ones of spear cor test the mets with VIP > 1 and p-val below 0.07

swarm_signMets <- rbind.data.frame(statMetrics[statMetrics$metabolites %in% names(spearCorSign), ], 
                                   vipSign2Add[!vipSign2Add$metabolites %in% names(spearCorSign), ])

save(swarm_signMets, file = "swarm_signMets.RData")
save(statMetrics, file = "statMetricsSwarm.RData")

quantVn <- qnorm(1 - spearCor / 2)
rmsQuantN <- sqrt(mean(quantVn^2))

plot(spearCor, VIP,
     col = "red",
     pch = 16,
     xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")

curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)

abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")


##############################################################################################################################################
# OPLS-DA Rhamn (2 categories)
##############################################################################################################################################

# OPLS-DA using rhamnolipid production as response variable (2 Categories).
pdf(file = "OPLSDA_rhamn2Cats.pdf", width = 10, height = 10)
rhamnOPLSDA2Cats <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 4)],
                         ccmn_norm_mets_good_old$rhamn2Cats, 
                         predI = 1, 
                         orthoI = 3,
                         permI = 2000)
dev.off()

predictions <- list()
for(i in 1:200){
        train <- sample(1:nrow(ccmn_norm_mets_good_old), floor(nrow(ccmn_norm_mets_good_old)/2))
        test <- 1:nrow(ccmn_norm_mets_good_old)
        test <- test[-train]
        mod <- opls(ccmn_norm_mets_good_old[train, 1:(ncol(ccmn_norm_mets_good_old) - 4)],
                    ccmn_norm_mets_good_old$rhamn2Cats[train], 
                    predI = 1, 
                    orthoI = 3,
                    permI = 1)
        pred <- predict(mod, ccmn_norm_mets_good_old[test, 1:(ncol(ccmn_norm_mets_good_old) - 4)])
        namesPred <- names(pred)
        pred <- as.numeric(pred) - 1
        names(pred) <- namesPred
        predictions[[i]] <- pred
}

library(ROCR)

labels <- lapply(predictions, function(x) ccmn_norm_mets_good_old$rhamn2Cats[match(names(x), rownames(ccmn_norm_mets_good_old))])

predictObjkt <- prediction(predictions, labels)
perf <- performance(predictObjkt, "tpr", "fpr")
plot(perf, 
     avg = "threshold",
     colorize = T#,
     #spread.estimate = "boxplot"
     )

plot(perf,
     lty=3,
     col="grey78",
     add=TRUE)

auc <- performance(predictObjkt, measure = "auc")

auc@y.values

# Histogram of permuted models

hist(rhamnOPLSDA2Cats@suppLs$permMN[, 3])

rhamnOPLSDALoads2Cats <- getLoadingMN(rhamnOPLSDA2Cats)

# Predictions

predOPLSDA2Cats <- predict(rhamnOPLSDA2Cats, ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)])
realValsOPLSDA2Cats <- ccmn_norm_mets_good_old$rhamn2Cats
names(realValsOPLSDA2Cats) <- rownames(ccmn_norm_mets_good_old)

confusion <- table(ccmn_norm_mets_good_old$rhamn2Cats,
                   predOPLSDA2Cats)




# Remove ambiguous mets
if(length(grep("_?", rownames(rhamnOPLSDALoads2Cats), fixed = T)) == 0){
        namesRhamnOPLSDALoads2Cats <- rownames(rhamnOPLSDALoads2Cats)
        rhamnOPLSDALoads2Cats <- as.vector(rhamnOPLSDALoads2Cats)
        names(rhamnOPLSDALoads2Cats) <- namesRhamnOPLSDALoads2Cats
        
}else{
        rhamnOPLSDALoads2Cats <- rhamnOPLSDALoads2Cats[-grep("_?", rownames(rhamnOPLSDALoads2Cats), fixed = T), ]
}

defNames <- dictionary$definitiveNames[match(names(rhamnOPLSDALoads2Cats), dictionary$Consensus)]

names(rhamnOPLSDALoads2Cats)[!is.na(defNames)] <- defNames[!is.na(defNames)]

rhamnOPLSDALoads2Cats <- cbind.data.frame(names(rhamnOPLSDALoads2Cats), rhamnOPLSDALoads2Cats)
colnames(rhamnOPLSDALoads2Cats) <- c("metabolite", "loading")

ggplot(data = rhamnOPLSDALoads2Cats,
       aes(x = reorder(metabolite, loading), y = loading)) +
        labs(x = "Metabolites", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity") +
        coord_flip()
ggsave(filename = "rhamnOPLSDALoads2CatsBarplot.pdf", limitsize = F, height = 15, width = 7)




rhamnOPLSDAScores2Cats <- cbind(getScoreMN(rhamnOPLSDA2Cats), getScoreMN(rhamnOPLSDA2Cats, orthoL = T)[, 1])

colnames(rhamnOPLSDAScores2Cats)[2] <- "o1"

rhamnOPLSDAScores2Cats <- as.data.frame(rhamnOPLSDAScores2Cats)

rhamnOPLSDAScores2Cats$rhamn <- as.numeric(ccmn_norm_mets_good_old$rhamn2Cats)-1


#low = "#003CFF", 
#mid = "#66FF00",
#high = "#FF0000",
#space = "Lab")

# get means of scores for plotting purposes

rhamnOPLSDAScores2CatsMeans <- c()
for(i in seq_along(strains)){
        strMean <- apply(rhamnOPLSDAScores2Cats[gsub("\\_.*|(PA14).*", 
                                                     rownames(rhamnOPLSDAScores2Cats), 
                                                     rep = "\\1") %in% strains[i], ], 2, mean)
        rhamnOPLSDAScores2CatsMeans <- rbind(rhamnOPLSDAScores2CatsMeans, strMean)
}
rownames(rhamnOPLSDAScores2CatsMeans) <- strains
rhamnOPLSDAScores2CatsMeans <- as.data.frame(rhamnOPLSDAScores2CatsMeans)


rhamnOPLSDAScorePlot2Cats <- ggplot(rhamnOPLSDAScores2Cats, aes(x=p1, y=o1, color=rhamn)) + geom_point()


rhamnOPLSDAScorePlot2Cats+scale_color_gradient2(midpoint=0.5, 
                                                low="blue",
                                                mid="white",
                                                high="red",
                                                space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScores2Cats), nudge_y = 0.3) + 
        labs(col="Rhamnolipid Production") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")

ggsave("OPLSDAFancy_rhamn2Cats.pdf", height = 10, width = 10)

rhamnOPLSDAScore2CatsMeansPlot <- ggplot(rhamnOPLSDAScores2CatsMeans, aes(x=p1, y=o1, color=rhamn)) + 
        geom_point() +
        scale_color_gradient2(midpoint=0.5, 
                              low="blue",
                              mid="white",
                              high="red",
                              space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScores2CatsMeans), nudge_y = 0.3) + 
        #labs(col="Rhamnolipid Production") + 
        theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("t1 (5%)") +
        ylab("to1")

rhamnOPLSDAScore2CatsMeansPlot

ggsave("OPLSDAFancy_rhamn2Cats_means.pdf", height = 6, width = 6)


mannWhitRhamn2Cats <- p.adjust(apply(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)],
                                     2,
                                     function(x) wilcox.test(x=x[ccmn_norm_mets_good_old$rhamn2Cats == 0], 
                                                             y=x[ccmn_norm_mets_good_old$rhamn2Cats == 1])$p.value),
                               method = "BH")
mannWhitRhamn2CatsSign <- mannWhitRhamn2Cats[mannWhitRhamn2Cats < 0.05]

# get VIP values

VIPRhamn2Cats <- getVipVn(rhamnOPLSDA2Cats)

VIPRhamnsign2Cats <- VIPRhamn2Cats[VIPRhamn2Cats > 1]

statMetricsRhamn2Cats <- data.frame(metabolites = names(mannWhitRhamn2Cats), 
                                    mannWhit = mannWhitRhamn2Cats, 
                                    VIP = VIPRhamn2Cats)

vipSign2Add_rhamn2Cats <- statMetricsRhamn2Cats[statMetricsRhamn2Cats$mannWhit < 0.07 & statMetricsRhamn2Cats$VIP > 1, ]
rhamn2Cats_signMets <- rbind.data.frame(statMetricsRhamn2Cats[statMetricsRhamn2Cats$metabolites %in% names(mannWhitRhamn2CatsSign), ], 
                                        vipSign2Add_rhamn2Cats[!vipSign2Add_rhamn2Cats$metabolites %in% names(mannWhitRhamn2CatsSign), ])

save(rhamn2Cats_signMets, file = "rhamn2Cats_signMets.RData")
save(statMetricsRhamn2Cats, file = "statMetricsRhamn2Cats.RData")

statMetricsRhamn2Cats[statMetricsRhamn2Cats$mannWhit < 0.05, ]

quantVn <- qnorm(1 - mannWhitRhamn2Cats / 2)
rmsQuantN <- sqrt(mean(quantVn^2))

plot(mannWhitRhamn2Cats, VIPRhamn2Cats,
     col = "red",
     pch = 16,
     xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")

curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)

abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")


##############################################################################################################################################
# OPLS-DA Rhamn (3 categories)
##############################################################################################################################################


# OPLS-DA using rhamnolipid production as response variable (3 Categories).
pdf(file = "OPLSDA_rhamn3Cats.pdf", width = 10, height = 10)
rhamnOPLSDA3Cats <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 4)],
                         ccmn_norm_mets_good_old$rhamn3Cats, 
                         predI = 1, 
                         orthoI = 3,
                         permI = 2000)
dev.off()

rhamnOPLSDALoads3Cats <- getLoadingMN(rhamnOPLSDA3Cats)

# Remove ambiguous mets
if(length(grep("_?", rownames(rhamnOPLSDALoads3Cats), fixed = T)) == 0){
        namesRhamnOPLSDALoads3Cats <- rownames(rhamnOPLSDALoads3Cats)
        rhamnOPLSDALoads3Cats <- as.vector(rhamnOPLSDALoads3Cats)
        names(rhamnOPLSDALoads3Cats) <- namesRhamnOPLSDALoads3Cats
        
}else{
        rhamnOPLSDALoads3Cats <- rhamnOPLSDALoads3Cats[-grep("_?", rownames(rhamnOPLSDALoads3Cats), fixed = T), ]
}

defNames <- dictionary$definitiveNames[match(names(rhamnOPLSDALoads3Cats), dictionary$Consensus)]

names(rhamnOPLSDALoads3Cats)[!is.na(defNames)] <- defNames[!is.na(defNames)]

rhamnOPLSDALoads3Cats <- cbind.data.frame(names(rhamnOPLSDALoads3Cats), rhamnOPLSDALoads3Cats)
colnames(rhamnOPLSDALoads3Cats) <- c("metabolite", "loading")

ggplot(data = rhamnOPLSDALoads3Cats,
       aes(x = reorder(metabolite, loading), y = loading)) +
        labs(x = "Metabolites", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity") +
        coord_flip()
ggsave(filename = "rhamnOPLSDALoads3CatsBarplot.pdf", limitsize = F, height = 15, width = 7)




rhamnOPLSDAScores3Cats <- cbind(getScoreMN(rhamnOPLSDA3Cats), getScoreMN(rhamnOPLSDA3Cats, orthoL = T)[, 1])

colnames(rhamnOPLSDAScores3Cats)[2] <- "o1"

rhamnOPLSDAScores3Cats <- as.data.frame(rhamnOPLSDAScores3Cats)

rhamnOPLSDAScores3Cats$rhamn <- ccmn_norm_mets_good_old$rhamn3Cats


#low = "#003CFF", 
#mid = "#66FF00",
#high = "#FF0000",
#space = "Lab")

# get means of scores for plotting purposes

rhamnOPLSDAScores3CatsMeans <- c()
for(i in seq_along(strains)){
        strMean <- apply(rhamnOPLSDAScores3Cats[gsub("\\_.*|(PA14).*", 
                                                     rownames(rhamnOPLSDAScores3Cats), 
                                                     rep = "\\1") %in% strains[i], ], 2, mean)
        rhamnOPLSDAScores3CatsMeans <- rbind(rhamnOPLSDAScores3CatsMeans, strMean)
}
rownames(rhamnOPLSDAScores3CatsMeans) <- strains
rhamnOPLSDAScores3CatsMeans <- as.data.frame(rhamnOPLSDAScores3CatsMeans)


rhamnOPLSDAScorePlot3Cats <- ggplot(rhamnOPLSDAScores3Cats, aes(x=p1, y=o1, color=rhamn)) + geom_point()


rhamnOPLSDAScorePlot3Cats+scale_color_gradient2(midpoint=1, 
                                                low="blue",
                                                mid="green",
                                                high="red",
                                                space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScores3Cats), nudge_y = 0.3) + 
        labs(col="Rhamnolipid Production") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")

ggsave("OPLSDAFancy_rhamn3Cats.pdf", height = 10, width = 10)

rhamnOPLSDAScore3CatsMeansPlot <- ggplot(rhamnOPLSDAScores3CatsMeans, aes(x=p1, y=o1, color=rhamn)) + 
        geom_point() +
        scale_color_gradient2(midpoint=1, 
                              low="blue",
                              mid="green",
                              high="red",
                              space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScores3CatsMeans), nudge_y = 0.3) + 
        labs(col="Rhamnolipid Production") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("t1 (10%)") +
        ylab("to1")

rhamnOPLSDAScore3CatsMeansPlot

ggsave("OPLSDAFancy_rhamn3Cats_means.pdf", height = 10, width = 10)


kruskalRhamn3Cats <- p.adjust(apply(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)],
                                    2,
                                    function(x) kruskal.test(x=list(x[ccmn_norm_mets_good_old$rhamn3Cats == 0],
                                                                    x[ccmn_norm_mets_good_old$rhamn3Cats == 1],
                                                                    x[ccmn_norm_mets_good_old$rhamn3Cats == 2]))$p.value),
                              method = "BH")
kruskalRhamn3CatsSign <- kruskalRhamn3Cats[kruskalRhamn3Cats < 0.05]



# Try also spearman correlation to see the number of significant compounds

spearCorRhamn3Cats <- p.adjust(apply(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-4)],
                                     2,
                                     function(x) cor.test(x=x, 
                                                          y = ccmn_norm_mets_good_old$rhamn3Cats,
                                                          method = "spearman")$p.value),
                               method = "BH")

spearCorRhamn3CatsSign <- spearCorRhamn3Cats[spearCorRhamn3Cats < 0.05]






# get VIP values

VIPRhamn3Cats <- getVipVn(rhamnOPLSDA3Cats)

VIPRhamnsign3Cats <- VIPRhamn3Cats[VIPRhamn3Cats > 1]

statMetricsRhamn3Cats <- data.frame(metabolites = names(kruskalRhamn3Cats), 
                                    kruskal = kruskalRhamn3Cats, 
                                    spearCor = spearCorRhamn3Cats,
                                    VIP = VIPRhamn3Cats)


vipSign2Add_rhamn3Cats <- statMetricsRhamn3Cats[statMetricsRhamn3Cats$kruskal < 0.07 & statMetricsRhamn3Cats$VIP > 1, ]
rhamn3Cats_signMets <- rbind.data.frame(statMetricsRhamn3Cats[statMetricsRhamn3Cats$metabolites %in% names(kruskalRhamn3CatsSign), ], 
                                        vipSign2Add_rhamn3Cats[!vipSign2Add_rhamn3Cats$metabolites %in% names(kruskalRhamn3CatsSign), ])


save(rhamn3Cats_signMets, file = "rhamn3Cats_signMets.RData")
save(statMetricsRhamn3Cats, file = "statMetricsRhamn3Cats.RData")

statMetricsRhamn3Cats[statMetricsRhamn3Cats$kruskal < 0.05, ]

quantVn <- qnorm(1 - kruskalRhamn3Cats / 2)
rmsQuantN <- sqrt(mean(quantVn^2))

plot(kruskalRhamn3Cats, VIPRhamn3Cats,
     col = "red",
     pch = 16,
     xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")

curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)

abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")

# Spear curve 

quantVn <- qnorm(1 - spearCorRhamn3Cats / 2)
rmsQuantN <- sqrt(mean(quantVn^2))

plot(spearCorRhamn3Cats, VIPRhamn3Cats,
     col = "red",
     pch = 16,
     xlab = "p-value", ylab = "VIP", xaxs = "i", yaxs = "i")

curve(qnorm(1 - x / 2) / rmsQuantN, 0, 1, add = TRUE, col = "red", lwd = 3)

abline(h = 1, col = "blue")
abline(v = 0.05, col = "blue")


##############################################################################################################################################
# OPLS-DA models with significant metabolites 
##############################################################################################################################################


# Let's do a model using only the metabolites that had a VIP over 1 in previous model()
if(!require(gplots)) install.packages("gplots")
library(gplots)

# Do color scales

gradient <- colorRampPalette(colors = c("#003CFF", "#66FF00", "#FF0000"))
swarmCols <- gradient(nrow(ccmn_norm_mets_good_old))[cut(ccmn_norm_mets_good_old$swarmData + abs(min(ccmn_norm_mets_good_old$swarmData)), 83)] 

rhamn2CatsCols <- rep("red", nrow(ccmn_norm_mets_good_old))
rhamn2CatsCols[ccmn_norm_mets_good_old$rhamn2Cats == 0] <- "blue"

rhamn3CatsCols <- rep("red", nrow(ccmn_norm_mets_good_old))
rhamn3CatsCols[ccmn_norm_mets_good_old$rhamn3Cats == 0] <- "blue"
rhamn3CatsCols[ccmn_norm_mets_good_old$rhamn3Cats == 1] <- "green"

# With metabolites significative in swarming model

if(length(grep("_?", 
               swarm_signMets$metabolites, 
               fixed = T)) == 0){
        ccmnSwarmSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% swarm_signMets$metabolites]
}else{
        ccmnSwarmSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% swarm_signMets$metabolites[-grep("_?", 
                                                                                                                               swarm_signMets$metabolites, 
                                                                                                                               fixed = T)]]
}


tiff("OPLSDA_swarm_signMets.tiff")
OPLSDA_swarm_signMets <- opls(ccmnSwarmSignMets,
                              ccmn_norm_mets_good_old$swarmData,
                              #predI = 1, 
                              #orthoI = 1,
                              #permI = 500
                              )
dev.off()

plot(predict(OPLSDA_swarm_signMets, ccmnSwarmSignMets), ccmn_norm_mets_good_old$swarmData)

scoreYDFSign <- data.frame("scores" = getScoreMN(OPLSDA_swarm_signMets),
                           "swarmScore" = ccmn_norm_mets_good_old$swarmData)



ggscatter(scoreYDFSign, x = "p1", y = "swarmScore", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10, label.y = 3)) +
        stat_regline_equation(label.x = 3)


plot(hclust(dist(ccmnSwarmSignMets, method = "euclidean"), method = "ward.D"))

tiff("heatSwarmSignMets.tiff", height = 500, width = 500)
heatSwarmSignMets <- heatmap.2(as.matrix(t(ccmnSwarmSignMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                               density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                               col = redgreen(75), breaks = 76, ColSideColors = swarmCols, notecol = NULL, trace = "none", xlab = "Strains", 
                               ylab = "Metabolites", margins = c(8, 8), 
                               cex.main = 20,
                               keysize = 0.7,
                               cexRow = 0.7,
                               cexCol = 1.2,
                               scale = "row",
                               #colCol = colCols,
                               colCol = swarmCols,
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
# With metabolites significative in rhamn (2 categories) model

if(length(grep("_?", 
               rhamn2Cats_signMets$metabolites, 
               fixed = T)) == 0){
        ccmnRhamn2CatsSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% rhamn2Cats_signMets$metabolites]
}else{
        ccmnRhamn2CatsSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% rhamn2Cats_signMets$metabolites[-grep("_?", 
                                                                                                                                         rhamn2Cats_signMets$metabolites, 
                                                                                                                                         fixed = T)]]
}

tiff("OPLSDA_rhamn2Cats_signMets.tiff")
OPLSDA_rhamn2Cats_signMets <- opls(ccmnRhamn2CatsSignMets,
                                   ccmn_norm_mets_good_old$rhamn2Cats,
                                   predI = 1, 
                                   orthoI = NA,
                                   permI = 500)
dev.off()

scoreYDFRhamn2CatsSign <- data.frame("scores" = getScoreMN(OPLSDA_rhamn2Cats_signMets),
                                     "rhamn" = ccmn_norm_mets_good_old$rhamn2Cats)

ggscatter(scoreYDFRhamn2CatsSign, x = "p1", y = "rhamn", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10, label.y = 3)) +
        stat_regline_equation(label.x = 3)

plot(hclust(dist(ccmnRhamn2CatsSignMets, method = "euclidean"), method = "ward.D"))

tiff("heatRhamn2CatsSignMets.tiff", height = 500, width = 500)
heatRhamn2CatsSignMets <- heatmap.2(as.matrix(t(ccmnRhamn2CatsSignMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                                    density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                                    col = redgreen(75), breaks = 76, ColSideColors = rhamn2CatsCols, notecol = NULL, trace = "none", xlab = "Strains", 
                                    ylab = "Metabolites", margins = c(8, 8), 
                                    cex.main = 20,
                                    keysize = 0.7,
                                    cexRow = 0.7,
                                    cexCol = 1.2,
                                    scale = "row",
                                    #colCol = colCols,
                                    colCol = rhamn2CatsCols,
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

# With metabolites significative in rhamn (3 categories) model

if(length(grep("_?", 
               rhamn3Cats_signMets$metabolites, 
               fixed = T)) == 0){
        ccmnRhamn3CatsSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% rhamn3Cats_signMets$metabolites]
}else{
        ccmnRhamn3CatsSignMets <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% rhamn3Cats_signMets$metabolites[-grep("_?", 
                                                                                                                                         rhamn3Cats_signMets$metabolites, 
                                                                                                                                         fixed = T)]]
}


tiff("OPLSDA_rhamn3Cats_signMets.tiff")
OPLSDA_rhamn3Cats_signMets <- opls(ccmnRhamn3CatsSignMets,
                                   ccmn_norm_mets_good_old$rhamn3Cats,
                                   predI = 1, 
                                   orthoI = NA,
                                   permI = 500)
dev.off()

scoreYDFRhamn3CatsSign <- data.frame("scores" = getScoreMN(OPLSDA_rhamn3Cats_signMets),
                                     "rhamn" = ccmn_norm_mets_good_old$rhamn3Cats)

ggscatter(scoreYDFRhamn3CatsSign, x = "p1", y = "rhamn", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10, label.y = 3)) +
        stat_regline_equation(label.x = 3)

plot(hclust(dist(ccmnRhamn3CatsSignMets, method = "euclidean"), method = "ward.D"))


tiff("heatRhamn3CatsSignMets.tiff", height = 500, width = 500)
heatRhamn3CatsSignMets <- heatmap.2(as.matrix(t(ccmnRhamn3CatsSignMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                                    density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                                    col = redgreen(75), breaks = 76, ColSideColors = rhamn3CatsCols, notecol = NULL, trace = "none", xlab = "Strains", 
                                    ylab = "Metabolites", margins = c(8, 8), 
                                    cex.main = 20,
                                    keysize = 0.7,
                                    cexRow = 0.7,
                                    cexCol = 1.2,
                                    scale = "row",
                                    #colCol = colCols,
                                    colCol = rhamn3CatsCols,
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


##############################################################################################################################################
# FELLA
##############################################################################################################################################



# Do FELLA for mets with highest loading:

if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)

# First with loadings of swarming OPLSDA

# Get top swarming predictors compounds: the significative ones according to Spear Cor test

swarmSignNames <- colnames(ccmnSwarmSignMets)
swarmSignNames[!is.na(dictionary$definitiveNames[match(swarmSignNames, 
                                                       dictionary$Consensus)])] <- dictionary$definitiveNames[match(swarmSignNames, 
                                                                                                                    dictionary$Consensus)][!is.na(dictionary$definitiveNames[match(swarmSignNames, 
                                                                                                                                                                                   dictionary$Consensus)])]

swarm_signMetsAmbigRem <- swarm_signMets[-grep("_?", swarm_signMets$metabolites, fixed = T), ]
swarm_signMetsAmbigRem$starTag <- rep("", nrow(swarm_signMetsAmbigRem))
swarm_signMetsAmbigRem$starTag[swarm_signMetsAmbigRem$spearCor < 0.07] <- "+"
swarm_signMetsAmbigRem$starTag[swarm_signMetsAmbigRem$spearCor < 0.05] <- "*"
swarm_signMetsAmbigRem$starTag[swarm_signMetsAmbigRem$spearCor < 0.01] <- "**"
swarm_signMetsAmbigRem$starTag[swarm_signMetsAmbigRem$spearCor < 0.001] <- "***"



KEGGIDs <- dictionary$`KEGG IDs`[match(swarmSignNames, dictionary$definitiveNames)]

graph <- buildGraphFromKEGGREST(
        organism = "pau",
        filter.path = c("01100", "01200", "01210", "01212", "01230")
)

buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = NULL,
        internalDir = TRUE,
        matrices = "none",
        normality = "diffusion",
        niter = 100)

loadKEGGdata()

entrez2ec <- KEGGREST::keggLink("enzyme", "pau")
entrez2path <- KEGGREST::keggLink("pathway", "pau")
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
        compounds = KEGGIDs,
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

pdf("FELLA_clusts_rfeRes_swarm.pdf", width = 15, height = 15)
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

KEGGIDs_all <- dictionary$`KEGG IDs`[match(rownames(OPLSDALoads), dictionary$definitiveNames)]
KEGGIDs_all <- KEGGIDs_all[order(OPLSDALoads$loading)]
KEGGIDs_all <- KEGGIDs_all[!is.na(KEGGIDs_all)]

pathsPerMet <- list()
for(i in seq_along(KEGGIDs_all)){
        print(i)
        paths <- keggLink("pathway", 
                          KEGGIDs_all[i])[gsub("path:map", 
                                               "pau", 
                                               keggLink("pathway", 
                                                        KEGGIDs_all[i])) %in% tab_rfeSwarm$KEGG.id]
        pathsPerMet[[i]] <- paths
}
names(pathsPerMet) <- KEGGIDs_all
pathsPerMet

allFELLAPaths <- unique(unlist(pathsPerMet))

metsInPaths <- list()

for(i in seq_along(allFELLAPaths)){
        allMetsinPath <- keggLink("compound", allFELLAPaths[i])
        allMetsinPath <- gsub("cpd:", "", allMetsinPath)
        mets <- allMetsinPath[allMetsinPath %in% KEGGIDs_all]
        metsInPaths[[i]] <- mets
}
names(metsInPaths) <- allFELLAPaths
metsInPaths$anyRelatedPath <- KEGGIDs_all[!KEGGIDs_all %in% unique(unlist(metsInPaths))]

loadsPaths <- data.frame()
for(i in seq_along(names(metsInPaths))){
        mets <- metsInPaths[[i]]
        metsDefNames <- dictionary$definitiveNames[match(mets, dictionary$`KEGG IDs`)]
        loadsSubMat <- OPLSDALoads[match(metsDefNames, OPLSDALoads$metabolite), ]
        loadsSubMat <- cbind.data.frame(loadsSubMat, rep(names(metsInPaths)[i], nrow(loadsSubMat)))
        loadsPaths <- rbind.data.frame(loadsPaths, loadsSubMat)
}

colnames(loadsPaths)[3] <- "Pathway"

allPathways <- keggList("pathway")

# Create a vector of colors to be assigned to each pathway
if(!require(randomcoloR)) install.packages("randomcoloR")
library(randomcoloR)

ncolors <- length(allPathways) + 1

pathCols <- distinctColorPalette(ncolors)

names(pathCols) <- c(allPathways, "Not Mapped")

loadsPaths$Pathway <- allPathways[match(loadsPaths$Pathway, names(allPathways))]
loadsPaths$Pathway[is.na(loadsPaths$Pathway)] <- "Not Mapped"

loadsPaths$metabolite <- make.unique(as.character(loadsPaths$metabolite))
loadsPaths$metabolite <- factor(loadsPaths$metabolite, levels = loadsPaths$metabolite)

# Create color code for pathways



pathCols[match(loadsPaths$Pathway, names(pathCols))]

# Save is commented to avoid overwriting this color code of pathways, which was nice

#save(pathCols, file = "pathCols.RData")
load("pathCols.RData")

ggplot(data = loadsPaths,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPaths$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size = 14)) #+ 
        #geom_text(data = swarm_signMetsAmbigRem, label = swarm_signMetsAmbigRem$starTag)

ggsave(filename = "loadsOPLSDAPaths.pdf", height = 20, width = 12)



# Then with rhamn (2 categories) significative mets

# Get top swarming predictors compounds

rhamn2CatsSignNames <- colnames(ccmnRhamn2CatsSignMets)
rhamn2CatsSignNames[!is.na(dictionary$definitiveNames[match(rhamn2CatsSignNames, 
                                                            dictionary$Consensus)])] <- dictionary$definitiveNames[match(rhamn2CatsSignNames, 
                                                                                                                         dictionary$Consensus)][!is.na(dictionary$definitiveNames[match(rhamn2CatsSignNames, 
                                                                                                                                                                                        dictionary$Consensus)])]


KEGGIDs <- dictionary$`KEGG IDs`[match(rhamn2CatsSignNames, dictionary$definitiveNames)]

KEGGIDs <- KEGGIDs[!is.na(KEGGIDs)]


analysis.rfeRes<- enrich(
        compounds = KEGGIDs,
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

pdf("FELLA_clusts_rfeRes_rhamn2Cats.pdf", width = 15, height = 15)
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

write.csv(tab_rfeSwarm, file = "tab_rfeRhamn2Cats.csv")
save(tab_rfeSwarm, file = "tab_rfeRhamn2Cats.RData")

KEGGIDs_all <- dictionary$`KEGG IDs`[match(rownames(rhamnOPLSDALoads2Cats), dictionary$definitiveNames)]
KEGGIDs_all <- KEGGIDs_all[order(rhamnOPLSDALoads2Cats$loading)]
KEGGIDs_all <- KEGGIDs_all[!is.na(KEGGIDs_all)]

pathsPerMet <- list()
for(i in seq_along(KEGGIDs_all)){
        print(i)
        paths <- keggLink("pathway", 
                          KEGGIDs_all[i])[gsub("path:map", 
                                               "pau", 
                                               keggLink("pathway", 
                                                        KEGGIDs_all[i])) %in% tab_rfeSwarm$KEGG.id]
        pathsPerMet[[i]] <- paths
}
names(pathsPerMet) <- KEGGIDs_all
pathsPerMet

allFELLAPaths <- unique(unlist(pathsPerMet))

metsInPaths <- list()

for(i in seq_along(allFELLAPaths)){
        allMetsinPath <- keggLink("compound", allFELLAPaths[i])
        allMetsinPath <- gsub("cpd:", "", allMetsinPath)
        mets <- allMetsinPath[allMetsinPath %in% KEGGIDs_all]
        metsInPaths[[i]] <- mets
}
names(metsInPaths) <- allFELLAPaths
metsInPaths$anyRelatedPath <- KEGGIDs_all[!KEGGIDs_all %in% unique(unlist(metsInPaths))]

loadsPaths <- data.frame()
for(i in seq_along(names(metsInPaths))){
        mets <- metsInPaths[[i]]
        metsDefNames <- dictionary$definitiveNames[match(mets, dictionary$`KEGG IDs`)]
        loadsSubMat <- rhamnOPLSDALoads2Cats[match(metsDefNames, rhamnOPLSDALoads2Cats$metabolite), ]
        loadsSubMat <- cbind.data.frame(loadsSubMat, rep(names(metsInPaths)[i], nrow(loadsSubMat)))
        loadsPaths <- rbind.data.frame(loadsPaths, loadsSubMat)
}

colnames(loadsPaths)[3] <- "Pathway"

allPathways <- keggList("pathway")

loadsPaths$Pathway <- allPathways[match(loadsPaths$Pathway, names(allPathways))]
loadsPaths$Pathway[is.na(loadsPaths$Pathway)] <- "Not Mapped"

loadsPaths$metabolite <- make.unique(as.character(loadsPaths$metabolite))
loadsPaths$metabolite <- factor(loadsPaths$metabolite, levels = loadsPaths$metabolite)

#Change names to the ones sugested by Kyu Rhee

KR_metNames <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/KRNewMetNames/metabolitesName-Xavier lab-KR.xlsx"))

loadsPaths$metabolite <- make.unique(KR_metNames$newName[match(gsub("\\..*", "", loadsPaths$metabolite), KR_metNames$oldName)])

pathsBarplot <- unique(loadsPaths$Pathway)


loadsPathsOrdered <- c()
for(i in 1:length(pathsBarplot)){
        metPath <- loadsPaths[loadsPaths$Pathway == pathsBarplot[i], ]
        metPath <- metPath[order(metPath$loading), ]
        if(is.null(nrow(loadsPathsOrdered))){
                loadsPathsOrdered <- metPath
        }else{
                loadsPathsOrdered <- rbind.data.frame(loadsPathsOrdered, metPath)
        }
}

loadsPathsOrdered$metabolite <- factor(loadsPathsOrdered$metabolite, levels = loadsPathsOrdered$metabolite)

ggplot(data = loadsPaths,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPaths$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size = 14))

ggsave(filename = "loadsOPLSDAPaths_rhamn2Cats.pdf", height = 20, width = 12)

ggplot(data = loadsPathsOrdered,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPaths$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size = 14))

ggsave(filename = "loadsOPLSDAPaths_rhamn2Cats_ordered_KRmetNames.pdf", height = 20, width = 12)

# Filter to show only pathways of interest
loadsPathsFilt <- loadsPaths[loadsPaths$Pathway %in% c("Valine, leucine and isoleucine biosynthesis", 
                                                       "Alanine, aspartate and glutamate metabolism", 
                                                       "Cysteine and methionine metabolism",
                                                       "Pyruvate metabolism", 
                                                       "Pentose phosphate pathway",
                                                       "Citrate cycle (TCA cycle)"), ]

ordVec <- unlist(sapply(c("Valine, leucine and isoleucine biosynthesis", 
                          "Alanine, aspartate and glutamate metabolism", 
                          "Cysteine and methionine metabolism",
                          "Pyruvate metabolism", 
                          "Pentose phosphate pathway",
                          "Citrate cycle (TCA cycle)"), function(x) which(loadsPathsFilt$Pathway %in% x)))

loadsPathsFilt <- loadsPathsFilt[ordVec, ]
loadsPathsFilt$metabolite <- factor(as.character(loadsPathsFilt$metabolite), 
                                    levels = as.character(loadsPathsFilt$metabolite))


ggplot(data = loadsPathsFilt,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPathsFilt$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size=20))

ggsave(filename = "loadsOPLSDAPathsFilt_rhamn2Cats.pdf", height = 20, width = 12)


loadsPathsFilt_custCategory <- loadsPathsFilt[!loadsPathsFilt$Pathway == "Pyruvate metabolism", ]

loadsPathsFilt_custCategory$Pathway[loadsPathsFilt_custCategory$Pathway == "Valine, leucine and isoleucine biosynthesis"] <- "Amino acid metabolism"
loadsPathsFilt_custCategory$Pathway[loadsPathsFilt_custCategory$Pathway == "Alanine, aspartate and glutamate metabolism"] <- "Amino acid metabolism"
loadsPathsFilt_custCategory$Pathway[loadsPathsFilt_custCategory$Pathway == "Cysteine and methionine metabolism"] <- "Amino acid metabolism"

loadsPathsFilt_custCategory <- rbind.data.frame(loadsPathsFilt_custCategory[order(loadsPathsFilt_custCategory$loading[loadsPathsFilt_custCategory$Pathway == "Amino acid metabolism"]), ],
                                                loadsPathsFilt_custCategory[order(loadsPathsFilt_custCategory$loading[loadsPathsFilt_custCategory$Pathway == "Pentose phosphate pathway"]) + 28, ],
                                                loadsPathsFilt_custCategory[order(loadsPathsFilt_custCategory$loading[loadsPathsFilt_custCategory$Pathway == "Citrate cycle (TCA cycle)"])+ 34, ])


loadsPathsFilt_custCategory <- loadsPathsFilt_custCategory[loadsPathsFilt_custCategory$metabolite %in% unique(gsub("\\..*", "", loadsPathsFilt_custCategory$metabolite)), ]

loadsPathsFilt_custCategory$metabolite <- factor(as.character(loadsPathsFilt_custCategory$metabolite), 
                                                 levels = as.character(loadsPathsFilt_custCategory$metabolite))

pathCols <- c(pathCols, "#FF4F33")
names(pathCols)[length(pathCols)] <- "Amino acid metabolism"


ggplot(data = loadsPathsFilt_custCategory,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPathsFilt_custCategory$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size=20))

ggsave(filename = "loadsOPLSDAPathsFilt_rhamn2Cats_custCategory.pdf", height = 10, width = 10.5)

# Then with rhamn (3 categories) significative mets

# Get top swarming predictors compounds

rhamn3CatsSignNames <- colnames(ccmnRhamn3CatsSignMets)
rhamn3CatsSignNames[!is.na(dictionary$definitiveNames[match(rhamn3CatsSignNames, 
                                                            dictionary$Consensus)])] <- dictionary$definitiveNames[match(rhamn3CatsSignNames, 
                                                                                                                         dictionary$Consensus)][!is.na(dictionary$definitiveNames[match(rhamn3CatsSignNames, 
                                                                                                                                                                                        dictionary$Consensus)])]


KEGGIDs <- dictionary$`KEGG IDs`[match(rhamn3CatsSignNames, dictionary$definitiveNames)]

KEGGIDs <- KEGGIDs[!is.na(KEGGIDs)]


analysis.rfeRes<- enrich(
        compounds = KEGGIDs,
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

pdf("FELLA_clusts_rfeRes_rhamn3Cats.pdf", width = 15, height = 15)
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

write.csv(tab_rfeSwarm, file = "tab_rfeRhamn3Cats.csv")
save(tab_rfeSwarm, file = "tab_rfeRhamn3Cats.RData")

KEGGIDs_all <- dictionary$`KEGG IDs`[match(rownames(rhamnOPLSDALoads3Cats), dictionary$definitiveNames)]
KEGGIDs_all <- KEGGIDs_all[order(rhamnOPLSDALoads3Cats$loading)]
KEGGIDs_all <- KEGGIDs_all[!is.na(KEGGIDs_all)]

pathsPerMet <- list()
for(i in seq_along(KEGGIDs_all)){
        print(i)
        paths <- keggLink("pathway", 
                          KEGGIDs_all[i])[gsub("path:map", 
                                               "pau", 
                                               keggLink("pathway", 
                                                        KEGGIDs_all[i])) %in% tab_rfeSwarm$KEGG.id]
        pathsPerMet[[i]] <- paths
}
names(pathsPerMet) <- KEGGIDs_all
pathsPerMet

allFELLAPaths <- unique(unlist(pathsPerMet))

metsInPaths <- list()

for(i in seq_along(allFELLAPaths)){
        allMetsinPath <- keggLink("compound", allFELLAPaths[i])
        allMetsinPath <- gsub("cpd:", "", allMetsinPath)
        mets <- allMetsinPath[allMetsinPath %in% KEGGIDs_all]
        metsInPaths[[i]] <- mets
}
names(metsInPaths) <- allFELLAPaths
metsInPaths$anyRelatedPath <- KEGGIDs_all[!KEGGIDs_all %in% unique(unlist(metsInPaths))]

loadsPaths <- data.frame()
for(i in seq_along(names(metsInPaths))){
        mets <- metsInPaths[[i]]
        metsDefNames <- dictionary$definitiveNames[match(mets, dictionary$`KEGG IDs`)]
        loadsSubMat <- rhamnOPLSDALoads3Cats[match(metsDefNames, rhamnOPLSDALoads3Cats$metabolite), ]
        loadsSubMat <- cbind.data.frame(loadsSubMat, rep(names(metsInPaths)[i], nrow(loadsSubMat)))
        loadsPaths <- rbind.data.frame(loadsPaths, loadsSubMat)
}

colnames(loadsPaths)[3] <- "Pathway"

allPathways <- keggList("pathway")

loadsPaths$Pathway <- allPathways[match(loadsPaths$Pathway, names(allPathways))]
loadsPaths$Pathway[is.na(loadsPaths$Pathway)] <- "Not Mapped"

loadsPaths$metabolite <- make.unique(as.character(loadsPaths$metabolite))
loadsPaths$metabolite <- factor(loadsPaths$metabolite, levels = loadsPaths$metabolite)


ggplot(data = loadsPaths,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPaths$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size = 14))

ggsave(filename = "loadsOPLSDAPaths_rhamn3Cats.pdf", height = 20, width = 12)

save(loadsPaths, file = "loadsPaths.RData")

# Filter to show only pathways of interest
loadsPathsFilt <- loadsPaths[loadsPaths$Pathway %in% c("Valine, leucine and isoleucine biosynthesis", 
                                                       "Pentose phosphate pathway", 
                                                       "Pyruvate metabolism", 
                                                       "Citrate cycle (TCA cycle)",
                                                       "Alanine, aspartate and glutamate metabolism", 
                                                       "Cysteine and methionine metabolism", 
                                                       "Arginine biosynthesis"), ]

ggplot(data = loadsPathsFilt,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() + 
        scale_fill_manual(values = pathCols[match(loadsPathsFilt$Pathway, names(pathCols))]) +
        scale_x_discrete(position = "top") +
        theme(axis.text = element_text(size=20))

ggsave(filename = "loadsOPLSDAPathsFilt_rhamn3Cats.pdf", height = 20, width = 12)

save(loadsPaths, file = "loadsPaths.RData")



if(!require(graph)) BiocManager::install("graph")
if(!require(RBGL)) BiocManager::install("RBGL")
if(!require(Vennerable)) install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library(Vennerable)

overlapSignMets <- list("OPLS Swarm Data" = swarm_signMets$metabolites, 
                        "OPLS-DA rhamnolipid production (2 categories)" = rhamn2Cats_signMets$metabolites, 
                        "OPLS-DA rhamnolipid production (3 categories)" = rhamn3Cats_signMets$metabolites)

vennSignMets <- Venn(Sets = overlapSignMets)

tiff(filename = "overlapSignMets.tiff", height = 1400, width = 1800, res = 300)
plot(vennSignMets, doWeights = T, type = "circles")
dev.off()
