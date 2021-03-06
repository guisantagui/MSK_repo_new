setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis")

swarmDat <- read.csv("C:/Users/Guillem/Documents/PhD/comput/pics/allSwarms.csv")
swarmDat$Strains <- gsub("H5707", "H5708", swarmDat$Strains)
swarmDat$Strains <- gsub("M6057", "M6075", swarmDat$Strains)

# Normalize each experiment to the mean value of PA14 (for each measure)

for(i in 1:length(exps)){
        expDat <- swarmDat[swarmDat$exp == exps[i], ]
        PA14 <- expDat[gsub("_.*", "", expDat$Strains) == "PA14", 2:8]
        PA14_mean <- apply(PA14, 2, mean)
        for(j in 2:8){
                expDat[, j] <- expDat[, j]/PA14_mean[j - 1]
        }
        swarmDat[swarmDat$exp == exps[i], ] <- expDat
}

# PCA for looking what variable among the ones we used for measuring the swarming is most separated in metabolites PCA.
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(readxl)) install.packages("readxl")
library(readxl)

# Load metabolomics data and dictionary
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]

ccmnNormMets <- ccmn_norm_mets_good_old
colnames(ccmnNormMets) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
ccmnNormMets <- ccmnNormMets[, !is.na(colnames(ccmnNormMets))]

groups <- unique(gsub("_.*", replacement = "", rownames(ccmnNormMets)))

getSwarmDataMeans <- function(swarmMat){
        groups <- unique(gsub("\\_.*", replacement = "", swarmMat$Strains))
        meanDF <- data.frame()
        for(i in seq_along(groups)){
                subMat <- swarmMat[grep(groups[i], swarmMat$Strains), ]
                means <- apply(subMat[, 2:8], 2, mean)
                meanDF <- rbind.data.frame(meanDF, means)
        }
        meanDF <- cbind.data.frame(groups, meanDF)
        colnames(meanDF) <- colnames(swarmMat)[1:8]
        return(meanDF)
}
swarmDatMeans <- getSwarmDataMeans(swarmDat)

save(swarmDatMeans, file = "swarmDatMeans.RData")

swarmDatMeansFilt <- swarmDatMeans[match(gsub("\\_.*|(PA14).*", 
                                              rownames(ccmnNormMets), 
                                              rep = "\\1"), 
                                         swarmDatMeans$Strains), ]
swarmDatMeansFilt

pcaMets <- prcomp(ccmnNormMets)
fviz_eig(pcaMets)

tiff("pcaSwarmMetsArea.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmDatMeansFilt$Area,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Area")
dev.off()

tiff("pcaSwarmMetsAreaPercentage.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmDatMeansFilt$AreaPercentage,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Area Percentage")
dev.off()

tiff("pcaSwarmMetsCircularity.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmDatMeansFilt$Circularity,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Circularity")
dev.off()

tiff("pcaSwarmMetsPerimeter.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmDatMeansFilt$Perimeter,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Perimeter")
dev.off()
tiff("pcaSwarmMetsMaxLength.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmDatMeansFilt$MaxLength,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Maximum Length")
dev.off()

# Obtain bidirectional correlations between swarming variables.
corMat <- matrix(ncol = ncol(swarmDat) - 2, nrow = ncol(swarmDat) - 2)
rownames(corMat) <- colnames(swarmDat[, 2:(ncol(swarmDat) - 1)])
colnames(corMat) <- colnames(swarmDat[, 2:(ncol(swarmDat) - 1)])
for(i in 1:ncol(corMat)){
        for(j in 1:nrow(corMat)){
                corMat[i, j] <- cor(swarmDat[[rownames(corMat)[j]]],
                                    swarmDat[[colnames(corMat)[i]]])
        }
}
corMat




# Do PCA with the different measures of swarming. We first normalize them by dividing each column by its
# highest value. That way all values range from 0 to 1. THe idea is to get a linear combination of the
# measures that spreads most of the variability of the strains and that we could use as the response 
# variable for the supervised analysis. 

save(swarmDat, file = "swarmDat.RData")

m <- swarmDatMeans[, 2:8]

m <- apply(m, 2, function(x) x/max(x))
rownames(m) <- make.unique(gsub(".tif", replacement = "", swarmDatMeans$Strains))


pcaSwarms <- prcomp(m)


tiff("pcaSwarmMeasures.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaSwarms#, 
             #col.ind = swarmDatMeansFilt$MaxLength,
             #gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             #legend.title = "Maximum Length"
             )
dev.off()

tiff("biplotSwarmMeasures.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_biplot(pcaSwarms)
dev.off()


m <- swarmDatMeans[, c(4, 6)]

m <- apply(m, 2, function(x) x/max(x))
rownames(m) <- make.unique(gsub(".tif", replacement = "", swarmDatMeans$Strains))


pcaSwarms <- prcomp(m)


tiff("biplotSwarmMeasures_CircAreaPerc.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_biplot(pcaSwarms)
dev.off()


distManh <- dist(m, method = "canberra")

c <- cmdscale(distManh,
              eig=TRUE,
              x.ret=TRUE)

mds.var.per <- round(c$eig/sum(c$eig)*100, 1)
mds.var.per

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- c$points
mds.data <- data.frame(Sample=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data

ggplot(data=mds.data, aes(x=X, y=Y, label=Sample)) +
        geom_text() +
        theme_bw() +
        xlab(paste("MDS1 - ", mds.var.per[1], "%", sep="")) +
        ylab(paste("MDS2 - ", mds.var.per[2], "%", sep="")) +
        ggtitle("PCOA Swarm Data")

library(ggplot2)

pcomp1 <- pcaSwarms$x[, 1]
groups <- unique(gsub("\\_.*", replacement = "", names(pcomp1)))
pcomp1Means <- c()
for(i in seq_along(groups)){
        stMean <- mean(pcomp1[grep(groups[i], names(pcomp1))])
        pcomp1Means <- c(pcomp1Means, stMean)
}
names(pcomp1Means) <- groups



pcomp1MeansFilt <- pcomp1Means[names(pcomp1Means) %in% rownames(ccmnNormMets)]

match(gsub("\\_.*|(PA14).*", rownames(ccmnNormMets), rep = "\\1"), names(pcomp1Means))


pcomp1MeansFilt <- pcomp1Means[match(gsub("\\_.*|(PA14).*", 
                                              rownames(ccmnNormMets), 
                                              rep = "\\1"), 
                                         names(pcomp1Means))]

tiff("pcaSwarmMetsPrComp1.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = pcomp1MeansFilt,
             gradient.cols = c("#FF0000", "#66FF00", "#003CFF"),
             legend.title = "PrComp1")
dev.off()

pcomp1Means

save(pcomp1Means, file = "pcomp1AreaPercCircularity.RData")

# The metabolites that are more correlated with swarming, according to OPLS-DA, are succinate and formyl-methionine.

plot(log(swarmDatMeansFilt$AreaPercentage), ccmnNormMets$`Formyl-methionine`)
plot(log(swarmDatMeansFilt$AreaPercentage), ccmnNormMets$succinate)

metSwarm <- cbind.data.frame(ccmnNormMets, log(swarmDatMeansFilt$AreaPercentage))
colnames(metSwarm)[ncol(metSwarm)] <- "logAreaPct"
colnames(metSwarm) <- make.names(colnames(metSwarm))

# Let's do a PCA of the five metabolites that have highest loadings in OPLS-DA model to get the linear combination of these in 
# x-axis. 
pcaSelMets <- prcomp(ccmnNormMets[, c("Formyl-methionine", 
                                      "2-aminoadipate-6-semialdehyde", 
                                      "Methylcitrate 1", 
                                      "succinate", 
                                      "Acetylhomoserine 2"
)])

fviz_pca_biplot(pcaSelMets)

# Obtain values of PC1 and add to metSwarm dataframe

metSwarm <- cbind.data.frame(metSwarm, pcaSelMets$x[, 1])
colnames(metSwarm)[ncol(metSwarm)] <- "prComp1SelMets"

if(!require(ggpubr)) install.packages("ggpubr")
library(ggplot2)
if(!require(ggpmisc)) BiocManager::install("ggpmisc")
library(ggpmisc)

ggscatter(metSwarm, x = "logAreaPct", y = "Formyl.methionine", add = "reg.line") +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 10)) +
        stat_regline_equation(label.x = 3, label.y = 13)

ggscatter(metSwarm, x = "logAreaPct", y = "succinate", add = "reg.line", label = rownames(metSwarm)) +
        stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                     label.x = 9)) +
        stat_regline_equation(label.x = 3, label.y = 10)

ggscatter(metSwarm, x = "logAreaPct", y = "prComp1SelMets", add = "reg.line") +
       stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                    label.x = 10, label.y = 3)) +
       stat_regline_equation(label.x = 3)


cor(log(swarmDatMeansFilt$AreaPercentage), pcaSelMets$x[, 1])
plot(log(swarmDatMeansFilt$AreaPercentage), pcaSelMets$x[, 1])

cor(log(swarmDatMeansFilt$AreaPercentage), ccmnNormMets$`Formyl-methionine`)
cor(log(swarmDatMeansFilt$AreaPercentage), ccmnNormMets$succinate)

