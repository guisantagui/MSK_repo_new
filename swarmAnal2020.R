setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")

swarmDat <- read.csv("/Users/santamag/Desktop/GUILLEM/pics/allSwarms.csv")
swarmDat$Strains <- gsub("H5707", "H5708", swarmDat$Strains)
swarmDat$Strains <- gsub("M6057", "M6075", swarmDat$Strains)

# PCA for looking what variable among the ones we used for measuring the swarming is most separated in metabolites PCA.
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(readxl)) install.packages("readxl")
library(readxl)

# Load metabolomics data and dictionary
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
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
                means <- apply(subMat[, 2:6], 2, mean)
                meanDF <- rbind.data.frame(meanDF, means)
        }
        meanDF <- cbind.data.frame(groups, meanDF)
        colnames(meanDF) <- colnames(swarmMat)[1:6]
        return(meanDF)
}
swarmDatMeans <- getSwarmDataMeans(swarmDat)

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
corMat <- matrix(ncol = 4, nrow = 4)
rownames(corMat) <- colnames(swarmDat[, 2:5])
colnames(corMat) <- colnames(swarmDat[, 2:5])
for(i in 1:ncol(corMat)){
        for(j in 1:nrow(corMat)){
                corMat[i, j] <- cor(swarmDat[[rownames(corMat)[j]]],
                                    swarmDat[[colnames(corMat)[i]]])
        }
}
corMat










