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
corMat <- matrix(ncol = 5, nrow = 5)
rownames(corMat) <- colnames(swarmDat[, 2:6])
colnames(corMat) <- colnames(swarmDat[, 2:6])
for(i in 1:ncol(corMat)){
        for(j in 1:nrow(corMat)){
                corMat[i, j] <- cor(swarmDat[[rownames(corMat)[j]]],
                                    swarmDat[[colnames(corMat)[i]]])
        }
}
corMat




swarmDat

m <- swarmDat[, c(3, 4)]

m <- apply(m, 2, function(x) x/max(x))
rownames(m) <- make.unique(gsub(".tif", replacement = "", swarmDat$Strains))


pcaSwarms <- prcomp(m)

fviz_pca_ind(pcaSwarms, 
             #col.ind = swarmDatMeansFilt$MaxLength,
             #gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             #legend.title = "Maximum Length"
             )

fviz_pca_biplot(pcaSwarms)

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



pcomp1MeansFilt <- pcomp1Means[names(pcomp1Means) %in% rownames(ccmnNormMeds)]

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


