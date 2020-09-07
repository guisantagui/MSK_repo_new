####################################################################################################################
##### Figure S2 paper 
####################################################################################################################

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figureS2")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/logdatatable.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure2/solvedAmbigMets.RData")

metNames <- dictionary$definitiveNames[match(logdatatable$metabolite, dictionary$`Old Data Names`)]


logdatatable$metabolite <- metNames


ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=function(x) -median(x)), y = peak, color=significant)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Abundance",
                x= "Metabolite"
        )

ggsave("metabolite_boxplot_oldData.pdf",width = 7, height = 11)  

logdatatable$metabolite[match(dictionary$definitiveNames[grep("?", dictionary$Consensus, fixed = T)], logdatatable$metabolite)]


# Biplots swarm measures

if(!require(factoextra)) install.packages("factoextra")
library(factoextra)

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/swarmDatMeans.RData")

swarmDat <- swarmDatMeans

m <- swarmDat[, 3:8]

m <- apply(m, 2, function(x) x/max(x))
rownames(m) <- make.unique(gsub(".tif", replacement = "", swarmDat$Strains))


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

pdf("biplotSwarmMeasures.pdf", height = 10, width = 10)
fviz_pca_biplot(pcaSwarms)
dev.off()

m <- swarmDat[, c(4, 6)]

m <- apply(m, 2, function(x) x/max(x))
rownames(m) <- make.unique(gsub(".tif", replacement = "", swarmDat$Strains))


pcaSwarms <- prcomp(m)


tiff("biplotSwarmMeasures_CircAreaPerc.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_biplot(pcaSwarms)
dev.off()

pdf("biplotSwarmMeasures_CircAreaPerc.pdf", height = 10, width = 10)
fviz_pca_biplot(pcaSwarms)
dev.off()


# Biplot of means of swarming measures

expers <- unique(swarmDat$exp)


# Do PCA by hand for obtaining eigenvector mat, and obtain formula to get swarm 
# score from max length and circularity.

pcaSwarms2 <- prcomp(m)

mCent <- apply(m, 2, function(x) x - mean(x))

covMat <- cov(apply(m, 2, function(x) x - mean(x)))
covMat <- cov(m)

eigenvalues <- eigen(covMat)$values
eigenvectors <- eigen(covMat)$vectors

pcaSwarms$rotation[2, 2] == eigenvectors[2, 2]
        
PC <- mCent %*% eigenvectors   

eigenvectors

cov(PC)      
        
eigenvalues

# This is the formula to get swarming score
apply(m, 1, function(x) (x[1] - pcaSwarms$center[1])*-pcaSwarms$rotation[1, 1] + (x[2] - pcaSwarms$center[2])*-pcaSwarms$rotation[2, 1])
