####################################################################################################################
##### Figure S2 paper 
####################################################################################################################

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figureS2")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/logdatatable.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure2/solvedAmbigMets.RData")

solvedAmbigMets <- solvedAmbigMets[-4]

dictionary$Consensus[match(names(solvedAmbigMets), dictionary$Consensus)] <- solvedAmbigMets

metNames <- dictionary$definitiveNames[match(logdatatable$metabolite, dictionary$Consensus)]

metNames[is.na(metNames)] <- logdatatable$metabolite[is.na(metNames)]

logdatatable$metabolite <- metNames

logdatatable$ambiguous <- rep("Non-Ambiguous", length(logdatatable$metabolite))

logdatatable$ambiguous[grep("_?", dictionary$Consensus[match(logdatatable$metabolite, dictionary$definitiveNames)], fixed = T)] <- "Ambiguous"

logdatatable$legend <- logdatatable$ambiguous

logdatatable$legend[logdatatable$significant == 0] <- "Internal Standards"

logdatatable$color <- as.factor(logdatatable$color)


logdatatable$legend[logdatatable$legend == "Non-Ambiguous"] <- "Variable Metabolites"
logdatatable$legend[logdatatable$legend == "Internal Standards"] <- "Non-changing Metabolites (IS)"
logdatatable$legend[logdatatable$legend == "Ambiguous"] <- "Putative Metabolites"
logdatatable$legend <- as.factor(logdatatable$legend)

#Change names to the ones suggested by Kyu Rhee
KR_metNames <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/KRNewMetNames/metabolitesName-Xavier lab-KR.xlsx"))

KR_metNames$oldName[grep("Hexose phosphates", KR_metNames$oldName)] <- "Hexose phosphate"

logdatatable$metabolite <- KR_metNames$newName[match(logdatatable$metabolite, KR_metNames$oldName)]

metBarplot <- ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=function(x) median(x[!is.na(x)])), y = peak, color=legend)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Log Abundance",
                x= "Metabolites"
        ) +
        scale_color_manual(values = c("#000000", "red", "#33B0FF")) +
        scale_y_reverse() +
        theme(axis.text.y = element_text(angle = 180, hjust = 0, size = 11),
              axis.text.x = element_text(angle = 90, hjust = 1, size =  11)) 

metBarplot

ggsave("metabolite_boxplot_oldData.pdf",width = 6.5, height = 10.2)  

metBarplot <- ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=function(x) median(x[!is.na(x)])), y = peak, color=legend)) +
        geom_boxplot() +
        #coord_flip() +
        labs(
                y= "Log Abundance",
                x= "Metabolites"
        ) +
        scale_color_manual(values = c("#000000", "red", "#33B0FF")) +
        #scale_y_reverse() +
        theme(axis.text.y = element_text(angle = 0, hjust = 0, size = 11),
              axis.text.x = element_text(angle = 90, hjust = 1, size =  11, vjust = 0)) 

metBarplot

ggsave("metabolite_boxplot_oldData.pdf",width = 13, height = 4)  

write.csv(logdatatable, file = "logdatatable.csv")

logdatatable$metabolite[match(dictionary$definitiveNames[grep("?", dictionary$Consensus, fixed = T)], logdatatable$metabolite)]

logdatatable_noAmbig <- logdatatable[logdatatable$legend != "Ambiguous", ]

logdatatable_noAmbig$legend[logdatatable_noAmbig$legend == "Non-Ambiguous"] <- "Variable Metabolites"
logdatatable_noAmbig$legend[logdatatable_noAmbig$legend == "Internal Standards"] <- "Non-changing Metabolites (IS)"

metBarplot_noAmbig <- ggplot(data = logdatatable_noAmbig, mapping = aes(x = reorder(metabolite, peak, FUN=function(x) median(x[!is.na(x)])), y = peak, color=legend)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Log Abundance",
                x= "Metabolites"
        ) +
        scale_color_manual(values = c("#000000", "#33B0FF")) +
        scale_y_reverse() +
        theme(axis.text.y = element_text(angle = 180, hjust = 0, size = 11),
              axis.text.x = element_text(angle = 180, hjust = 1, size = 11)) 

metBarplot_noAmbig

ggsave("metabolite_boxplot_oldData_noAmbig.pdf",width = 7.5, height = 13)

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

pdf("biplotSwarmMeasures.pdf", height = 8, width = 8)
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
