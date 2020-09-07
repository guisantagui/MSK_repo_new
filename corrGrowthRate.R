setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/corGrowthRate")

curveFeatures <- read.csv("C:/Users/Guillem/Documents/PhD/comput/scriptOrhologues/Source_code_for_Pseudomonas_Metabolomics_Paper-master/nnmf/tblgcfeatures.csv")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")

library(KEGGREST)


TCACpds <- gsub("cpd:", "", keggLink("compound", "map00020"))

rownames(ccmnNormMets) <- gsub("W70322", "W70332", rownames(ccmnNormMets))

specMaxGrowthRt <- curveFeatures[curveFeatures$phase == 1, ]$specific_growth_rate_max

names(specMaxGrowthRt) <- curveFeatures$strain[curveFeatures$phase == 1]

names(specMaxGrowthRt) <- gsub("M6057", "M6075", names(specMaxGrowthRt))

specMaxGrowthRtOrdered <- specMaxGrowthRt[match(gsub("\\_.*|(PA14).*", 
                                                     "\\1", rownames(ccmnNormMets)),
                                                names(specMaxGrowthRt))]

cTests <- c()
corrs <- c()
for(i in 1:ncol(ccmnNormMets)){
        cTest <- cor.test(ccmnNormMets[, i],
                          specMaxGrowthRtOrdered,
                          method = "spearman")
        cTests <- c(cTests, cTest$p.value)
        corr <- cor(ccmnNormMets[, i], 
                    specMaxGrowthRtOrdered,
                    method = "spearman")
        corrs <- c(corrs, corr)
}
names(cTests) <- colnames(ccmnNormMets)

cTestsPadj <- p.adjust(cTests, method = "BH")

corTestMat <- data.frame(cor = corrs,
                         pVals = cTests, 
                         pAdj = cTestsPadj)

rownames(corTestMat) <- names(cTests)

corTestMatSign <- corTestMat[corTestMat$pAdj < .05, ]

# Remove ambiguous compounds
corTestMat <- corTestMat[-grep("_?", rownames(corTestMat), fixed = T), ]

newNames <- dictionary$definitiveNames[match(rownames(corTestMat), dictionary$Consensus)]

rownames(corTestMat)[!is.na(newNames)] <- newNames[!is.na(newNames)]

# See correlation of compounds in TCA cycle with max growth rate in phase I

corTestMat[dictionary$`KEGG IDs`[match(rownames(corTestMat), dictionary$definitiveNames)] %in% TCACpds, ]
