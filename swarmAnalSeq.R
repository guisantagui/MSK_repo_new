####################################################################################################################
####################################################################################################################
#   Swarm-Sequence analysis
####################################################################################################################
####################################################################################################################
setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/swarmFunctions.R")
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/variantAnalysis_functions.R")
if(!require(DECIPHER)) install.packages("DECIPHER")
library(DECIPHER)
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/ParsedSubGraphs_allStrains_named.RData")
load("swarmMeans.RData")
load("swarmMeansRearr.RData")

# PA14 is not in the alignment, and M55212 yes (position eleven, so let's add it's value to the swarmMeans)

swarmMeansGood <- swarmMeansRearr[-1]
swarmMeansGood <- c(swarmMeansGood[1:10], swarmMeans[grep("M55212", names(swarmMeans))], swarmMeansGood[11:length(swarmMeansGood)])
swarmMeansGoodBin <- binarizeSwarm(swarmMeansGood, threshold = -1.7)

load("ORA_rfeResults.RData")
load("ORA_OPLSDAQual.RData")
load("ORA_OPLSDAQuant.RData")

load("tab_rfeSwarm.RData")
load("ORA_OPLSDAQual.RData")
load("tab_OPLSDAQuant.RData")

load("overlapORA.RData")
load("overlapFELLA.RData")
#BrowseSeqs(ParsedSubGraphs_allStrains_named$`Chromosomal replication initiator protein DnaA`$DNAAlignment)

diffGenesSwarmBin <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, 
                                       strainGroups = swarmMeansGoodBin,
                                       p_adjust = T, 
                                       method = "BY", 
                                       alpha = 0.05)

overlapFELLA

enzymesFELLAOPLSDAQuant <- paste(paste("\\(EC", tab_OPLSDAQuant[tab_OPLSDAQuant$Entry.type == "enzyme", ]$KEGG.id), "\\)", sep = "")

fellaEnzymesSignGenes <- data.frame()
for(i in seq_along(enzymesFELLAOPLSDAQuant)){
        fellaEnzymesSignGenes <- rbind.data.frame(fellaEnzymesSignGenes, diffGenesSwarmBin$significativeGenes_DNA[grep(enzymesFELLAOPLSDAQuant[i], 
                                                                                                                       diffGenesSwarmBin$significativeGenes_DNA$Annotation), ])
}

write.csv(fellaEnzymesSignGenes, file = "fellaEnzymesSignGenes.csv")

# Keep only the blocks that are not a mix of different annotations
fellaEnzymesSignGenesInteresting <- fellaEnzymesSignGenes[sapply(as.character(fellaEnzymesSignGenes$Annotation), nchar) < 100, ]

BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[1]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[2]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[3]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[4]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[5]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[6]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaEnzymesSignGenesInteresting$Annotation[7]]]$DNAAlignment)

variantsDiffGenesFellaOPLSDAQuant <- variantsPerGene(ParsedSubGraphsNamed_objct = ParsedSubGraphs_allStrains_named, 
                                                     IDDiffGenes = diffGenesSwarmBin, 
                                                     strainNames = names(swarmMeansGoodBin))
save(variantsDiffGenesFellaOPLSDAQuant, file = "variantsDiffGenesFellaOPLSDAQuant.RData")

 
variantsDiffGenesFellaOPLSDAQuantInteresting <- lapply(variantsDiffGenesFellaOPLSDAQuant, function(x) x[names(x) %in% fellaEnzymesSignGenesInteresting$Annotation])


variantsDiffGenesFellaOPLSDAQuant[[1]][names(variantsDiffGenesFellaOPLSDAQuant[[1]]) %in% fellaEnzymesSignGenesInteresting$Annotation]
