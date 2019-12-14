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
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/geneCallsAllStrains4Heron.RData")
# PA14 is not in the alignment, and M55212 yes (position eleven, so let's add it's value to the swarmMeans)

swarmMeansGood <- swarmMeansRearr[-1]
swarmMeansGood <- c(swarmMeansGood[1:10], swarmMeans[grep("M55212", names(swarmMeans))], swarmMeansGood[11:length(swarmMeansGood)])
swarmMeansGoodBin <- binarizeSwarm(swarmMeansGood, threshold = -1.7)
save(swarmMeansGoodBin, file = "swarmMeansGoodBin.RData")

load("ORA_rfeResults.RData")
load("ORA_OPLSDAQual.RData")
load("ORA_OPLSDAQuant.RData")

load("tab_rfeSwarm.RData")
load("tab_OPLSDAQual.RData")
load("tab_OPLSDAQuant.RData")

load("overlapORA.RData")
load("overlapFELLA.RData")

load("rfeResultKEGGIDs.RData")
load("OPLSDAQualResultKEGGIDs.RData")
load("OPLSDAQuantResultKEGGIDs.RData")
#BrowseSeqs(ParsedSubGraphs_allStrains_named$`Chromosomal replication initiator protein DnaA`$DNAAlignment)

diffGenesSwarmBin <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, 
                                       strainGroups = swarmMeansGoodBin,
                                       p_adjust = T, 
                                       method = "BY", 
                                       alpha = 0.05)

overlapFELLA

getEnzymeIDFELLA <- function(fellaTab){
        enzymeIDs <- paste(paste("\\(EC", fellaTab[fellaTab$Entry.type == "enzyme", ]$KEGG.id), "\\)", sep = "")
}
getEnzymeIDORA <- function(oraTab){
        if(!require(KEGGREST)) install.packages("KEGGREST")
        library(KEGGREST)
        subPaths <- gsub(x = rownames(ORA_rfeResults), pattern = "path:pae", replacement = "map")
        allEnzs <- unique(unlist(sapply(subPaths, keggLink, target = "enzyme")))
        allEnzs <- paste(paste("\\(EC", gsub(allEnzs, pattern = "ec:", replacement = "EC ")), "\\)", sep = "")
        return(allEnzs)
}
enzymesFELLAOPLSDAQuant <- getEnzymeIDFELLA(tab_OPLSDAQuant)
enzymesFELLAOPLSDAQual <- getEnzymeIDFELLA(tab_OPLSDAQual)
enzymesORAOPLSDAQuant <- getEnzymeIDORA(ORA_OPLSDAQuant)

getRelEnzs <- function(IDDiffGenes, enz){
        relEnz <- data.frame()
        for(i in seq_along(enz)){
                relEnz <- rbind.data.frame(relEnz, IDDiffGenes$significativeGenes_DNA[grep(enz[i],
                                                                                           IDDiffGenes$significativeGenes_DNA$Annotation), ])
        }
        return(relEnz)
}

fellaSignEnzymesOPLSDAQuant <- getRelEnzs(diffGenesSwarmBin, enzymesFELLAOPLSDAQuant)

write.csv(fellaEnzymesSignGenes, file = "fellaEnzymesSignGenes.csv")

# Keep only the blocks that are not a mix of different annotations
getInteresting <- function(relEnzs, threshold = 100){
        relEnzsInter <- relEnzs[sapply(as.character(relEnzs$Annotation), nchar) < 100, ]
        KEGGIDs <- sapply(relEnzsInter$Annotation, function(x) sub("\\).*", "", sub(".*\\(", "", x)))
        KEGGIDs <- gsub("EC ", "ec:", KEGGIDs)
        relEnzsInter <- cbind.data.frame(KEGGIDs, relEnzsInter)
        colnames(relEnzsInter)[1] <- "KEGGIDs"
        return(relEnzsInter)
}

fellaSignEnzymesOPLSDAQuantInteresting <- getInteresting(fellaSignEnzymesOPLSDAQuant)
write.csv(fellaSignEnzymesOPLSDAQuantInteresting, file = "fellaSignEnzymesOPLSDAQuantInteresting.csv")



BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[1]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[2]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[3]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[4]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[5]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[6]]]$DNAAlignment)
BrowseSeqs(ParsedSubGraphs_allStrains_named[[fellaSignEnzymesOPLSDAQuantInteresting$Annotation[7]]]$DNAAlignment)

variantsDiffGenesFellaOPLSDAQuant <- variantsPerGene(ParsedSubGraphsNamed_objct = ParsedSubGraphs_allStrains_named, 
                                                     IDDiffGenes = diffGenesSwarmBin, 
                                                     strainNames = names(swarmMeansGoodBin))
save(variantsDiffGenesFellaOPLSDAQuant, file = "variantsDiffGenesFellaOPLSDAQuant.RData")
load("variantsDiffGenesFellaOPLSDAQuant.RData")
 
variantsDiffGenesFellaOPLSDAQuantInteresting <- lapply(variantsDiffGenesFellaOPLSDAQuant, function(x) x[names(x) %in% fellaSignEnzymesOPLSDAQuantInteresting$Annotation])

cpdFellaSignEnzymesOPLSDAQuant <- sapply(fellaSignEnzymesOPLSDAQuantInteresting$KEGGIDs, function(x) gsub("cpd:", replacement = "", keggLink(x, target = "compound")))
lapply(cpdFellaSignEnzymesOPLSDAQuant, function(x) x[x %in% OPLSDAQuantResultKEGGIDs])
cpdFellaSignEnzymesOPLSDAQuant[[1]][cpdFellaSignEnzymesOPLSDAQuant[[1]] %in% OPLSDAQuantResultKEGGIDs]

variantsDiffGenesFellaOPLSDAQuantInteresting
