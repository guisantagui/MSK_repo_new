####################################################################################################################
####################################################################################################################
#   Swarm-Sequence analysis
####################################################################################################################
####################################################################################################################
setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
source("swarmFunctions.R")
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/variantAnalysis_functions.R")
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
BrowseSeqs(ParsedSubGraphs_allStrains_named$`Chromosomal replication initiator protein DnaA`$DNAAlignment)


diffGenesSwarmBin <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, 
                                       strainGroups = swarmBin,
                                       p_adjust = T, 
                                       method = "BY", 
                                       alpha = 0.05)
