setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/solveAmbig")

# Load datasets
load(file = "C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/metabolitePeaks_oldData.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

# Obtain the molecular weights of the ambiguous compounds

molWght <- metabolitePeaks_preCur$Mass[match(colnames(ccmn_norm_mets_good_old), metabolitePeaks_preCur$`Compound Name`)]
names(molWght) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
molWght


ambigMolWght <- molWght[grep("_?", names(molWght), fixed = T)]

# Obtain all the compounds with that molecular weight from KEGG

possCpds <- sapply(ambigMolWght, keggFind, database = "compound", option = "exact_mass")
keggFind("compound", 129.04, "exact_mass")

possCpds <- possCpds[sapply(possCpds, function(x) length(x) != 0)]

# Obtain all the enzymes that are in reactions where each compound with that molecular appears, and see what of these 
# enzymes are in the annotation. Count for each compound the total number of times a enzyme that produces or consumes
# the candidate compound is in the genome.

enzListsAmbigComp <- list()
for(i in seq_along(possCpds)){
        enzsPossComps <- list()
        for(j in seq_along(possCpds[[i]])){
                enzs <- keggLink(target = "enzyme", source = gsub("cpd:", replacement = "", names(possCpds[[i]])[j]))
                enzs <- gsub("ec:", replacement = "", enzs)
                inAnnot <- enzs %in% dictEnzymes$ECnums
                names(inAnnot) <- enzs
                inAnnot <- c(inAnnot, sum(inAnnot))
                names(inAnnot)[length(inAnnot)] <- "sums"
                enzsPossComps[[j]] <- inAnnot
        }
        names(enzsPossComps) <- names(possCpds[[i]])
        enzListsAmbigComp[[i]] <- enzsPossComps
}
names(enzListsAmbigComp) <- names(possCpds)

# Identify, for each of the ambiguous compounds, the candidates that we have genetic evidence of their presence.  
# Keep the ones that are unique (only one of the candidates has genetic evidence)

uniqueCandidates <- unlist(lapply(enzListsAmbigComp, function(x){
        sums <- sapply(x, function(y) y[length(y)]) > 0
        if(sum(sums) <= 1){
                return(names(x)[sums])
                }
        })) 

save(uniqueCandidates, file = "uniqueCandidates.RData")



