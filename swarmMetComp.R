# This script takes the output of the sequence analysis done to look for correlations between differences in 
# sequence between strains in cluster 1 and in cluster 2, and compares it to swarming predictors according to 
# OPLS-DA


setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")
load("OPLSDAQuantResultKEGGIDs.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/outPut1_2_old.RData")
signSeqGenesOPLSDApreds <- outPut1_2_old$DNA_enz[-1, colnames(outPut1_2_old$DNA_enz) %in% OPLSDAQuantResultKEGGIDs]
# Keep only the ones that are not empty
signSeqGenesOPLSDApreds <- signSeqGenesOPLSDApreds[, apply(signSeqGenesOPLSDApreds, 2, function(x) !all(is.na(x)))]
#View(signSeqGenesOPLSDApreds)
sapply(signSeqGenesOPLSDApreds)

signSeqGenesOPLSDApredsInteresting <- list()
for(i in 1:ncol(signSeqGenesOPLSDApreds)){
        interesting <- signSeqGenesOPLSDApreds[which(nchar(signSeqGenesOPLSDApreds[, i]) < 150), i]
        signSeqGenesOPLSDApredsInteresting[[i]] <- interesting
}
names(signSeqGenesOPLSDApredsInteresting) <- colnames(signSeqGenesOPLSDApreds)
signSeqGenesOPLSDApredsInterestingCounts <- table(unlist(signSeqGenesOPLSDApredsInteresting))
dim(signSeqGenesOPLSDApredsInterestingCounts)
write.csv(signSeqGenesOPLSDApredsInterestingCounts, file = "signSeqGenesOPLSDApredsInterestingCounts.csv")

load("tab_OPLSDAQuant.RData")
tab_OPLSDAQuant

names(signSeqGenesOPLSDApredsInterestingCounts)

OPLSDAQuantEnzymes <- tab_OPLSDAQuant$KEGG.id[tab_OPLSDAQuant$Entry.type == "enzyme"]

for(i in seq_along(OPLSDAQuantEnzymes)){
        print(grep(OPLSDAQuantEnzymes[i], names(signSeqGenesOPLSDApredsInterestingCounts)))
}
grep(OPLSDAQuantEnzymes[34], names(signSeqGenesOPLSDApredsInterestingCounts))

signSeqGenesOPLSDApredsInterestingCounts
