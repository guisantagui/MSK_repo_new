setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary")

load("../normMetAnal/newData/metabolitePeaks_newData.RData")
load("../normMetAnal/oldDataGood/metabolitePeaks_oldData.RData")

metabolitePeaks_newData$`Compound Name`[!metabolitePeaks_newData$`Compound Name` %in% metabolitePeaks_preCur$`Compound Name`]

metabolitePeaks <- metabolitePeaks_newData

metPeaksOld_notInNew <- metabolitePeaks_preCur[!metabolitePeaks_preCur$`Compound Name` %in% metabolitePeaks_newData$`Compound Name`, ]

##################################################################################################################################################################
# METABOLITES THAT APPEAR IN NEW DATA BUT NOT IN OLD #############################################################################################################
##################################################################################################################################################################

# 2-Dehydro-3-deoxy-L-rhamnoate 3 ---------> Doesn't appear. Added in the new data
# 2-Oxobutanoate 2 ------------------------> Doesn't appear Redundant: same as succinate semialdehyde
# Acetyl-2, 6-diaminopimelate -------------> Appears, but the space makes it not visible
# alanine 1 -------------------------------> Doesn't appear. The closest one is b-alanine 2: Let's check if the numbers are the same
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "alanine 1"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "b-alanine 2"), 5:10] # Similar numbers, but not the same. Maybe they changed the name...
# Aminobutyraldehyde ----------------------> Doesn't appear. Added in the new data
# Aspartate 4-semialdehyde ----------------> Doesn't appear. Added in the new data
# Deoxyribose 2 ---------------------------> Doesn't appear. Maybe Dihydroxy-isovalerate 1? Similar values and same trend. Mann Whitney says are the same.
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Deoxyribose 2"), 2:ncol(metabolitePeaks)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Dihydroxy-isovalerate 1"), 5:(ncol(metabolitePeaks_preCur) - 1)]

wilcox.test(as.numeric(metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Deoxyribose 2"), 2:(ncol(metabolitePeaks) - 1)]),
            as.numeric(metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Dihydroxy-isovalerate 1"), 
                                              5:(ncol(metabolitePeaks_preCur) - 1)]))
# Galactono-1,4-lactone -------------------> Appears, but is mispelled
# gliceraldehyde3p ------------------------> Doesn't appear. --> New in the new file
# hexose diphosphate ----------------------> Doesn't appear. The closest one is hexose phosphates. Let's check the numbers. 
#                                            --> Added new in the new file.
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "hexose diphosphate"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Hexose phosphates"), 5:10] # Completely different numbers
# Methylcitrate 1 -------------------------> Doesn't appear. Maybe Homocitrate?
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Methylcitrate 1"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Homocitrate 2"), 5:(ncol(metabolitePeaks_preCur) - 1)]

wilcox.test(as.numeric(metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Methylcitrate 1"), 2:(ncol(metabolitePeaks) - 1)]),
            as.numeric(metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Homocitrate 2"), 
                                              5:(ncol(metabolitePeaks_preCur) - 1)]))
# pentose 5-phosphate ---------------------> Doesn't appear Added new in the new data
# Sedoheptulose  7-phosphate 1 ------------> Doesn't appear Added in the new data


metPeaksOld_notInNew$`Compound Name`

##################################################################################################################################################################
# METABOLITES THAT APPEAR IN OLD DATA BUT NOT IN NEW #############################################################################################################
##################################################################################################################################################################

# -Acetyl-L-glutamate 5-semialdehyde ------> Removed from old data
# 1-Pyrroline-4-hydroxy-2-carboxylate 3 ---> In new data It's "Pyrroline-3-hydroxy-5-carboxylate 2", which also appears in old data, 
#                                            but not the same --> Removed from old
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "1-Pyrroline-4-hydroxy-2-carboxylate 3"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Pyrroline-3-hydroxy-5-carboxylate 2"), 2:7]
# 2-Acetyl-aminoadipate -------------------> Removed from old data
# 2-Aminoadipate 1 ------------------------> Same formula than "Acetylhomoserine", but not the same---> Removed from old data.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2-Acetyl-aminoadipate"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Acetylhomoserine 2"), 2:7]
# 2,6-Diaminopelate 2 ---------------------> Same formula than "2,6-Diaminoheptanedioate 1". ---> Same, values, but in Old data 2,6-Diaminoheptanedioate 1 
#                                            also appears, with a lot of missing Values. According to KEGG those compounds are stereoisomers
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2,6-Diaminopelate 2"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "2,6-Diaminoheptanedioate 1"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2,6-Diaminoheptanedioate 1"), ]
# 2'5-Dioxopentanoate 2 -------------------> Same formula than "methylmaleate 1". --> Here is weird, because they are not the same but values in 20
#                                            samples are the same. They are structural isomers. Maybe the peak values got mixed?
#                                            In old data 2'5-Dioxopentanoate 2 and methylmaleate 1 are redundant.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2'5-Dioxopentanoate 2"), 5:(ncol(metabolitePeaks_preCur)- 1)] == 
        metabolitePeaks[which(metabolitePeaks$`Compound Name` == "methylmaleate 1"), 2:(ncol(metabolitePeaks) - 1)]

metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2'5-Dioxopentanoate 2"), 5:(ncol(metabolitePeaks_preCur) - 1)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "methylmaleate 1"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "methylmaleate 1"), 5:(ncol(metabolitePeaks_preCur) - 1)]
# 3-Methyl-2-oxopentanoate 1 --------------> Removed from old data.
# 4-Acetylaminobutanal --------------------> Removed from old data.
# 4-Amino-4-cyanobutanoic acid ------------> Removed from old data.
# 4-Aminobutanoate 3 ----------------------> Same formula than "3-Amino-isobutanoate 1". Same values, but "3-Amino-isobutanoate 1" also appears in 
#                                            old data, and with different values.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "4-Aminobutanoate 3"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "3-Amino-isobutanoate 1"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "3-Amino-isobutanoate 1"), 5:10]
# 5-Oxo-proline 1 -------------------------> Same formula than "Pyrroline-3-hydroxy-5-carboxylate 2". Slightly different values, but same trend, 
#                                            but "Pyrroline-3-hydroxy-5-carboxylate 2" also appears in the old dataset, with different values 
#                                            than the new dataset has.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "5-Oxo-proline 1"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Pyrroline-3-hydroxy-5-carboxylate 2"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Pyrroline-3-hydroxy-5-carboxylate 2"), 5:10]

wilcox.test(as.numeric(metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "5-Oxo-proline 1"), 5:(ncol(metabolitePeaks_preCur) - 1)]),
            as.numeric(metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Pyrroline-3-hydroxy-5-carboxylate 2"), 2:(ncol(metabolitePeaks) -1)]))
# Acetolactate ----------------------------> Removed from old data.              
# Acetyl-2,6-diaminopimelate --------------> Mispelling, new dataset has a space between the 2 and the 6.
# Acetylhomoserine ------------------------> Same formula than "Acetylhomoserine 2". --> Reundant: Acetylhomoserine 2 is also in old data, and 
#                                            the values between the 3 vectors are only slightly different.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Acetylhomoserine"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Acetylhomoserine 2"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Acetylhomoserine 2"), 5:10]
# Adenosine 2 -----------------------------> Same formula than "deoxyguanosine 1". Values are the same between two datasets, but "deoxyguanosine 1" 
#                                            also appears in the old dataset, but with different values. The 2 compounds are structural isomers.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Adenosine 2"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "deoxyguanosine 1"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "deoxyguanosine 1"), 5:10]
# b-alanine 2 -----------------------------> Same formula than "alanine 1". ---> Same compound: slightly different values, and "alanine 1" doesn't
#                                            appear in old dataset.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "b-alanine 2"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "alanine 1"), 2:7]
# cis-2-Methylaconitate 1 -----------------> Same formula than "Homo-cisaconitate 2". ---> Values between the 2 compounds in the 2 datasets are the same, 
#                                            but "Homo-cisaconitate 2" also appears in the old dataset with different values.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "cis-2-Methylaconitate 1"), 5:10]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Homo-cisaconitate 2"), 2:7]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Homo-cisaconitate 2"), 5:10]
# Dihydroxy-isovalerate 1 -----------------> Same formula than "Deoxyribose 2". ---> Probably same compound, slightly different values, but
#                                            following the same trend.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Dihydroxy-isovalerate 1"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Deoxyribose 2"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Deoxyribose 2"), 5:10] # ---> Deoxyribose 2 doesn't appear in 2nd dataset
# Dimethylmalate 1 ------------------------> Same formula than "2-Dehydro-3-deoxy-L-rhamnonate 3". ---> Different compounds, 
#                                            Dimethylmalate 1 removed from old dataset.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Dimethylmalate 1"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "2-Dehydro-3-deoxy-L-rhamnonate 3"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2-Dehydro-3-deoxy-L-rhamnonate 3"), 5:10]
# Ethyl malate 2 --------------------------> Same formula than "2-Dehydro-3-deoxy-L-rhamnonate 3". --> Different values, not so slightly, but the trend
#                                            is the same. Those 2 compounds are structural isomers. We can apply an statistical test to see if the vectors 
#                                            are or not the same.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Ethyl malate 2"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "2-Dehydro-3-deoxy-L-rhamnonate 3"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "2-Dehydro-3-deoxy-L-rhamnonate 3"), 5:10] # --> "2-Dehydro-3-deoxy-L-rhamnonate 3"
#                                            doesn't appear in old data.
wilcox.test(as.numeric(metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Ethyl malate 2"), 5:(ncol(metabolitePeaks_preCur) - 1)]),
            as.numeric(metabolitePeaks[which(metabolitePeaks$`Compound Name` == "2-Dehydro-3-deoxy-L-rhamnonate 3"), 2:(ncol(metabolitePeaks))])) #--> Different

        
# fumarate --------------------------------> Removed from old data.
# Galactono-1'4-lactone -------------------> Mismpelling: the old data has an apostrophe separating 1 and 4 instead of a comma.
# Glucose ---------------------------------> Removed from old data.
# Glutamylcysteine ------------------------> Removed from old data.
# Guanosine -------------------------------> Removed from old data.
# Homocitrate 2 ---------------------------> Same formula than "Methylcitrate 1". The values are not the same, but the trend is similar, though not the same.
#                                            Apply an statistical test to see if they are the same.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Homocitrate 2"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Methylcitrate 1"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Methylcitrate 1"), 5:10]
# homoserine 2 ----------------------------> Same formula than "Threonine 1". --> Different compounds: "homoserine 2" removed from old dataset
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "homoserine 2"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Threonine 1"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Threonine 1"), 5:10] # --> Appears in old dataset
# Homovanillate ---------------------------> Removed from old data.    
# Hydroxyectoine --------------------------> Removed from old data.
# Methylglyoxal ---------------------------> Removed from old data.
# N-acetyl-2,4-diaminobutanoate 1 ---------> Same formula than "d-ala-d-ala 2". Same values in both metabolites and datasets, but "d-ala-d-ala 2"
#                                            also appears in old dataset, with different values. Those compounds are structural isomers. 
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "N-acetyl-2,4-diaminobutanoate 1"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "d-ala-d-ala 2"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "d-ala-d-ala 2"), 5:10]
# PABA 2 ----------------------------------> Same formula than "Anthranilate 1". Different values: PABA2 removed from old dataset.
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "PABA 2"), 5:ncol(metabolitePeaks_preCur)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Anthranilate 1"), 2:(ncol(metabolitePeaks) - 1)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Anthranilate 1"), 5:10]
# phophoenolpyruvate ----------------------> Removed from old dataset.
# serine ----------------------------------> Removed from old dataset.
# Tartronate semialdehyde -----------------> Removed from old dataset.
# Tyramine --------------------------------> Same formula than "Norepinephrine", but different values
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Tyramine"), 5:(ncol(metabolitePeaks_preCur) - 1)]
metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Norepinephrine"), 2:ncol(metabolitePeaks)]
metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Norepinephrine"), 5:10]

wilcox.test(as.numeric(metabolitePeaks_preCur[which(metabolitePeaks_preCur$`Compound Name` == "Tyramine"), 5:(ncol(metabolitePeaks_preCur) - 1)]),
            as.numeric(metabolitePeaks[which(metabolitePeaks$`Compound Name` == "Norepinephrine"), 2:ncol(metabolitePeaks)])) #--> Different



# Let's build dictionary

dictionary <- cbind(as.character(metabolitePeaks$`Compound Name`[metabolitePeaks$`Compound Name` %in% metabolitePeaks_preCur$`Compound Name`]),
                    as.character(metabolitePeaks$`Compound Name`[metabolitePeaks$`Compound Name` %in% metabolitePeaks_preCur$`Compound Name`]))

colnames(dictionary) <- c("Old Data Names", "New Data Names")

presNew <- cbind(rep(NA, length(metabolitePeaks$`Compound Name`[!metabolitePeaks$`Compound Name` %in% metabolitePeaks_preCur$`Compound Name`])),
                 as.character(metabolitePeaks$`Compound Name`[!metabolitePeaks$`Compound Name` %in% metabolitePeaks_preCur$`Compound Name`]))


dictionary <- rbind(dictionary, presNew)
dictionary <- cbind(dictionary, rep(NA, nrow(dictionary)))
colnames(dictionary)[3] <- "State"
dictionary[1:60, 3] <- rep("Same", 60)
dictionary[65, 3] <- "Added"
dictionary[66, 1] <- "succinate semialdehyde 1"
dictionary[66, 3] <- "Redundant"
dictionary[61, 1] <- "Acetyl-2,6-diaminopimelate"
dictionary[61, 3] <- "Misspelled"
dictionary[62, 1] <- "b-alanine 2"
dictionary[62, 3] <- "Renamed"
dictionary[63, 3] <- "Added"
dictionary[64, 3] <- "Added"
dictionary[67, 1] <- "Dihydroxy-isovalerate 1"
dictionary[67, 3] <- "Renamed"
dictionary[68, 1] <- "Galactono-1'4-lactone"
dictionary[68, 3] <- "Misspelled"
dictionary[69, 3] <- "Added"
dictionary[70, 3] <- "Added"
dictionary[71, 1] <- "Homocitrate 2"
dictionary[71, 3] <- "Renamed"
dictionary[72, 3] <- "Added"
dictionary[73, 3] <- "Added"

presOld <- cbind(as.character(metabolitePeaks_preCur$`Compound Name`[!metabolitePeaks_preCur$`Compound Name` %in% metabolitePeaks$`Compound Name`]),
                 rep(NA, length(metabolitePeaks_preCur$`Compound Name`[!metabolitePeaks_preCur$`Compound Name` %in% metabolitePeaks$`Compound Name`])),
                 rep(NA, length(metabolitePeaks_preCur$`Compound Name`[!metabolitePeaks_preCur$`Compound Name` %in% metabolitePeaks$`Compound Name`])))

dictionary <- rbind(dictionary, presOld)

dictionary[74, 3] <- "Removed"
dictionary[75, 3] <- "Removed"
dictionary[76, 3] <- "Removed"
dictionary[77, 3] <- "Removed"
dictionary[78, 2] <- "2,6-Diaminoheptanedioate 1"
dictionary[78, 3] <- "Renamed, and duplicated removed"
dictionary[79, 2] <- "methylmaleate 1"
dictionary[79, 3] <- "Redundant"
dictionary[80, 3] <- "Removed"
dictionary[81, 3] <- "Removed"
dictionary[82, 3] <- "Removed"
dictionary[83, 2] <- "3-Amino-isobutanoate 1"
dictionary[83, 3] <- "Renamed, and duplicated removed"
dictionary[84, 2] <- "Pyrroline-3-hydroxy-5-carboxylate 2"
dictionary[84, 3] <- "Renamed, and duplicated removed"
dictionary[85, 3] <- "Removed"
dictionary <- dictionary[-86, ]
dictionary[86, 2] <- "Acetylhomoserine 2"
dictionary[86, 3] <- "Redundant"
dictionary[87, 2] <- "deoxyguanosine 1"
dictionary[87, 3] <- "Renamed, and duplicated removed"
dictionary <- dictionary[-88, ]
dictionary[88, 2] <- "Homo-cisaconitate 2"
dictionary[88, 3] <- "Renamed, and duplicated removed"
dictionary <- dictionary[-89, ]
dictionary[89, 3] <- "Removed"
dictionary[90, 3] <- "Removed"
dictionary[91, 3] <- "Removed"
dictionary <- dictionary[-92, ]
dictionary[92, 3] <- "Removed"
dictionary[93, 3] <- "Removed"
dictionary[94, 3] <- "Removed"
dictionary <- dictionary[-95, ]
dictionary[95, 3] <- "Removed"
dictionary[96:98, 3] <- "Removed"
dictionary[99, 2] <- "d-ala-d-ala 2"
dictionary[99, 3] <- "Renamed, and duplicated removed"
dictionary[100:104, 3] <- "Removed"


curate <- dictionary[!ifelse(dictionary[, 3] %in% c("Redundant", "Renamed, and duplicated removed"), TRUE, FALSE), ]
# This object has a total number f 95 metabolites, but should have 96. Let's check what happened: I was not accounting
# for the redundant metabolites

curateOld <- curate[, 1] 
for(i in seq_along(curateOld)){
        if(is.na(curateOld[i])){
                curateOld[i] <- curate[i, 2]
        }
}
metabolitePeaks_preCur$`Compound Name`[!metabolitePeaks_preCur$`Compound Name` %in% curateOld]

dictionary <- cbind(dictionary, rep(NA, nrow(dictionary)))
colnames(dictionary)[4] <- "Consensus"


dictionary[, 4][ifelse(dictionary[, 3] %in% c("Same", "Added", "Misspelled", "Renamed", "Renamed, and duplicated removed"), 
                       T, F)] <- dictionary[, 2][ifelse(dictionary[, 3] %in% c("Same", 
                                                                               "Added", 
                                                                               "Misspelled",
                                                                               "Renamed",
                                                                               "Renamed, and duplicated removed"), T, F)]

dictionary[, 4][which(dictionary[, 3] == "Removed")] <- paste(dictionary[, 1][which(dictionary[, 3] == "Removed")],
                                                              "?", sep = "_")

dictionary[, 4][dictionary[, 1] %in% 
                        dictionary[, 2][dictionary[, 3] == 
                                                "Renamed, and duplicated removed"]] <- paste(dictionary[, 1][dictionary[, 1] %in% dictionary[, 2][dictionary[, 3] == 
                                                                                                                                                          "Renamed, and duplicated removed"]], 
                                                                                             "?", sep = "_")
dictionary <- cbind(dictionary, rep(NA, nrow(dictionary)))
colnames(dictionary)[5] <- "Formula"

for(i in 1:nrow(dictionary)){
        if(gsub("_?", "", dictionary[i, 4], fixed = T) %in% metabolitePeaks$`Compound Name`){
                dictionary[i, 5] <- as.character(metabolitePeaks$Formula[gsub("_?", "", 
                                                                              dictionary[i, 4], 
                                                                              fixed = T) == metabolitePeaks$`Compound Name`])
        }else if(gsub("_?", "", dictionary[i, 4], fixed = T) %in% metabolitePeaks_preCur$`Compound Name`){
                dictionary[i, 5] <- gsub("[[:blank:]]", "", metabolitePeaks_preCur$Formula[gsub("_?", "", 
                                                                                                dictionary[i, 4], 
                                                                                                fixed = T) == metabolitePeaks_preCur$`Compound Name`])
        }
}

dictionary <- cbind(dictionary, rep(NA, nrow(dictionary)))
colnames(dictionary)[6] <- "KEGG IDs"

dictionary[1, 6] <- "C01234"
dictionary[2, 6] <- "C03871"
dictionary[3, 6] <- "C02218"
dictionary[4, 6] <- "C04076"
dictionary[5, 6] <- "C00666"
dictionary[6, 6] <- "C03284*"
dictionary[7, 6] <- "C00170"
dictionary[8, 6] <- "C01077"
dictionary[9, 6] <- "C02713"
dictionary[10, 6] <- "C00437"
dictionary[11, 6] <- "C00147"
dictionary[12, 6] <- "C04677"
dictionary[13, 6] <- "C05714"
dictionary[14, 6] <- "C00108"
dictionary[15, 6] <- "C00062"
dictionary[16, 6] <- "C00049"
dictionary[17, 6] <- "C00438"
dictionary[18, 6] <- "C00327"
dictionary[19, 6] <- "C00097"
dictionary[20, 6] <- "C00993*"
dictionary[21, 6] <- "C00330*"
dictionary[22, 6] <- "C05947"
dictionary[23, 6] <- "C03145"
dictionary[24, 6] <- "C00037"
dictionary[25, 6] <- "C00135"
dictionary[26, 6] <- "C00262"
dictionary[27, 6] <- "C00123"
dictionary[28, 6] <- "C00047"
dictionary[29, 6] <- "C00547"
dictionary[30, 6] <- "C00077"
dictionary[31, 6] <- "C00079"
dictionary[32, 6] <- "C00148"
dictionary[33, 6] <- "C05942"
dictionary[34, 6] <- "C04281*"
dictionary[35, 6] <- "C00019"
dictionary[36, 6] <- "C04421"
dictionary[37, 6] <- "C00188"
dictionary[38, 6] <- "C00078"
dictionary[39, 6] <- "C00082"
dictionary[40, 6] <- "C00183"
dictionary[41, 6] <- "C00672"
dictionary[42, 6] <- "C06007"
dictionary[43, 6] <- "C00417"
dictionary[44, 6] <- "C00158"
dictionary[45, 6] <- "C06032"
dictionary[46, 6] <- "C00025"
dictionary[47, 6] <- "C00064"
dictionary[48, 6] <- "C00258"
# dictionary[49, 6] <-
dictionary[50, 6] <- "C04002*"
dictionary[51, 6] <- "C02504"
dictionary[52, 6] <- "C02631"
dictionary[53, 6] <- "C00026"
dictionary[54, 6] <- "C00186"
dictionary[55, 6] <- "C00149"
dictionary[56, 6] <- "C02226"
dictionary[57, 6] <- "C00022"
dictionary[58, 6] <- "C00042"
dictionary[59, 6] <- "C00232"
dictionary[60, 6] <- "C01050"
dictionary[61, 6] <- "C04390"
dictionary[62, 6] <- "C00041"
dictionary[63, 6] <- "C00555"
dictionary[64, 6] <- "C00441"
dictionary[65, 6] <- "C03979"
# dictionary[66, 6] <- "C00109"
dictionary[67, 6] <- "C01801"
dictionary[68, 6] <- "C03383"
dictionary[69, 6] <- "C00118"
# dictionary[70, 6] <-
dictionary[71, 6] <- "C02225"
dictionary[72, 6] <- "C00117"
dictionary[73, 6] <- "C05382"
dictionary[74, 6] <- "C01250"
dictionary[75, 6] <- "C04282"
dictionary[76, 6] <- "C12986"
dictionary[77, 6] <- "C00956"
dictionary[78, 6] <- "C00680"
# dictionary[79, 6] <-
dictionary[80, 6] <- "C00671"
dictionary[81, 6] <- "C05936"
dictionary[82, 6] <- "C05715"
dictionary[83, 6] <- "C03284"
dictionary[84, 6] <- "C04281"
dictionary[85, 6] <- "C06010"
# dictionary[86, 6] <-
dictionary[87, 6] <- "C00330"
dictionary[88, 6] <- "C04002"
dictionary[89, 6] <- "C03652"
dictionary[90, 6] <- "C02488"
dictionary[91, 6] <- "C00122"
dictionary[92, 6] <- "C00031"
dictionary[93, 6] <- "C00669"
dictionary[94, 6] <- "C00387"
dictionary[95, 6] <- "C00263"
dictionary[96, 6] <- "C05582"
dictionary[97, 6] <- "C16432"
dictionary[98, 6] <- "C00546"
dictionary[99, 6] <- "C00993"
dictionary[100, 6] <- "C00568"
dictionary[101, 6] <- "C00074"
dictionary[102, 6] <- "C00065"
dictionary[103, 6] <- "C01146"
dictionary[104, 6] <- "C00483"

dictionary <- as.data.frame(dictionary)
dictionary$`Old Data Names` <- as.character(dictionary$`Old Data Names`)
dictionary$`New Data Names` <- as.character(dictionary$`New Data Names`)
dictionary$Consensus <- as.character(dictionary$Consensus)
dictionary$Formula <- as.character(dictionary$Formula)
dictionary$`KEGG IDs` <- as.character(dictionary$`KEGG IDs`)

save(dictionary, file = "dictionary.RData")

