setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes as the ParsedSubGraphs_allStrains_named.RData object and the CCMN-normalized data and looks for    #################################
# correlations between the different variants found in each syntenic block and the differences in                      #################################
# metabolite abundance.                                                                                                #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

#######################################################################################################################################################
#
# Compare with metabolomics data
#
#######################################################################################################################################################

# Load ParsedSubGraphs_allStrains_named.RData object

load("ParsedSubGraphs_allStrains_named_new.RData")
ParsedSubGraphs_allStrains_named <- ParsedSubGraphs_allStrains_named_new


# Load CCMN normalized metabolomics data and dictionary
load("../normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("../normMetAnal/newData/ccmn_norm_mets_newData.RData")
load("../normMetAnal/hybrid/ccmn_norm_mets_hybrid.RData")
load("../dictionary/dictionary.RData")

# Change names of M6075 to M55212 (as in sequence data M6075 doesn't appear, and in metabolite data M55212 neither)

match(gsub("\\_.*|(PA14).*", 
           rownames(ccmn_norm_mets_good_old), 
           rep = "\\1"), 
      namesSeqs)



rownames(ccmn_norm_mets_good_old) <- gsub(pattern = "W70322", replacement = "W70332", rownames(ccmn_norm_mets_good_old))
rownames(ccmn_norm_mets_newData) <- gsub(pattern = "W70322", replacement = "W70332", rownames(ccmn_norm_mets_newData))
rownames(ccmn_norm_mets_hybrid) <- gsub(pattern = "W70322", replacement = "W70332", rownames(ccmn_norm_mets_hybrid))

metKEGGIDs_old <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
metKEGGIDs_new <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_newData), dictionary$`New Data Names`)]
metKEGGIDs_hybrid <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_hybrid), dictionary$Consensus)]

colnames(ccmn_norm_mets_good_old)[!is.na(metKEGGIDs_old)] <- metKEGGIDs_old[!is.na(metKEGGIDs_old)]
colnames(ccmn_norm_mets_newData)[!is.na(metKEGGIDs_new)] <- metKEGGIDs_new[!is.na(metKEGGIDs_new)]
colnames(ccmn_norm_mets_hybrid)[!is.na(metKEGGIDs_hybrid)] <- metKEGGIDs_hybrid[!is.na(metKEGGIDs_hybrid)]


# Obtain dist matrix and divide in 2 major clusters, and each major cluster in its 2 subclusters.
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/genePresAbs_functions.R")
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/variantAnalysis_functions.R")
clustsOld <- getMetClusts(ccmn_norm_mets_good_old)
clustsNew <- getMetClusts(ccmn_norm_mets_newData)
clustsHyb <- getMetClusts(ccmn_norm_mets_hybrid)

metDists1_2 <- dist(ccmn_norm_mets_good_old, "euclidean")
HCA_mets1_2 <- hclust(metDists1_2, method = "ward.D")
#HCA_mets1_2_cut <- cutree(HCA_mets1_2, k = 2)

#metDists2.1_2.2 <- dist(ccmn_norm_mets_good_old[which(HCA_mets1_2_cut == 2), ], "euclidean")
#HCA_mets2.1_2.2 <- hclust(metDists2.1_2.2, method = "ward.D")
#HCA_mets2.1_2.2_cut <- cutree(HCA_mets2.1_2.2, k = 2)

#strainsC1_2 <- HCA_mets1_2_cut
#names(strainsC1_2) <- gsub("_.*", names(HCA_mets1_2_cut), replacement = "")
#strainsC1_2 <- strainsC1_2[unique(names(strainsC1_2))]
#strainsC1_2 <- strainsC1_2[-c(1, 2)] #--> This line is for removing PA14A & PA14B

#strainsC1.1_1.2 <- HCA_mets1.1_1.2_cut
#names(strainsC1.1_1.2) <- gsub("_.*", names(HCA_mets1.1_1.2_cut), replacement = "")
#strainsC1.1_1.2 <- strainsC1.1_1.2[unique(names(strainsC1.1_1.2))]
#strainsC1.1_1.2 <- strainsC1.1_1.2[-c(1, 2)] #--> This line is for removing PA14A & PA14B


#strainsC2.1_2.2 <- HCA_mets2.1_2.2_cut
#names(strainsC2.1_2.2) <- gsub("_.*", names(HCA_mets2.1_2.2_cut), replacement = "")
#strainsC2.1_2.2 <- strainsC2.1_2.2[unique(names(strainsC2.1_2.2))]

# Obtain median of abundance of replicates

getStrainMedian <- function(normMets){
        normMetsMedians <- matrix(nrow = nrow(normMets)/3, ncol = ncol(normMets))
        for(j in 1:nrow(normMetsMedians)){
                for(i in 1:ncol(normMetsMedians)){
                        normMetsMedians[j, i] <- median(normMets[(1+3*(j-1)):(3*j), i])
                }
        }
        rownames(normMetsMedians) <- unique(gsub("\\_.*", "", rownames(normMets)))
        colnames(normMetsMedians) <-colnames(normMets)
        return(normMetsMedians)
}

goodOldMeds <- getStrainMedian(ccmn_norm_mets_good_old)
newMeds <- getStrainMedian(ccmn_norm_mets_newData)
hybMeds <- getStrainMedian(ccmn_norm_mets_hybrid)

metDists_median <- dist(ccmn_quant_norm_median, "euclidean")
HCA_mets_median <- hclust(metDists_median, method = "ward.D")
HCA_mets_1_2_median <- cutree(HCA_mets_median, k = 2)
strainsC1_2_median <- HCA_mets_1_2_median[-c(1, 2)]
# As M6075 is probable that is M55212 let's change the name
#names(strainsC1_2_median)[11] <- "M55212"

strainNames <- gsub(".fasta|.fna", 
                    replacement = "", 
                    list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive"))

strainNames <- gsub("PA14_jbx", "PA14", strainNames)

#######################################################################################################################################################
#                                                                                                                     #################################\
# Obtain then the vector of distances between pairs of strains belonging to the same cluster by one side, and         #################################
# to different cluster by the other side, and compare vectors gene a gene with a mann Whitney test, to then           #################################
# isolate those with alpha < 0.05. Do that using DNA distances and AA distances.                                      #################################
#                                                                                                                     #################################
#######################################################################################################################################################

# Here we are coding the cluster groups. As in the alignment are some strains that don't appear in our metabolic data we're adding them in the 
# vector of the group and we are ordering it according to the order of the sequences in the alignment to keep them in order inside of the function.

C1_2_old <- levels(clustsOld)[clustsOld]
names(C1_2_old) <- names(clustsOld)
C1_2_old[grep("1.", C1_2_old)] <- "1"
C1_2_old[grep("2.", C1_2_old)] <- "2"
C1_2_old <- as.factor(C1_2_old)
C1_2_old <- C1_2_old[match(strainNames, names(C1_2_old))]
names(C1_2_old) <- strainNames

C1_2_new <- levels(clustsNew)[clustsNew]
names(C1_2_new) <- names(clustsNew)
C1_2_new[grep("1.", C1_2_new)] <- "1"
C1_2_new[grep("2.", C1_2_new)] <- "2"
C1_2_new <- as.factor(C1_2_new)
C1_2_new <- C1_2_new[match(strainNames, names(C1_2_new))]
names(C1_2_new) <- strainNames


C1_2_hyb <- levels(clustsHyb)[clustsHyb]
names(C1_2_hyb) <- names(clustsHyb)
C1_2_hyb[grep("1.", C1_2_hyb)] <- "1"
C1_2_hyb[grep("2.", C1_2_hyb)] <- "2"
C1_2_hyb <- as.factor(C1_2_hyb)
C1_2_hyb <- C1_2_hyb[match(strainNames, names(C1_2_hyb))]
names(C1_2_hyb) <- strainNames

IDDiffGenes1_2_Old <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_old, p_adjust = T, alpha = 0.05, method = "BY")
IDDiffGenes1_2_New <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_new, p_adjust = T, alpha = 0.05, method = "BY")
IDDiffGenes1_2_Hyb <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_hyb, p_adjust = T, alpha = 0.05, method = "BY")
#IDDiffGenes1.1_1.2 <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, strainsC1.1_1.2, p_adjust = F)#, alpha = 0.05, method = "BY")
#IDDiffGenes2.1_2.2 <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, strainsC2.1_2.2, p_adjust = F)#, alpha = 0.05, method = "BY")

# Get the genes that are significatively different in terms of distance between the 2 clusters and that don't have a composed name

signFiltered1_2OldDNA <- IDDiffGenes1_2_Old$significativeGenes_DNA[sapply(as.character(IDDiffGenes1_2_Old$significativeGenes_DNA$Annotation), nchar) < 10, ]
signFiltered1_2OldAA <- IDDiffGenes1_2_Old$significativeGenes_AA[sapply(as.character(IDDiffGenes1_2_Old$significativeGenes_AA$Annotation), nchar) < 10, ]
# Load enzyme dictionary and see what enzymatic genes are in there. 
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")
signFiltered1_2OldDNAEnzymatic <- signFiltered1_2OldDNA[signFiltered1_2OldDNA$Annotation %in% dictEnzymes$Gene, ]
signFiltered1_2OldDNAEnzymatic <- cbind.data.frame(signFiltered1_2OldDNAEnzymatic, 
                                                   dictEnzymes$ECnums[match(signFiltered1_2OldDNAEnzymatic$Annotation, 
                                                                            dictEnzymes$Gene)])
colnames(signFiltered1_2OldDNAEnzymatic)[ncol(signFiltered1_2OldDNAEnzymatic)] <- "EC_numbers"

signFiltered1_2OldAAEnzymatic <- signFiltered1_2OldAA[signFiltered1_2OldAA$Annotation %in% dictEnzymes$Gene, ]
signFiltered1_2OldAAEnzymatic <- cbind.data.frame(signFiltered1_2OldAAEnzymatic, 
                                                   dictEnzymes$ECnums[match(signFiltered1_2OldAAEnzymatic$Annotation, 
                                                                            dictEnzymes$Gene)])
colnames(signFiltered1_2OldAAEnzymatic)[ncol(signFiltered1_2OldAAEnzymatic)] <- "EC_numbers"

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis/tab_OPLSDAQuant.RData")

signFiltered1_2OldDNAEnzymatic[which(signFiltered1_2OldDNAEnzymatic$EC_numbers %in% tab_OPLSDAQuant$KEGG.id), ] # EC 2.6.1.84 overlaps (aruH gene)
signFiltered1_2OldAAEnzymatic[which(signFiltered1_2OldAAEnzymatic$EC_numbers %in% tab_OPLSDAQuant$KEGG.id), ] # EC 2.7.1.221, EC 3.5.1.28 and EC 4.3.1.19.


aruH_dnaAlign <- ParsedSubGraphs_allStrains_named_new$aruH$DNAAlignment
names(aruH_dnaAlign) <- strainNames

BrowseSeqs(aruH_dnaAlign)

# Boxplot aruH
getGroupDists <- function(parsed, gene, typeSeq, strainGroups){
        if(typeSeq == "DNA"){
                mat1 <- parsed[[gene]]$DNADist[as.integer(gsub("_.*", 
                                                               "",
                                                               colnames(parsed[[gene]]$DNADist)))
                                               %in% which(strainGroups == 1),
                                               as.integer(gsub("_.*",
                                                               "",
                                                               colnames(parsed[[gene]]$DNADist)))
                                               %in% which(strainGroups == 1)]
                mat2 <- parsed[[gene]]$DNADist[as.integer(gsub("_.*", 
                                                               "",
                                                               colnames(parsed[[gene]]$DNADist)))
                                               %in% which(strainGroups == 2),
                                               as.integer(gsub("_.*",
                                                               "",
                                                               colnames(parsed[[gene]]$DNADist)))
                                               %in% which(strainGroups == 2)]
                matDiff <- parsed[[gene]]$DNADist[as.integer(gsub("_.*", 
                                                                  "",
                                                                  colnames(parsed[[gene]]$DNADist)))
                                                  %in% which(strainGroups == 1),
                                                  as.integer(gsub("_.*",
                                                                  "",
                                                                  colnames(parsed[[gene]]$DNADist)))
                                                  %in% which(strainGroups == 2)] 
        }
        if(typeSeq == "AA"){
                mat1 <- parsed[[gene]]$AADist[as.integer(gsub("_.*", 
                                                              "",
                                                              colnames(parsed[[gene]]$AADist)))
                                              %in% which(strainGroups == 1),
                                              as.integer(gsub("_.*",
                                                              "",
                                                              colnames(parsed[[gene]]$AADist)))
                                              %in% which(strainGroups == 1)]
                mat2 <- parsed[[gene]]$AADist[as.integer(gsub("_.*", 
                                                              "",
                                                              colnames(parsed[[gene]]$AADist)))
                                              %in% which(strainGroups == 2),
                                              as.integer(gsub("_.*",
                                                              "",
                                                              colnames(parsed[[gene]]$AADist)))
                                              %in% which(strainGroups == 2)]
                matDiff <- parsed[[gene]]$AADist[as.integer(gsub("_.*", 
                                                                 "",
                                                                 colnames(parsed[[gene]]$AADist)))
                                                 %in% which(strainGroups == 1),
                                                 as.integer(gsub("_.*",
                                                                 "",
                                                                 colnames(parsed[[gene]]$AADist)))
                                                 %in% which(strainGroups == 2)] 
        }
        sameGroupDists <- c(mat1[upper.tri(mat1)], mat2[upper.tri(mat2)])
        diffGroupDists <- as.vector(matDiff)
        group <- c(rep("withinGroup", length(sameGroupDists)),
                   rep("acrossGroup", length(diffGroupDists)))
        vals <- c(sameGroupDists, diffGroupDists)
        res <- cbind.data.frame(vals, group)
        colnames(res) <- c("distance", "measure")
        return(res)
}        

aruH <- getGroupDists(ParsedSubGraphs_allStrains_named_new, "aruH", "DNA", C1_2_old)
amgK <- getGroupDists(ParsedSubGraphs_allStrains_named_new, "amgK", "AA", C1_2_old)
amiC_3 <- getGroupDists(ParsedSubGraphs_allStrains_named_new, "amiC_3", "AA", C1_2_old)
ilvA_2 <- getGroupDists(ParsedSubGraphs_allStrains_named_new, "ilvA_2", "AA", C1_2_old)

aruH_boxplot <- ggplot(aruH, aes(x=measure, y=distance)) + 
        geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
aruH_boxplot
wilcox.test(aruH$distance[aruH$measure == "withinGroup"],
            aruH$distance[aruH$measure == "acrossGroup"])

amgK_boxplot <- ggplot(amgK, aes(x=measure, y=distance)) + 
        geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
amgK_boxplot
wilcox.test(amgK$distance[aruH$measure == "withinGroup"],
            amgK$distance[aruH$measure == "acrossGroup"])

amiC_3_boxplot <- ggplot(amiC_3, aes(x=measure, y=distance)) + 
        geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
amiC_3_boxplot
wilcox.test(amiC_3$distance[aruH$measure == "withinGroup"],
            amiC_3$distance[aruH$measure == "acrossGroup"])

ilvA_2_boxplot <- ggplot(ilvA_2, aes(x=measure, y=distance)) + 
        geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
ilvA_2_boxplot
wilcox.test(ilvA_2$distance[aruH$measure == "withinGroup"],
            ilvA_2$distance[aruH$measure == "acrossGroup"])

plotListDNA <- list()
for(i in 1:nrow(signFiltered1_2OldDNAEnzymatic)){
        dists <- getGroupDists(ParsedSubGraphs_allStrains_named_new, 
                               signFiltered1_2OldDNAEnzymatic$Annotation[i], 
                               "DNA", 
                               C1_2_old)
        p <- ggplot(dists, aes(x=measure, y=distance)) + 
                geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
        plotListDNA[[i]] <- p
}
names(plotListDNA) <- signFiltered1_2OldDNAEnzymatic$Annotation

plotListAA <- list()
for(i in 1:nrow(signFiltered1_2OldAAEnzymatic)){
        dists <- getGroupDists(ParsedSubGraphs_allStrains_named_new, 
                               signFiltered1_2OldAAEnzymatic$Annotation[i], 
                               "AA", 
                               C1_2_old)
        p <- ggplot(dists, aes(x=measure, y=distance)) + 
                geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2))
        plotListAA[[i]] <- p
}
names(plotListAA) <- signFiltered1_2OldAAEnzymatic$Annotation

#######################################################################################################################################################
#                                                                                                                     #################################
# Classify the strains according to the variants that exist for each one of the genes.                                #################################
#                                                                                                                     #################################
#######################################################################################################################################################

variants1_2_old <- variantsPerGene(ParsedSubGraphs_allStrains_named, IDDiffGenes1_2_Old, names(C1_2_old))
variants1_2_new <- variantsPerGene(ParsedSubGraphs_allStrains_named, IDDiffGenes1_2_New, names(C1_2_new))
variants1_2_hyb <- variantsPerGene(ParsedSubGraphs_allStrains_named, IDDiffGenes1_2_Hyb, names(C1_2_hyb))
#variants1.1_1.2 <- variantsPerGene(ParsedSubGraphs_allStrains_named, IDDiffGenes1.1_1.2, strainNames)
#variants2.1_2.2 <- variantsPerGene(ParsedSubGraphs_allStrains_named, IDDiffGenes2.1_2.2, strainNames)
save(variants1_2_old, file = "variants1_2_old.RData")
save(variants1_2_new, file = "variants1_2_new.RData")
save(variants1_2_hyb, file = "variants1_2_hyb.RData")

#######################################################################################################################################################
#                                                                                                                     #################################
# Divide the metabolite abundance matrix in submatrixes according to the common sequence variants between             #################################
# strains.                                                                                                            #################################
#                                                                                                                     #################################
#######################################################################################################################################################

divMetMats1_2_old <- divStrainMetsVars(variants1_2_old, ccmn_norm_mets_good_old)
divMetMats1_2_new <- divStrainMetsVars(variants1_2_new, ccmn_norm_mets_newData)
divMetMats1_2_hyb <- divStrainMetsVars(variants1_2_hyb, ccmn_norm_mets_hybrid)
#divMetMats1.1_1.2 <- divStrainMetsVars(variants1.1_1.2, ccmn_quant_norm)
#divMetMats2.1_2.2 <- divStrainMetsVars(variants2.1_2.2, ccmn_quant_norm)

#######################################################################################################################################################
#                                                                                                                     #################################
# Kruskal-Wallis test between the vectors corresponding to each metabolite for each group of strains that has         #################################
# a common variant, for each gene, with alpha = 0.05.                                                                 #################################
#                                                                                                                     #################################
#######################################################################################################################################################

kruskResults1_2_old <- kruskTestVars(divMetMats1_2_old)
kruskResults1_2_new <- kruskTestVars(divMetMats1_2_new)
kruskResults1_2_hyb <- kruskTestVars(divMetMats1_2_hyb)

kruskSign1_2_old <- kruskSignGenes(kruskResults1_2_old, alpha = 0.05)
kruskSign1_2_new <- kruskSignGenes(kruskResults1_2_new, alpha = 0.05)
kruskSign1_2_hyb <- kruskSignGenes(kruskResults1_2_hyb, alpha = 0.05)

#######################################################################################################################################################
#                                                                                                                     #################################
# Generate CSVs with all the significative genes for each one of the reportedly differential metabolites              #################################
# between clusters 1 and 2.                                                                                           #################################
#                                                                                                                     #################################
#######################################################################################################################################################


# Load the table corresponding to the differential metabolites between each pair of clusters.

load("../diffMetAnal/oldDataGood/diffMets_oldGood.RData")
load("../diffMetAnal/newData/diffMets_newData.RData")
load("../diffMetAnal/hybrid/diffMets_hybrid.RData")

csvOutSignDiff <- function(kruskSign, diffMets = NULL, metNames, spitOut = F){
        if(is.null(diffMets) == T){
                kruskSignFiltDNA <- kruskSign$DNA
                kruskSignFiltAA <- kruskSign$AA
        }else{
                kruskSignFiltDNA <- kruskSign$DNA[diffMets]
                kruskSignFiltAA <- kruskSign$AA[diffMets]
        }
        signDiffDNA <- list()
        signDiffDNA_enz <- list()
        signDiffAA <- list()
        signDiffAA_enz <- list()
        for(i in seq_along(kruskSignFiltDNA)){
                signDiffDNA[[i]] <- t(as.matrix(names(kruskSignFiltDNA[[i]])))
                signDiffDNA_enz[[i]] <- t(as.matrix(names(kruskSignFiltDNA[[i]])[grep("EC", names(kruskSignFiltDNA[[i]]))]))
                signDiffAA[[i]] <- t(as.matrix(names(kruskSignFiltAA[[i]])))
                signDiffAA_enz[[i]] <- t(as.matrix(names(kruskSignFiltAA[[i]])[grep("EC", names(kruskSignFiltAA[[i]]))]))
        }
        if(!require(plyr)) install.packages('plyr')
        library(plyr)
        MatDNA <- t(rbind.fill.matrix(signDiffDNA))
        MatDNA_enz <- t(rbind.fill.matrix(signDiffDNA_enz))
        MatAA <- t(rbind.fill.matrix(signDiffAA))
        MatAA_enz <- t(rbind.fill.matrix(signDiffAA_enz))
        colnames(MatDNA) <- diffMets
        colnames(MatDNA_enz) <- diffMets
        colnames(MatAA) <- diffMets
        colnames(MatAA_enz) <- diffMets
        MatDNA <- rbind(metNames[match(diffMets, names(kruskSign$DNA))], MatDNA)
        MatDNA_enz <- rbind(metNames[match(diffMets, names(kruskSign$DNA))], MatDNA_enz)
        MatAA <- rbind(metNames[match(diffMets, names(kruskSign$AA))], MatAA)
        MatAA_enz <- rbind(metNames[match(diffMets, names(kruskSign$AA))], MatAA_enz)
        rownames(MatDNA)[1] <- "gene_annot"
        rownames(MatDNA_enz)[1] <- "gene_annot"
        rownames(MatAA)[1] <- "gene_annot"
        rownames(MatAA_enz)[1] <- "gene_annot"
        if(spitOut == T){
                write.csv(MatDNA, file = paste(deparse(substitute(kruskSign)), "DNA.csv"))
                write.csv(MatDNA_enz, file = paste(deparse(substitute(kruskSign)), "DNA_enz.csv"))
                write.csv(MatAA, file = paste(deparse(substitute(kruskSign)), "AA.csv"))
                write.csv(MatAA_enz, file = paste(deparse(substitute(kruskSign)), "AA_enz.csv"))
        }
        return(list("DNA" = MatDNA, "DNA_enz" = MatDNA_enz, "AA" = MatAA, "AA_enz" = MatAA_enz))
}
outPut1_2_old <- csvOutSignDiff(kruskSign = kruskSign1_2_old, 
                                diffMets =  diffMets_oldGood[, 1], 
                                metNames = colnames(ccmn_norm_mets_good_old), 
                                spitOut = T)
save(outPut1_2_old, file = "outPut1_2_old.RData")
outPut1_2_new <- csvOutSignDiff(kruskSign = kruskSign1_2_new, 
                                diffMets =  diffMets_newData[, 1], 
                                metNames = colnames(ccmn_norm_mets_newData), 
                                spitOut = T)
save(outPut1_2_old, file = "outPut1_2_new.RData")
outPut1_2_hyb <- csvOutSignDiff(kruskSign = kruskSign1_2_hyb, 
                                diffMets =  diffMets_hybrid[, 1], 
                                metNames = colnames(ccmn_norm_mets_hybrid), 
                                spitOut = T)
save(outPut1_2_old, file = "outPut1_2_hyb.RData")

outPut1_2_oldfiltOPLSDAEnzDNA <- outPut1_2_old$DNA_enz[, colnames(outPut1_2_old$DNA_enz) %in% OPLSDAQuantResultKEGGIDs]
outPut1_2_oldfiltOPLSDAEnzAA <- outPut1_2_old$AA_enz[, colnames(outPut1_2_old$AA_enz) %in% OPLSDAQuantResultKEGGIDs]
write.csv(outPut1_2_oldfiltOPLSDAEnzDNA, file = "outPut1_2_oldfiltOPLSDAEnzDNA.csv")
write.csv(outPut1_2_oldfiltOPLSDAEnzAA, file = "outPut1_2_oldfiltOPLSDAEnzAA.csv")

t(as.matrix(names(kruskSign1_2_new$DNA[[29]])))
t(as.matrix(names(kruskSign1_2_new$DNA[29])[grep("EC", names(kruskSign1_2_new$DNA[[29]]))]))
signDiffAA[[i]] <- t(as.matrix(names(kruskSignFiltAA[[i]])))
signDiffAA_enz[[i]] <- t(as.matrix(names(kruskSignFiltAA[[i]])[grep("EC", names(kruskSignFiltAA[[i]]))]))


#######################################################################################################################################################
#                                                                                                                     #################################
# Build matrix with the variant each strain has for each gene, and paste it to gene presence/absence matrix,          #################################
# for then obtain Hamming Distance between strains and do a HCA with it, to see if with this we can reproduce         #################################
# the HCA we obtained from metabolite abundance data.                                                                 #################################
#                                                                                                                     #################################
#######################################################################################################################################################
if(!require(clustertend)) install.packages("clustertend")
library(clustertend)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(e1071)) install.packages("e1071")
library(e1071)
if(!require(ggdendro)) install.packages("ggdendro")
library(ggdendro)

genesMets_old <- c()
for(i in 1:ncol(outPut1_2_old$AA)){
        genesMets_old <- c(genesMets_old, outPut1_2_old$AA[2:nrow(outPut1_2_old$AA), i])
}
genesMets_old <- unique(genesMets_old[!is.na(genesMets_old)])

varMat_old <- matrix(nrow = 26, ncol = length(genesMets_old), dimnames = list(names(C1_2_old), genesMets_old))
for(i in seq_along(genesMets_old)){
        for(j in 1:length(variants1_2_old$AA[[genesMets_old[i]]])){
                varMat_old[names(variants1_2_old$AA[[genesMets_old[i]]])[j], i] <- variants1_2_old$AA[[i]][j]
        }
}
varMat_old[is.na(varMat_old)] <- 0
write.csv(varMat_old, file = "varMat_old.csv")
# Load gene_tab_filt  and gene_enz_tab_filt

load("../decTreesLogReg/mixedSign_old.RData")
mixedSign_old <- t(mixedSign_old)
mixedSign_old <- mixedSign_old[order(rownames(mixedSign_old)),]

load("../Diff_metabolites_analysis/gene_tab_filt.RData")
load("../Diff_metabolites_analysis/gene_enz_tab_filt.RData")

varMat_geneTab <- cbind(as.matrix(gene_tab_filt), varMat)

varMat_mixedSign_old <- cbind(mixedSign_old, varMat_old)

# Change the name of "M6075" by "M55212", as we're assuming they are the same.
rownames(varMat_geneTab)[which(rownames(varMat_geneTab) == "M6075")] <- "M55212"
save(varMat_geneTab, file = "varMat_geneTab.RData")

rownames(varMat_mixedSign_old)[which(rownames(varMat_mixedSign_old) == "M6075")] <- "M55212"
save(varMat_mixedSign_old, file = "varMat_mixedSign_old.RData")

load("varMat_geneTab.RData")
write.csv(varMat_geneTab, file = "varMat_geneTab.csv")

hopkins(varMat_geneTab, n = 24, header = T)
get_clust_tendency(varMat_geneTab, n = 24)


hamDist <- hamming.distance(varMat_geneTab)
save(hamDist, file = "hamDist.RData")
HCA_hamDist <- hclust(as.dist(hamDist), method = "complete")

hamDist_mixOld <- hamming.distance(varMat_mixedSign_old)
save(hamDist_mixOld, file = "hamDist_mixOld.RData")
HCA_hamDist_mixOld <- hclust(as.dist(hamDist_mixOld), method = "complete")

# Compare with the dendrogram obtained when clustering the metabolomic data
if(!require(dendextend)) install.packages("dendextend")
library(dendextend)
ccmn_quant_norm_median_woPA14 <- ccmn_quant_norm_median[-c(1, 2), ]

metDists_median_woPA14 <- dist(ccmn_quant_norm_median_woPA14, "euclidean")
HCA_mets_median_woPA14 <- hclust(metDists_median_woPA14, method = "ward.D")

metDists_goodOldMeds <- dist(goodOldMeds, "euclidean")
HCA_goodOldMeds <- hclust(metDists_goodOldMeds, method = "ward.D")

goodOldMeds_noPA14 <- goodOldMeds[-c(1, 2), ]
metDists_goodOldMeds_noPA14 <- dist(goodOldMeds_noPA14, "euclidean")
HCA_goodOldMeds_noPA14 <- hclust(metDists_goodOldMeds_noPA14, method = "ward.D")

plot(HCA_hamDist)
plot(HCA_mets_median_woPA14)
plot(HCA_hamDist_mixOld)
plot(HCA_goodOldMeds_noPA14)

tissues <- c("wild type", 
             "wild type", 
             "pubic bone", 
             "blood", 
             "urine", 
             "tissue", 
             "body fluid", 
             "stool", 
             "urine", 
             "tissue", 
             "blood", 
             "body fluid", 
             "blood",
             "urine",
             "sputum",
             "urine",
             "blood",
             "blood",
             "blood",
             "blood",
             "sputum",
             "urine",
             "blood",
             "body fluid",
             "blood",
             "sputum",
             "blood",
             "blood")

names(tissues) <- rownames(goodOldMeds)
tissues <- as.factor(tissues)


library(dendextend)

small_iris <- iris[c(1, 51, 101, 2, 52, 102), ]
dend <- as.dendrogram(hclust(dist(small_iris[,-5])))
# Like: 
# dend <- small_iris[,-5] %>% dist %>% hclust %>% as.dendrogram

# By default, the dend has no colors to the labels
labels_colors(dend)
par(mfrow = c(1,2))
plot(dend, main = "Original dend")

# let's add some color:
colors_to_use <- as.numeric(small_iris[,5])
colors_to_use
# But sort them based on their order in dend:
colors_to_use <- colors_to_use[order.dendrogram(dend)]
colors_to_use
# Now we can use them
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 
plot(dend, main = "A color for every Species")





df <- goodOldMeds


#if(!require(magittr)) BiocManager::install('magittr')
require(ggplot)
require(dendextend)

dend <- df %>% dist(method = "euclidean") %>%
        hclust(method = "ward.D") %>% as.dendrogram %>%
        set("branches_k_color", k = 2) %>% set("branches_lwd", 0.7) %>%
        set("labels_cex", 0.8) %>% set("leaves_pch", 19) %>% 
        set("leaves_cex", 0.5) 

unique(tissues)
tissuesColors <- c("firebrick",
                   "blue3",
                   "burly wood", 
                   "cadetblue3",
                   "chocolate4",
                   "darkgreen",
                   "darkgoldenrod1",
                   "black")
colorsToUse <- as.numeric(tissues)
colorsToUse <- tissuesColors[colorsToUse]
colorsToUse <- colorsToUse[order.dendrogram(dend)]
labels_colors(dend) <- colorsToUse

ggd1 <- as.ggdend(dend)
ggplot(ggd1)
ggsave(filename = "tst.pdf", plot = last_plot(), device = "pdf", 
       path = NULL,
       scale = 1, width = 40, height = 35, units = "cm",
       dpi = 300, limitsize = F)







ggdendrogram(data = HCA_goodOldMeds, main = "HCA medians Mets  (Old Data)")
ggsave(filename = "HCA_CCMN_goodOld.pdf", plot = last_plot(), device = "pdf", 
       path = NULL,
       scale = 1, width = 40, height = 20, units = "cm",
       dpi = 300, limitsize = F)
cophMat <- cophenetic(HCA_CCMN)

ggdendrogram(data = HCA_hamDist_mixOld, main = "HCA HamDist Mic (Old Data)")
ggdendrogram(data = HCA_goodOldMeds_noPA14, main = "HCA medians Mets without PA14 (Old Data)")
ggsave(filename = "HCA_CCMN_goodOld_noPA14.pdf", plot = last_plot(), device = "pdf", 
       path = NULL,
       scale = 1, width = 40, height = 20, units = "cm",
       dpi = 300, limitsize = F)
cophMat <- cophenetic(HCA_CCMN)
cor(distMatCCMN, cophMat)

dend_hamDist <- as.dendrogram(HCA_hamDist)
dend_mets_median_woPA14 <- as.dendrogram(HCA_mets_median_woPA14)
dend_hamDist_mixOld <- as.dendrogram(HCA_hamDist_mixOld)
dend_goodOldMeds_noPA14 <- as.dendrogram(HCA_goodOldMeds_noPA14)


dend_compList <- dendlist(dend_hamDist, dend_mets_median_woPA14) %>% untangle(method = "step1side") %>%
        tanglegram()  

dendlist(dend_hamDist, dend_mets_median_woPA14) %>% untangle(method = "step1side") %>%
        entanglement() #--> 0.17 --> Close to cero, so it means that are related.

dend_compList_old <- dendlist(dend_hamDist_mixOld, dend_goodOldMeds_noPA14) %>% untangle(method = "step1side") %>%
        tanglegram()  

dendlist(dend_hamDist_mixOld, dend_goodOldMeds_noPA14) %>% untangle(method = "step1side") %>%
        entanglement() 

cor.dendlist(dend_compList, method = "cophenetic")
cor.dendlist(dend_compList, method = "baker")

cor.dendlist(dend_compList_old, method = "cophenetic")
cor.dendlist(dend_compList_old, method = "baker")

cor_cophenetic(dend_compList_old) # Some correlation, but not very high.

