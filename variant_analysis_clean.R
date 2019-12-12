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

load("ParsedSubGraphs_allStrains_named.RData")

# Load CCMN normalized metabolomics data and dictionary
load("../normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("../normMetAnal/newData/ccmn_norm_mets_newData.RData")
load("../normMetAnal/hybrid/ccmn_norm_mets_hybrid.RData")
load("../dictionary/dictionary.RData")

# Change names of M6075 to M55212 (as in sequence data M6075 doesn't appear, and in metabolite data M55212 neither)
rownames(ccmn_norm_mets_good_old) <- gsub(pattern = "M6075", replacement = "M55212", rownames(ccmn_norm_mets_good_old))
rownames(ccmn_norm_mets_newData) <- gsub(pattern = "M6075", replacement = "M55212", rownames(ccmn_norm_mets_newData))
rownames(ccmn_norm_mets_hybrid) <- gsub(pattern = "M6075", replacement = "M55212", rownames(ccmn_norm_mets_hybrid))

metKEGGIDs_old <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
metKEGGIDs_new <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_newData), dictionary$`New Data Names`)]
metKEGGIDs_hybrid <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_hybrid), dictionary$Consensus)]

colnames(ccmn_norm_mets_good_old)[!is.na(metKEGGIDs_old)] <- metKEGGIDs_old[!is.na(metKEGGIDs_old)]
colnames(ccmn_norm_mets_newData)[!is.na(metKEGGIDs_new)] <- metKEGGIDs_new[!is.na(metKEGGIDs_new)]
colnames(ccmn_norm_mets_hybrid)[!is.na(metKEGGIDs_hybrid)] <- metKEGGIDs_hybrid[!is.na(metKEGGIDs_hybrid)]


# Obtain dist matrix and divide in 2 major clusters, and each major cluster in its 2 subclusters.
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/genePresAbs_functions.R")
clustsOld <- getMetClusts(ccmn_norm_mets_good_old)
clustsNew <- getMetClusts(ccmn_norm_mets_newData)
clustsHyb <- getMetClusts(ccmn_norm_mets_hybrid)

metDists1_2 <- dist(ccmn_norm_mets_good_old, "euclidean")
HCA_mets1_2 <- hclust(metDists1_2, method = "ward.D")
HCA_mets1_2_cut <- cutree(HCA_mets1_2, k = 2)

metDists2.1_2.2 <- dist(ccmn_quant_norm[which(HCA_mets_1_2 == 2), ], "euclidean")
HCA_mets2.1_2.2 <- hclust(metDists2.1_2.2, method = "ward.D")
HCA_mets2.1_2.2_cut <- cutree(HCA_mets2.1_2.2, k = 2)

strainsC1_2 <- HCA_mets_1_2
names(strainsC1_2) <- gsub("_.*", names(HCA_mets_1_2), replacement = "")
strainsC1_2 <- strainsC1_2[unique(names(strainsC1_2))]
strainsC1_2 <- strainsC1_2[-c(1, 2)] #--> This line is for removing PA14A & PA14B

strainsC1.1_1.2 <- HCA_mets1.1_1.2_cut
names(strainsC1.1_1.2) <- gsub("_.*", names(HCA_mets1.1_1.2_cut), replacement = "")
strainsC1.1_1.2 <- strainsC1.1_1.2[unique(names(strainsC1.1_1.2))]
strainsC1.1_1.2 <- strainsC1.1_1.2[-c(1, 2)] #--> This line is for removing PA14A & PA14B


strainsC2.1_2.2 <- HCA_mets2.1_2.2_cut
names(strainsC2.1_2.2) <- gsub("_.*", names(HCA_mets2.1_2.2_cut), replacement = "")
strainsC2.1_2.2 <- strainsC2.1_2.2[unique(names(strainsC2.1_2.2))]

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
names(strainsC1_2_median)[11] <- "M55212"

strainNames <- gsub(".fasta", "", list.files(path = "../../Sequences/FASTA_full_genomes/my_strains/Seqs")[1:26])

#######################################################################################################################################################
#                                                                                                                     #################################
# Obtain by one side median distances for each gene between all pairs of strains belonging to same cluster and        #################################
# belonging to the other cluster (Clusters 1 and 2) and compare those two vectors with a Mann-Whitney test.           #################################
# Obtain then the vector of distances between pairs of strains belonging to the same cluster by one side, and         #################################
# to different cluster by the other side, and compare vectors gene a gene with a mann Whitney test, to then           #################################
# isolate those with alpha < 0.05. Do that using DNA distances and AA distances.                                      #################################
#                                                                                                                     #################################
#######################################################################################################################################################
C1_2_old <- levels(clustsOld)[clustsOld]
names(C1_2_old) <- names(clustsOld)
C1_2_old[grep("1.", C1_2_old)] <- "1"
C1_2_old[grep("2.", C1_2_old)] <- "2"
C1_2_old <- C1_2_old[-1]
C1_2_old <- as.factor(C1_2_old)

C1_2_new <- levels(clustsNew)[clustsNew]
names(C1_2_new) <- names(clustsNew)
C1_2_new[grep("1.", C1_2_new)] <- "1"
C1_2_new[grep("2.", C1_2_new)] <- "2"
C1_2_new <- C1_2_new[-1]
C1_2_new <- as.factor(C1_2_new)

C1_2_hyb <- levels(clustsHyb)[clustsHyb]
names(C1_2_hyb) <- names(clustsHyb)
C1_2_hyb[grep("1.", C1_2_hyb)] <- "1"
C1_2_hyb[grep("2.", C1_2_hyb)] <- "2"
C1_2_new <- C1_2_new[-1]
C1_2_hyb <- as.factor(C1_2_hyb)

identifyDiffGenes <- function(ParsedSubGraphsNamed_objt, strainGroups, p_adjust = T, method = NA, alpha = 0.05){
        medDistSameDNACl1 <- c()
        medDistSameDNACl2 <- c()
        medDistSameDNA <- c()
        medDistDiffDNA <- c()
        distsDiffCDNA <- list()
        medDistSameAACl1 <- c()
        medDistSameAACl2 <- c()
        medDistSameAA <- c()
        medDistDiffAA <- c()
        distsDiffCAA <- list()
        mannWhitResDNA <- c()
        mannWhitResAA <- c()
        namesGenes <- names(ParsedSubGraphsNamed_objt)
        pb = txtProgressBar(min = 0, max = length(ParsedSubGraphsNamed_objt), initial = 0)
        for(i in seq_along(ParsedSubGraphs_allStrains_named)){
                subMat1DNA <- ParsedSubGraphsNamed_objt[[i]]$DNADist[as.integer(gsub("_.*", 
                                                                                     replacement = "", 
                                                                                     colnames(ParsedSubGraphs_allStrains_named[[i]]$DNADist)))
                                                                     %in% which(strainGroups == 1), 
                                                                     as.integer(gsub("_.*", 
                                                                                     replacement = "", 
                                                                                     colnames(ParsedSubGraphsNamed_objt[[i]]$DNADist))) 
                                                                     %in% which(strainGroups == 1)]
                medDistSameDNACl1 <- c(medDistSameDNACl1, median(subMat1DNA[upper.tri(subMat1DNA)]))
                subMat2DNA <- ParsedSubGraphsNamed_objt[[i]]$DNADist[as.integer(gsub("_.*", 
                                                                                     replacement = "", 
                                                                                     colnames(ParsedSubGraphsNamed_objt[[i]]$DNADist)))
                                                                     %in% which(strainGroups == 2), 
                                                                     as.integer(gsub("_.*", 
                                                                                     replacement = "", 
                                                                                     colnames(ParsedSubGraphsNamed_objt[[i]]$DNADist))) 
                                                                     %in% which(strainGroups == 2)]
                medDistSameDNACl2 <- c(medDistSameDNACl2, median(subMat2DNA[upper.tri(subMat2DNA)]))
                subMatDiffDNA <- ParsedSubGraphsNamed_objt[[i]]$DNADist[as.integer(gsub("_.*", 
                                                                                        replacement = "", 
                                                                                        colnames(ParsedSubGraphsNamed_objt[[i]]$DNADist)))
                                                                        %in% which(strainGroups == 1), 
                                                                        as.integer(gsub("_.*", 
                                                                                        replacement = "", 
                                                                                        colnames(ParsedSubGraphsNamed_objt[[i]]$DNADist))) 
                                                                        %in% which(strainGroups == 2)]
                medDistDiffDNA <- c(medDistDiffDNA, median(as.vector(subMatDiffDNA)))
                distsDiffCDNA[[i]] <- as.vector(subMatDiffDNA)
                if(length(distsDiffCDNA[[i]]) > 6){
                        mannWhitResDNA <- c(mannWhitResDNA, wilcox.test(c(subMat1DNA[upper.tri(subMat1DNA)], 
                                                                          subMat2DNA[upper.tri(subMat2DNA)]), 
                                                                        as.vector(subMatDiffDNA), alternative = "two.sided")[[3]])
                }else{
                        mannWhitResDNA <- c(mannWhitResDNA, NA)
                }
                medDistSameDNA <- c(medDistSameDNA, median(c(subMat1DNA[upper.tri(subMat1DNA)],
                                                             subMat2DNA[upper.tri(subMat2DNA)])))
                subMat1AA <- ParsedSubGraphsNamed_objt[[i]]$AADist[as.integer(gsub("_.*", 
                                                                                   replacement = "", 
                                                                                   colnames(ParsedSubGraphsNamed_objt[[i]]$AADist)))
                                                                   %in% which(strainGroups == 1), 
                                                                   as.integer(gsub("_.*", 
                                                                                   replacement = "", 
                                                                                   colnames(ParsedSubGraphsNamed_objt[[i]]$AADist))) 
                                                                   %in% which(strainGroups == 1)]
                medDistSameAACl1 <- c(medDistSameAACl1, median(subMat1AA[upper.tri(subMat1AA)]))
                subMat2AA <- ParsedSubGraphsNamed_objt[[i]]$AADist[as.integer(gsub("_.*", 
                                                                                   replacement = "", 
                                                                                   colnames(ParsedSubGraphsNamed_objt[[i]]$AADist)))
                                                                   %in% which(strainGroups == 2), 
                                                                   as.integer(gsub("_.*", 
                                                                                   replacement = "", 
                                                                                   colnames(ParsedSubGraphsNamed_objt[[i]]$AADist))) 
                                                                   %in% which(strainGroups == 2)]
                medDistSameAACl2 <- c(medDistSameAACl2, median(subMat2AA[upper.tri(subMat2AA)]))
                subMatDiffAA <- ParsedSubGraphsNamed_objt[[i]]$AADist[as.integer(gsub("_.*", 
                                                                                      replacement = "", 
                                                                                      colnames(ParsedSubGraphsNamed_objt[[i]]$AADist)))
                                                                      %in% which(strainGroups == 1), 
                                                                      as.integer(gsub("_.*", 
                                                                                      replacement = "", 
                                                                                      colnames(ParsedSubGraphsNamed_objt[[i]]$AADist))) 
                                                                      %in% which(strainGroups == 2)]
                medDistDiffAA <- c(medDistDiffAA, median(as.vector(subMatDiffAA)))
                distsDiffCAA[[i]] <- as.vector(subMatDiffAA)
                if(length(distsDiffCAA[[i]]) > 6){
                        mannWhitResAA <- c(mannWhitResAA, wilcox.test(c(subMat1AA[upper.tri(subMat1AA)], 
                                                                        subMat2AA[upper.tri(subMat2AA)]), 
                                                                      as.vector(subMatDiffAA), alternative = "two.sided")[[3]])
                }else{
                        mannWhitResAA <- c(mannWhitResAA, NA)
                }
                medDistSameAA <- c(medDistSameAA, median(c(subMat1AA[upper.tri(subMat1AA)],
                                                           subMat2AA[upper.tri(subMat2AA)])))
                setTxtProgressBar(pb, i)
        }
        mannWhitMeds <- list("DiffDNA" = wilcox.test(medDistSameDNA, medDistDiffDNA, alternative = "two.sided"),
                             "SameDNA" = wilcox.test(medDistSameDNACl1, medDistSameDNACl2, alternative = "two.sided"),
                             "DiffAA" = wilcox.test(medDistSameAA, medDistDiffAA, alternative = "two.sided"),
                             "SameAA" = wilcox.test(medDistSameAACl1, medDistSameAACl2, alternative = "two.sided"))
        distDifferDNA <- medDistDiffDNA - medDistSameDNA
        distDifferAA <- medDistDiffAA - medDistSameAA
        names(mannWhitResDNA) <- names(ParsedSubGraphsNamed_objt)
        names(mannWhitResAA) <- names(ParsedSubGraphsNamed_objt)
        if(p_adjust == T){
                mannWhitResDNA <- p.adjust(mannWhitResDNA, method = method)
                mannWhitResAA <- p.adjust(mannWhitResAA, method = method) 
        }
        signGenesClusts1_2DNA <- cbind(namesGenes[which(mannWhitResDNA <= alpha)], mannWhitResDNA[which(mannWhitResDNA <= alpha)])
        signGenesClusts1_2DNA <- cbind(signGenesClusts1_2DNA, distDifferDNA[which(mannWhitResDNA <= alpha)])
        colnames(signGenesClusts1_2DNA) <- c("Annotation", "p-value", "medians_difference")
        signGenesClusts1_2DNA <- signGenesClusts1_2DNA[order(signGenesClusts1_2DNA[, 3], decreasing = T), ]
        signGenesClusts1_2AA <- cbind(namesGenes[which(mannWhitResAA <= alpha)], mannWhitResAA[which(mannWhitResAA <= alpha)])
        signGenesClusts1_2AA <- cbind(signGenesClusts1_2AA, distDifferAA[which(mannWhitResAA <= alpha)])
        colnames(signGenesClusts1_2AA) <- c("Annotation", "p-value", "medians_difference")
        signGenesClusts1_2AA <- signGenesClusts1_2AA[order(signGenesClusts1_2AA[, 3], decreasing = T), ]
        results <- list("MannWhitney_meds" = mannWhitMeds, 
                        "significativeGenes_DNA" = signGenesClusts1_2DNA,
                        "significativeGenes_AA" = signGenesClusts1_2AA)
        return(results)
}

IDDiffGenes1_2_Old <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_old, p_adjust = T, alpha = 0.05, method = "BY")
IDDiffGenes1_2_New <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_new, p_adjust = T, alpha = 0.05, method = "BY")
IDDiffGenes1_2_Hyb <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, C1_2_hyb, p_adjust = T, alpha = 0.05, method = "BY")
#IDDiffGenes1.1_1.2 <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, strainsC1.1_1.2, p_adjust = F)#, alpha = 0.05, method = "BY")
#IDDiffGenes2.1_2.2 <- identifyDiffGenes(ParsedSubGraphs_allStrains_named, strainsC2.1_2.2, p_adjust = F)#, alpha = 0.05, method = "BY")

#######################################################################################################################################################
#                                                                                                                     #################################
# Classify the strains according to the variants that exist for each one of the genes.                                #################################
#                                                                                                                     #################################
#######################################################################################################################################################

variantsPerGene <- function(ParsedSubGraphsNamed_objct, IDDiffGenes, strainNames){
        variantsDNA <- list()
        pbDNA = txtProgressBar(min = 0, max = nrow(IDDiffGenes$significativeGenes_DNA), initial = 0)
        for(i in 1:nrow(IDDiffGenes$significativeGenes_DNA)){
                DNAAlign <- ParsedSubGraphsNamed_objct[[IDDiffGenes$significativeGenes_DNA[i, 1]]]$DNAAlignment
                strainVarVecDNA <- rep(NA, length(strainNames))
                names(strainVarVecDNA) <- strainNames
                for(j in 1:length(unique(DNAAlign))){
                        for(k in 1:length(DNAAlign)){
                                if(DNAAlign[k] == unique(DNAAlign)[j]){
                                        strainVarVecDNA[as.integer(gsub("_.*",
                                                                        replacement = "",
                                                                        names(DNAAlign)))[k]] <- j
                                }
                        }
                }
                strainVarVecDNA <- strainVarVecDNA[!is.na(strainVarVecDNA)]
                variantsDNA[[i]] <- strainVarVecDNA
                setTxtProgressBar(pbDNA, i)
        }
        names(variantsDNA) <- IDDiffGenes$significativeGenes_DNA[, 1]
        variantsAA <- list()
        pbAA = txtProgressBar(min = 0, max = nrow(IDDiffGenes$significativeGenes_AA), initial = 0)
        for(l in 1:nrow(IDDiffGenes$significativeGenes_AA)){
                AAAlign <- ParsedSubGraphsNamed_objct[[IDDiffGenes$significativeGenes_AA[l, 1]]]$AAAlignment
                strainVarVecAA <- rep(NA, length(strainNames))
                names(strainVarVecAA) <- strainNames
                for(m in 1:length(unique(AAAlign))){
                        for(n in 1:length(AAAlign)){
                                if(AAAlign[n] == unique(AAAlign)[m]){
                                        strainVarVecAA[as.integer(gsub("_.*",
                                                                       replacement = "",
                                                                       names(AAAlign)))[n]] <- m
                                }
                        }
                }
                strainVarVecAA <- strainVarVecAA[!is.na(strainVarVecAA)]
                variantsAA[[l]] <- strainVarVecAA
                setTxtProgressBar(pbAA, l)
        }
        names(variantsAA) <- IDDiffGenes$significativeGenes_AA[, 1]
        variants <- list("DNA" = variantsDNA, "AA" = variantsAA)
        return(variants)
}

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

divStrainMetsVars <- function(variants, metsAbundance){
        divMetMatsDNA <- list()
        for(i in seq_along(variants$DNA)){
                subMatsDNA <- list()
                for(j in 1:max(variants$DNA[[i]])){
                        subMatsDNA[[j]] <- metsAbundance[grepl(pattern = paste(names(which(variants$DNA[[i]] == j)), 
                                                                                     collapse = "|"),
                                                               rownames(metsAbundance)), ]
                }
                divMetMatsDNA[[i]] <- subMatsDNA
        }
        names(divMetMatsDNA) <- names(variants$DNA)
        divMetMatsAA <- list()
        for(k in seq_along(variants$AA)){
                subMatsAA <- list()
                for(l in 1:max(variants$AA[[k]])){
                        subMatsAA[[l]] <- metsAbundance[grepl(pattern = paste(names(which(variants$AA[[k]] == l)), 
                                                                              collapse = "|"),
                                                              rownames(metsAbundance)), ]
                }
                divMetMatsAA[[k]] <- subMatsAA
        }
        names(divMetMatsAA) <- names(variants$AA)
        divMetsMats <- list("DNA" = divMetMatsDNA, "AA" = divMetMatsAA)
        return(divMetsMats)
}

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

kruskTestVars <- function(divMet_objct){
        resultDNA <- list()
        resultAA <- list()
        for(i in 1:ncol(divMet_objct$DNA[[1]][[1]])){
                kruskResDNA <- c()
                for(j in seq_along(divMet_objct$DNA)){
                        kruskVecsDNA <- list()
                        for(k in seq_along(divMet_objct$DNA[[j]])){
                                kruskVecsDNA[[k]] <- divMet_objct$DNA[[j]][[k]][, i]
                        }
                        kruskResDNA <- c(kruskResDNA, kruskal.test(kruskVecsDNA)[[3]])
                }
                names(kruskResDNA) <- names(divMet_objct$DNA)
                kruskResDNA <- p.adjust(kruskResDNA, method = "BY")
                resultDNA[[i]] <- kruskResDNA
                
                kruskResAA <- c()
                for(l in seq_along(divMet_objct$AA)){
                        kruskVecsAA <- list()
                        for(m in seq_along(divMet_objct$AA[[l]])){
                                kruskVecsAA[[m]] <- divMet_objct$AA[[l]][[m]][, i]
                        }
                        kruskResAA <- c(kruskResAA, kruskal.test(kruskVecsAA)[[3]])
                }
                names(kruskResAA) <- names(divMet_objct$AA)
                kruskResAA <- p.adjust(kruskResAA, method = "BY")
                resultAA[[i]] <- kruskResAA
        }
        names(resultDNA) <- colnames(divMet_objct$DNA[[1]][[1]])
        names(resultAA) <- colnames(divMet_objct$AA[[1]][[1]])
        result <- list("DNA" = resultDNA, "AA" = resultAA)
        return(result)
}

kruskResults1_2_old <- kruskTestVars(divMetMats1_2_old)
kruskResults1_2_new <- kruskTestVars(divMetMats1_2_new)
kruskResults1_2_hyb <- kruskTestVars(divMetMats1_2_hyb)

kruskSignGenes <- function(kruskRes, alpha = 0.05){
        kruskSignDNA <- list()
        kruskSignAA <- list()
        for(i in seq_along(kruskRes$DNA)){
                kruskSignDNA[[i]] <- kruskRes$DNA[[i]][which(kruskRes$DNA[[i]] <= alpha)]
                kruskSignAA[[i]] <- kruskRes$AA[[i]][which(kruskRes$AA[[i]] <= alpha)]
        }
        names(kruskSignDNA) <- names(kruskRes$DNA)
        names(kruskSignAA) <- names(kruskRes$AA)
        kruskSign <- list("DNA" = kruskSignDNA, "AA" = kruskSignAA)
        return(kruskSign)
}

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
                print(i)
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
outPut1_2_new <- csvOutSignDiff(kruskSign = kruskSign1_2_new, 
                                diffMets =  diffMets_newData[, 1], 
                                metNames = colnames(ccmn_norm_mets_newData), 
                                spitOut = T)
outPut1_2_hyb <- csvOutSignDiff(kruskSign = kruskSign1_2_hyb, 
                                diffMets =  diffMets_hybrid[, 1], 
                                metNames = colnames(ccmn_norm_mets_hybrid), 
                                spitOut = T)




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


if(!require(magittr)) BiocManager::install('magittr')
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

