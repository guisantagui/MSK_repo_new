setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/oldDataGood")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes as input CCMN normalized data and selects the differential metabolites between the 4 major groups  #################################
# of strains found according to metabolite abundance.                                                                  #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

#Load data
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
ccmn_norm_mets <- ccmn_norm_mets_good_old

########################################################################################################################################################
#
# Apply Hierarchical Clustering Analysis
#
########################################################################################################################################################

if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(ggdendro)) install.packages("ggdendro")
library(ggdendro)
if(!require(clustertend)) install.packages("clustertend")
library(clustertend)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)

#Evaluate cluster tendency of the data
hopkins(ccmn_norm_mets, n = 24, header = T) #--> ~0.3 -->Regularly spaced data?!?

get_clust_tendency(ccmn_norm_mets, n = 24)

distMatCCMN <- dist(ccmn_norm_mets, method = "euclidean")

HCA_CCMN <- hclust(distMatCCMN, method = "ward.D")
plot(HCA_CCMN, main = "HCA of all the strains", sub = "Ward.D distance used", xlab = NULL)

ggdendrogram(data = HCA_CCMN, main = "HCA of all the strains")
ggsave(filename = "HCA_CCMN_goodOld.pdf", plot = last_plot(), device = "pdf", 
       path = NULL,
       scale = 1, width = 40, height = 20, units = "cm",
       dpi = 300, limitsize = F)
cophMat <- cophenetic(HCA_CCMN)
cor(distMatCCMN, cophMat) #--> 0.5104988 --> Low correlation 


#With HeatMaps function from metabolomics package we can see the clustering of the strains over
#a heatmap of the abundance of metabolites
if(!require(remotes)) install.packages("remotes")
library(remotes)
if(!require(pcaMethods)) BiocManager::install("pcaMethods")
library(pcaMethods)
if(!require(metabolomics)) install_version("metabolomics", "0.1.4")
library(metabolomics)
if(!require(IRanges)) install.packages("IRanges")
library(IRanges)


strain_names <- unique(gsub("\\_.*", "", rownames(ccmn_norm_mets)))
cols <- topo.colors(length(strain_names))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(ccmn_norm_mets))))
for(ii in 1:length(ccmn_norm_mets[, 1])) {
        selected <- which(strain_names == gsub("\\_.*", "", rownames(ccmn_norm_mets))[ii])
        Cols[ii] <- cols[selected]
}

HeatMap(ccmn_norm_mets, 
        scale = 'column',
        ColSideColors = Cols, 
        colramp = bluered(75),
        distmethod = "euclidean",
        aggmethod = "ward.D",
        cexRow = 0.65)

########################################################################################################################################################
#
# Quantile-Normalize CCMN normalized data in order to make differences visualizable (all replicates are clustering together, so it
# must be differences)
#
########################################################################################################################################################
if(!require(tibble)) install.packages("tibble")
library(tibble)

if(!require(tibble)) install.packages("tibble")
library(tibble)

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")



compKEGGIDs <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets), dictionary$`Old Data Names`)]

newMetNames <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
colnames(ccmn_norm_mets) <- newMetNames
ccmn_norm_mets <- ccmn_norm_mets[, !is.na(colnames(ccmn_norm_mets))]

quantNorm <- function(m){
        meanrank <- function(x){rank(x)/length(x)}
        apply(m, 2, meanrank)
}

ccmn_norm_mets_quant <- quantNorm(ccmn_norm_mets)

# Plot in HeatMap

HeatMap(add_column(as.data.frame(ccmn_norm_mets_quant), Group = rownames(ccmn_norm_mets), .before = colnames(ccmn_norm_mets)[1]),
        scale = 'column',
        ColSideColors = Cols, 
        colramp = redgreen(75),
        distmethod = "euclidean",
        aggmethod = "ward.D",
        cexRow = 0.65,
        dendrogram = "both")

tiff("heatmapCCMN_quant_oldGood.tiff", width = 5000, height = 5000, units = "px", pointsize = 50)

heatmap.2(as.matrix(t(ccmn_norm_mets_quant)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized, quantile-normalized", margins = c(6, 13), 
          keysize = 1,
          key.xtickfun=function() {
                  cex <- par("cex")*par("cex.axis")
                  side <- 1
                  line <- 0
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=FALSE, tick=FALSE))
          })

dev.off()

# After that almost all replicates still cluster together.

########################################################################################################################################################
#
# Divide the dendrogram obtained after Quantile Normalization in the 2 biggest branches and apply a Wilcoxon test to see what metabolites
# are different between those branches (alpha = 0.05)
#
########################################################################################################################################################

HCA_ccmn_quant <- hclust(dist(ccmn_norm_mets_quant, "euclidean"), "ward.D")

ggdendrogram(data = HCA_ccmn_quant, main = "HCA of all the strains after CCMN normalization and Quantile normalization")
ggsave(filename = "HCA_CCMN_quant_oldGood.pdf", plot = last_plot(), device = "pdf", 
       path = NULL,
       scale = 1, width = 40, height = 20, units = "cm",
       dpi = 300, limitsize = F)

HCA_ccmn_quant_split <- cutree(HCA_ccmn_quant, k = 2)

ccmn_quant_clust1 <- ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 1), ]
ccmn_quant_clust2 <- ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 2), ]


HCA_ccmn_quant_median_group1 <- apply(ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 1), ], 2, median)
HCA_ccmn_quant_median_group2 <- apply(ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 2), ], 2, median)

wilcres <- list()
for(i in 1:ncol(ccmn_norm_mets_quant)){
        wilcres[[i]] <- wilcox.test(ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 1), ][, i], 
                                    ccmn_norm_mets_quant[which(HCA_ccmn_quant_split == 2), ][, i])
        names(wilcres)[i] <- colnames(ccmn_norm_mets_quant)[i]
}

is_different <- c()
for(i in 1:length(wilcres)){
        if(wilcres[[i]][3] < 0.05){
                is_different[i] <- T
        }else{is_different[i] <- F}
        names(is_different)[i] <- names(wilcres)[i]
}

# Sort differential metabolites between Clusters 1 and 2 according to how different are, in a decrescent way.

topDiffMets_branches1_2 <- rev(sort(abs(HCA_ccmn_quant_median_group1[is_different] - HCA_ccmn_quant_median_group2[is_different])))

# Repeat with ccmn normalized data to check if quantile normalization has an effect on the differential metabolite 
# determination, first splitting the tree according to the pre-quantile normalization data, and after splitting it according to the
# post quantile normalization data.

HCA_CCMN_split <- cutree(HCA_CCMN, k = 2)

ccmn_clust1 <- ccmn_norm_mets[which(HCA_CCMN_split == 1), ]
ccmn_clust2 <- ccmn_norm_mets[which(HCA_CCMN_split == 2), ]


HCA_ccmn_median_group1 <- apply(ccmn_norm_mets[which(HCA_CCMN_split == 1), ], 2, median)
HCA_ccmn_median_group2 <- apply(ccmn_norm_mets[which(HCA_CCMN_split == 2), ], 2, median)

wilcres_undisc <- list()
for(i in 1:ncol(ccmn_norm_mets)){
        wilcres_undisc[[i]] <- wilcox.test(ccmn_norm_mets[which(HCA_CCMN_split == 1), ][, i], 
                                           ccmn_norm_mets[which(HCA_CCMN_split == 2), ][, i])
        names(wilcres_undisc)[i] <- colnames(ccmn_norm_mets)[i]
}

is_different_undisc <- c()
pvals_undisc1_2 <- c()
for(i in 1:length(wilcres_undisc)){
        if(wilcres_undisc[[i]][3] < 0.05){
                is_different_undisc[i] <- T
                pvals_undisc1_2 <- c(pvals_undisc1_2, wilcres_undisc[[i]][3])
        }else{is_different_undisc[i] <- F}
        names(is_different_undisc)[i] <- names(wilcres_undisc)[i]
}
names(pvals_undisc1_2) <- names(wilcres_undisc[is_different_undisc])

# Remove those with _? ---> Those are the ones that we're not sure about their ID.

rmAmbig <- function(topDiffMets){
        removed <- topDiffMets[-grep("_?", names(topDiffMets), fixed = T)] 
        return(removed)
}

pvals_undisc1_2 <- rmAmbig(pvals_undisc1_2)

topDiffMets_branches1_2_undisc <- rev(sort(abs(HCA_ccmn_median_group1[is_different_undisc] - HCA_ccmn_median_group2[is_different_undisc])))



topDiffMets_branches1_2_undisc <- rmAmbig(topDiffMets_branches1_2_undisc)

pvals_undisc1_2 <- pvals_undisc1_2[match(names(topDiffMets_branches1_2_undisc), names(pvals_undisc1_2))]

topDiffMets_branches1_2_undisc_KEGGIDs <- dictionary$`KEGG IDs`[match(names(topDiffMets_branches1_2_undisc), dictionary$Consensus)]



sum(names(topDiffMets_branches1_2) %in% names(topDiffMets_branches1_2_undisc), na.rm = T)

# 55 out of 64 of the reported differential metabolites are the same in pre and post quantile normalization

wilcres_undisc2 <- list()
for(i in 1:ncol(ccmn_norm_mets)){
        wilcres_undisc2[[i]] <- wilcox.test(ccmn_norm_mets[which(HCA_ccmn_quant_split == 1), ][, i], 
                                            ccmn_norm_mets[which(HCA_ccmn_quant_split == 2), ][, i])
        names(wilcres_undisc2)[i] <- colnames(ccmn_norm_mets)[i]
}

is_different_undisc2 <- c()
for(i in 1:length(wilcres_undisc2)){
        if(wilcres_undisc2[[i]][3] < 0.05){
                is_different_undisc2[i] <- T
        }else{is_different_undisc2[i] <- F}
        names(is_different_undisc2)[i] <- names(wilcres_undisc2)[i]
}

topDiffMets_branches1_2_undisc2 <- rev(sort(abs(HCA_ccmn_quant_median_group1[is_different_undisc2] - HCA_ccmn_quant_median_group2[is_different_undisc2])))

length(topDiffMets_branches1_2_undisc2)

# 62 out of 62 of the reported differential metabolites are the same in pre and post quantile normalization

########################################################################################################################################################
#
# Change names of the metabolites by the KEGG indexes, and obtain the pathways where each metabolite has a role. Then get the most represented
# pathways among the metabolites reported as differential: rank, among all the pathways that are related to the observed differential metabolites,
# according to the percentage of the total metabolites of the pathway we have among the differential metabolites.
#
########################################################################################################################################################

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
metlist <- keggList("compound")
reactionlist <- keggList("reaction")
pathwaylist <- keggList("pathway")
if(!require(limma)) install.packages("limma")
library(limma)



getRelatedPaths <- function(keggidxs, org = NULL, a = NULL, b = NULL){
        library("KEGGREST")
        mets <- keggList("compound")
        keggindexes <- list()
        
        if(is.null(org) == F){
                pathlist <- keggList("pathway", org)
        }else{pathlist <- keggList("pathway")}
        pathways <- list()
        path_counts <- c()
        for(j in 1:length(keggidxs)){
                paths <- keggLink("pathway", keggidxs[j])
                if(length(paths) > 0 & is.null(org) == F){
                        paths <- sapply(paths, gsub, pattern = "map", replacement = org)
                        pathways[[j]] <- paths[paths %in% names(pathlist)]
                }
                else if(length(paths) > 0 & is.null(org) ==T){
                        pathways[[j]] <- paths
                }else{pathways[[j]] <- paste("Any pathway found for", keggidxs[j])}
                names(pathways)[j] <- keggidxs[j]
                for(k in 1:length(pathways[[j]])){
                        if(is.na(strsplit(pathways[[j]], ":")[[k]][2])){
                                path_counts <- path_counts
                        }else{path_counts <- c(path_counts, strsplit(pathways[[j]], ":")[[k]][2])}
                        
                }
        }
        path_counts <- as.data.frame(table(path_counts))
        path_counts[, "Pathways"] <- NA
        for(l in 1:length(path_counts[, 1])){
                path_counts[l, 3] <- pathlist[[paste("path:", path_counts[l,1], sep = "")]]
        }
        path_counts <- path_counts[order(path_counts$Freq, decreasing = T), ]
        rownames(path_counts) <- c()
        
        if(!is.null(a) & !is.null(b) == T){
                pathdif <- cbind("Cluster 1" = a, "Cluster 2" = b)
                pathdif <- pathdif[match(keggidxs, rownames(pathdif)), ]
                numbPaths <- c()
                for(i in 1:length(keggidxs)){
                        numbPaths[i] <- length(pathways[[i]])
                }
                pathdif <- cbind(pathdif, "Number of Related Pathways" = numbPaths)
                
                pathdif <- pathdif[rep(1:nrow(pathdif), times = pathdif[, 3]), ]
                path_indexes <- c()
                path_names <- c()
                for(i in 1:length(pathways)){
                        for(j in 1:length(pathways[[i]])){
                                path_indexes <- c(path_indexes, pathways[[i]][j])
                                path_names <- c(path_names, pathlist[pathways[[i]][j]])
                        }
                }
                path_names[which(is.na(path_names))] <- path_indexes[which(is.na(path_names))]
                
                pathdif <- cbind.data.frame(pathdif, "Pathway Indexes" = path_indexes)
                pathdif <- cbind.data.frame(pathdif, "Pathway Names" = path_names)
                
                compPerPath <- list()
                if(!is.null(org) == T){
                        for(ii in 1:length(unique(pathdif[, 4]))){
                                compPerPath[[ii]] <- keggLink("compound", as.character(unique(sapply(pathdif[, 4], gsub, pattern = org, replacement = "map"))[ii]))
                        } 
                }
                if(is.null(org) == T){
                        for(ii in 1:length(unique(pathdif[, 4]))){
                                compPerPath[[ii]] <- keggLink("compound", as.character(unique(pathdif[, 4])[ii]))
                        }
                }
                names(compPerPath) <- as.character(unique(pathdif[, 5]))
                reprPath <- c()
                comp_in_path <- c()
                for(ii in 1:length(compPerPath)){
                        reprPath[ii] <- (length(compPerPath[[ii]][compPerPath[[ii]] %in% paste("cpd:", names(pathways), sep = "")])/length(compPerPath[[ii]]))*100
                        comp_in_path[ii] <- paste(compPerPath[[ii]][compPerPath[[ii]] %in% paste("cpd:", names(pathways), sep = "")], collapse = " ")
                }
                names(reprPath) <- names(compPerPath)
                
                pathdif_ord <- pathdif[order(pathdif[, 5]),]
                
                pos <- c()
                neg <- c()
                medPathDif <- c()
                for(ii in 1:length(unique(path_indexes))){
                        submat <- pathdif_ord[which(pathdif_ord[, 4] == unique(path_indexes)[ii]), c(1, 2)]
                        pos[ii] <- length(which((submat[, 1] - submat[, 2]) > 0))
                        neg[ii] <- length(which((submat[, 1] - submat[, 2]) < 0))
                        medPathDif[ii] <- median(submat[, 1] - submat[, 2])
                }
                names(pos) <- unique(path_indexes)
                names(neg) <- unique(path_indexes)
                names(medPathDif) <- unique(path_indexes)
                
                reprPath <- cbind.data.frame(reprPath, pos, neg, medPathDif, comp_in_path)
                
                reprPath <- reprPath[order(reprPath[, 1], decreasing = T), ]
                
                colnames(reprPath) <- c("% of Representation", "Higher in Clust 1", "Higher in Clust 2", "Median of Difference", "Comp. Indx.")
                
                return(list("Pathways per Metabolite" = pathways, "Pathway Counts" = path_counts, 
                            "Mets & Rel. Pathways" = pathdif, "Representation of Pathways"= reprPath))
        }else{
                return(list("Pathways per Metabolite" = pathways, 
                            "Pathway Counts" = path_counts))
        }
        
}

medsG1KEGGIDs_undisc <- HCA_ccmn_quant_median_group1[is_different_undisc]
medsG2KEGGIDs_undisc <- HCA_ccmn_quant_median_group2[is_different_undisc]

names(medsG1KEGGIDs_undisc) <- dictionary$`KEGG IDs`[match(names(medsG1KEGGIDs_undisc), dictionary$Consensus)]
names(medsG2KEGGIDs_undisc) <- dictionary$`KEGG IDs`[match(names(medsG2KEGGIDs_undisc), dictionary$Consensus)]

topDiffMets_branches1_2_Paths <- getRelatedPaths(keggidxs = topDiffMets_branches1_2_undisc_KEGGIDs[!is.na(topDiffMets_branches1_2_undisc_KEGGIDs)], 
                                                 org = "pae",
                                                 a = medsG1KEGGIDs_undisc[!is.na(names(medsG1KEGGIDs_undisc))], 
                                                 b = medsG2KEGGIDs_undisc[!is.na(names(medsG2KEGGIDs_undisc))])

doORA <- function(diffMetObjkt, allMetsObjkt, org = NULL, alpha = 0.05){
        if(!require(KEGGREST)) install.packages("KEGGREST")
        library(KEGGREST)
        if(is.null(org)){
                paths <- keggList("pathway")
        }else{
                paths <- keggList("pathway", org)
        }
        diffMet <- diffMetObjkt[!is.na(diffMetObjkt)]
        diffMet <- sapply("cpd:", diffMet, FUN = paste, sep = "")[, 1]
        totPaths <- unique(unlist(sapply(allMetsObjkt, keggLink, target = "pathway")))
        if(!is.null(org)){
                totPaths <- totPaths[gsub("map", replacement = org, totPaths) %in% names(paths)]
        }
        compsPerPath <- sapply(totPaths, keggLink, target = "compound")
        allComps <- unique(unlist(compsPerPath))
        allCompsLen <- length(allComps)
        contMat <- function(x) {
                compsInPath <- length(x)
                mat <- matrix(c(compsInPath, allCompsLen - compsInPath, sum(diffMet %in% x), sum(!diffMet %in% x)), 
                              ncol = 2, 
                              nrow = 2,
                              dimnames = list(c("in_path", "not_in_path"),
                                              c("met_not_interest", "met_in_interest")))
                return(mat)
        }
        contMatsPaths <- lapply(compsPerPath, contMat)
        fishRes <- lapply(contMatsPaths, fisher.test)
        filt <- function(x) x$p.value <= alpha
        vecTrue <- unlist(lapply(fishRes, filt))
        sign <- fishRes[vecTrue]
        pVals <- sapply(sign, function(f) f$p.value)
        signMat <- cbind.data.frame(pathwaylist[names(pathwaylist) %in% names(pVals)], pVals)
        colnames(signMat) <- c("Pathways", "p.values")
        return(signMat)
}

allMets <- dictionary$`KEGG IDs`[!grepl("_?", dictionary$Consensus, fixed = T)][!is.na(dictionary$`KEGG IDs`[!grepl("_?", dictionary$Consensus, fixed = T)])]

ORA_C1_2 <- doORA(diffMetObjkt = topDiffMets_branches1_2_undisc_KEGGIDs,
                  allMetsObjkt = allMets,
                  org = "pae")

ORA_C1_2[match(gsub(" -.*", rownames(topDiffMets_branches1_2_Paths$`Representation of Pathways`), replacement = ""), ORA_C1_2$Pathways), ]

########################################################################################################################################################
#
# Divide clusters 1 and 2 in its 2 major subclusters and determine differential metabolites between the subclusters (Mann Whitney test, alpha = 0.05), 
# and do the same "pathway enrichment" we did before.
#
########################################################################################################################################################

# Cluster 1

HCA_ccmn_clust1 <- hclust(dist(ccmn_clust1, "euclidean"), "ward.D")
HCA_ccmn_clust1_split <- cutree(HCA_ccmn_clust1, k = 2)

HCA_ccmn_median_group1.1 <- apply(ccmn_clust1[which(HCA_ccmn_clust1_split == 1), ], 2, median)
HCA_ccmn_median_group1.2 <- apply(ccmn_clust1[which(HCA_ccmn_clust1_split == 2), ], 2, median)

wilcres_clust1 <- list()
for(i in 1:ncol(ccmn_clust1)){
        wilcres_clust1[[i]] <- wilcox.test(ccmn_clust1[which(HCA_ccmn_clust1_split == 1), ][, i], 
                                           ccmn_clust1[which(HCA_ccmn_clust1_split == 2), ][, i])
        names(wilcres_clust1)[i] <- colnames(ccmn_norm_mets)[i]
}

is_different_clust1 <- c()
pvals_1.1_1.2 <- c()
for(i in 1:length(wilcres_clust1)){
        if(wilcres_clust1[[i]][3] < 0.05){
                is_different_clust1[i] <- T
                pvals_1.1_1.2 <- c(pvals_1.1_1.2, wilcres_clust1[[i]][3])
        }else{is_different_clust1[i] <- F}
        names(is_different_clust1)[i] <- names(wilcres_clust1)[i]
}
names(pvals_1.1_1.2) <- names(wilcres_clust1[is_different_clust1])

pvals_1.1_1.2 <- rmAmbig(pvals_1.1_1.2)

topDiffMets_branches1.1_1.2 <- rev(sort(abs(HCA_ccmn_median_group1.1[is_different_clust1] - HCA_ccmn_median_group1.2[is_different_clust1])))

topDiffMets_branches1.1_1.2<- rmAmbig(topDiffMets_branches1.1_1.2)

pvals_1.1_1.2 <- pvals_1.1_1.2[match(names(topDiffMets_branches1.1_1.2), names(pvals_1.1_1.2))]

topDiffMets_branches1.1_1.2_KEGGIDs <- dictionary$`KEGG IDs`[match(names(topDiffMets_branches1.1_1.2), dictionary$Consensus)]

topDiffMets_branches1.1_1.2

medsG1.1KEGGIDs <- HCA_ccmn_median_group1.1[is_different_clust1]
medsG1.2KEGGIDs <- HCA_ccmn_median_group1.2[is_different_clust1]

names(medsG1.1KEGGIDs) <- dictionary$`KEGG IDs`[match(names(medsG1.1KEGGIDs), dictionary$Consensus)]
names(medsG1.2KEGGIDs) <- dictionary$`KEGG IDs`[match(names(medsG1.2KEGGIDs), dictionary$Consensus)]

topDiffMets_branches1.1_1.2_Paths <- getRelatedPaths(keggidxs = topDiffMets_branches1.1_1.2_KEGGIDs[!is.na(topDiffMets_branches1.1_1.2_KEGGIDs)], 
                                                     org = "pae",
                                                     a = medsG1.1KEGGIDs[!is.na(names(medsG1.1KEGGIDs))], 
                                                     b = medsG1.2KEGGIDs[!is.na(names(medsG1.2KEGGIDs))])

ORA_C1.1_1.2 <- doORA(diffMetObjkt = topDiffMets_branches1.1_1.2_KEGGIDs,
                      allMetsObjkt = allMets,
                      org = "pae")

ORA_C1.1_1.2[match(gsub(" -.*", rownames(topDiffMets_branches1.1_1.2_Paths$`Representation of Pathways`), replacement = ""), ORA_C1.1_1.2$Pathways), ]

# Cluster 2

HCA_ccmn_clust2 <- hclust(dist(ccmn_clust2, "euclidean"), "ward.D")
HCA_ccmn_clust2_split <- cutree(HCA_ccmn_clust2, k = 2)

HCA_ccmn_median_group2.1 <- apply(ccmn_clust2[which(HCA_ccmn_clust2_split == 1), ], 2, median)
HCA_ccmn_median_group2.2 <- apply(ccmn_clust2[which(HCA_ccmn_clust2_split == 2), ], 2, median)

wilcres_clust2 <- list()
for(i in 1:ncol(ccmn_clust2)){
        wilcres_clust2[[i]] <- wilcox.test(ccmn_clust2[which(HCA_ccmn_clust2_split == 1), ][, i], 
                                           ccmn_clust2[which(HCA_ccmn_clust2_split == 2), ][, i])
        names(wilcres_clust2)[i] <- colnames(ccmn_norm_mets)[i]
}

is_different_clust2 <- c()
pvals_2.1_2.2 <- c()
for(i in 1:length(wilcres_clust2)){
        if(wilcres_clust2[[i]][3] < 0.05){
                is_different_clust2[i] <- T
                pvals_2.1_2.2 <- c(pvals_2.1_2.2, wilcres_clust2[[i]][3])
        }else{is_different_clust2[i] <- F}
        names(is_different_clust2)[i] <- names(wilcres_clust2)[i]
}
names(pvals_2.1_2.2) <- names(wilcres_clust2[is_different_clust2])

pvals_2.1_2.2 <- rmAmbig(pvals_2.1_2.2)
pvals_2.1_2.2 <- pvals_2.1_2.2[match(names(topDiffMets_branches2.1_2.2), names(pvals_2.1_2.2))]

topDiffMets_branches2.1_2.2 <- rev(sort(abs(HCA_ccmn_median_group2.1[is_different_clust2] - HCA_ccmn_median_group2.2[is_different_clust2])))

topDiffMets_branches2.1_2.2 <- rmAmbig(topDiffMets_branches2.1_2.2)

topDiffMets_branches2.1_2.2_KEGGIDs <- dictionary$`KEGG IDs`[match(names(topDiffMets_branches2.1_2.2), dictionary$Consensus)]

topDiffMets_branches2.1_2.2

medsG2.1KEGGIDs <- HCA_ccmn_median_group2.1[is_different_clust2]
medsG2.2KEGGIDs <- HCA_ccmn_median_group2.2[is_different_clust2]

names(medsG2.1KEGGIDs) <- dictionary$`KEGG IDs`[match(names(medsG2.1KEGGIDs), dictionary$Consensus)]
names(medsG2.2KEGGIDs) <- dictionary$`KEGG IDs`[match(names(medsG2.2KEGGIDs), dictionary$Consensus)]

topDiffMets_branches2.1_2.2_Paths <- getRelatedPaths(keggidxs = topDiffMets_branches2.1_2.2_KEGGIDs[!is.na(topDiffMets_branches2.1_2.2_KEGGIDs)], 
                                                     org = "pae",
                                                     a = medsG2.1KEGGIDs[!is.na(names(medsG2.1KEGGIDs))], 
                                                     b = medsG2.2KEGGIDs[!is.na(names(medsG2.2KEGGIDs))])

ORA_C2.1_2.2 <- doORA(diffMetObjkt = topDiffMets_branches2.1_2.2_KEGGIDs, 
                      allMetsObjkt = allMets,
                      org = "pae")

ORA_C2.1_2.2[match(gsub(" -.*", rownames(topDiffMets_branches2.1_2.2_Paths$`Representation of Pathways`), replacement = ""), ORA_C2.1_2.2$Pathways), ]


########################################################################################################################################################
#
# Substitute names of differential metabolites by the KEGG indexes and build a table with differential metabolites between clusters 1 and 2, 
# clusters 1.1 and 1.2 and 2.1 and 2.2.
#
########################################################################################################################################################


diffMets_1_2 <- topDiffMets_branches1_2_undisc_KEGGIDs[!is.na(topDiffMets_branches1_2_undisc_KEGGIDs)]
diffMets_1_2 <- unique(diffMets_1_2)

diffMets_1.1_1.2 <- topDiffMets_branches1.1_1.2_KEGGIDs[!is.na(topDiffMets_branches1.1_1.2_KEGGIDs)]
diffMets_1.1_1.2 <- unique(diffMets_1.1_1.2)

diffMets_2.1_2.2 <- topDiffMets_branches2.1_2.2_KEGGIDs[!is.na(topDiffMets_branches2.1_2.2_KEGGIDs)]
diffMets_2.1_2.2 <- unique(diffMets_2.1_2.2)

diffMets <- list(t(as.matrix(diffMets_1_2)), t(as.matrix(diffMets_1.1_1.2)), t(as.matrix(diffMets_2.1_2.2)))

if(!require(plyr)) install.packages("plyr")
library(plyr)
diffMets <- t(rbind.fill.matrix(diffMets))
colnames(diffMets) <- c("1 and 2", "1.1 and 1.2", "2.1 and 2.2")

diffMets_oldGood <- diffMets
save(diffMets_oldGood, file = "diffMets_oldGood.RData")
write.csv(diffMets, file = "diffMets_allClusts_oldGood.csv")

write.csv(topDiffMets_branches1_2_Paths$`Representation of Pathways`, file = "pathRep1_2_old.csv")
write.csv(topDiffMets_branches1.1_1.2_Paths$`Representation of Pathways`, file = "pathRep1.1_1.2_old.csv")
write.csv(topDiffMets_branches2.1_2.2_Paths$`Representation of Pathways`, file = "pathRep2.1_2.2_old.csv")

ccmn_norm_mets_good_old_diffMets <- ccmn_norm_mets[, dictionary$Consensus[match(unique(as.character(diffMets_oldGood)[!is.na(as.character(diffMets_oldGood))]), 
                                                                                dictionary$`KEGG IDs`)]]

source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs/Diff_metabolites_analysis/metabo_functions.R")

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

ccmn_norm_mets_good_old_diffMets_quant <- quantNorm(ccmn_norm_mets_good_old_diffMets)

ccmn_norm_mets_good_old_diffMets_quantMeds <- getStrainMedian(ccmn_norm_mets_good_old_diffMets_quant)



tiff("heatmapCCMN_diffMets_quant_oldGood.tiff", width = 3000, height = 3000, units = "px", pointsize = 50)

heatmap.2(as.matrix(t(ccmn_norm_mets_good_old_diffMets_quantMeds)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = unique(Cols), notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized, quantile-normalized, only differential metabolites", margins = c(7, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 1,
          cexCol = 1.2,
          cellnote = round(as.matrix(t(ccmn_norm_mets_good_old_diffMets_quantMeds)), 2),
          notecex = 0.7,
          key.xtickfun=function() {
                  cex <- par("cex")*par("cex.axis")
                  side <- 1
                  line <- 0
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=FALSE, tick=FALSE))
          })

dev.off()

