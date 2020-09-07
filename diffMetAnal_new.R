setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/newData")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes as input CCMN normalized data and selects the differential metabolites between the 4 major groups  #################################
# of strains found according to metabolite abundance.                                                                  #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

#Load data
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/newData/ccmn_norm_mets_newData.RData")
ccmn_norm_mets <- ccmn_norm_mets_newData

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
ggsave(filename = "HCA_CCMN_newData.pdf", plot = last_plot(), device = "pdf", 
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



compKEGGIDs <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets), dictionary$`New Data Names`)]



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

tiff("heatmapCCMN_quant_newData.tiff", width = 5000, height = 5000, units = "px", pointsize = 50)

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
ggsave(filename = "HCA_CCMN_quant_newData.pdf", plot = last_plot(), device = "pdf", 
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
for(i in 1:length(wilcres_undisc)){
        if(wilcres_undisc[[i]][3] < 0.05){
                is_different_undisc[i] <- T
        }else{is_different_undisc[i] <- F}
        names(is_different_undisc)[i] <- names(wilcres_undisc)[i]
}

topDiffMets_branches1_2_undisc <- rev(sort(abs(HCA_ccmn_median_group1[is_different_undisc] - HCA_ccmn_median_group2[is_different_undisc])))

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
for(i in 1:length(wilcres_clust1)){
        if(wilcres_clust1[[i]][3] < 0.05){
                is_different_clust1[i] <- T
        }else{is_different_clust1[i] <- F}
        names(is_different_clust1)[i] <- names(wilcres_clust1)[i]
}

topDiffMets_branches1.1_1.2 <- rev(sort(abs(HCA_ccmn_median_group1.1[is_different_clust1] - HCA_ccmn_median_group1.2[is_different_clust1])))

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
for(i in 1:length(wilcres_clust2)){
        if(wilcres_clust2[[i]][3] < 0.05){
                is_different_clust2[i] <- T
        }else{is_different_clust2[i] <- F}
        names(is_different_clust2)[i] <- names(wilcres_clust2)[i]
}

topDiffMets_branches2.1_2.2 <- rev(sort(abs(HCA_ccmn_median_group2.1[is_different_clust2] - HCA_ccmn_median_group2.2[is_different_clust2])))

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

diffMets_newData <- diffMets
save(diffMets_newData, file = "diffMets_newData.RData")
write.csv(diffMets, file = "diffMets_allClusts_newData.csv")

