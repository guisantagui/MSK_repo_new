###################################################################################################################################
# Functions for doing the gene presence?absence analysis
###################################################################################################################################

# Given a metabolite abundance normalized matrix, gets the cluster each strain belong among 4 major clusters
getMetClusts <- function(normMets, methodDist = "euclidean", methodAgg = "ward.D"){
        distMat <- dist(normMets, method = methodDist)
        HCA <- hclust(distMat, method = methodAgg)
        split <- cutree(HCA, k = 2)
        distMat1 <- dist(normMets[split == 1, ], method = methodDist)
        HCA1 <- hclust(distMat1, method = methodAgg)
        split1 <- cutree(HCA1, k = 2)
        distMat2 <- dist(normMets[split == 2, ], method = methodDist)
        HCA2 <- hclust(distMat2, method = methodAgg)
        split2 <- cutree(HCA2, k = 2)
        metClusts <- c()
        clustVec <- c()
        for(i in 1:nrow(normMets)){
                if(split[i] == 1){
                        if(split1[names(split[i])] == 1){clustVec <- c(clustVec, "1.1")}
                        if(split1[names(split[i])] == 2){clustVec <- c(clustVec, "1.2")}
                }
                if(split[i] == 2){
                        if(split2[names(split[i])] == 1){clustVec <- c(clustVec, "2.1")}
                        if(split2[names(split[i])] == 2){clustVec <- c(clustVec, "2.2")}
                }
        }
        clustVec <- as.factor(clustVec)
        names(clustVec) <- names(split)
        strainNames <- unique(gsub("\\_.*|(PA14).*", rep = "\\1", names(clustVec)))
        metClusts <- c()
        for(i in seq_along(strainNames)){
                metClusts <- c(metClusts, 
                               names(which.max(table(clustVec[gsub("\\_.*|(PA14).*", 
                                                                   rep = "\\1", 
                                                                   names(clustVec)) == strainNames[i]]))))
        }
        metClusts <- as.factor(metClusts)
        names(metClusts) <- strainNames
        return(metClusts)
}

# Given the gene presence/absence matrix and the metClust vector, yields a contingency table for each gene and cluster
getContTab <- function(genePresAbsObjkt, metClustObjkt){
        cont_tab <- c()
        groups <- sort(unique(metClustObjkt))
        for(i in seq_along(groups)){
                cont_tab <- rbind(cont_tab, colSums(genePresAbsObjkt[metClustObjkt == groups[i], ]))
        }
        rownames(cont_tab) <- groups
        return(cont_tab)
}

# Does a Fisher test between the 4 major clusters for each gene
doFisher <- function(cont_tab, 
                     metClustObjkt,
                     alpha = 0.05){
        pb = txtProgressBar(min = 0, max = ncol(cont_tab), initial = 0)
        clust_total <- table(metClustObjkt)
        fishResult <- matrix(nrow = ncol(cont_tab), ncol = 2, dimnames = list(NULL, c("Gene", "p.value")))
        for(i in 1:ncol(cont_tab)){
                geneContTab <- rbind(cont_tab[, i], clust_total - cont_tab[, i])
                fishResult[i, 1] <- colnames(cont_tab)[i]
                fishResult[i, 2] <- fisher.test(geneContTab)[[1]]
                setTxtProgressBar(pb, i)
        }
        fishSign <- fishResult[fishResult[, 2] <= 0.05, ]
        return(list("all.genes" = fishResult,
                    "sign.genes" = fishSign))
        
}

# Given a metabolite abundance normalized matrix and a gene presence/absence matrix, performs a Mann Whitney test for 
# each metabolite versus each gene, according to if the strains have or not the gene: determines if the strains that 
# have or not have a gene are different in the terms of the abundance of each metabolite. The output is a matrix of 
# p-values and a matrix of z-scores.

getMannWhitPerGene <- function(genePresAbsObjkt, 
                               metMatObjkt, 
                               p_adjust = T, 
                               method = NULL, 
                               remove.underrepresented = T,
                               threshold = 1){
        if(!require(uwIntroStats)) install.packages('uwIntroStats')
        library(uwIntroStats)
        pb = txtProgressBar(min = 0, max = ncol(genePresAbsObjkt), initial = 0)
        mannWhitPerGene <- matrix(nrow = ncol(genePresAbsObjkt),
                                  ncol = ncol(metMatObjkt))
        zScore <- matrix(nrow = ncol(genePresAbsObjkt),
                         ncol = ncol(metMatObjkt))
        for(j in 1:ncol(genePresAbsObjkt)){
                genePosi <- rownames(genePresAbsObjkt)[genePresAbsObjkt[, j] == 1]
                geneNega <- rownames(genePresAbsObjkt)[genePresAbsObjkt[, j] == 0]
                for(h in 1:ncol(metMatObjkt)){
                        if(remove.underrepresented == T){
                                if(length(genePosi) > threshold && length(geneNega) > threshold){
                                        p_val <- wilcoxon(metMatObjkt[grepl(paste(genePosi, collapse = "|"), 
                                                                            rownames(metMatObjkt)), h], 
                                                          metMatObjkt[grepl(paste(geneNega, collapse = "|"),
                                                                            rownames(metMatObjkt)), h])$p.value
                                        mannWhitPerGene[j, h] <- p_val
                                        z_score <- as.numeric(wilcoxon(metMatObjkt[grepl(paste(genePosi, collapse = "|"), 
                                                                                         rownames(metMatObjkt)), h],
                                                                       metMatObjkt[grepl(paste(geneNega, collapse = "|"),
                                                                                         rownames(metMatObjkt)), h])$inf[2, 1])
                                        zScore[j, h ] <- z_score
                                }  
                        }else{
                                p_val <- wilcoxon(metMatObjkt[grepl(paste(genePosi, collapse = "|"), 
                                                                    rownames(metMatObjkt)), h], 
                                                  metMatObjkt[grepl(paste(geneNega, collapse = "|"),
                                                                    rownames(metMatObjkt)), h])$p.value
                                mannWhitPerGene[j, h] <- p_val
                                z_score <- as.numeric(wilcoxon(metMatObjkt[grepl(paste(genePosi, collapse = "|"), 
                                                                                 rownames(metMatObjkt)), h],
                                                               metMatObjkt[grepl(paste(geneNega, collapse = "|"),
                                                                                 rownames(metMatObjkt)), h])$inf[2, 1])
                                zScore[j, h ] <- z_score
                        }
                        
                }
                if(p_adjust == T && is.null(method)){
                        stop( "Please select p-adjusting method")
                }else if(p_adjust == T && !is.null(method)){
                        if(anyNA(mannWhitPerGene[j, ]) == F){
                                mannWhitPerGene[j, ] <- p.adjust(mannWhitPerGene[j, ], method)
                        }
                }
                setTxtProgressBar(pb, j)
        }
        rownames(mannWhitPerGene) <- colnames(genePresAbsObjkt)
        colnames(mannWhitPerGene) <- colnames(metMatObjkt)
        mannWhitPerGene <- mannWhitPerGene[!is.na(rowSums(mannWhitPerGene)), ]
        
        rownames(zScore) <- colnames(genePresAbsObjkt)
        colnames(zScore) <- colnames(metMatObjkt)
        zScore <- zScore[!is.na(rowSums(zScore)), ]
        return(list("p.value" = mannWhitPerGene,
                    "z.score" = zScore))
}






# Filters the p-value matrix obtained in previous step (mannWhitPerGeneObjkt) to a given alpha (default: alpha = 0.05)
filtMannWhitPerGene <- function(mannWhitPerGeneObjkt, alpha = 0.05){
        if(!require(plyr)) install.packages('plyr')
        library(plyr)
        mannWhitPerGeneObjkt <- mannWhitPerGeneObjkt$p.value
        FUN <- function(x) t(as.matrix(names(x)[(x <= alpha)]))
        difMetsPerGene <- apply(mannWhitPerGeneObjkt, 1, FUN)
        difMetsPerGeneMat <- t(rbind.fill.matrix(difMetsPerGene))
        colnames(difMetsPerGeneMat) <- names(difMetsPerGene)
        return(list("as_list" = difMetsPerGene,
                    "as_mat" = difMetsPerGeneMat))
}

# Given the diffMetObjkt (differential metabolites between each pair of clusters) and the object generated in past 
# function (mannWhitPerGeneFiltObjkt) gives the names of all the genes whose presence/absence is correlated to 
# differences in abundance of each one of the differential metabolites between clusters (the ones that appear in
# diffMetObjkt)
getDiffGenePerDiffMet <- function(mannWhitPerGeneFiltObjkt, diffMetObjkt){
        pb = txtProgressBar(min = 0, 
                            max = length(diffMetObjkt[!is.na(diffMetObjkt)]), 
                            initial = 0)
        pbAv <- 0
        mannWhit <- mannWhitPerGeneFiltObjkt$as_list
        diffGenesPerDiffMet <- list()
        for(i in 1:ncol(diffMetObjkt)){
                genesRelated2Met <- list()
                for(j in 1:length(diffMetObjkt[, i][!is.na(diffMetObjkt[, i])])){
                        FUN <- function(x) any(x %in% diffMetObjkt[j, i])
                        genes <- names(which(lapply(mannWhit, FUN) == T))
                        if(length(genes) > 0){
                                genesRelated2Met[[j]] <- genes
                        }else{
                                genesRelated2Met[[j]] <- "Any relationship between the gene and the metabolite"
                        }
                        pbAv <- pbAv + 1
                        setTxtProgressBar(pb, pbAv)
                }
                names(genesRelated2Met) <- diffMetObjkt[, i][!is.na(diffMetObjkt[, i])]
                diffGenesPerDiffMet[[i]] <- genesRelated2Met
        }
        names(diffGenesPerDiffMet) <- colnames(diffMetObjkt)
        return(diffGenesPerDiffMet)
}

# Given the presence/absence matrix (genePresAbsObjkt), the output of the past function (DGenesPerDMetsObjkt) and the 
# metClustObjkt (cluster each strain belongs to), gives a matrix for each one of the differential metabolites between 
# each clusters. In that matrix appear the genes that appect that metabolite, indicating which strain have it and 
# dividing the strains according to the cluster they belong to, as well as a column of totals per strain in each matrix.
getGeneMatsPerMet <- function(DGenesPerDMetsObjkt, 
                              genePresAbsObjkt,
                              metClustObjkt){
        geneMatsPerMet <- list()
        for(i in seq_along(DGenesPerDMetsObjkt)){
                groups <- strsplit(names(DGenesPerDMetsObjkt)[i], "&")[[1]]
                metMats <- list()
                notEmpty <- DGenesPerDMetsObjkt[[i]][DGenesPerDMetsObjkt[[i]] != "Any relationship between the gene and the metabolite"]
                for(h in seq_along(notEmpty)){
                        if(i == 1){
                                group1 <- names(metClustObjkt[gsub("\\..*", rep = "", metClustObjkt) == groups[1]])
                                group2 <- names(metClustObjkt[gsub("\\..*", rep = "", metClustObjkt) == groups[2]])
                        }else{
                                group1 <- names(metClustObjkt[metClustObjkt == groups[1]])
                                group2 <- names(metClustObjkt[metClustObjkt == groups[2]])
                        }
                        metMat1 <- as.matrix(genePresAbsObjkt[group1,
                                                              notEmpty[[h]]])
                        metMat2 <- as.matrix(genePresAbsObjkt[group2,
                                                              notEmpty[[h]]])
                        metMat <- rbind(metMat1, 
                                        rep(NA, length(notEmpty[[h]])),
                                        metMat2)
                        colnames(metMat) <- notEmpty[[h]]
                        metMat <- cbind(metMat, Total = rowSums(metMat))
                        metMats[[h]] <- metMat
                }
                names(metMats) <- names(notEmpty)
                geneMatsPerMet[[i]] <- metMats
        }
        names(geneMatsPerMet) <- names(DGenesPerDMetsObjkt)
        return(geneMatsPerMet)
}

# Takes as an input the matrixes generated in past step (geneMatMetObjkt), and does a Mann Whitney test between the vectors of strains corresponding to 
# each pair of clusters gene by gene by one side, and between totals, by the other side. 
mannWhit_tabs <- function(geneMatMetObjkt){
        mannWhitTabs <- list()
        allDiffMets <- lapply(geneMatMetObjkt, names)
        allDiffMets <- unique(unlist(allDiffMets))
        totalMannWhit <- matrix(nrow =length(allDiffMets), 
                                ncol = length(geneMatMetObjkt), 
                                dimnames = list(allDiffMets, 
                                                names(geneMatMetObjkt)))
        for(ii in seq_along(geneMatMetObjkt)){
                allDiffGenes <- lapply(geneMatMetObjkt[[ii]], colnames)
                allDiffGenes <- unique(unlist(allDiffGenes))
                allDiffGenes <- allDiffGenes[-which(allDiffGenes == "Total")]
                
                equator <- which(is.na(geneMatMetObjkt[[ii]][[1]][, 1]))
                upper <- length(geneMatMetObjkt[[ii]][[1]][, 1])
                
                metGenes <- matrix(nrow = length(geneMatMetObjkt[[ii]]), ncol = length(allDiffGenes))
                rownames(metGenes) <- names(geneMatMetObjkt[[ii]])
                colnames(metGenes) <- allDiffGenes
                for(j in 1:length(geneMatMetObjkt[[ii]])){
                        for(i in 1:(ncol(geneMatMetObjkt[[ii]][[j]]) - 1)){
                                metGenes[names(geneMatMetObjkt[[ii]][j]), 
                                         colnames(geneMatMetObjkt[[ii]][[j]])[i]] <- wilcox.test(as.numeric(geneMatMetObjkt[[ii]][[j]][1:(equator - 1), i]), 
                                                                                                 as.numeric(geneMatMetObjkt[[ii]][[j]][(equator + 1):upper, i]))[[3]]
                        }
                        
                        totalMannWhit[names(geneMatMetObjkt[[ii]])[j], ii] <- wilcox.test(as.numeric(geneMatMetObjkt[[ii]][[j]][1:(equator - 1), "Total"]), 
                                                                                          as.numeric(geneMatMetObjkt[[ii]][[j]][(equator + 1):upper, "Total"]))[[3]]        
                }
                mannWhitTabs[[ii]] <- metGenes
        }
        names(mannWhitTabs) <- names(geneMatMetObjkt)
        return(list("Individual_genes" = mannWhitTabs,
                    "Totals" = totalMannWhit))
}

# Filters the result of mannWhitTabsObjkt (the p-values generated in past function) to a given alpha (default: alpha = 0.05)
filtMannWhitTabs <- function(manWhitTabsObjkt, alpha = 0.05){
        single <- manWhitTabsObjkt$Individual_genes
        totals <- manWhitTabsObjkt$Totals
        filteredTabs <- list()
        filteredTabs_Total <- list()
        for(ii in seq_along(single)){
                diffGenesPerMet <- list()
                for(i in 1:nrow(single[[ii]])){
                        diffGenes <- c()
                        for(j in 1:ncol(single[[ii]])){
                                if(!is.na(single[[ii]][i, j]) && single[[ii]][i, j] <= alpha){
                                        diffGenes <- c(diffGenes, colnames(single[[ii]])[j])
                                }
                        }
                        if(is.null(diffGenes)){
                                diffGenes <- "Any correlated gene with the metabolite found"
                        }
                        diffGenesPerMet[[i]] <- diffGenes 
                }
                names(diffGenesPerMet) <- rownames(single[[ii]])
                filteredTabs[[ii]] <- diffGenesPerMet
                
                genesMet <- list()
                totals_sign <- names(which(totals[, ii] <= alpha))
                for(i in seq_along(totals_sign)){
                        trueVec <- !is.na(single[[ii]][totals_sign[i], ])
                        genesMet[[i]] <- names(single[[ii]][totals_sign[i], ])[trueVec]
                }
                names(genesMet) <- totals_sign
                filteredTabs_Total[[ii]] <- genesMet
        }
        names(filteredTabs) <- names(single)
        names(filteredTabs_Total) <- names(single)
        return(list("Individual_genes" = filteredTabs,
                    "Total" = filteredTabs_Total))
}