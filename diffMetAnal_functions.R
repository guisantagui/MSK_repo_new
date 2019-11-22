
# Quantile normalizes metabolomic matrix
quantNorm <- function(m){
        meanrank <- function(x){rank(x)/length(x)}
        apply(m, 2, meanrank)
}

# Removes the metabolites with "_?" (the tag we used to label ambiguous metabolites)
rmAmbig <- function(topDiffMets){
        removed <- topDiffMets[-grep("_?", names(topDiffMets), fixed = T)] 
        return(removed)
}

# Gets the pathways per differential metabolite, the counts each pathway appears and the percentage of representation, 
# for a certain organism.
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

# Does Over representation analysis given a vector of differential metabolites, the total metabolites detected, and a 
# level of significance, for a certain organism.
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

# Gets the median of the replicates of the metabolomic matrix.
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