####################################################################################################################
# FUNCTIONS USED IN SWARMING ANALYSIS                 
####################################################################################################################

# Get means of swarming replicates 

getSwarmMeans <- function(swarmDat){
        swarmMeans <- c()
        for(i in seq_along(levels(swarmDat[, 1]))){
                vec <- swarmDat[which(swarmDat[, 1] == levels(swarmDat[, 1])[i]), 2]
                swarmMeans <- c(swarmMeans, mean(vec))
        }
        names(swarmMeans) <- levels(swarmDat[, 1])
        return(swarmMeans)
}

# Binarize swarming data according to a given threshold
binarizeSwarm <- function(x, threshold = -2){
        swarmMeans <- x
        binVec <- rep(NA, length(swarmMeans))
        names(binVec) <- names(swarmMeans)
        binVec[swarmMeans > threshold] <- "Swarmer"
        binVec[swarmMeans < threshold] <- "nonSwarmer"
        binVec <- as.factor(binVec)
        return(binVec)
}

# Modified version of doORA function to save time in multiple iterations (for check bias)
doORAMod <- function(diffMetObjkt, alpha = 0.05){
        if(!require(KEGGREST)) install.packages("KEGGREST")
        library(KEGGREST)
        paths <- paePaths
        diffMet <- diffMetObjkt[!is.na(diffMetObjkt)]
        diffMet <- sapply("cpd:", diffMet, FUN = paste, sep = "")[, 1]
        totPaths <- totPathsInMets
        
        
        
        compsPerPath <- compsPerPathAllPaths
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
        names(pVals) <- gsub("map", replacement = "pae", names(pVals))
        signMat <- cbind.data.frame(paths[match(names(pVals), names(paths))], pVals)
        colnames(signMat) <- c("Pathways", "p.values")
        return(signMat)
}

# Given a matrix of the loadings of each one of the components result of a OPLS-DA model, outputs the n variables 
# with most extreme values (positive and negative) of each component.
getExtremVals <- function(loadsMat, n = 5){
        extrVals <- list()
        for(i in 1:ncol(loadsMat)){
                extrVals[[i]] <- c(head(sort(loadsMat[, i]), n),
                                   tail(sort(loadsMat[, i]), n)) 
        }
        names(extrVals) <- colnames(loadsMat)
        uniqueExtrVars <- unique(unlist(lapply(extrVals, names)))
        return(list("extremeVals" = extrVals, 
                    "uniqueExtremeVars" = uniqueExtrVars))
}

# Group genes according to a given Jaccard similarity
jaccGroup <- function(genePresAbsObjkt, threshold = 0.05){
        if(!require(proxy)) install.packages('proxy')
        library(proxy)
        geneTab <- sapply(genePresAbsObjkt, as.integer)
        rownames(geneTab) <- rownames(genePresAbsObjkt)
        colnames(geneTab) <- make.names(colnames(geneTab), unique = T)
        jaccDist <- proxy::dist(geneTab, by_rows = F, method = "Jaccard", diag = T)
        jaccHCA <- hclust(jaccDist, method = "complete")
        jaccGroups <- cutree(jaccHCA, h = threshold)
        groups <- c()
        for(i in 1:max(jaccGroups)){
                groups <- c(groups, names(sort(jaccGroups)[which(sort(jaccGroups) == i)])[1])
        }
        geneTabGrouped <- geneTab[, groups]
        return(list("GroupedMat" = geneTabGrouped, "Groups" = jaccGroups))
}

# From the list of result matrixes of multiple iteration of random ORAs or FELLAs, obtains null distribution of
# p-values for each pathway.
getPathNullDistrs <- function(randEnrichList){
        allPaths <- as.character(unique(unlist(lapply(randEnrichList, function(x) x$Pathways))))
        pathNullDistrs <- list()
        for(i in seq_along(allPaths)){
                nullVec <- c()
                for(j in seq_along(randEnrichList)){
                        if(allPaths[i] %in% randEnrichList[[j]]$Pathways){
                                nullVec <- c(nullVec, 
                                             randEnrichList[[j]][randEnrichList[[j]]$Pathways == allPaths[i], ]$p.values)
                        }
                }
                pathNullDistrs[[i]] <- nullVec
        }
        names(pathNullDistrs) <- allPaths
        return(pathNullDistrs)
}