# given a vector of compounds' KEGG IDs and a pathway ID retrieves all the enzymes related with these compounds that are in
# this pathway.
getEnzsRelWmetsInPath <- function(metIDVec, pathID){
        pathEnzs <- keggLink("enzyme", pathID)
        enzsPerMet <- sapply(metIDVec,
                             function(x) gsub("ec:",
                                              "",
                                              keggLink(x, target = "enzyme")[keggLink(x, target = "enzyme") %in% pathEnzs]))
        return(enzsPerMet)
} 

# given a list of enzymes per metabolite gives the core enzymes

getCoreEnzsXMet <- function(enzsPerMetObjkt){
        coreEnzsXMet <- list()
        for(i in seq_along(enzsPerMetObjkt)){
                enzsMetTab <- dictEnzymes[which(dictEnzymes$ECnums %in% enzsPerMetObjkt[[i]]), ]
                m <- enzsMetTab[, 4:(ncol(enzsMetTab)-1)]
                coreEnzs <- enzsMetTab$Gene[apply(m, 1, function(x) sum(nchar(x) > 0) == ncol(m))]
                coreEnzsXMet[[i]] <- coreEnzs
        }
        names(coreEnzsXMet) <- names(enzsPerMetObjkt)
        coreEnzsXMet <- coreEnzsXMet[sapply(coreEnzsXMet, function(x) length(x) > 0)]
        return(coreEnzsXMet)
}

# Given an alignment and the rhamnolipid production matrix (1 & 2 categories) outputs a dataframe
# consisting in the SNPs for the given alignment, with the 2 columns appended correspondent to rhamn
# phenotype. This matriz is intended to be the input of doRForest function. 

buildSNPsMat <- function(alignment, rhamnMat){
        alignMat <- as.matrix(alignment)
        SNPsLoc <-which(apply(alignMat, 2, function(x) length(unique(x)) > 1))
        if(length(SNPsLoc) == 0){
                return("Any SNPs in this alignment")
        }else{
                alignMatSNPS <- as.data.frame(alignMat[, SNPsLoc])
                colnames(alignMatSNPS) <- paste("pos", sep = "_", as.character(SNPsLoc))
                alignMatSNPS$rhamn2cats <- as.factor(rhamnMat$rhamn2cats[match(rownames(alignMatSNPS), 
                                                                               rhamnMat$strains)])
                alignMatSNPS$rhamn3cats <- as.factor(rhamnMat$rhamn3cats[match(rownames(alignMatSNPS), 
                                                                               rhamnMat$strains)])
                return(alignMatSNPS)
        }
        
}

# Transforms the previous matrix into a 0/1 one.

binaryCodeSNPs <- function(SNPsMat){
        SNPsMatBin <- data.frame()
        for(i in 1:(ncol(SNPsMat)-2)){
                uniqueVars <- sort(unique(SNPsMat[, i]))
                binarized <- data.frame()
                for(j in seq_along(uniqueVars)){
                        binVec <- rep(0, nrow(SNPsMat))
                        binVec[SNPsMat[, i] == uniqueVars[j]] <- 1
                        binVec <- as.factor(binVec)
                        if(length(as.vector(binarized)) == 0){
                                binarized <- binVec
                        }else{
                                binarized <- cbind.data.frame(binarized, binVec)
                        }
                }
                colNames <- paste(colnames(SNPsMat)[i],
                                  uniqueVars, sep = "_")
                colNames <- make.names(colNames)
                colnames(binarized) <- colNames
                if(ncol(SNPsMatBin) == 0){
                        SNPsMatBin <- binarized
                }else{
                        SNPsMatBin <- cbind.data.frame(SNPsMatBin, binarized)
                }
                
        }
        SNPsMatBin <- cbind.data.frame(SNPsMatBin, SNPsMat[, (ncol(SNPsMat)-1):ncol(SNPsMat)])
        return(SNPsMatBin)
}

# Given a SNPs matrix (binarized or non binarized) does a random forest classification, using as Y
# the 2 categories and the 3, outputting accuracy and p-value of both.
doRForest <- function(SNPsMat, ntree = 500, nCV = 10){
        trControl <- trainControl(method = "cv",
                                  number = nCV,
                                  search = "grid",
                                  allowParallel = T)
        w2catsFit <- train(rhamn2cats ~ . ,
                           data = SNPsMat[,-(ncol(SNPsMat))],
                           method = "rf",
                           metric = "Accuracy",
                           trControl = trControl)
        w2catsPred <- predict(w2catsFit, SNPsMat[,-(ncol(SNPsMat))])
        w2confMat <- confusionMatrix(w2catsPred, SNPsMat$rhamn2cats)
        w2catsAcc <- w2confMat$overall[1]
        w2catsAccPVal <- w2confMat$overall[6]
        
        w3catsFit <- train(rhamn3cats ~ . ,
                           data = SNPsMat[,-(ncol(SNPsMat)-1)],
                           method = "rf",
                           metric = "Accuracy",
                           trControl = trControl)
        w3catsPred <- predict(w3catsFit, SNPsMat[,-(ncol(SNPsMat)-1)])
        w3confMat <- confusionMatrix(w3catsPred, SNPsMat$rhamn3cats)
        w3catsAcc <- w3confMat$overall[1]
        w3catsAccPVal <- w3confMat$overall[6]
        output <- c(w2catsAcc, w2catsAccPVal, w3catsAcc, w3catsAccPVal)
        names(output) <- c("2_cats_Acc", "2_cats_Acc_pValue", "3_cats_Acc", "3_cats_Acc_pValue")
        return(output)
}

# This monster does, given a list of core enzymes in a pathway (output of getCoreEnzsXMet), dendrograms for each 
# alignment and rforests. Does a rforest consisting of all the snps of the inputed enzymes (global)
retrAlignAndDends <- function(coreEnzsPath, aggrMethod = "complete", doDendrograms = F){
        if(!require(caret)) install.packages("caret")
        library(caret)
        unlisted <- unique(unlist(coreEnzsPath))
        coreInParsedObjkt <- unlisted[unlisted %in% names(ParsedSubGraphs_allStrains_named_new)]
        output <- list()
        allSNPsInPathDNA <- data.frame()
        allSNPsInPathAA <- data.frame()
        DNALabels <- list()
        AALabels <- list()
        rForestDNA <- c()
        rForestAA <- c()
        for(i in seq_along(coreInParsedObjkt)){
                geneDat <- ParsedSubGraphs_allStrains_named_new[names(ParsedSubGraphs_allStrains_named_new) == coreInParsedObjkt[i]][[1]]
                DNADist <- geneDat$DNADist
                strainNames4Gene <- allSeqs[as.numeric(gsub("\\_.*", "", rownames(DNADist)))]
                colnames(DNADist) <- strainNames4Gene
                rownames(DNADist) <- strainNames4Gene
                AADist <- geneDat$AADist
                colnames(AADist) <- strainNames4Gene
                rownames(AADist) <- strainNames4Gene
                DNADist <- as.dist(DNADist)
                AADist <- as.dist(AADist)
                DNAAlign <- geneDat$DNAAlignment
                names(DNAAlign) <- strainNames4Gene
                AAAlign <- geneDat$AAAlignment
                names(AAAlign) <- strainNames4Gene
                if(doDendrograms == T){
                        DNAHCA <- IdClusters(myDistMatrix = DNADist, myXStringSet = DNAAlign, method = aggrMethod,
                                             type = "dendrogram")
                        AAHCA <- IdClusters(myDistMatrix = AADist, myXStringSet = AAAlign, method = aggrMethod,
                                            type = "dendrogram")
                        DNADendDat <- dendro_data(DNAHCA)
                        AADendDat <- dendro_data(AAHCA)
                        DNALabs <- label(DNADendDat)
                        DNALabs$rhamn <- as.factor(rhamnMat$rhamn3cats[match(DNALabs$label, rhamnMat$strains)])
                        AALabs <- label(AADendDat)
                        AALabs$rhamn <- as.factor(rhamnMat$rhamn3cats[match(AALabs$label, rhamnMat$strains)])
                        DNALabels[[i]] <- DNALabs
                        AALabels[[i]] <- AALabs
                        dendDNA <- ggplot(segment(DNADendDat)) +
                                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                                geom_text(data=label(DNADendDat),
                                          aes(label=label, x=x, y=-0.0001, colour=DNALabs$rhamn, angle = 90, hjust = "top"),
                                          nudge_y = -0.0001) +
                                ylim(-0.001, 0.008) + 
                                scale_color_manual(values = c("blue", "green","red"))
                        dendAA <- ggplot(segment(AADendDat)) +
                                geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
                                geom_text(data=label(AADendDat),
                                          aes(label=label, x=x, y=-0.0001, colour=AALabs$rhamn, angle = 90, hjust = "top"),
                                          nudge_y = -0.0001) +
                                ylim(-0.001, 0.008) + 
                                scale_color_manual(values = c("blue", "green","red"))
                        
                        dendrograms <- list("DNA" = dendDNA,
                                            "AA" = dendAA)
                        geneOutObjkt <- dendrograms
                        output[[i]] <- geneOutObjkt
                }
                
                
                
                alignSNPsDNA <- buildSNPsMat(DNAAlign, rhamnMat)
                
                alignSNPsAA <- buildSNPsMat(AAAlign, rhamnMat)
                
                if(alignSNPsDNA == "Any SNPs in this alignment"){
                        forestDNA <- rep(NA, 4)
                }else{
                        alignSNPsDNA <- binaryCodeSNPs(alignSNPsDNA)
                        forestDNA <- doRForest(alignSNPsDNA)
                }
                if(alignSNPsAA == "Any SNPs in this alignment"){
                        forestAA <- rep(NA, 4)
                }else{
                        alignSNPsAA <- binaryCodeSNPs(alignSNPsAA)
                        forestAA <- doRForest(alignSNPsAA)
                }
                colNamesForests <- c("2_cats_Acc",
                                     "2_cats_Acc_pValue",
                                     "3_cats_Acc",
                                     "3_cats_Acc_pValue")
                forestDNA <- matrix(forestDNA, ncol = 4)
                colnames(forestDNA) <- colNamesForests
                forestAA <- matrix(forestAA, ncol = 4)
                colnames(forestAA) <- colNamesForests
                
                if(length(as.vector(rForestDNA)) == 0){
                        rForestDNA <- forestDNA
                }else{
                        rForestDNA <- rbind(rForestDNA, forestDNA)
                }
                if(length(as.vector(rForestAA)) == 0){
                        rForestAA <- forestAA
                }else{
                        rForestAA <- rbind(rForestAA, forestAA)
                }
                
                if(alignSNPsDNA != "Any SNPs in this alignment"){
                        SNPsDNA2bind <- alignSNPsDNA[, 1:(ncol(alignSNPsDNA)-2), drop = F]
                        colnames(SNPsDNA2bind) <- paste(unlisted[i], colnames(SNPsDNA2bind), sep = "_")
                        if(ncol(allSNPsInPathDNA) == 0){
                                allSNPsInPathDNA <- SNPsDNA2bind  
                        }else{
                                allSNPsInPathDNA <- cbind.data.frame(allSNPsInPathDNA,
                                                                     SNPsDNA2bind)
                        }
                }
                if(alignSNPsAA != "Any SNPs in this alignment"){
                        print(ncol(alignSNPsAA))
                        SNPsAA2bind <- alignSNPsAA[, 1:(ncol(alignSNPsAA)-2), drop = F]
                        colnames(SNPsAA2bind) <- paste(unlisted[i], colnames(SNPsAA2bind), sep = "_")
                        if(ncol(allSNPsInPathAA) == 0){
                                allSNPsInPathAA <- SNPsAA2bind  
                        }else{
                                allSNPsInPathAA <- cbind.data.frame(allSNPsInPathAA,
                                                                    SNPsAA2bind)
                        }
                }
        }
        allSNPsInPathDNA$rhamn2cats <- as.factor(rhamnMat$rhamn2cats[match(rownames(allSNPsInPathDNA), 
                                                                           rhamnMat$strains)])
        allSNPsInPathDNA$rhamn3cats <- as.factor(rhamnMat$rhamn3cats[match(rownames(allSNPsInPathDNA), 
                                                                           rhamnMat$strains)])
        allSNPsInPathAA$rhamn2cats <- as.factor(rhamnMat$rhamn2cats[match(rownames(allSNPsInPathAA), 
                                                                          rhamnMat$strains)])
        allSNPsInPathAA$rhamn3cats <- as.factor(rhamnMat$rhamn3cats[match(rownames(allSNPsInPathAA), 
                                                                          rhamnMat$strains)])
        rForestAllDNA <- doRForest(allSNPsInPathDNA)
        rForestAllAA <- doRForest(allSNPsInPathAA)
        
        globalRForest <- list("DNA" = rForestAllDNA,
                              "AA" = rForestAllAA)
        allSNPs <- list("DNA" = allSNPsInPathDNA,
                        "AA" = allSNPsInPathAA)
        global <- list("SNPs" = allSNPs,
                       "rForest" = globalRForest)
        
        print(rForestDNA)
        print(rForestAA)
        
        rForestMat <- cbind(rForestDNA, rForestAA)
        rownames(rForestMat) <- coreInParsedObjkt
        colnames(rForestMat)[1:4] <- paste("DNA", colnames(rForestMat)[1:4])
        colnames(rForestMat)[5:8] <- paste("AA", colnames(rForestMat)[5:8])
        
        if(doDendrograms == T){
                names(output) <- coreInParsedObjkt
                names(DNALabels) <- coreInParsedObjkt
                names(AALabels) <- coreInParsedObjkt
                labels <- list("DNA"= DNALabels,
                               "AA" = AALabels)
                return(list(dendrograms = output,
                            individualRforest = rForestMat,
                            global = global,
                            labels = labels))
        }else{
                return(list(individualRforest = rForestMat,
                            global = global))
        }
}

# Given a vector of genes performes a rForest of the SNPs, both for DNA and AA sequence, 
# and with 2 categories and 3 categories as Y. 
doMultRForests <- function(geneVec, nCV = 10, maxWaitTime = 30){
        DNArForests <- c()
        AArForests <- c()
        pb = txtProgressBar(min = 0, max = length(geneVec), initial = 0)
        for(i in seq_along(geneVec)){
                print(i)
                posGene <- match(geneVec[i], 
                                 names(ParsedSubGraphs_allStrains_named_new))
                DNAAlignment <- ParsedSubGraphs_allStrains_named_new[[posGene]]$DNAAlignment
                AAAlignment <- ParsedSubGraphs_allStrains_named_new[[posGene]]$AAAlignment
                
                namesStrains <- allSeqs[as.numeric(gsub("\\_.*",
                                                        "", 
                                                        names(DNAAlignment)))]
                
                names(DNAAlignment) <- namesStrains
                names(AAAlignment) <- namesStrains
                
                DNA_SNPs <- buildSNPsMat(DNAAlignment, rhamnMat)
                AA_SNPs <- buildSNPsMat(AAAlignment, rhamnMat)
                
                if(DNA_SNPs == "Any SNPs in this alignment"){
                        rForestDNA <- rep(NA, 4)
                }else if(ncol(DNA_SNPs)-2 == 1 && length(table(DNA_SNPs[, 1])) == 2 && min(table(DNA_SNPs[, 1])) == 1){
                        rForestDNA <- rep(NA, 4)
                }else{
                        rForestDNA <- withTimeout({
                                doRForest(DNA_SNPs, nCV = nCV)}, 
                                timeout = maxWaitTime, 
                                onTimeout = "silent")
                }
                if(is.null(rForestDNA)){
                        rForestDNA <- rep(NA, 4)
                }
                if(AA_SNPs == "Any SNPs in this alignment"){
                        rForestAA <- rep(NA, 4)
                }else if(ncol(AA_SNPs)-2 == 1 && length(table(AA_SNPs[, 1])) == 2 && min(table(AA_SNPs[, 1])) == 1){
                        rForestAA <- rep(NA, 4)
                }else{
                        rForestAA <- withTimeout({
                                doRForest(AA_SNPs, nCV = nCV)}, 
                                timeout = maxWaitTime, 
                                onTimeout = "silent")
                }
                if(is.null(rForestAA)){
                        rForestAA <- rep(NA, 4)
                }
                colNamesForests <- c("2_cats_Acc",
                                     "2_cats_Acc_pValue",
                                     "3_cats_Acc",
                                     "3_cats_Acc_pValue")
                rForestDNA <- matrix(rForestDNA, ncol = 4)
                colnames(rForestDNA) <- colNamesForests
                rForestAA <- matrix(rForestAA, ncol = 4)
                colnames(rForestAA) <- colNamesForests
                
                if(length(as.vector(DNArForests)) == 0){
                        DNArForests <- rForestDNA
                }else{
                        DNArForests <- rbind(DNArForests, rForestDNA)
                }
                if(length(as.vector(AArForests)) == 0){
                        AArForests <- rForestAA
                }else{
                        AArForests <- rbind(AArForests, rForestAA)
                }
                setTxtProgressBar(pb, i)
        }
        rForestMat <- cbind(DNArForests, AArForests)
        rownames(rForestMat) <- geneVec
        colnames(rForestMat)[1:4] <- paste("DNA", colnames(rForestMat)[1:4])
        colnames(rForestMat)[5:8] <- paste("AA", colnames(rForestMat)[5:8])
        return(rForestMat)
}


# as in parsedSubGraph object some genes were tossed out with this gene we can do the alignment retrieving it from
# the positions in GFF file. Given a gene name it outputs the DNA and AA alignment. 
doAlignFromGeneCalls <- function(gene){
        geneCallGeneInfo <- lapply(geneCallsAllStrains4HeronProkka,
                                   function(x) x[grep(gene, x$Annotation), ])
        geneCallGeneInfo <- geneCallGeneInfo[sapply(geneCallGeneInfo, function(x) nrow(x) > 0)]
        geneSeqs <- c()
        for(i in seq_along(geneCallGeneInfo)){
                if(geneCallGeneInfo[[i]]$Strand == 0){
                        seq <- as.character(allGenomes[[grep(names(geneCallGeneInfo)[i], 
                                                             names(allGenomes))]][[1]][geneCallGeneInfo[[i]]$Start:geneCallGeneInfo[[i]]$Stop])
                }else{
                        seq <- allGenomes[[grep(names(geneCallGeneInfo)[i], 
                                                names(allGenomes))]][[1]][geneCallGeneInfo[[i]]$Stop:geneCallGeneInfo[[i]]$Start]
                        seq <- complement(seq)
                        seq <- as.character(seq)
                }
                
                geneSeqs <- c(geneSeqs, seq)
        }
        strSet <- DNAStringSet(geneSeqs)
        AAAlign <- AlignTranslation(strSet, type = "AAStringSet")
        names(AAAlign) <- names(geneCallGeneInfo)
        DNAAlign <- AlignTranslation(strSet, type = "DNAStringSet")
        names(DNAAlign) <- names(geneCallGeneInfo)
        return(list("DNA" = DNAAlign,
                    "AA" = AAAlign))
}
