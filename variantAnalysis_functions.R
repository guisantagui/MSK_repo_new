# GIven the ParsedSubGraphsNamed object, compares hamming sequence distances for each gene intra and inter groups, and applies a Mann Whitney Test. 
# The output is given as the genes that are different at sequence level using AA sequence by one side, and by other theones that are different 
# using DNA sequence.

identifyDiffGenes <- function(ParsedSubGraphsNamed_objt, strainGroups, p_adjust = T, method = NA, alpha = 0.05){
        namesGroups <- names(strainGroups)
        strainGroups <- as.factor(as.integer(strainGroups))
        names(strainGroups) <- namesGroups
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
        signGenesClusts1_2DNA <- cbind.data.frame(namesGenes[which(mannWhitResDNA <= alpha)], mannWhitResDNA[which(mannWhitResDNA <= alpha)])
        signGenesClusts1_2DNA <- cbind.data.frame(signGenesClusts1_2DNA, distDifferDNA[which(mannWhitResDNA <= alpha)])
        colnames(signGenesClusts1_2DNA) <- c("Annotation", "p-value", "medians_difference")
        signGenesClusts1_2DNA <- signGenesClusts1_2DNA[order(signGenesClusts1_2DNA[, 3], decreasing = T), ]
        signGenesClusts1_2AA <- cbind.data.frame(namesGenes[which(mannWhitResAA <= alpha)], mannWhitResAA[which(mannWhitResAA <= alpha)])
        signGenesClusts1_2AA <- cbind.data.frame(signGenesClusts1_2AA, distDifferAA[which(mannWhitResAA <= alpha)])
        colnames(signGenesClusts1_2AA) <- c("Annotation", "p-value", "medians_difference")
        signGenesClusts1_2AA <- signGenesClusts1_2AA[order(signGenesClusts1_2AA[, 3], decreasing = T), ]
        results <- list("MannWhitney_meds" = mannWhitMeds, 
                        "significativeGenes_DNA" = signGenesClusts1_2DNA,
                        "significativeGenes_AA" = signGenesClusts1_2AA)
        return(results)
}

# Given the output from past function, the strain names and the parsedSubGraphNamed object, classifies the strains according to 
# what unique variant they have (for each gene), assigning one number to the unique variant
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

# Given the output from past function and the metabolite abundance matrix, divides that last one according to 
# the unique variant the strains have for each gene.
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

#G iven the output from past function, does a Kruskal-Wallis test between the vectors corresponding to each 
# metabolite for each group of strains that has a common variant, for each gene.
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
# Filter the output of the past function to a given alpha value (0.05 is the default).
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

# Given the significative genes from past function and the metabolite names of metabolite abundance matrix, outputs all 
# the genes that are found to be correlated with each metabolite. Optionally the output can be filtered to show only 
# the ones that are reported to be differential. Also can be choosen to generate csv files.
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