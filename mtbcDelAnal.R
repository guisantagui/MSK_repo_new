if(!require(UniProt.ws)) BiocManager::install("UniProt.ws")
library(UniProt.ws)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)
if(!require(ropls)) install.packages("ropls")
library(ropls)
if(!require(vcfR)) install.packages("vcfR")
library(vcfR)
if(!require(ape)) install.packages("ape")
library(ape)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal")

fileNames <- list.files("C:/Users/Guillem/Documents/PhD/comput/data/mtbcLineagesDels")
H37Rv_annot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/genes.gff")

allDels <- data.frame()
for(i in 1:length(fileNames)){
        delsLin <- read.csv(paste("C:/Users/Guillem/Documents/PhD/comput/data/mtbcLineagesDels", fileNames[i], sep = "/"))
        rownames(delsLin) <- delsLin[, 1]
        delsLin <- delsLin[, 2:ncol(delsLin)]
        lin <- gsub("nr", "", strsplit(fileNames[i], "_|\\.|-")[[1]][2])
        colnames(delsLin) <- paste(colnames(delsLin), lin, sep = "_")
        if(nrow(allDels) == 0){
                allDels <- delsLin
        }else{
                allDels <- cbind.data.frame(allDels, delsLin)
        }
}

getGeneDBIDs <- function(RvCodes){
        ort <- c()
        ncbi_gID <- c()
        ncbi_pID <- c()
        uniProt <- c()
        for(i in 1:length(RvCodes)){
                kGet <- keggGet(paste("mtu", RvCodes[i], sep = ":"))[[1]]
                a <- kGet$ORTHOLOGY[1]
                dbs <- kGet$DBLINKS
                if(!is.null(a)){
                        ort <- c(ort, a)
                }else{
                        ort <- c(ort, NA)
                }
                gID <- dbs[grep("GeneID", dbs)]
                pID <- dbs[grep("ProteinID", dbs)]
                up <- dbs[grep("UniProt", dbs)]
                if(length(gID) != 0){
                        ncbi_gID <- c(ncbi_gID, gsub("NCBI-GeneID: ", "", gID))
                }else{
                        ncbi_gID <- c(ncbi_gID, NA)
                }
                if(length(pID) != 0){
                        ncbi_pID <- c(ncbi_pID, gsub("NCBI-ProteinID: ", "", pID))
                }else{
                        ncbi_pID <- c(ncbi_pID, NA)
                }
                if(length(up) != 0){
                        uniProt <- c(uniProt, gsub("UniProt: ", "", up))
                }else{
                        uniProt <- c(uniProt, NA)
                }
        }
        annotDF <- data.frame(ncbi_geneID = ncbi_gID,
                              ncbi_proteinID = ncbi_pID,
                              uniProt = uniProt,
                              orthology = ort)
        
        rownames(annotDF) <- RvCodes
        return(annotDF)
}

annotDF <- getGeneDBIDs(rownames(allDels))

annotDF$product <- H37Rv_annot$product[match(rownames(annotDF), H37Rv_annot$locus_tag) + 1]

ECnum <- sub("\\].*", "", sub(".*\\[", "", as.character(annotDF$orthology))) 

ECnum[!grepl("EC:", ECnum)] <- NA

#strsplit(ECnum, " ")

annotDF$EC_number <- ECnum 

allDels <- cbind.data.frame(annotDF, allDels)

save(allDels, file = "allDels.RData")
load("allDels.RData")
write.csv(allDels, file = "allDels.csv")

lineages <- sapply(colnames(allDels)[7:ncol(allDels)], function(x) strsplit(x, "_")[[1]][2])
linGroup <- gsub("[[:digit:]]+", "", lineages)
linCols <- c("#c4bf62",
             "#87e0a0",
             "#6db3d6",
             "#f279ce",
             "#ff30c1",
             "#001aff",
             "#8826b5",
             "#ff0000",
             "#871414",
             "#24ad37",
             "#fbff00",
             "#ff9d00",
             "#37ff30")
names(linCols) <- unique(lineages)

allDelsPCA <- prcomp(t(allDels[, 7:ncol(allDels)]))
save(allDelsPCA, file = "mtbcAllDelsPCA.RData")
load("mtbcAllDelsPCA.RData")

pdf(file = "allDelsPCA.pdf")
fviz_pca_ind(allDelsPCA, 
             col.ind = lineages,
             habillage = lineages,
             geom = "point") + 
        scale_color_manual(name = "Lineages", 
                           labels = lineages,
                           values = linCols) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           labels = sapply(linGroup, function(x) if(x == "A") x <- "Animal" else x <- "Human"))
dev.off()

allDelsDist <- dist(t(allDels[, 7:ncol(allDels)]), method = "euclidean")
allDelsHCA <- hclust(allDelsDist, method = "ward.D")
save(allDelsHCA, file = "allDelsHCA.RData")
load("allDelsHCA.RData")

enzDels <- allDels[!is.na(allDels$EC_number), ]

enzDelsPCA <- prcomp(t(enzDels[, 7:ncol(enzDels)]))
save(enzDelsPCA, file = "enzDelsPCA.RData")
load("enzDelsPCA.RData")

pdf(file = "enzDelsPCA.pdf")
fviz_pca_ind(enzDelsPCA, 
             col.ind = lineages,
             habillage = lineages,
             geom = "point") + 
        scale_color_manual(name = "Lineages", 
                           labels = lineages,
                           values = linCols) +
        scale_shape_manual(name = "Lineages", 
                           values = c(rep(2, 4),
                                      rep(19, 9)),
                           labels = sapply(linGroup, function(x) if(x == "A") x <- "Animal" else x <- "Human"))
dev.off()

enzDelsDist <- dist(t(enzDels[, 7:ncol(allDels)]), method = "euclidean")
enzDelsHCA <- hclust(enzDelsDist, method = "ward.D")
save(enzDelsHCA, file = "enzDelsHCA.RData")
load("enzDelsHCA.RData")

enzDelsDend <- as.dendrogram(enzDelsHCA)
enzDelsGGdend <- as.ggdend(enzDelsDend)

# See if samples of same lineage are together in tree

dendCut <- sort(cutree(enzDelsDend, k = 13))

uniqueLinsInClust <- list()
for(i in 1:13){
        sampsInLin <- names(dendCut[dendCut == i])
        lin <- unique(sapply(sampsInLin, function(x) strsplit(x, "_")[[1]][2]))
        uniqueLinsInClust[[i]] <- lin
}
uniqueLinsInClust


ggplot(enzDelsGGdend, labels = FALSE) + 
        scale_y_reverse(expand = c(0.2, 0)) +
        coord_polar(theta="x")

dendCols <- linCols[match(sapply(labels(enzDelsDend), function(x) strsplit(x, "_")[[1]][2]), names(linCols))]
dendCols <- linCols[match(sapply(colnames(allDels[, 7:ncol(allDels)]), function(x) strsplit(x, "_")[[1]][2]), names(linCols))]

pdf("enzDend.pdf")
plot(as.phylo(enzDelsHCA), 
     type = "fan",
     tip.color = dendCols,
     cex = 0.05,
     label.offset = 1)
#legend(1, 
#       95, 
#       legend = unique(lineages),
#       col = linCols)
dev.off()

pdf("enzDelsOPLSDA.pdf")
enzDelsOPLSDA <- opls(t(enzDels[, 7:ncol(enzDels)]), 
                      linGroup,
                      predI = 1, 
                      orthoI = 3)
dev.off()

enzDelsOPLSDA_loads <- getLoadingMN(enzDelsOPLSDA,)

enzDelsOPLSDA_loads <- cbind.data.frame(enzDels[, 1:6],
                                        enzDelsOPLSDA_loads)

colnames(enzDelsOPLSDA_loads)[ncol(enzDelsOPLSDA_loads)] <- "loadings_p1"

enzDelsOPLSDA_loadsOrdered <- enzDelsOPLSDA_loads[order(enzDelsOPLSDA_loads$loadings_p1), ]

enzDelsPCAROPLS <- opls(t(enzDels[, 7:ncol(enzDels)]), 
                        typeVc = "x-score", 
                        parAsColFcVn = linGroup)
save(enzDelsPCAROPLS, file = "enzDelsPCAROPLS.RData")

enzDelsMannWhit <- apply(enzDels[, 7:ncol(enzDels)], 1, function(x) wilcox.test(x[linGroup == "A"], x[linGroup == "L"])$p.value)
enzDelsMannWhitAdjust <- p.adjust(enzDelsMannWhit, "BH")

enzDelsSignifPvals <- enzDelsMannWhitAdjust[enzDelsMannWhitAdjust < 0.05]

enzDelsSignif <- allDels[match(names(enzDelsSignifPvals), rownames(allDels)), ]
save(enzDelsSignif, file = "enzDelsSignif.RData")
write.csv(enzDelsSignif, file = "enzDelsSignif_AvsL.csv")

thrshld <- 10

enzDels[enzDels > thrshld]
class(enzDels$G01493_A1)

binarizeDels <- function(DM, thrshld){
        nums <- DM[, sapply(DM, class) == "numeric"]
        gInf <- DM[, sapply(DM, class) != "numeric"]
        bin <- nums
        bin[bin < thrshld] <- 0
        bin[bin >= thrshld] <- 1
        binOut <- cbind.data.frame(gInf, bin)
        return(binOut)
}

enzDelsBin <- binarizeDels(DM = enzDels, thrshld = thrshld)

chiPVals <- c()
fishPVals <- c()
enzDelsBinNumCls <- enzDelsBin[, sapply(enzDelsBin, class) == "numeric"]
for(i in 1:nrow(enzDelsBinNumCls)){
        cTab <- table(factor(enzDelsBinNumCls[i, ]),
                      factor(linGroup))
        chi <- chisq.test(cTab)
        fish <- fisher.test(cTab)
        chiPVal <- chi$p.value
        fishPVal <- fish$p.value
        chiPVals <- c(chiPVals, chiPVal)
        fishPVals <- c(fishPVals, fishPVal)
}

names(chiPVals) <- rownames(enzDelsBin)
chiPValsAdj <- p.adjust(chiPVals, "BH")
chiPValsAdjSign <- chiPValsAdj[chiPValsAdj < 0.05]

names(fishPVals) <- rownames(enzDelsBin)
fishPValsAdj <- p.adjust(fishPVals, "BH")
fishPValsAdjSign <- fishPValsAdj[fishPValsAdj < 0.05]

enzDelsBinFishSignAvsL <- enzDelsBin[match(names(fishPValsAdjSign), 
                                           rownames(enzDelsBin)), ]

write.csv(enzDelsBinFishSignAvsL, file = "enzDelsBinFishSignAvsL.csv")

ApropMut <- apply(enzDelsBinFishSignAvsL[, sapply(enzDelsBin, class) == "numeric"][, linGroup == "A"], 1, function(x) sum(x)/length(x))
LpropMut <- apply(enzDelsBinFishSignAvsL[, sapply(enzDelsBin, class) == "numeric"][, linGroup == "L"], 1, function(x) sum(x)/length(x))

propMutAL <- data.frame(A = ApropMut, L = LpropMut)

View(propMutAL)

enzDelsLoadsAvsLsignif <- enzDelsOPLSDA_loadsOrdered[rownames(enzDelsOPLSDA_loadsOrdered) %in% rownames(enzDelsSignif), ]

save(enzDelsLoadsAvsLsignif, file = "enzDelsLoadsAvsLsignif.RData")

write.csv(enzDelsLoadsAvsLsignif, file = "enzDelsLoadsAvsLsignif.csv")

View(enzDelsLoadsAvsLsignif)

View(enzDelsOPLSDA_loadsOrdered)
