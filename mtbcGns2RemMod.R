setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcGns2RemMod")

enzBinFishSign <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/enzDelsBinFishSignAvsL.csv")
enzBinChiSign <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/enzDelsBinChiSignAvsL.csv")
allPVallsAvsL <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/allPvalsEnzDelsAvsL.csv")
allDels <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/allDels.csv")
colnames(allDels)[1] <- "gene"
enzDels <- allDels[!is.na(allDels$EC_number), ]


#propMutAL <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/propMutAL.csv")
propMutAL <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/propMutAL_2.csv")


rownames(propMutAL) <- enzDels$gene
propMutAL$X <- enzDels$gene
colnames(propMutAL)[1] <- "gene"

rownames(allPVallsAvsL) <- enzDels$gene
allPVallsAvsL$X <- enzDels$gene
colnames(allPVallsAvsL)[1] <- "gene"

thrshld <- c(0.5, 0.35)
propMutChiSignAL <- propMutAL[((propMutAL$A > thrshld[1] & propMutAL$L < thrshld[2]) | (propMutAL$L > thrshld[1] & propMutAL$A < thrshld[2])) & propMutAL$gene %in% enzBinChiSign$gene, ]
propMutFishSignAL <- propMutAL[((propMutAL$A > thrshld[1] & propMutAL$L < thrshld[2]) | (propMutAL$L > thrshld[1] & propMutAL$A < thrshld[2])) & propMutAL$gene %in% enzBinFishSign$gene, ]

toRemOrt <- enzDels$orthology[match(propMutChiSignAL$gene, enzDels$gene)]
toRemECn <- enzDels$EC_number[match(propMutChiSignAL$gene, enzDels$gene)]

toDeleteFromMods <- data.frame(gene = propMutChiSignAL$gene, 
                               EC_number = toRemECn, 
                               orthology = toRemOrt, 
                               propA = propMutChiSignAL$A, 
                               propL = propMutChiSignAL$L)

toDeleteFromMods <- cbind.data.frame(toDeleteFromMods, 
                                     allPVallsAvsL[match(toDeleteFromMods$gene, allPVallsAvsL$gene), ])

write.csv(toDeleteFromMods, file = "toDeleteFromMods.csv")

# Generate files of genes to remove based on differences between animal lineages and L5, L6 and L9
propMutAL569 <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/propMutALpVals.csv")

colnames(propMutAL569)[1] <- "gene"
rownames(propMutAL569) <- propMutAL569$gene

propMutAL569 <- propMutAL569[!is.na(propMutAL569$chiSqPValAdj), ]
propMutAL569 <- propMutAL569[((propMutAL569$A > thrshld[1] & propMutAL569$L > thrshld[1]) | (propMutAL569$A > thrshld[1] & propMutAL569$L < thrshld[2]) | (propMutAL569$L > thrshld[1] & propMutAL569$A < thrshld[2])), ]

propMutAL569$EC_number <- enzDels$EC_number[match(propMutAL569$gene, enzDels$gene)]
propMutAL569$orthology <- enzDels$orthology[match(propMutAL569$gene, enzDels$gene)]
propMutAL569 <- propMutAL569[, c("gene", "EC_number", "orthology", "A", "L", "mWhitPVal", "mWhitPValAdj", "chiSqPVal", "chiSqPValAdj", "fishPVal", "fishPValAdj")]

toDeleteFromL569 <- propMutAL569[(propMutAL569$A > thrshld[1] & propMutAL569$L > thrshld[1]), ]
toDeleteFromA <- propMutAL569

write.csv(toDeleteFromL569, file = "toDeleteFromL569.csv")
write.csv(toDeleteFromA, file = "toDeleteFromA.csv")

# This function accepts a vector of genes and lineages and a threshold and outputs the proportion 
# of sequences in each lineage that have each gene deleted
getDelProp <- function(delObjkt, geneVec, linVec, thrhld = 5){
        #geneVec <- factor(geneVec)
        propVecList <- list()
        for(i in 1:length(linVec)){
                lin <- linVec[i]
                genePropVec <- c()
                for(j in 1:length(geneVec)){
                        gene <- geneVec[j]
                        geneProp <- sum(delObjkt[gene==delObjkt$gene, grep(lin, colnames(delObjkt))] > thrhld)/length(delObjkt[gene==delObjkt$gene, grep(lin, colnames(delObjkt))])
                        genePropVec <- c(genePropVec, geneProp)
                }
                names(genePropVec) <- geneVec
                propVecList[[i]] <-genePropVec
        }
        names(propVecList) <- linVec
        propVecDF <- data.frame(propVecList)
        rownames(propVecDF) <- geneVec
        return(propVecDF)
}


getDelProp(allDels, toDeleteFromA$gene, linVec = c("A1", "A2", "A3", "A4"), thrhld = 5)

getDelProp(allDels, c("Rv2349c", "Rv2350c", "Rv2351c", "Rv3802c"), linVec = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "A1", "A2", "A3", "A4"), thrhld = 99.9)
getDelProp(allDels, c("Rv2349c", "Rv2350c", "Rv2351c", "Rv3802c"), linVec = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "A1", "A2", "A3", "A4"), thrhld = 2)

"Rv3802c" %in% enzDels$gene

#getDelProp(enzDels, c("Rv2349c", "Rv2350c", "Rv2351c", "Rv1755c"), linVec = c("A1", "A2", "A3", "A4"), thrhld = 5)