setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcGns2RemMod")

enzBinFishSign <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/enzDelsBinFishSignAvsL.csv")
enzBinChiSign <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/enzDelsBinChiSignAvsL.csv")
allPVallsAvsL <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/allPvalsEnzDelsAvsL.csv")
allDels <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/allDels.csv")
colnames(allDels)[1] <- "gene"
enzDels <- allDels[!is.na(allDels$EC_number), ]


propMutAL <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/propMutAL.csv")

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