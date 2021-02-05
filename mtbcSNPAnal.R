setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcSNPAnal")

SNPdir <- "C:/Users/Guillem/Documents/PhD/comput/data/mtbcSNPs/"

btwSNPs <- read.table(paste(SNPdir, "btwlin_positions.annot", sep = ""), sep = "\t", header = T)
SNPsA1 <- read.table(paste(SNPdir, "A1_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsA2 <- read.table(paste(SNPdir, "A2_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsA3 <- read.table(paste(SNPdir, "A3_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsA4 <- read.table(paste(SNPdir, "A4_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL1 <- read.table(paste(SNPdir, "L1_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL2 <- read.table(paste(SNPdir, "L2_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL3 <- read.table(paste(SNPdir, "L3_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL4 <- read.table(paste(SNPdir, "L4_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL5 <- read.table(paste(SNPdir, "L5_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL6 <- read.table(paste(SNPdir, "L6_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL7 <- read.table(paste(SNPdir, "L7_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL8 <- read.table(paste(SNPdir, "L8_intralineage_DR.annot", sep = ""), sep = "\t", header = T)
SNPsL9 <- read.table(paste(SNPdir, "L9_intralineage_DR.annot", sep = ""), sep = "\t", header = T)

SNPsA1$annotation <- gsub(" ", "", SNPsA1$annotation)
SNPsA2$annotation <- gsub(" ", "", SNPsA2$annotation)
SNPsA3$annotation <- gsub(" ", "", SNPsA3$annotation)
SNPsA4$annotation <- gsub(" ", "", SNPsA4$annotation)
SNPsL1$annotation <- gsub(" ", "", SNPsL1$annotation)
SNPsL2$annotation <- gsub(" ", "", SNPsL2$annotation)
SNPsL3$annotation <- gsub(" ", "", SNPsL3$annotation)
SNPsL4$annotation <- gsub(" ", "", SNPsL4$annotation)
SNPsL5$annotation <- gsub(" ", "", SNPsL5$annotation)
SNPsL6$annotation <- gsub(" ", "", SNPsL6$annotation)
SNPsL7$annotation <- gsub(" ", "", SNPsL7$annotation)
SNPsL8$annotation <- gsub(" ", "", SNPsL8$annotation)
SNPsL9$annotation <- gsub(" ", "", SNPsL9$annotation)

SNPsA1$Gene <- gsub(" ", "", SNPsA1$Gene)
SNPsA2$Gene <- gsub(" ", "", SNPsA2$Gene)
SNPsA3$Gene <- gsub(" ", "", SNPsA3$Gene)
SNPsA4$Gene <- gsub(" ", "", SNPsA4$Gene)
SNPsL1$Gene <- gsub(" ", "", SNPsL1$Gene)
SNPsL2$Gene <- gsub(" ", "", SNPsL2$Gene)
SNPsL3$Gene <- gsub(" ", "", SNPsL3$Gene)
SNPsL4$Gene <- gsub(" ", "", SNPsL4$Gene)
SNPsL5$Gene <- gsub(" ", "", SNPsL5$Gene)
SNPsL6$Gene <- gsub(" ", "", SNPsL6$Gene)
SNPsL7$Gene <- gsub(" ", "", SNPsL7$Gene)
SNPsL8$Gene <- gsub(" ", "", SNPsL8$Gene)
SNPsL9$Gene <- gsub(" ", "", SNPsL9$Gene)

SNPsA1$annotation[grep("stopgain", SNPsA1$annotation)] <- "stop_gained"
SNPsA2$annotation[grep("stopgain", SNPsA2$annotation)] <- "stop_gained"
SNPsA3$annotation[grep("stopgain", SNPsA3$annotation)] <- "stop_gained"
SNPsA4$annotation[grep("stopgain", SNPsA4$annotation)] <- "stop_gained"
SNPsL1$annotation[grep("stopgain", SNPsL1$annotation)] <- "stop_gained"
SNPsL2$annotation[grep("stopgain", SNPsL2$annotation)] <- "stop_gained"
SNPsL3$annotation[grep("stopgain", SNPsL3$annotation)] <- "stop_gained"
SNPsL4$annotation[grep("stopgain", SNPsL4$annotation)] <- "stop_gained"
SNPsL5$annotation[grep("stopgain", SNPsL5$annotation)] <- "stop_gained"
SNPsL6$annotation[grep("stopgain", SNPsL6$annotation)] <- "stop_gained"
SNPsL7$annotation[grep("stopgain", SNPsL7$annotation)] <- "stop_gained"
SNPsL8$annotation[grep("stopgain", SNPsL8$annotation)] <- "stop_gained"
SNPsL9$annotation[grep("stopgain", SNPsL9$annotation)] <- "stop_gained"

btwSNPs[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c", btwSNPs$GENE), ]

SNPsIntLst <- list(A1 = SNPsA1, 
                   A2 = SNPsA2,
                   A3 = SNPsA3, 
                   A4 = SNPsA4,
                   L1 = SNPsL1,
                   L2 = SNPsL2,
                   L3 = SNPsL3,
                   L4 = SNPsL4,
                   L5 = SNPsL5,
                   L6 = SNPsL6,
                   L7 = SNPsL7,
                   L8 = SNPsL8,
                   L9 = SNPsL9)
                   

plcSNPsLst <- lapply(SNPsIntLst, function(x) x[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c", x$Gene), ])


humanPlcSNPs <- plcSNPsLst[grep("L", names(plcSNPsLst))]

intraSNPAnnotProp <- lapply(humanPlcSNPs, function(x) table(x$annotation)/length(x$annotation))

unique(unlist(lapply(humanPlcSNPs, function(x) x$annotation)))


intraSNPAnnotPropDF <- data.frame("missense_variant" = rep(0, 9),
                                  "synonymous_variant" = rep(0, 9),
                                  "stop_gained" = rep(0, 9),
                                  "synonymous_variant,synonymous_variant" = rep(0, 9),
                                  "stop_lost&splice_region_variant" = rep(0, 9),
                                  "intergenic_region" = rep(0, 9),
                                  "synonymous_variant,missense_variant" = rep(0, 9),
                                  "missense_variant,synonymous_variant" = rep(0, 9),
                                  "missense_variant,missense_variant" = rep(0, 9),
                                  "splice_region_variant&stop_retained_variant" = rep(0, 9),
                                  "missense_variant,missense_variant,synonymous_variant" = rep(0, 9), 
                                  row.names = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9"))

for(i in 1:length(intraSNPAnnotProp)){
        colPos <- match(make.names(names(intraSNPAnnotProp[[i]])), colnames(intraSNPAnnotPropDF))
        intraSNPAnnotPropDF[i, colPos] <- intraSNPAnnotProp[[i]]
}

write.csv(intraSNPAnnotPropDF, file = "intraSNPAnnotPropDF.csv")

lapply(humanPlcSNPs, function(x) x$Gene[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c", names(table(x$Gene)))])

for(i in 1:length(humanPlcSNPs)){
        print(names(humanPlcSNPs)[i])
        print(table(humanPlcSNPs[[i]]$Gene)[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c", names(table(humanPlcSNPs[[i]]$Gene)))]/sum(table(humanPlcSNPs[[i]]$Gene)))
}


intraSNPAnnotPropDFPerGene <- data.frame("gene" = rep(c("Rv1755c", "plcA_Rv2351c", "plcB_Rv2350c", "plcC_Rv2349c", "IG2388_Rv2349c-Rv2350c"), 9),
                                         "missense_variant" = rep(0, 9*5),
                                         "synonymous_variant" = rep(0, 9*5),
                                         "stop_gained" = rep(0, 9),
                                         "synonymous_variant,synonymous_variant" = rep(0, 9*5),
                                         "stop_lost&splice_region_variant" = rep(0, 9*5),
                                         "intergenic_region" = rep(0, 9*5),
                                         "synonymous_variant,missense_variant" = rep(0, 9*5),
                                         "missense_variant,synonymous_variant" = rep(0, 9*5),
                                         "missense_variant,missense_variant" = rep(0, 9*5),
                                         "splice_region_variant&stop_retained_variant" = rep(0, 9*5),
                                         "missense_variant,missense_variant,synonymous_variant" = rep(0, 9*5), 
                                         row.names = make.unique(c(rep("L1", 5), 
                                                                   rep("L2", 5), 
                                                                   rep("L3", 5), 
                                                                   rep("L4", 5), 
                                                                   rep("L5", 5), 
                                                                   rep("L6", 5), 
                                                                   rep("L7", 5), 
                                                                   rep("L8", 5), 
                                                                   rep("L9", 5))))


hLin <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9")
for(i in 1:9){
        print(i)
        tab <- table(humanPlcSNPs[[i]]$Gene, humanPlcSNPs[[i]]$annotation)[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c|IG2388_Rv2349c-Rv2350c", rownames(table(humanPlcSNPs[[i]]$Gene, humanPlcSNPs[[i]]$annotation))), ]
        print(tab)
        colPos <- match(make.names(colnames(tab)), colnames(intraSNPAnnotPropDFPerGene))
        rowPos <- which(as.character(intraSNPAnnotPropDFPerGene$gene) %in% rownames(tab) & gsub("\\..*", "",  rownames(intraSNPAnnotPropDFPerGene)) == hLin[i])
        print(rowPos)
        intraSNPAnnotPropDFPerGene[rowPos, colPos] <- tab
        #rowPos <- 
        #intraSNPAnnotPropDF[i, colPos] <- intraSNPAnnotProp[[i]]
}

# Filter to keep just the genes and lineages with SNPs
intraSNPAnnotPropDFPerGene <- intraSNPAnnotPropDFPerGene[apply(intraSNPAnnotPropDFPerGene[, 2:ncol(intraSNPAnnotPropDFPerGene)], 1, sum) != 0 , ]

write.csv(intraSNPAnnotPropDFPerGene, file = "intraSNPAnnotPropDFPerGene.csv")

table(humanPlcSNPs$L1$Gene, humanPlcSNPs$L1$annotation)[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c|IG2388_Rv2349c-Rv2350c", rownames(table(humanPlcSNPs$L1$Gene, humanPlcSNPs$L1$annotation))), ]

which(as.character(intraSNPAnnotPropDFPerGene$gene) %in% rownames(tab) & gsub("\\..*", "",  rownames(intraSNPAnnotPropDFPerGene)) == "L1")

table(humanPlcSNPs$L1$Gene)[grep("Rv2349c|Rv2350c|Rv2351c|Rv1755c|IG2388_Rv2349c-Rv2350c", names(table(humanPlcSNPs$L1$Gene)))]/sum(table(humanPlcSNPs$L1$Gene))
