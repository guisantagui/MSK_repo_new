# This script compares the four annotation files we have, builds a version with 
# all the genes that appear in the four and outputs a tab separated txt file to 
# be used as input for the deletion detection script. 
if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)
if(!require(Vennerable)) install.packages("Vennerable")
library(Vennerable)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcGenesGFF2TXT")

grntAnnot <- as.data.frame(readGFF("../../data/MTBC_ancestralAnnot/genes.gff"))
ncbiAnnot <- as.data.frame(readGFF("../../data/MTBC_ancestralAnnot/h37RvNCBIAnnot.gff3"))
myc1Annot <- as.data.frame(readGFF("../../data/MTBC_ancestralAnnot/Mycobacterium_tuberculosis_H37Rv_gff_v1.gff"))
myc2Annot <- as.data.frame(readGFF("../../data/MTBC_ancestralAnnot/Mycobacterium_tuberculosis_H37Rv_gff_v2.gff"))

grntAnnot_gns <- grntAnnot$locus_tag[!is.na(grntAnnot$locus_tag)]
ncbiAnnot_gns <- unique(ncbiAnnot$locus_tag)[!is.na(unique(ncbiAnnot$locus_tag))]
myc1Annot_gns <- myc1Annot$Locus[grep("Rv", myc1Annot$Locus)]
myc2Annot_gns <- myc2Annot$Locus[grep("Rv", myc2Annot$Locus)]

annotGenLst <- list(grntxa = grntAnnot_gns,
                    ncbi = ncbiAnnot_gns,
                    mycbrwsrV1 = myc1Annot_gns,
                    mycbrwsrV2 = myc2Annot_gns)


vennDg <- Venn(annotGenLst)

plot(vennDg, doWeights = F)

ncbiAnnot_gns[!ncbiAnnot_gns %in% grntAnnot_gns]
grntAnnot_gns[!grntAnnot_gns %in% ncbiAnnot_gns]


grntAnnot$locus_tag[!is.na(grntAnnot$locus_tag)]

## The annotation in garnatxa has 92 genes that are not in the other annotations (corresponding to 
## noncoding genes), while the other 3 have 115 genes that are not in the garnatxa one. So let's add 
## this 115 genes to the output.

genes <- grntAnnot

genes4txt <- genes[!is.na(genes$locus_tag), c("start", "end", "strand", "locus_tag")]

ncbiAnnot <- ncbiAnnot[1:nrow(ncbiAnnot) %% 2 == 1 & !is.na(ncbiAnnot$locus_tag), ]

ncbiAnnot2Add <- ncbiAnnot[!ncbiAnnot$locus_tag %in% genes4txt$locus_tag, c("start", "end", "strand", "locus_tag")]

genes4txt <- rbind.data.frame(genes4txt, ncbiAnnot2Add)

genes4txt <- genes4txt[order(genes4txt$start), ]

genes4txt$strand <- as.factor(sapply(genes4txt$strand, function(x) if(x == "+") x <- "F" else x <- "R"))

write.table(genes4txt, file = "allGenesMTBC.txt", sep = "\t", quote = F, row.names = F, col.names = F)
