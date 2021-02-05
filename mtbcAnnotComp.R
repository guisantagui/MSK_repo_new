library(rtracklayer)
library(Vennerable)

grntAnnot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/genes.gff")
ncbiAnnot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/h37RvNCBIAnnot.gff3")
myc1Annot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/Mycobacterium_tuberculosis_H37Rv_gff_v1.gff")
myc2Annot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/Mycobacterium_tuberculosis_H37Rv_gff_v2.gff")

View(as.data.frame(ncbiAnnot))

class(ncbiAnnot)

geneList



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

View(as.data.frame(myc1Annot)[grep("MT", myc1Annot$Locus), ])

ncbiAnnot_gns[!ncbiAnnot_gns %in% grntAnnot_gns]
grntAnnot_gns[!grntAnnot_gns %in% ncbiAnnot_gns]


grntAnnot$locus_tag[!is.na(grntAnnot$locus_tag)]
