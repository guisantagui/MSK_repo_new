# This script takes as input the GFF file of the annotation of the ancestral MTBC strain and outputs a 
# tab separated txt file to be used as input for the deletion detection script. 

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcGenesGFF2TXT")

if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)

genes <- as.data.frame(readGFF("../../data/MTBC_ancestralAnnot/genes.gff"))

genes4txt <- genes[!is.na(genes$locus_tag), c("start", "end", "strand", "locus_tag")]

genes4txt$strand <- as.factor(sapply(genes4txt$strand, function(x) if(x == "+") x <- "F" else x <- "R"))

write.table(genes4txt, file = "allGenesMTBC.txt", sep = "\t", quote = F, row.names = F, col.names = F)