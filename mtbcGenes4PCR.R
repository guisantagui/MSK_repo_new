setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcGenes4PCR")
annotDir <- "C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot"
annot1 <- as.data.frame(rtracklayer::readGFF(paste(annotDir, "genes.gff", sep = "/")))
annot2 <- as.data.frame(rtracklayer::readGFF(paste(annotDir, "h37RvNCBIAnnot.gff3", sep = "/")))
if(!require(Biostrings)) install.packages("Biostrings")
if(!require(rtracklayer)) install.packages("rtracklayer")
allNTSeqs <- Biostrings::readDNAStringSet("C:/Users/Guillem/Documents/PhD/comput/data/mtbcSeqs/Mycobacterium_tuberculosis_H37Rv_genes_v2.fasta")

genes4PCR <- c("Rv1755c",
               "Rv2349c",
               "Rv2350c",
               "Rv2351c",
               "Rv3802c",
               "Rv0217c",
               "Rv0220",
               "Rv0646c",
               "Rv1076",
               "Rv1104",
               "Rv1105",
               "Rv1399c",
               "Rv1400c",
               "Rv1426c",
               "Rv1497",
               "Rv1900c",
               "Rv1923",
               "Rv2045c",
               "Rv2284",
               "Rv2385",
               "Rv2463",
               "Rv2485c",
               "Rv2970c",
               "Rv3084",
               "Rv3176c",
               "Rv3203",
               "Rv3487c",
               "Rv3775")

plcRelatedGenes <- annot2[match(genes4PCR, annot2$locus_tag), c("seqid", "start", "end", "strand", "gene", "gene_biotype", "locus_tag")]

plcRelatedGenes$product <- unlist(annot3$Product[match(genes4PCR, annot3$Locus)])



plcRelatedGenes$ntSeq <- sapply(genes4PCR, function(x) as.character(allNTSeqs[[grep(x, names(allNTSeqs))[1]]]))

write.csv(plcRelatedGenes, file = "plcRelatedGenes.csv")