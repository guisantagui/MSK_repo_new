setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis")
# We are going to align separatedly gene ilvB because it could be causing differences in abundance of compounds related with swarming 
# that appear in Valine, Leucine and Isoleucine Biosynthesis. It's only present in 4 strains. In the MSA this gene doesn't appear, 
# and is not because the names were not assigned correctly, so we're gonna align this gene separatedly and "manually" using the 
# coordinates that appear in the GFF file.
if(!require(DECIPHER)) install.packages("DECIPHER")
library(DECIPHER)
if(!require(seqinr)) install.packages("seqinr")
library(seqinr)

allSeqs <- list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive")
load("geneCallsAllStrains4HeronProkka.RData")
F22031 <- readDNAStringSet(filepath = paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive",
                                            allSeqs[1], sep = "/"))
F34365 <- readDNAStringSet(filepath = paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive",
                                            allSeqs[4], sep = "/"))
T6313 <- readDNAStringSet(filepath = paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive",
                                           allSeqs[19], sep = "/"))
W16407 <- readDNAStringSet(filepath = paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive",
                                            allSeqs[22], sep = "/"))

ilvB_F22031Dat <- geneCallsAllStrains4HeronProkka$F22031[which(geneCallsAllStrains4HeronProkka$F22031$Annotation == "ilvB"), ]
ilvB_F34365Dat <- geneCallsAllStrains4HeronProkka$F34365[which(geneCallsAllStrains4HeronProkka$F34365$Annotation == "ilvB"), ]
ilvB_T6313Dat <- geneCallsAllStrains4HeronProkka$T6313[which(geneCallsAllStrains4HeronProkka$T6313$Annotation == "ilvB"), ]
ilvB_W16407Dat <- geneCallsAllStrains4HeronProkka$W16407[which(geneCallsAllStrains4HeronProkka$W16407$Annotation == "ilvB"), ]

ilvB_F22031 <- as.character(F22031$`12 F22031 | closed genome `[ilvB_F22031Dat$Start:ilvB_F22031Dat$Stop])
ilvB_F34365 <- as.character(F34365$`14 F34365 | 4 contig(s) `[ilvB_F34365Dat$Start:ilvB_F34365Dat$Stop])
ilvB_T6313 <- as.character(T6313$`41 T6313 | 5 contig(s) `[ilvB_T6313Dat$Start:ilvB_T6313Dat$Stop])
ilvB_W16407 <- as.character(W16407$`43 W16407 | closed genome `[ilvB_W16407Dat$Start:ilvB_W16407Dat$Stop])

set <- c(ilvB_F22031, ilvB_F34365, ilvB_T6313, ilvB_W16407)
ilvB_set <- DNAStringSet(set)
ilvB_AAAlign <- AlignTranslation(ilvB_set, type = "AAStringSet")
ilvB_DNAAlign <- AlignTranslation(ilvB_set, type = "DNAStringSet")
BrowseSeqs(ilvB_AAAlign)
BrowseSeqs(ilvB_DNAAlign)

BrowseSeqs(ParsedSubGraphs_allStrains_named_new$ilvE$AAAlignment) # Same in all strains
BrowseSeqs(ParsedSubGraphs_allStrains_named_new$`ilvG /*/ ilvG_2`$AAAlignment) # Shorter and different in T63266. Many point mutations in different genes. 


# Formym-Methionine is the compound that is most correlated with swarming. fmt is the gene that codes for the enzyme that formylates the methionil-tRNA
fmtAA <- ParsedSubGraphs_allStrains_named_new$fmt$AAAlignment
names(fmtAA) <- gsub(".fasta",  replacement = "", allSeqs)

BrowseSeqs(fmtAA)

write.fasta(ilvB_F22031, names = "ilvB", file.out = "ilvB_F22031.fasta")

