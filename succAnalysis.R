setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/succAnalysis")

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/ParsedSubGraphs_allStrains_named_new.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/prokkaAnnot.RData")




source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(DECIPHER)) install.packages("DECIPHER")
library(DECIPHER)
if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)

prokkAnot <- read.csv("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/untitled_folder/gene_presence_absence.csv")
orderAlignSeqs <- gsub(".fasta", "", list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive"))

GFFfileNames <- list.files("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3")

allGFFs <- list()
for(i in seq_along(GFFfileNames)){
        gff <- readGFF(paste("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3", GFFfileNames[i], sep = "/"))
        allGFFs[[i]] <- gff
}
names(allGFFs) <- make.names(gsub(".gff", replacement = "", GFFfileNames))
 
# Succinate is negatively correlated with swarming according OPLS-DA. We did a swarming assay adding succinate 
# to the plates using PA14, F30658 and F34365 strains, and observed that in PA14 and F34365 strains succinate 
# addition increased swarming (their swarming is the highest), while in F30658 swarming decreased (this strain 
# desn't swarm the same than the other two). Other studies reported that mutations in sdh genes (it has 
# four genes, one for each subunit of the enzyme) (succinate DH, EC 1.3.5.1) decreased swarming. This enzyme 
# catalyzes the reaction of succinate--> Fumarate. Fumarate is positively correlated with swarming (swarmer 
# strains have a high abundance of intracellular fumarate). When aligning sdh genes there are no differences
# in terms of AA sequence. Maybe the nt differences could cause structural differences in the mRNA that 
# might produce different transcription rates. These compounds have important roles in the TCA cycle, 
# so we're gonna check the levels of the compounds that have a role here and the sequence of the enzymes.

TCACpds <- gsub("cpd:", "", keggLink("compound", "map00020"))

dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
colnames(ccmn_norm_mets_good_old) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]

# When checking compounds with same molecular weights than the ambiguous ones and comparing the results with
# annotation data we found that fumarate should't be ambiguous. Let's untag it. 

# Phosphoenolpyruvate also is tagged as ambiguous, but is the only compound with that mol weight that 
# we have genetic evidence that could be in the cell. 
#colnames(ccmn_norm_mets_good_old) <- gsub("fumarate_?", "fumarate", colnames(ccmn_norm_mets_good_old), fixed = T)
#colnames(ccmn_norm_mets_good_old) <- gsub("phophoenolpyruvate_?", "phophoenolpyruvate", colnames(ccmn_norm_mets_good_old), fixed = T)


TCACpdMetTab <- ccmn_norm_mets_good_old[, match(dictionary$Consensus[dictionary$`KEGG IDs` %in% TCACpds], 
                                                colnames(ccmn_norm_mets_good_old))[!is.na(match(dictionary$Consensus[dictionary$`KEGG IDs` %in% TCACpds], 
                                                                                                colnames(ccmn_norm_mets_good_old)))]]

groups <- unique(gsub("_.*", replacement = "", rownames(ccmn_norm_mets_good_old)))
cols <- topo.colors(length(groups))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(ccmn_norm_mets_good_old))))
for(ii in 1:nrow(ccmn_norm_mets_good_old)) {
        selected <- which(groups == gsub("\\_.*", "", rownames(ccmn_norm_mets_good_old))[ii])
        Cols[ii] <- cols[selected]
}
#ccmnNormMets <- ccmnNormMets[-25, ]
#Cols <- Cols[-25]

heatMapNoRowNorm <- heatmap.2(as.matrix(t(TCACpdMetTab)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                              density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                              col = redgreen(75), 
                              #breaks = 76, 
                              ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
                              ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
                              cex.main = 20,
                              keysize = 0.7,
                              cexRow = 0.7,
                              cexCol = 1.2,
                              scale = "row",
                              #colCol = colCols,
                              cellnote = round(as.matrix(t(TCACpdMetTab)), 2),
                              notecex = 0.7,
                              key.xtickfun=function() {
                                      cex <- par("cex")*par("cex.axis")
                                      side <- 1
                                      line <- 0
                                      col <- par("col.axis")
                                      font <- par("font.axis")
                                      mtext("low", side=side, at=0, adj=0,
                                            line=line, cex=cex, col=col, font=font)
                                      mtext("high", side=side, at=1, adj=1,
                                            line=line, cex=cex, col=col, font=font)
                                      return(list(labels=FALSE, tick=FALSE))
                              })


TCACpdMetTab_meds <- getStrainMedian(TCACpdMetTab)

heatMapNoRowNorm <- heatmap.2(as.matrix(t(TCACpdMetTab_meds)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                              density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                              col = redgreen(75), 
                              #breaks = 76, 
                              ColSideColors = unique(Cols), notecol = NULL, trace = "none", xlab = "Strains", 
                              ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
                              cex.main = 20,
                              keysize = 0.7,
                              cexRow = 0.7,
                              cexCol = 1.2,
                              scale = "row",
                              #colCol = colCols,
                              cellnote = round(as.matrix(t(TCACpdMetTab)), 2),
                              notecex = 0.7,
                              key.xtickfun=function() {
                                      cex <- par("cex")*par("cex.axis")
                                      side <- 1
                                      line <- 0
                                      col <- par("col.axis")
                                      font <- par("font.axis")
                                      mtext("low", side=side, at=0, adj=0,
                                            line=line, cex=cex, col=col, font=font)
                                      mtext("high", side=side, at=1, adj=1,
                                            line=line, cex=cex, col=col, font=font)
                                      return(list(labels=FALSE, tick=FALSE))
                              })


dictEnzymes[which(dictEnzymes$ECnums == "6.2.1.5"), ]
which(dictEnzymes$ECnums == "2.8.3.18")

dictEnzymes[grep("succi", dictEnzymes$Annotation),]

# scoA and scoB are the genes that code for ec 2.8.3.5, which turns succinyl-CoA into succinate. 
# None of them are in M1608 (number 11 of the seqs)
scoA_AA <- ParsedSubGraphs_allStrains_named_new$scoA$AAAlignment
names(scoA_AA) <- orderAlignSeqs[-11]
scoB_AA <- ParsedSubGraphs_allStrains_named_new$scoB$AAAlignment
names(scoB_AA) <- orderAlignSeqs[-11]

BrowseSeqs(scoA_AA)
BrowseSeqs(scoB_AA)

scoA_DNA <- ParsedSubGraphs_allStrains_named_new$scoA$DNAAlignment
names(scoA_DNA) <- orderAlignSeqs[-11]
scoB_DNA <- ParsedSubGraphs_allStrains_named_new$scoB$DNAAlignment
names(scoB_DNA) <- orderAlignSeqs[-11]

BrowseSeqs(scoA_DNA)
BrowseSeqs(scoB_DNA)

# Are the same in terms of AA seq, but not nt.

# Look all enzymes in TCA

# Load all genomes

allSeqs <- list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive")

allGenomes <- list()
for(i in seq_along(allSeqs)){
        fasta <- readDNAStringSet(filepath = paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive",
                                                   allSeqs[i], sep = "/"))
        allGenomes[[i]] <- fasta
}
names(allGenomes) <- gsub(".fasta|.fna", replacement = "", allSeqs)

TCAEnzs <- gsub("ec:", "", keggLink("enzyme", "map00020"))
TCAEnzs <- TCAEnzs[TCAEnzs %in% dictEnzymes$ECnums]

TCAEnzsTab <- dictEnzymes[dictEnzymes$ECnums %in% TCAEnzs, ]


TCAAlign <- list()
for(i in seq_along(TCAEnzsTab$Gene)){
        print(paste("i:", i))
        if(TCAEnzsTab$Gene[i] %in% names(ParsedSubGraphs_allStrains_named_new)){
                DNA <- ParsedSubGraphs_allStrains_named_new[names(ParsedSubGraphs_allStrains_named_new) == TCAEnzsTab$Gene[i]][[1]]$DNAAlignment
                names(DNA) <- orderAlignSeqs[as.numeric(gsub("\\_.*", "", names(DNA)))]
                AA <- ParsedSubGraphs_allStrains_named_new[names(ParsedSubGraphs_allStrains_named_new) == TCAEnzsTab$Gene[i]][[1]]$AAAlignment
                names(AA) <- orderAlignSeqs[as.numeric(gsub("\\_.*", "", names(AA)))]
        }else{
                strWGene <- colnames(TCAEnzsTab[, 5:ncol(TCAEnzsTab)-1])[nchar(TCAEnzsTab[i, 5:ncol(TCAEnzsTab)-1]) > 1]
                geneSet <- c()
                for(j in seq_along(strWGene)){
                        print(paste("j:",j))
                        seqDat <- allGFFs[[strWGene[j]]][which(allGFFs[[strWGene[j]]]$ID == TCAEnzsTab[i, strWGene[j]]), ][1, ]
                        gene <- as.character(allGenomes[[strWGene[j]]][[1]][seqDat$start:seqDat$end])
                        geneSet <- c(geneSet, gene)
                }
                geneSet <- DNAStringSet(geneSet)
                names(geneSet) <- strWGene
                if(length(geneSet) > 1){
                        DNA <- AlignTranslation(geneSet, type = "DNAStringSet")
                        AA <- AlignTranslation(geneSet, type = "AAStringSet")
                }else{
                        DNA <- geneSet
                        AA <- geneSet
                }
        }
        both <- list(DNA, AA)
        names(both) <- c("DNAAlignment", "AAAlignment")
        TCAAlign[[i]] <- both
}
names(TCAAlign) <- TCAEnzsTab$Gene
strWGene <- colnames(TCAEnzsTab[, 5:ncol(TCAEnzsTab)-1])[nchar(TCAEnzsTab[9, 5:ncol(TCAEnzsTab)-1]) > 1]
seqDat <- allGFFs[[strWGene[1]]][which(allGFFs[[strWGene[21]]]$ID == TCAEnzsTab[9, strWGene[21]]), ][1, ]

TCAEnzsTab[27, strWGene[1]] 
BrowseSeqs(TCAAlign[[1]]$AAAlignment[names(TCAAlign[[1]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.1.5.4. mqo_2 Mut in F23197, F9670, H47921, M1608, M37351, PA14, S869968, T38079, T6313 & T63266
BrowseSeqs(TCAAlign[[2]]$AAAlignment[names(TCAAlign[[2]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:2.3.1.12. aceF
BrowseSeqs(TCAAlign[[3]]$AAAlignment[names(TCAAlign[[3]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.8.1.4. lpdG 2 muts
BrowseSeqs(TCAAlign[[4]]$AAAlignment[names(TCAAlign[[4]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:2.3.1.12. acoC
BrowseSeqs(TCAAlign[[5]]$AAAlignment[names(TCAAlign[[5]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:6.2.1.5. sucC
BrowseSeqs(TCAAlign[[6]]$AAAlignment[names(TCAAlign[[6]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.1.1.42 icd_2 # Mut in W91453, W24637, M55212 & F30658
BrowseSeqs(TCAAlign[[7]]$AAAlignment[names(TCAAlign[[7]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.1.1.49. pckA
BrowseSeqs(TCAAlign[[8]]$AAAlignment[names(TCAAlign[[8]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.2.4.2. sucA No muts 
BrowseSeqs(TCAAlign[[9]]$AAAlignment[names(TCAAlign[[9]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx")]) # ec:6.2.1.5. sucD
BrowseSeqs(TCAAlign[[10]]$AAAlignment[names(TCAAlign[[10]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.2.1.3. acnB
BrowseSeqs(TCAAlign[[11]]$AAAlignment[names(TCAAlign[[11]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.2.1.3. acnA
BrowseSeqs(TCAAlign[[12]]$AAAlignment[names(TCAAlign[[12]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.3.5.1. sdhA No muts 
BrowseSeqs(TCAAlign[[13]]$AAAlignment[names(TCAAlign[[13]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec: 4.2.1.2. fumC_1
BrowseSeqs(TCAAlign[[14]]$AAAlignment[names(TCAAlign[[14]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:2.3.1.61. sucB
BrowseSeqs(TCAAlign[[15]]$AAAlignment[names(TCAAlign[[15]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.2.1.3. acnM 
BrowseSeqs(TCAAlign[[16]]$AAAlignment[names(TCAAlign[[16]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.2.1.2. fumC_2
BrowseSeqs(TCAAlign[[17]]$AAAlignment[names(TCAAlign[[17]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.1.1.42. icd_1
BrowseSeqs(TCAAlign[[18]]$AAAlignment[names(TCAAlign[[18]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.2.4.1. pdhA 3 Muts, many of them occur in paralel.
BrowseSeqs(TCAAlign[[19]]$AAAlignment[names(TCAAlign[[19]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:2.3.1.12. pdhC
BrowseSeqs(TCAAlign[[20]]$AAAlignment[names(TCAAlign[[20]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:4.2.1.2. fumA
BrowseSeqs(TCAAlign[[21]]$AAAlignment[names(TCAAlign[[21]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.1.5.4. mqo_1
BrowseSeqs(TCAAlign[[22]]$AAAlignment[names(TCAAlign[[22]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.3.5.1. sdhB
BrowseSeqs(TCAAlign[[23]]$AAAlignment[names(TCAAlign[[23]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.2.4.1. pdhB
BrowseSeqs(TCAAlign[[24]]$AAAlignment[names(TCAAlign[[24]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.8.1.4. lpd3
BrowseSeqs(TCAAlign[[25]]$AAAlignment[names(TCAAlign[[25]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.2.4.1. aceE
BrowseSeqs(TCAAlign[[26]]$AAAlignment[names(TCAAlign[[26]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.8.1.4. lpdV
BrowseSeqs(TCAAlign[[27]]$AAAlignment[names(TCAAlign[[27]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:1.1.5.4. mqo_3
BrowseSeqs(TCAAlign[[28]]$AAAlignment[names(TCAAlign[[28]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:6.2.1.5. sucD_1
BrowseSeqs(TCAAlign[[29]]$AAAlignment[names(TCAAlign[[29]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:6.2.1.5. sucD_2
BrowseSeqs(TCAAlign[[30]]$AAAlignment[names(TCAAlign[[30]]$AAAlignment) %in% c("F30658", "F34365", "PA14_jbx.fna")]) # ec:6.2.1.5. sucC_1

dictEnzymes$ECnums[match(names(TCAAlign), dictEnzymes$Gene)]

ParsedSubGraphs_allStrains_named_new$icd
BrowseSeqs(TCAAlign[[12]]$AAAlignment)

plot()