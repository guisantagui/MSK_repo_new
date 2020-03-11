setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/succAnalysis")

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis/ParsedSubGraphs_allStrains_named_new.RData")
orderAlignSeqs <- gsub(".fasta", "", list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive"))

source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

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
                              #cellnote = round(as.matrix(t(TCACpdMetTab)), 2),
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
TCAEnzs <- gsub("ec:", "", keggLink("enzyme", "map00020"))[TCAEnzs %in% dictEnzymes$ECnums]

TCAEnzsTab <- dictEnzymes[match(TCAEnzs, dictEnzymes$ECnums), ]

TCAAlign <- list()
for(i in seq_along(TCAEnzsTab$Gene)){
        DNA <- ParsedSubGraphs_allStrains_named_new[names(ParsedSubGraphs_allStrains_named_new) == TCAEnzsTab$Gene[i]][[1]]$DNAAlignment
        names(DNA) <- orderAlignSeqs[as.numeric(gsub("\\_.*", "", names(DNA)))]
        AA <- ParsedSubGraphs_allStrains_named_new[names(ParsedSubGraphs_allStrains_named_new) == TCAEnzsTab$Gene[i]][[1]]$AAAlignment
        names(AA) <- orderAlignSeqs[as.numeric(gsub("\\_.*", "", names(AA)))]
        both <- list(DNA, AA)
        names(both) <- c("DNAAlignment", "AAAlignment")
        TCAAlign[[i]] <- both
}
names(TCAAlign) <- TCAEnzsTab$Gene

BrowseSeqs(TCAAlign[[1]]$AAAlignment) # Mut in W91453, W24637, M55212 & F30658
BrowseSeqs(TCAAlign[[2]]$AAAlignment) # Mut in F23197, F9670, H47921, M1608, M37351, PA14, S869968, T38079, T6313 & T63266
BrowseSeqs(TCAAlign[[3]]$AAAlignment) # 3 Muts, many of them occur in paralel. 
BrowseSeqs(TCAAlign[[4]]$AAAlignment) # No muts
BrowseSeqs(TCAAlign[[5]]$AAAlignment) # No muts
BrowseSeqs(TCAAlign[[6]]$AAAlignment) # 2 Muts
BrowseSeqs(TCAAlign[[7]]$AAAlignment)
BrowseSeqs(TCAAlign[[8]]$AAAlignment)
BrowseSeqs(TCAAlign[[9]]$AAAlignment)
BrowseSeqs(TCAAlign[[10]]$AAAlignment)
BrowseSeqs(TCAAlign[[11]]$AAAlignment)
BrowseSeqs(TCAAlign[[12]]$AAAlignment)