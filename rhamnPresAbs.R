if(!require(rtracklayer)) BiocManager::install("rtracklayer")
library(rtracklayer)
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
library(KEGGREST)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/genePresAbs")


load("gene_tab_filt.RData")
load("gene_enz_tab.RData")
load("gene_enz_tab_filt.RData")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("dictEnzymes.RData")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamn2Cats_signMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamn3Cats_signMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/swarm_signMets.RData")

source("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/MSK_repo_new/genePresAbs_functions.R")

rownames(ccmnNormMets) <- gsub("W70322", "W70332", rownames(ccmnNormMets))

ccmnNormMets <- ccmnNormMets[-grep("M6075", rownames(ccmnNormMets)), ]

rownames(ccmnNormMets)[1:6] <- paste(c(rep("PA14.jbx", 3),
                                       rep("UCBPP.PA14", 3)),
                                     as.character(1:6),
                                     sep = "_")

gene_tab_filt <- gene_tab_filt[rownames(gene_tab_filt) %in% unique(gsub("\\_.*", rownames(ccmnNormMets), replacement = "")), ]
gene_enz_tab <- gene_enz_tab[rownames(gene_enz_tab) %in% unique(gsub("\\_.*", rownames(ccmnNormMets), replacement = "")), ]
gene_enz_tab_filt <- gene_enz_tab_filt[rownames(gene_enz_tab_filt) %in% unique(gsub("\\_.*", rownames(ccmnNormMets), replacement = "")), ]

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat.RData")


rhamn2CatsSort <- rhamnSwarmMat$rhamn2cats[match(unique(gsub("\\_.*|(PA14).*", 
                                                             "\\1", 
                                                             rownames(gene_tab_filt))),  rhamnSwarmMat$strain)]+1

rhamn2CatsSort[is.na(rhamn2CatsSort)] <- 2

names(rhamn2CatsSort) <- rownames(gene_tab_filt)

rhamn2CatsSort <- as.factor(rhamn2CatsSort)

contTab_oldGood_filt <- getContTab(gene_tab_filt, rhamn2CatsSort)

# Only enzymatic genes (the ones with EC number)
contTab_oldGood_enz <- getContTab(gene_enz_tab, rhamn2CatsSort)

# Enzimatic genes without PA14 and that are not shared across all strains
contTab_oldGood_filtEnz <- getContTab(gene_enz_tab_filt, rhamn2CatsSort)



fishGeneFilt <- doFisher(contTab_oldGood_filt, rhamn2CatsSort)

fishGeneEnz <- doFisher(contTab_oldGood_enz, rhamn2CatsSort)

fishGeneEnzFilt <- doFisher(contTab_oldGood_filtEnz, rhamn2CatsSort)

load("mannWhitPerGene_oldGood_filt.RData")
load("mannWhitPerGene_oldGood_enz.RData")
load("mannWhitPerGene_oldGood_filtEnz.RData")


mannWhitPerGeneFilt_oldGood_filt <- filtMannWhitPerGene(mannWhitPerGene_oldGood_filt)
mannWhitPerGeneFilt_oldGood_enz <- filtMannWhitPerGene(mannWhitPerGene_oldGood_enz)
mannWhitPerGeneFilt_oldGood_filtEnz <- filtMannWhitPerGene(mannWhitPerGene_oldGood_filtEnz)

rhamn2Cats_signMets <- rhamn2Cats_signMets[- grep("_?", 
                                                  rhamn2Cats_signMets$metabolites,
                                                  fixed = T), ]

rhamn2CatsSignMets <- dictionary$definitiveNames[match(rhamn2Cats_signMets$metabolites, dictionary$Consensus)]
rhamn2CatsSignMets[is.na(rhamn2CatsSignMets)] <- as.character(rhamn2Cats_signMets$metabolites[is.na(rhamn2CatsSignMets)])

rhamn2CatsSignMetsKEGGIDs <- dictionary$`KEGG IDs`[match(rhamn2CatsSignMets, dictionary$definitiveNames)]

rhamn2CatsSignMetsKEGGIDs <- matrix(rhamn2CatsSignMetsKEGGIDs, ncol = 1, dimnames = list(rhamn2CatsSignMets, "1&2"))


diffGenesPerMetRhamn <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_filt,
                                              diffMetObjkt = rhamn2CatsSignMetsKEGGIDs)

diffGenesPerMetRhamn_enz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_enz,
                                                  diffMetObjkt = rhamn2CatsSignMetsKEGGIDs)

diffGenesPerMetRhamn_enzFilt <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_filtEnz,
                                                      diffMetObjkt = rhamn2CatsSignMetsKEGGIDs)

geneMatsRhamn2Cats <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerMetRhamn,
                                        genePresAbsObjkt = gene_tab_filt,
                                        metClustObjkt = rhamn2CatsSort)

geneMatsRhamn2Cats_enz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerMetRhamn_enz,
                                            genePresAbsObjkt = gene_enz_tab,
                                            metClustObjkt = rhamn2CatsSort)

geneMatsRhamn2Cats_enzFilt <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerMetRhamn_enzFilt,
                                                genePresAbsObjkt = gene_enz_tab_filt,
                                                metClustObjkt = rhamn2CatsSort)

mannWhitTabsRhamn2Cats <- mannWhit_tabs(geneMatsRhamn2Cats)

mannWhitTabsRhamn2Cats_enz <- mannWhit_tabs(geneMatsRhamn2Cats_enz)

mannWhitTabsRhamn2Cats_enzFilt <- mannWhit_tabs(geneMatsRhamn2Cats_enzFilt)

# C00097
C00097relEnzs <- colnames(geneMatsRhamn2Cats_enzFilt$`1&2`$C00097)[-ncol(geneMatsRhamn2Cats_enzFilt$`1&2`$C00097)]

dictEnzymes[match(C00097relEnzs, dictEnzymes$Gene), ]
