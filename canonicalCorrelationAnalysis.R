setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/canonicalCorrelationAnalysis")

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(GGally)) install.packages("GGally")
if(!require(CCA)) install.packages("CCA")

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]

ccmnNormMets <- ccmn_norm_mets_good_old
colnames(ccmnNormMets) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
ccmnNormMets <- ccmnNormMets[, !is.na(colnames(ccmnNormMets))]

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis/pcomp1AreaPercCircularity.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis/swarmDatMeans.RData")

ggpairs(ccmnNormMets)

ggpairs(swarmDatMeans[, 2:ncol(swarmDatMeans)])

swarmDatMeans <- swarmDatMeans[match(gsub("\\_.*|(PA14).*", 
                                          rownames(ccmnNormMets), 
                                          rep = "\\1"), 
                                     swarmDatMeans$Strains), ]

matcor(ccmnNormMets, swarmDatMeans[, 2:ncol(swarmDatMeans)])

cc1 <- cc(ccmnNormMets, swarmDatMeans[, 2:ncol(swarmDatMeans)])

cc2 <- comput(ccmnNormMets, swarmDatMeans[, 2:ncol(swarmDatMeans)], cc1)

cc1$cor
cc1$xcoef
cc1$ycoef

cc2
