#####################################################################################################################################
#######################################                                                   ###########################################
######################################           GENE PRESENCE/ABSENCE ANALYSIS            ##########################################
#######################################          Remodeled for new hybrid data            ###########################################
########################################                                                 ############################################
#####################################################################################################################################

setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs")

####################
#
# Load the data
#
####################################################################################################################################
PA14 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/PA14/PATRIC_genome_feature.csv")
F22031 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/F22031/PATRIC_genome_feature.csv")
F30658 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/F30658/PATRIC_genome_feature.csv")
F34365 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/F34365/PATRIC_genome_feature.csv")
F63912 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/F63912/PATRIC_genome_feature.csv")
F9670 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/F9670/PATRIC_genome_feature.csv")
H27930 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/H27930/PATRIC_genome_feature.csv")
H47921 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/H47921/PATRIC_genome_feature.csv")
H5708 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/H5708/PATRIC_genome_feature.csv")
M1608 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/M1608/PATRIC_genome_feature.csv")
M37351 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/M37351/PATRIC_genome_feature.csv")
M6075 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/M55212/PATRIC_genome_feature.csv")
M74707 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/M74707/PATRIC_genome_feature.csv")
S86968 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/S86968/PATRIC_genome_feature.csv")
T38079 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/T38079/PATRIC_genome_feature.csv")
T52373 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/T52373/PATRIC_genome_feature.csv")
T6313 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/T6313/PATRIC_genome_feature.csv")
T63266 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/T63266/PATRIC_genome_feature.csv")
W16407 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W16407/PATRIC_genome_feature.csv")
W25637 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W25637/PATRIC_genome_feature.csv")
W36662 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W36662/PATRIC_genome_feature.csv")
W45909 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W45909/PATRIC_genome_feature.csv")
W60856 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W60856/PATRIC_genome_feature.csv")
W70322 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W70332/PATRIC_genome_feature.csv")
W91453 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/W91453/PATRIC_genome_feature.csv")
X78812 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/X78812/PATRIC_genome_feature.csv")
X9820 <- read.csv("/Users/santamag/Desktop/GUILLEM/PA genomic info/PA annotated info/X9820/PATRIC_genome_feature.csv")

####################
#
# Create the table
#
####################################################################################################################################
strain_names<- c("PA14", "F22031", "F30658", "F34365", "F63912", "F9670", "H27930", "H47921", "H5708", "M1608", "M37351", "M6075",
                 "M74707", "S86968", "T38079", "T52373", "T6313", "T63266", "W16407", "W25637", "W36662", "W45909", "W60856", 
                 "W70322", "W91453", "X78812", "X9820")

all_genes <- unique(c(as.character(PA14[, 20]), as.character(F22031[, 20]), as.character(F30658[, 20]), as.character(F34365[, 20]), 
                      as.character(F63912[, 20]), as.character(F9670[, 20]), as.character(H27930[, 20]), as.character(H47921[, 20]),
                      as.character(H5708[, 20]), as.character(M1608[, 20]), as.character(M37351[, 20]), as.character(M6075[, 20]), 
                      as.character(M74707[, 20]), as.character(S86968[, 20]), as.character(T38079[, 20]), as.character(T52373[, 20]), 
                      as.character(T6313[, 20]), as.character(T63266[, 20]), as.character(W16407[, 20]), as.character(W25637[, 20]), 
                      as.character(W36662[, 20]), as.character(W45909[, 20]), as.character(W60856[, 20]), as.character(W70322[, 20]), 
                      as.character(W91453[, 20]), as.character(X78812[, 20]), as.character(X9820[, 20])))

gene_tab <- matrix(nrow = length(strain_names), ncol = length(all_genes), dimnames = list(strain_names, all_genes))

strain_genes <- list("PA14" = as.character(PA14[, 20]), "F22031" = as.character(F22031[, 20]), 
                     "F30658" = as.character(F30658[, 20]), "F34365" = as.character(F34365[, 20]), 
                     "F63912" = as.character(F63912[, 20]), "F9670" = as.character(F9670[, 20]), 
                     "H27930" = as.character(H27930[, 20]), "H47921" = as.character(H47921[, 20]),
                     "H5708" = as.character(H5708[, 20]), "M1608" = as.character(M1608[, 20]), 
                     "M37351" = as.character(M37351[, 20]), "M6075" = as.character(M6075[, 20]), 
                     "M74707" = as.character(M74707[, 20]), "S86968" = as.character(S86968[, 20]), 
                     "T38079" = as.character(T38079[, 20]), "T52373" = as.character(T52373[, 20]), 
                     "T6313" = as.character(T6313[, 20]), "T63266" = as.character(T63266[, 20]), 
                     "W16407" = as.character(W16407[, 20]), "W25637" = as.character(W25637[, 20]), 
                     "W36662" = as.character(W36662[, 20]), "W45909" = as.character(W45909[, 20]), 
                     "W60856" = as.character(W60856[, 20]), "W70322" = as.character(W70322[, 20]), 
                     "W91453" = as.character(W91453[, 20]), "X78812" = as.character(X78812[, 20]), 
                     "X9820" = as.character(X9820[, 20]))

for(i in 1:length(gene_tab[, 1])){
        for(j in 1:length(gene_tab[1, ])){
                if(colnames(gene_tab)[j] %in% strain_genes[[i]]){
                        gene_tab[i, j] <- 1
                }else{
                        gene_tab[i, j] <- 0
                }
        }
}
save(gene_tab, file = "gene_tab.RData")
load(file = "gene_tab.RData")

####################
#
# Plot it
#
####################################################################################################################################
library(gplots)
library(ggplot2)

#Make colorscale equivalent to the one that we've been using in previous heatmaps

cols <- topo.colors(length(strain_names) + 1)
cols <- cols[-1]
cols[1] <- '#000000'

tiff("gene_tab.tiff", width = 1000, height = 4000, units = "px", pointsize = 50)

heatmap.2(t(gene_tab), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "column", 
          col = redgreen(75), breaks = 76, ColSideColors = cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of genes", margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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

dev.off()

####################
#
# Obtain filtered versions of the set
#
####################################################################################################################################

# PA14 removed, and after that every gene that is present in all strains.
gene_tab_filt <- gene_tab[-1, ]


for(i in 1:dim(gene_tab)[2]){
        if(length(unique(gene_tab_filt[, i])) != 1){
                gene_tab_filt[, i] <- gene_tab_filt[, i] 
        }else{
                gene_tab_filt[, i] <- c(rep(NA, length(gene_tab_filt[, i])))
        }
}

if (!require(janitor)) install.packages('janitor')
library(janitor)

gene_tab_filt <- remove_empty(as.data.frame(gene_tab_filt), which = "cols")
save(gene_tab_filt, file = "gene_tab_filt.RData")

# Starting from all genes, keep only enzymatic ones (the ones with EC)
gene_enz_tab <- gene_tab[, sapply(colnames(gene_tab), grepl, pattern = "EC")]
save(gene_enz_tab, file = "gene_enz_tab.RData")

# Remove PA14 and the genes that are in all strains and keep only enzymatic ones.
gene_enz_tab_filt <- gene_tab_filt[, sapply(colnames(gene_tab_filt), grepl, pattern = "EC")]
gene_enz_tab_filt <- gene_enz_tab_filt[, colnames(gene_enz_tab_filt) != "Late competence protein ComEC, DNA transport"]
gene_enz_tab_filt <- gene_enz_tab_filt[, colnames(gene_enz_tab_filt) != "hypothetical protein APECO1_2271"]
gene_enz_tab_filt <- gene_enz_tab_filt[, colnames(gene_enz_tab_filt) != "RNA polymerase ECF-type sigma factor"]
save(gene_enz_tab_filt, file = "gene_enz_tab_filt.RData")

####################
#
# Obtain heatmaps previous sets
#
####################################################################################################################################
tiff("gene_tab_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(t(as.matrix(gene_tab_filt)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "column", 
          col = redgreen(75), breaks = 76, ColSideColors = cols[-1], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of genes (core filtered)", margins = c(6, 13), 
          keysize = 1,
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
dev.off()

tiff("gene_enz_tab.tiff", width = 5000, height = 25000, units = "px", pointsize = 50)
heatmap.2(t(gene_enz_tab), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "column", 
          col = redgreen(75), breaks = 76, ColSideColors = cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of enzymatic genes", margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          sepcolor = "black",
          colsep=1:ncol(t(gene_tab)),
          rowsep=1:nrow(t(gene_tab)),
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
dev.off()

tiff("gene_enz_tab_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(t(as.matrix(gene_enz_tab_filt)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "column", 
          col = redgreen(75), breaks = 76, ColSideColors = cols[-1], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of enzymatic genes (core filtered)", margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          sepcolor = "black",
          colsep=1:ncol(t(gene_tab)),
          rowsep=1:nrow(t(gene_tab)),
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
dev.off()

####################
#
# Load normalized metabolite datasets, diffMet objects and dictionary (to translate names to KEGG IDs).
#
####################################################################################################################################
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/newData/ccmn_norm_mets_newData.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/hybrid/ccmn_norm_mets_hybrid.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/oldDataGood/diffMets_oldGood.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/newData/diffMets_newData.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/diffMetAnal/hybrid/diffMets_hybrid.RData")
colnames(diffMets_oldGood) <- c("1&2", "1.1&1.2", "2.1&2.2")
colnames(diffMets_newData) <- c("1&2", "1.1&1.2", "2.1&2.2")
colnames(diffMets_hybrid) <- c("1&2", "1.1&1.2", "2.1&2.2")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")

metKEGGIDs_old <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
metKEGGIDs_new <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_newData), dictionary$`New Data Names`)]
metKEGGIDs_hybrid <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_hybrid), dictionary$Consensus)]

colnames(ccmn_norm_mets_good_old)[!is.na(metKEGGIDs_old)] <- metKEGGIDs_old[!is.na(metKEGGIDs_old)]
colnames(ccmn_norm_mets_newData)[!is.na(metKEGGIDs_new)] <- metKEGGIDs_new[!is.na(metKEGGIDs_new)]
colnames(ccmn_norm_mets_hybrid)[!is.na(metKEGGIDs_hybrid)] <- metKEGGIDs_hybrid[!is.na(metKEGGIDs_hybrid)]

####################################################################################################################################
###################################                      GOOD OLD DATA ANALYSIS                  ###################################
####################################################################################################################################
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs/gene_tab/genePresAbs_functions.R")

####################
#
# Get strain cluster classification and obtain contingency tabs for doing Fisher tests
#
####################################################################################################################################

metClusts_oldGood <- getMetClusts(ccmn_norm_mets_good_old)
save(metClusts_oldGood, file = "metClusts_oldGood.RData")

# All genes
contTab_oldGood_allGenes <- getContTab(gene_tab, metClusts_oldGood)

# Removing PA14 and, after that, all across-strain common genes
contTab_oldGood_filt <- getContTab(gene_tab_filt, metClusts_oldGood[2:length(metClusts_oldGood)])

# Only enzymatic genes (the ones with EC number)
contTab_oldGood_enz <- getContTab(gene_enz_tab, metClusts_oldGood)

# Enzimatic genes without PA14 and that are not shared across all strains
contTab_oldGood_filtEnz <- getContTab(gene_enz_tab_filt, metClusts_oldGood[2:length(metClusts_oldGood)])

####################
#
# Fisher test between 4 groups of the tables obtained in previous step
#
####################################################################################################################################

fish_oldGood_allGen <- doFisher(contTab_oldGood_allGenes, metClusts_oldGood)                                # 180 sign genes
fish_oldGood_filt <- doFisher(contTab_oldGood_filt, metClusts_oldGood[2:length(metClusts_oldGood)])         # 139 sign genes
fish_oldGood_enz <- doFisher(contTab_oldGood_enz, metClusts_oldGood)                                        # 12 sign genes
fish_oldGood_filtEnz <- doFisher(contTab_oldGood_filtEnz, metClusts_oldGood[2:length(metClusts_oldGood)])   # 9 sign genes
save(fish_oldGood_allGen, file = "fish_oldGood_allGen.RData")
save(fish_oldGood_filt, file = "fish_oldGood_filt.RData")
save(fish_oldGood_enz, file = "fish_oldGood_enz.RData")
save(fish_oldGood_filtEnz, file = "fish_oldGood_allGen.RData")

####################
#
# Obtain genes and groups of genes that produce differences in one metabolite of the reportedly differential ones
#
#####################################################################################################################################

#####################################################################################################################################
### Build perform a Mann-Whitney test for each gene against each metabolite, by dividing the strains ################################
### according to if they have or not a gene. Adjust p-value with Bonferroni-Yekuteli. Only get the   ################################
### p-value if the number of strains having the gene is greater than one. Pick genes with            ################################
### alpha < 0.05.                                                                                    ################################ 
#####################################################################################################################################

mannWhitPerGene_oldGood_allGenes <- getMannWhitPerGene(gene_tab, ccmn_norm_mets_good_old, p_adjust = T, method = "BY")
save(mannWhitPerGene_oldGood_allGenes, file = "mannWhitPerGene_oldGood_allGenes.RData")
mannWhitPerGene_oldGood_filt <- getMannWhitPerGene(gene_tab_filt, ccmn_norm_mets_good_old, p_adjust = T, method = "BY")
save(mannWhitPerGene_oldGood_filt, file = "mannWhitPerGene_oldGood_filt.RData")
mannWhitPerGene_oldGood_enz <- getMannWhitPerGene(gene_enz_tab, ccmn_norm_mets_good_old, p_adjust = T, method = "BY")
save(mannWhitPerGene_oldGood_enz, file = "mannWhitPerGene_oldGood_enz.RData")
mannWhitPerGene_oldGood_filtEnz <- getMannWhitPerGene(gene_enz_tab_filt, ccmn_norm_mets_good_old, p_adjust = T, method = "BY")
save(mannWhitPerGene_oldGood_filtEnz, file = "mannWhitPerGene_oldGood_filtEnz.RData")

# Plot z score for each one 
tiff("gene_met_corr_oldGood_allGen.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_oldGood_allGenes$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Metabolites", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_oldGood_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_oldGood_filt$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Metabolites", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_oldGood_enz.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_oldGood_enz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Metabolites", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_oldGood_filtEnz.tiff", width = 6000, height = 5000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_oldGood_filtEnz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Metabolites", 
          ylab = "Genes", 
          main = "Correlation accessory enzymatic genome and metabolome", 
          margins = c(8, 50), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

####################################################################################################################################
### Take only metabolites with alpha < 0.05 and build a matrix where in each column are all         ################################
### the metabolites affected individually by the presence/absence of an individual gene.            ################################
####################################################################################################################################

mannWhitPerGeneFilt_oldGood_allGenes <- filtMannWhitPerGene(mannWhitPerGene_oldGood_allGenes)
mannWhitPerGeneFilt_oldGood_filt <- filtMannWhitPerGene(mannWhitPerGene_oldGood_filt)
mannWhitPerGeneFilt_oldGood_enz <- filtMannWhitPerGene(mannWhitPerGene_oldGood_enz)
mannWhitPerGeneFilt_oldGood_filtEnz <- filtMannWhitPerGene(mannWhitPerGene_oldGood_filtEnz)

unique(as.vector(mannWhitPerGeneFilt_oldGood_allGenes$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_allGenes$as_mat))]) # 83/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_filt$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_filt$as_mat))])         # 82/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_enz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_enz$as_mat))])           # 59/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_filtEnz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_filtEnz$as_mat))])   # 52/86 mets

####################################################################################################################################
### Build a list with genes whose presence/absence individually was found to be correlated with     ################################
### differences in abundance of each one of the metabolites of the 3 sets of differential           ################################
### metabolites.                                                                                    ################################
####################################################################################################################################

diffGenesPerDiffMet_oldGood_allGenes <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_allGenes,
                                                              diffMetObjkt = diffMets_oldGood)
diffGenesPerDiffMet_oldGood_filt <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_filt,
                                                          diffMetObjkt = diffMets_oldGood)
diffGenesPerDiffMet_oldGood_enz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_enz,
                                                         diffMetObjkt = diffMets_oldGood)
diffGenesPerDiffMet_oldGood_filtEnz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_oldGood_filtEnz,
                                                             diffMetObjkt = diffMets_oldGood)

####################################################################################################################################
### Obtain, for each one of the metabolites, a matrix of presence/absence of the genes reported     ################################
### to affect that metabolite, including the strains that belong to booth of the clusters           ################################
### to whom the set of differental metabolites belong. Include a total of times those genes appear  ################################
### in each strain of the matrix.                                                                   ################################
####################################################################################################################################

geneMatMets_oldGood_allGenes <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_allGenes,
                                                  genePresAbsObjkt = gene_tab,
                                                  metClustObjkt = metClusts_oldGood)
geneMatMets_oldGood_filt <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_filt,
                                              genePresAbsObjkt = gene_tab_filt,
                                              metClustObjkt = metClusts_oldGood[2:length(metClusts_oldGood)])
geneMatMets_oldGood_enz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_enz,
                                             genePresAbsObjkt = gene_enz_tab,
                                             metClustObjkt = metClusts_oldGood)
geneMatMets_oldGood_filtEnz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_filtEnz,
                                                 genePresAbsObjkt = gene_enz_tab_filt,
                                                 metClustObjkt = metClusts_oldGood[2:length(metClusts_oldGood)])

####################################################################################################################################
### Apply Mann-Whitney test for the genes associated to the differential metabolites for each       ################################
### pair of clusters, with alpha = 0.05, by one side comparing vectors of individual genes,         ################################
### and by other comparing the vector of the total times the genes that affect one metabolite       ################################
### appears in each strain.                                                                         ################################
####################################################################################################################################

allDiffGenes_oldGood_allGenes <- mannWhit_tabs(geneMatMets_oldGood_allGenes)
allDiffGenes_oldGood_filt <- mannWhit_tabs(geneMatMets_oldGood_filt)
allDiffGenes_oldGood_enz <- mannWhit_tabs(geneMatMets_oldGood_enz)
allDiffGenes_oldGood_filtEnz <- mannWhit_tabs(geneMatMets_oldGood_filtEnz)

diffGenesPerMet_oldGood_allGenes <- filtMannWhitTabs(allDiffGenes_oldGood_allGenes, alpha = 0.05)
diffGenesPerMet_oldGood_filt <- filtMannWhitTabs(allDiffGenes_oldGood_filt, alpha = 0.05)
diffGenesPerMet_oldGood_enz <- filtMannWhitTabs(allDiffGenes_oldGood_enz, alpha = 0.05)
diffGenesPerMet_oldGood_filtEnz <- filtMannWhitTabs(allDiffGenes_oldGood_filtEnz, alpha = 0.05)

# Build heatmap with presence/absence of the genes that affect diffmets between 2 major groups, with strains ordered according 
# to metabolic clusters.
C1_2Genes_old_filtEnz <- unique(unlist(diffGenesPerMet_oldGood_filtEnz$Total$`1&2`))[unique(unlist(diffGenesPerMet_oldGood_filtEnz$Total$`1&2`)) != "Any relationship between the gene and the metabolite"]
presAbsC1_2Genes_old_filtEnz <- gene_enz_tab_filt[, C1_2Genes_old_filtEnz]
clustOrdC1_2_old <- names(sort(metClusts_oldGood[2:length(metClusts_oldGood)]))
presAbsC1_2Genes_old_filtEnz <- t(presAbsC1_2Genes_old_filtEnz[clustOrdC1_2_old, ])
save(presAbsC1_2Genes_old_filtEnz, file = "presAbsC1_2Genes_old_filtEnz.RData")
save(C1_2Genes_old_filtEnz, file = "C1_2Genes_old_filtEnz.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/presAbsC1_2Genes_old_filtEnz.RData")

tiff("minedDiffGenes_old_filtEnz.tiff", width = 7000, height = 5000, units = "px", pointsize = 100)
heatmap.2(presAbsC1_2Genes_old_filtEnz, Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(clustOrdC1_2_old, rownames(gene_tab))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of diff genes", margins = c(6, 50), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          #sepcolor = "black",
          #colsep=1:ncol(t(gene_tab)),
          #rowsep=1:nrow(t(gene_tab)),
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
dev.off()

####################################################################################################################################
###################################                        NEW DATA ANALYSIS                     ###################################
####################################################################################################################################

####################
#
# Get strain cluster classification and obtain contingency tabs for doing Fisher tests
#
####################################################################################################################################

metClusts_new <- getMetClusts(ccmn_norm_mets_newData)
save(metClusts_new, file = "metClusts_new.RData")

# All genes
contTab_new_allGenes <- getContTab(gene_tab, metClusts_new)

# Removing PA14 and, after that, all across-strain common genes
contTab_new_filt <- getContTab(gene_tab_filt, metClusts_new[2:length(metClusts_new)])

# Only enzymatic genes (the ones with EC number)
contTab_new_enz <- getContTab(gene_enz_tab, metClusts_new)

# Enzimatic genes without PA14 and that are not shared across all strains
contTab_new_filtEnz <- getContTab(gene_enz_tab_filt, metClusts_new[2:length(metClusts_new)])

####################
#
# Fisher test between 4 groups of the tables obtained in previous step
#
####################################################################################################################################

fish_new_allGen <- doFisher(contTab_new_allGenes, metClusts_new)                                # 52 sign genes
fish_new_filt <- doFisher(contTab_new_filt, metClusts_new[2:length(metClusts_new)])             # 52 sign genes
fish_new_enz <- doFisher(contTab_new_enz, metClusts_new)                                        # 1 sign genes
fish_new_filtEnz <- doFisher(contTab_new_filtEnz, metClusts_new[2:length(metClusts_new)])       # 1 sign genes

####################
#
# Obtain genes and groups of genes that produce differences in one metabolite of the reportedly differential ones
#
#####################################################################################################################################

#####################################################################################################################################
### Build perform a Mann-Whitney test for each gene against each metabolite, by dividing the strains ################################
### according to if they have or not a gene. Adjust p-value with Bonferroni-Yekuteli. Only get the   ################################
### p-value if the number of strains having the gene is greater than one. Pick genes with            ################################
### alpha < 0.05.                                                                                    ################################ 
#####################################################################################################################################

mannWhitPerGene_new_allGenes <- getMannWhitPerGene(gene_tab, ccmn_norm_mets_newData, p_adjust = T, method = "BY")
save(mannWhitPerGene_new_allGenes, file = "mannWhitPerGene_new_allGenes.RData")
mannWhitPerGene_new_filt <- getMannWhitPerGene(gene_tab_filt, ccmn_norm_mets_newData, p_adjust = T, method = "BY")
save(mannWhitPerGene_new_filt, file = "mannWhitPerGene_new_filt.RData")
mannWhitPerGene_new_enz <- getMannWhitPerGene(gene_enz_tab, ccmn_norm_mets_newData, p_adjust = T, method = "BY")
save(mannWhitPerGene_new_enz, file = "mannWhitPerGene_new_enz.RData")
mannWhitPerGene_new_filtEnz <- getMannWhitPerGene(gene_enz_tab_filt, ccmn_norm_mets_newData, p_adjust = T, method = "BY")
save(mannWhitPerGene_new_filtEnz, file = "mannWhitPerGene_new_filtEnz.RData")

# Plot z score for each one 
tiff("gene_met_corr_new_allGen.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_new_allGenes$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_new_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_new_filt$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_new_enz.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_new_enz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_new_filtEnz.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_new_filtEnz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

####################################################################################################################################
### Take only metabolites with alpha < 0.05 and build a matrix where in each column are all         ################################
### the metabolites affected individually by the presence/absence of an individual gene.            ################################
####################################################################################################################################

mannWhitPerGeneFilt_new_allGenes <- filtMannWhitPerGene(mannWhitPerGene_new_allGenes)
mannWhitPerGeneFilt_new_filt <- filtMannWhitPerGene(mannWhitPerGene_new_filt)
mannWhitPerGeneFilt_new_enz <- filtMannWhitPerGene(mannWhitPerGene_new_enz)
mannWhitPerGeneFilt_new_filtEnz <- filtMannWhitPerGene(mannWhitPerGene_new_filtEnz)

unique(as.vector(mannWhitPerGeneFilt_new_allGenes$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_new_allGenes$as_mat))]) # 70/70 mets
unique(as.vector(mannWhitPerGeneFilt_new_filt$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_new_filt$as_mat))])         # 69/70 mets
unique(as.vector(mannWhitPerGeneFilt_new_enz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_new_enz$as_mat))])           # 52/70 mets
unique(as.vector(mannWhitPerGeneFilt_new_filtEnz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_new_filtEnz$as_mat))])   # 54/70 mets

####################################################################################################################################
### Build a list with genes whose presence/absence individually was found to be correlated with     ################################
### differences in abundance of each one of the metabolites of the 3 sets of differential           ################################
### metabolites.                                                                                    ################################
####################################################################################################################################

diffGenesPerDiffMet_new_allGenes <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_new_allGenes,
                                                          diffMetObjkt = diffMets_newData)
diffGenesPerDiffMet_new_filt <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_new_filt,
                                                      diffMetObjkt = diffMets_newData)
diffGenesPerDiffMet_new_enz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_new_enz,
                                                     diffMetObjkt = diffMets_newData)
diffGenesPerDiffMet_new_filtEnz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_new_filtEnz,
                                                         diffMetObjkt = diffMets_newData)

####################################################################################################################################
### Obtain, for each one of the metabolites, a matrix of presence/absence of the genes reported     ################################
### to affect that metabolite, including the strains that belong to booth of the clusters           ################################
### to whom the set of differental metabolites belong. Include a total of times those genes appear  ################################
### in each strain of the matrix.                                                                   ################################
####################################################################################################################################

geneMatMets_new_allGenes <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_new_allGenes,
                                                  genePresAbsObjkt = gene_tab,
                                                  metClustObjkt = metClusts_new)
geneMatMets_new_filt <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_new_filt,
                                              genePresAbsObjkt = gene_tab_filt,
                                              metClustObjkt = metClusts_new[2:length(metClusts_new)])
geneMatMets_new_enz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_new_enz,
                                             genePresAbsObjkt = gene_enz_tab,
                                             metClustObjkt = metClusts_new)
geneMatMets_new_filtEnz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_new_filtEnz,
                                                 genePresAbsObjkt = gene_enz_tab_filt,
                                                 metClustObjkt = metClusts_new[2:length(metClusts_new)])

####################################################################################################################################
### Apply Mann-Whitney test for the genes associated to the differential metabolites for each       ################################
### pair of clusters, with alpha = 0.05, by one side comparing vectors of individual genes,         ################################
### and by other comparing the vector of the total times the genes that affect one metabolite       ################################
### appears in each strain.                                                                         ################################
####################################################################################################################################

allDiffGenes_new_allGenes <- mannWhit_tabs(geneMatMets_new_allGenes)
allDiffGenes_new_filt <- mannWhit_tabs(geneMatMets_new_filt)
allDiffGenes_new_enz <- mannWhit_tabs(geneMatMets_new_enz)
allDiffGenes_new_filtEnz <- mannWhit_tabs(geneMatMets_new_filtEnz)

diffGenesPerMet_new_allGenes <- filtMannWhitTabs(allDiffGenes_new_allGenes, alpha = 0.05)
diffGenesPerMet_new_filt <- filtMannWhitTabs(allDiffGenes_new_filt, alpha = 0.05)
diffGenesPerMet_new_enz <- filtMannWhitTabs(allDiffGenes_new_enz, alpha = 0.05)
diffGenesPerMet_new_filtEnz <- filtMannWhitTabs(allDiffGenes_new_filtEnz, alpha = 0.05)

# Build heatmap with presence/absence of the genes that affect diffmets between 2 major groups, with strains ordered according 
# to metabolic clusters.
C1_2Genes_new_filtEnz <- unique(unlist(diffGenesPerMet_new_filtEnz$Total$`1&2`))[unique(unlist(diffGenesPerMet_new_filtEnz$Total$`1&2`)) != "Any relationship between the gene and the metabolite"]
presAbsC1_2Genes_new_filtEnz <- gene_enz_tab_filt[, C1_2Genes_new_filtEnz]
clustOrdC1_2_new <- names(sort(metClusts_new[2:length(metClusts_new)]))
presAbsC1_2Genes_new_filtEnz <- t(presAbsC1_2Genes_new_filtEnz[clustOrdC1_2_new, ])
save(presAbsC1_2Genes_new_filtEnz, file = "presAbsC1_2Genes_new_filtEnz.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/presAbsC1_2Genes_new_filtEnz.RData")

tiff("minedDiffGenes_new_filtEnz.tiff", width = 6000, height = 5000, units = "px", pointsize = 140)
heatmap.2(presAbsC1_2Genes_new_filtEnz, Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(clustOrdC1_2_new, rownames(gene_tab))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of diff genes", margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
          #sepcolor = "black",
          #colsep=1:ncol(t(gene_tab)),
          #rowsep=1:nrow(t(gene_tab)),
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
dev.off()

####################################################################################################################################
###################################                        HYBRID DATA ANALYSIS                     ###################################
####################################################################################################################################

####################
#
# Get strain cluster classification and obtain contingency tabs for doing Fisher tests
#
####################################################################################################################################

metClusts_hyb <- getMetClusts(ccmn_norm_mets_hybrid)
save(metClusts_hyb, file = "metClusts_hyb.RData")

# All genes
contTab_hyb_allGenes <- getContTab(gene_tab, metClusts_hyb)

# Removing PA14 and, after that, all across-strain common genes
contTab_hyb_filt <- getContTab(gene_tab_filt, metClusts_hyb[2:length(metClusts_hyb)])

# Only enzymatic genes (the ones with EC number)
contTab_hyb_enz <- getContTab(gene_enz_tab, metClusts_hyb)

# Enzimatic genes without PA14 and that are not shared across all strains
contTab_hyb_filtEnz <- getContTab(gene_enz_tab_filt, metClusts_hyb[2:length(metClusts_hyb)])

####################
#
# Fisher test between 4 groups of the tables obtained in previous step
#
####################################################################################################################################

fish_hyb_allGen <- doFisher(contTab_hyb_allGenes, metClusts_hyb)                                # 48 sign genes
fish_hyb_filt <- doFisher(contTab_hyb_filt, metClusts_hyb[2:length(metClusts_hyb)])             # 69 sign genes
fish_hyb_enz <- doFisher(contTab_hyb_enz, metClusts_hyb)                                        # 4 sign genes
fish_hyb_filtEnz <- doFisher(contTab_hyb_filtEnz, metClusts_hyb[2:length(metClusts_hyb)])       # 5 sign genes

####################
#
# Obtain genes and groups of genes that produce differences in one metabolite of the reportedly differential ones
#
#####################################################################################################################################

#####################################################################################################################################
### Build perform a Mann-Whitney test for each gene against each metabolite, by dividing the strains ################################
### according to if they have or not a gene. Adjust p-value with Bonferroni-Yekuteli. Only get the   ################################
### p-value if the number of strains having the gene is greater than one. Pick genes with            ################################
### alpha < 0.05.                                                                                    ################################ 
#####################################################################################################################################

mannWhitPerGene_hyb_allGenes <- getMannWhitPerGene(gene_tab, ccmn_norm_mets_hybrid, p_adjust = T, method = "BY")
save(mannWhitPerGene_hyb_allGenes, file = "mannWhitPerGene_hyb_allGenes.RData")
mannWhitPerGene_hyb_filt <- getMannWhitPerGene(gene_tab_filt, ccmn_norm_mets_hybrid, p_adjust = T, method = "BY")
save(mannWhitPerGene_hyb_filt, file = "mannWhitPerGene_hyb_filt.RData")
mannWhitPerGene_hyb_enz <- getMannWhitPerGene(gene_enz_tab, ccmn_norm_mets_hybrid, p_adjust = T, method = "BY")
save(mannWhitPerGene_hyb_enz, file = "mannWhitPerGene_hyb_enz.RData")
mannWhitPerGene_hyb_filtEnz <- getMannWhitPerGene(gene_enz_tab_filt, ccmn_norm_mets_hybrid, p_adjust = T, method = "BY")
save(mannWhitPerGene_hyb_filtEnz, file = "mannWhitPerGene_hyb_filtEnz.RData")

# Plot z score for each one 
tiff("gene_met_corr_hyb_allGen.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_hyb_allGenes$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_hyb_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_hyb_filt$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_hyb_enz.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_hyb_enz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

tiff("gene_met_corr_hyb_filtEnz.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(mannWhitPerGene_hyb_filtEnz$z.score, Rowv = T, Colv = T, 
          distfun = function(x) dist(x, method = "euclidean"), 
          #density.info = "none", 
          hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(100), 
          breaks = 101, 
          notecol = NULL, 
          trace = "none", 
          xlab = "Strains", 
          ylab = "Genes", 
          main = "Correlation Accessory genome and metabolome", 
          margins = c(6, 13), 
          keysize = 1,
          sepwidth = c(0.1, 0.05),
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
dev.off()

####################################################################################################################################
### Take only metabolites with alpha < 0.05 and build a matrix where in each column are all         ################################
### the metabolites affected individually by the presence/absence of an individual gene.            ################################
####################################################################################################################################

mannWhitPerGeneFilt_hyb_allGenes <- filtMannWhitPerGene(mannWhitPerGene_hyb_allGenes)
mannWhitPerGeneFilt_hyb_filt <- filtMannWhitPerGene(mannWhitPerGene_hyb_filt)
mannWhitPerGeneFilt_hyb_enz <- filtMannWhitPerGene(mannWhitPerGene_hyb_enz)
mannWhitPerGeneFilt_hyb_filtEnz <- filtMannWhitPerGene(mannWhitPerGene_hyb_filtEnz)

unique(as.vector(mannWhitPerGeneFilt_hyb_allGenes$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_hyb_allGenes$as_mat))]) # 81/89 mets
unique(as.vector(mannWhitPerGeneFilt_hyb_filt$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_hyb_filt$as_mat))])         # 63/89 mets
unique(as.vector(mannWhitPerGeneFilt_hyb_enz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_hyb_enz$as_mat))])           # 62/89 mets
unique(as.vector(mannWhitPerGeneFilt_hyb_filtEnz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_hyb_filtEnz$as_mat))])   # 57/89 mets

####################################################################################################################################
### Build a list with genes whose presence/absence individually was found to be correlated with     ################################
### differences in abundance of each one of the metabolites of the 3 sets of differential           ################################
### metabolites.                                                                                    ################################
####################################################################################################################################

diffGenesPerDiffMet_hyb_allGenes <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_hyb_allGenes,
                                                          diffMetObjkt = diffMets_hybrid)
diffGenesPerDiffMet_hyb_filt <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_hyb_filt,
                                                      diffMetObjkt = diffMets_hybrid)
diffGenesPerDiffMet_hyb_enz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_hyb_enz,
                                                     diffMetObjkt = diffMets_hybrid)
diffGenesPerDiffMet_hyb_filtEnz <- getDiffGenePerDiffMet(mannWhitPerGeneFiltObjkt = mannWhitPerGeneFilt_hyb_filtEnz,
                                                         diffMetObjkt = diffMets_hybrid)

####################################################################################################################################
### Obtain, for each one of the metabolites, a matrix of presence/absence of the genes reported     ################################
### to affect that metabolite, including the strains that belong to booth of the clusters           ################################
### to whom the set of differental metabolites belong. Include a total of times those genes appear  ################################
### in each strain of the matrix.                                                                   ################################
####################################################################################################################################

geneMatMets_hyb_allGenes <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_hyb_allGenes,
                                              genePresAbsObjkt = gene_tab,
                                              metClustObjkt = metClusts_hyb)
geneMatMets_hyb_filt <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_hyb_filt,
                                          genePresAbsObjkt = gene_tab_filt,
                                          metClustObjkt = metClusts_hyb[2:length(metClusts_hyb)])
geneMatMets_hyb_enz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_hyb_enz,
                                         genePresAbsObjkt = gene_enz_tab,
                                         metClustObjkt = metClusts_hyb)
geneMatMets_hyb_filtEnz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_hyb_filtEnz,
                                             genePresAbsObjkt = gene_enz_tab_filt,
                                             metClustObjkt = metClusts_hyb[2:length(metClusts_hyb)])

####################################################################################################################################
### Apply Mann-Whitney test for the genes associated to the differential metabolites for each       ################################
### pair of clusters, with alpha = 0.05, by one side comparing vectors of individual genes,         ################################
### and by other comparing the vector of the total times the genes that affect one metabolite       ################################
### appears in each strain.                                                                         ################################
####################################################################################################################################

allDiffGenes_hyb_allGenes <- mannWhit_tabs(geneMatMets_hyb_allGenes)
allDiffGenes_hyb_filt <- mannWhit_tabs(geneMatMets_hyb_filt)
allDiffGenes_hyb_enz <- mannWhit_tabs(geneMatMets_hyb_enz)
allDiffGenes_hyb_filtEnz <- mannWhit_tabs(geneMatMets_hyb_filtEnz)

diffGenesPerMet_hyb_allGenes <- filtMannWhitTabs(allDiffGenes_hyb_allGenes, alpha = 0.05)
diffGenesPerMet_hyb_filt <- filtMannWhitTabs(allDiffGenes_hyb_filt, alpha = 0.05)
diffGenesPerMet_hyb_enz <- filtMannWhitTabs(allDiffGenes_hyb_enz, alpha = 0.05)
diffGenesPerMet_hyb_filtEnz <- filtMannWhitTabs(allDiffGenes_hyb_filtEnz, alpha = 0.05)

# Build heatmap with presence/absence of the genes that affect diffmets between 2 major groups, with strains ordered according 
# to metabolic clusters.
C1_2Genes_hyb_filtEnz <- unique(unlist(diffGenesPerMet_hyb_filtEnz$Total$`1&2`))[unique(unlist(diffGenesPerMet_hyb_filtEnz$Total$`1&2`)) != "Any relationship between the gene and the metabolite"]
presAbsC1_2Genes_hyb_filtEnz <- gene_enz_tab_filt[, C1_2Genes_hyb_filtEnz]
clustOrdC1_2_hyb <- names(sort(metClusts_hyb[2:length(metClusts_hyb)]))
presAbsC1_2Genes_hyb_filtEnz <- t(presAbsC1_2Genes_hyb_filtEnz[clustOrdC1_2_hyb, ])
save(presAbsC1_2Genes_hyb_filtEnz, file = "presAbsC1_2Genes_hyb_filtEnz.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/presAbsC1_2Genes_hyb_filtEnz.RData")


# Using hybrid data group 1&2 is empty