#####################################################################################################################################
#######################################                                                   ###########################################
######################################           GENE PRESENCE/ABSENCE ANALYSIS            ##########################################
#######################################          Remodeled for new hybrid data            ###########################################
########################################                                                 ############################################
#####################################################################################################################################
if(!require(rtracklayer)) BiocManager::install("rtracklayer")
library(rtracklayer)
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
library(KEGGREST)

setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs")

####################
#
# Load the data
#
####################################################################################################################################

prokkAnot <- read.csv("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/untitled_folder/gene_presence_absence.csv")

gffsAllStrains <- list()
for(i in 1:length(list.files(path = "../../prokkAnnotStuff/gff3"))){
        gff <- readGFF(filepath = paste("../../prokkAnnotStuff/gff3/",
                                        list.files(path = "../../prokkAnnotStuff/gff3/")[i],
                                        sep = ""))
        gffsAllStrains[[i]] <- gff
}
names(gffsAllStrains) <- gsub(".gff", "", list.files(path = "../../prokkAnnotStuff/gff3/"))

tail(apply(prokkAnot[, 15:ncol(prokkAnot)], 2, asa))

prokkaBinarize <- function(prokkaMat){
        bin <- function(x){
                pos <- nchar(x) > 0
                bin <- rep(0, length(x))
                bin[pos] <- 1
                return(bin)
        }
        binarized <- apply(prokkaMat[, 15:ncol(prokkaMat)], 2, bin)
        rownames(binarized) <- prokkaMat[, 1]
        return(binarized)
}


gene_tab <- t(prokkaBinarize(prokkAnot))

rownames(gene_tab) <- gsub("_", rownames(gene_tab), replacement =  ".")

####################
#
# Plot it
#
####################################################################################################################################
library(gplots)
library(ggplot2)

#Make colorscale equivalent to the one that we've been using in previous heatmaps

strain_names <- gsub(".gff", "", list.files(path = "../../prokkAnnotStuff/gff3/"))

cols <- topo.colors(length(strain_names))
cols[15] <- '#000000'
cols[21] <- '#555555'

tiff("gene_tab_prokka.tiff", width = 2000, height = 4000, units = "px", pointsize = 50)

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

# Remove every gene that is present in all strains.
gene_tab_filt <- gene_tab


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

enzIdxs <- c()
for(i in seq_along(strain_names)){
        idx <- which(prokkAnot[[strain_names[i]]] %in% gffsAllStrains[[strain_names[i]]]$ID[!is.na(gffsAllStrains[[strain_names[i]]]$eC_number)]) 
        enzIdxs <- c(enzIdxs, idx)
}
enzIdxs <- unique(enzIdxs)

gene_enz_tab <- gene_tab[, enzIdxs]

save(gene_enz_tab, file = "gene_enz_tab.RData")

# Remove PA14 and the genes that are in all strains and keep only enzymatic ones.
gene_enz_tab_filt <- gene_tab_filt[, colnames(gene_tab_filt) %in% colnames(gene_enz_tab)]

save(gene_enz_tab_filt, file = "gene_enz_tab_filt.RData")

dictEnzymes <- prokkAnot[enzIdxs, -(4:14)]
names(gffsAllStrains)
strain_names
for(i in seq_along(strain_names)){
        posInd <- which(nchar(as.character(dictEnzymes[[make.names(strain_names[i])]])) > 0)
        ecVec <- rep(NA, nrow(dictEnzymes))
        for(j in seq_along(posInd)){
                ecNum <- gffsAllStrains[[strain_names[i]]]$eC_number[which(gffsAllStrains[[strain_names[i]]]$ID == dictEnzymes[[make.names(strain_names[i])]][posInd[j]])]
                if(!is.null(ecNum)){
                        ecVec[posInd[j]] <- as.character(ecNum)
                }
        }
        dictEnzymes <- cbind.data.frame(dictEnzymes, ecVec, stringsAsFactors = F)
}


#GFF file of UCBPP-PA14 doesn't have ec numbers

# Look if there are any disrepancies regarding the EC numbers of the genes in the table.
discrepVec <- c()
for(i in 1:nrow(dictEnzymes)){
        discrep <- length(unique(unlist(dictEnzymes[i, 34:ncol(dictEnzymes)]))) > 2
        discrepVec <- c(discrepVec, discrep)
}
sum(discrepVec)

genesWDiscrep <- list()
for(i in seq_along(which(discrepVec))){
        dup <- unique(unlist(dictEnzymes[which(discrepVec)[i], 34:ncol(dictEnzymes)]))
        dup <- dup[!is.na(dup)]
        genesWDiscrep[[i]] <- dup
}
names(genesWDiscrep) <- dictEnzymes$Gene[which(discrepVec)]

ECnums <- apply(dictEnzymes[, 34:ncol(dictEnzymes)], 1, function(x) names(sort(table(x), decreasing = T))[1])

dictEnzymes <- cbind.data.frame(dictEnzymes[, 1:33], ECnums)

save(dictEnzymes, file = "dictEnzymes.RData")

# There are 56 genes with discrepancies according to the EC number. We're gonna align them to see if any of them looks very 
# different to the other ones.
if(!require(DECIPHER)) install.packages("DECIPHER")
library(DECIPHER)

F22031 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F22031/PROKKA_01212020.faa")
F23197 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F23197/PROKKA_01212020.faa")
F30658 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F30658/PROKKA_01212020.faa")
F34365 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F34365/PROKKA_01212020.faa")
F5677 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F5677/PROKKA_01212020.faa")
F63912 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F63912/PROKKA_01212020.faa")
F9670 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/F9670/PROKKA_01212020.faa")
H27930 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/H27930/PROKKA_01212020.faa")
H47921 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/H47921/PROKKA_01212020.faa")
H5708 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/H5708/PROKKA_01212020.faa")
M1608 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/M1608/PROKKA_01212020.faa")
M37351 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/M37351/PROKKA_01212020.faa")
M55212 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/M55212/PROKKA_01212020.faa")
M74707 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/M74707/PROKKA_01212020.faa")
PA14_jbx <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/PA14_jbx/PROKKA_01212020.faa")
S86968 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/S86968/PROKKA_01212020.faa")
T38079 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/T38079/PROKKA_01212020.faa")
T52373 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/T52373/PROKKA_01212020.faa")
T6313 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/T6313/PROKKA_01212020.faa")
T63266 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/T63266/PROKKA_01212020.faa")
W16407 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W16407/PROKKA_01212020.faa")
W25637 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W25637/PROKKA_01212020.faa")
W36662 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W36662/PROKKA_01212020.faa")
W45909 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W45909/PROKKA_01212020.faa")
W60856 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W60856/PROKKA_01212020.faa")
W70332 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W70332/PROKKA_01212020.faa")
W91453 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/W91453/PROKKA_01212020.faa")
X78812 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/X78812/PROKKA_01212020.faa")
X9820 <- readAAStringSet("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/prokkaOutput/X9820/PROKKA_01212020.faa")

allFAAs <- list("F22031" = F22031,
                "F23197" = F23197,
                "F30658" = F30658,
                "F34365" = F34365,
                "F5677" = F5677,
                "F63912" = F63912,
                "F9670" = F9670,
                "H27930" = H27930,
                "H47921" = H47921, 
                "H5708" = H5708,
                "M1608" = M1608,
                "M37351" = M37351,
                "M55212" = M55212,
                "M74707" = M74707,
                "PA14_jbx" = PA14_jbx,
                "S86968" = S86968,
                "T38079" = T38079,
                "T52373" = T52373,
                "T6313" = T6313,
                "T63266" = T63266,
                "W16407" = W16407,
                "W25637" = W25637,
                "W36662" = W36662,
                "W45909" = W45909,
                "W60856" = W60856,
                "W70332" = W70332,
                "W91453" = W91453,
                "X78812" = X78812,
                "X9820" = X9820)

discrepGenesAAStrSets <- list()
discrepGenesAlignments <- list()
for(i in seq_along(genesWDiscrep)){
        print(i)
        geneAllStrains <- list()
        strainsWgene <- colnames(dictEnzymes[, 4:33])[nchar(as.character(unlist(dictEnzymes[dictEnzymes$Gene == names(genesWDiscrep)[i], 4:33]))) > 0]
        strainsWgene <- strainsWgene[strainsWgene %in% names(allFAAs)]
        for(j in seq_along(strainsWgene)){
                locName <- as.character(dictEnzymes[[strainsWgene[j]]][dictEnzymes$Gene == names(genesWDiscrep)[i]])
                genSeq <- allFAAs[[strainsWgene[j]]][[grep(locName, names(allFAAs[[strainsWgene[j]]]))]]
                geneAllStrains[[j]] <- genSeq
        }
        names(geneAllStrains) <- strainsWgene
        geneAAStrSet <- AAStringSet(geneAllStrains)
        alignmnt <- AlignSeqs(geneAAStrSet)
        discrepGenesAAStrSets[[i]] <- geneAAStrSet
        discrepGenesAlignments[[i]] <- alignmnt
}
names(discrepGenesAAStrSets) <- names(genesWDiscrep)
names(discrepGenesAlignments) <- names(genesWDiscrep)

BrowseSeqs(discrepGenesAlignments[[1]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[2]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[3]])  # T52373 is shorter. Rest ovelaps, with some point mutations. T52373 is the one with a different EC number (1.10.3.-).
BrowseSeqs(discrepGenesAlignments[[4]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[5]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[6]])  # F9670 is shorter. Rest ovelaps, with some point mutations. F9670 is the one with a different EC number (2.7.1.17).
BrowseSeqs(discrepGenesAlignments[[7]])  # F9670 and T38079 are shorter. Rest ovelaps, with some point mutations. They have same EC than most (3.5.1.-). H47921 and S86968 are different (4.2.99.20, 1.11.1.-).
BrowseSeqs(discrepGenesAlignments[[8]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[9]])  # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[10]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[11]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[12]]) # T52373 is shorter. Rest ovelaps, with some point mutations. T52373 is the one with a different EC number (2.3.1.41).
BrowseSeqs(discrepGenesAlignments[[13]]) # F23197, F9670 and T52373 are shorter. W16407 is larger. W16407 is different (2.4.1.250)
BrowseSeqs(discrepGenesAlignments[[14]]) # F9670 and T38079 are larger.  They also have a different EC (7.6.2.8)
BrowseSeqs(discrepGenesAlignments[[15]]) # F9670 is shorter. But it's same than most of the others (1.7.1.17)
BrowseSeqs(discrepGenesAlignments[[16]]) # X9820 is shorter. Also has a different EC (7.6.2.8)
BrowseSeqs(discrepGenesAlignments[[17]]) # T63266 is shorter. It's also the different one (2.7.13.3)
BrowseSeqs(discrepGenesAlignments[[18]]) # M74707 is larger. It's also the different one (1.3.8.13)
BrowseSeqs(discrepGenesAlignments[[19]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[20]]) # F9670 is shorter. It's also the different one (4.4.-.-)
BrowseSeqs(discrepGenesAlignments[[21]]) # H47921 is shorter. It's also the different one (1.1.1.-)
BrowseSeqs(discrepGenesAlignments[[22]]) # F9670 is larger. It's also the different one (7.6.2.8)
BrowseSeqs(discrepGenesAlignments[[23]]) # Half of them have 10 bp more than the other half. The rest overlaps, with point mutations.
BrowseSeqs(discrepGenesAlignments[[24]]) # W16407 is shorter. It's also the different one (2.1.1.-)
BrowseSeqs(discrepGenesAlignments[[25]]) # F9670 doesn't overlap as well as the other ones. It's also the different one (1.14.12.15)
BrowseSeqs(discrepGenesAlignments[[26]]) # F5677 is way larger. It's also the different one (2.7.11.1).
BrowseSeqs(discrepGenesAlignments[[27]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[28]]) # Same. 
BrowseSeqs(discrepGenesAlignments[[29]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[30]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[31]]) # W16407 is shorter and a fragment doesn't match. It's also the different one (2.1.1.297)
BrowseSeqs(discrepGenesAlignments[[32]]) # F9670 is shorter. It's also the different one (1.20.4.1)
BrowseSeqs(discrepGenesAlignments[[33]]) # PA14_jbx is shorter. It's also the different one (1.3.8.10)
BrowseSeqs(discrepGenesAlignments[[34]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[35]]) # H27930 is larger. It's also the different one (3.5.1.-)
BrowseSeqs(discrepGenesAlignments[[36]]) # H27930 and T6313 are shorter. They're also the different ones (1.13.11.33)
BrowseSeqs(discrepGenesAlignments[[37]]) # W16407 and W36662 are shorter. They're also the different ones (2.3.1.41 and 1.6.5.5)
BrowseSeqs(discrepGenesAlignments[[38]]) # H47921 is larger and W45909 is shorter. W45909 is the different one (1.1.1.387)
BrowseSeqs(discrepGenesAlignments[[39]]) # X78812 is shorter and X9820 has nonhmologous regions. X9820 is the different one (1.4.1.1)
BrowseSeqs(discrepGenesAlignments[[40]]) # W16407 is way shorter. It's also the different one (1.14.13.29)  
BrowseSeqs(discrepGenesAlignments[[41]]) # W91453 is shorter. It's also the different one (7.6.2.8)
BrowseSeqs(discrepGenesAlignments[[42]]) # H27930 and M1608 are differend and shorter. H27930 is the different one (1.2.1.71)
BrowseSeqs(discrepGenesAlignments[[43]]) # W70332 and W91453 are shorter. They're also the different ones (1.3.1.95)
BrowseSeqs(discrepGenesAlignments[[44]]) # F22031 is shorter. It's also the different one (3.4.6.13) 
BrowseSeqs(discrepGenesAlignments[[45]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[46]]) # F30658, M55212, W25637 and W91453 are way shorter. They're also the different ones (2.8.3.16)
BrowseSeqs(discrepGenesAlignments[[47]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[48]]) # F9670 and H47921 are shorter. F9670 has no name neither ec number for this locus. H47921 is the different one (2.7.13.3)
BrowseSeqs(discrepGenesAlignments[[49]]) # Same overall. Some point mutations.
BrowseSeqs(discrepGenesAlignments[[50]]) # F63912 and M74707 shorter. F63912 is the different one (2.6.1.106)
BrowseSeqs(discrepGenesAlignments[[51]]) # T52373 is shorter. It's also the different one (1.-.-.-)
BrowseSeqs(discrepGenesAlignments[[52]]) # Only 2 strains have that. Very different sequences.
BrowseSeqs(discrepGenesAlignments[[53]]) # Only 2 strains have that. Very different sequences.
BrowseSeqs(discrepGenesAlignments[[54]]) # Only 2 strains have that. T52373 shorter.
BrowseSeqs(discrepGenesAlignments[[55]]) # H27930 larger. It's also the different one (1.1.1.87)
BrowseSeqs(discrepGenesAlignments[[56]]) # W91453 shorter. It's also the different one (1.9.6.1)

####################
#
# Obtain heatmaps previous sets
#
####################################################################################################################################
tiff("gene_tab_filt.tiff", width = 5000, height = 15000, units = "px", pointsize = 50)
heatmap.2(t(as.matrix(gene_tab_filt)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "column", 
          col = redgreen(75), breaks = 76, ColSideColors = cols, notecol = NULL, trace = "none", xlab = "Strains", 
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
          col = redgreen(75), breaks = 76, ColSideColors = cols, notecol = NULL, trace = "none", xlab = "Strains", 
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

#colnames(ccmn_norm_mets_good_old)[!is.na(metKEGGIDs_old)] <- metKEGGIDs_old[!is.na(metKEGGIDs_old)]
#colnames(ccmn_norm_mets_newData)[!is.na(metKEGGIDs_new)] <- metKEGGIDs_new[!is.na(metKEGGIDs_new)]
#colnames(ccmn_norm_mets_hybrid)[!is.na(metKEGGIDs_hybrid)] <- metKEGGIDs_hybrid[!is.na(metKEGGIDs_hybrid)]

####################################################################################################################################
###################################                      GOOD OLD DATA ANALYSIS                  ###################################
####################################################################################################################################
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/genePresAbs_functions.R")

####################
#
# Get strain cluster classification and obtain contingency tabs for doing Fisher tests
#
####################################################################################################################################

# Remove strains that don't match between genomic and metabolomic datasets. Change names of PA14 strains in 
# metabolomic dataset to make them coincide.

ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[-grep("M6075", rownames(ccmn_norm_mets_good_old)), ]
rownames(ccmn_norm_mets_good_old)[grep("W70322",rownames(ccmn_norm_mets_good_old))] <- gsub("W70322", 
                                                                                            "W70332", 
                                                                                            rownames(ccmn_norm_mets_good_old)[grep("W70322", 
                                                                                                                                   rownames(ccmn_norm_mets_good_old))])


rownames(ccmn_norm_mets_good_old)[1:3] <- paste("PA14.jbx", as.character(1:3), sep = "_")
rownames(ccmn_norm_mets_good_old)[4:6] <- paste("UCBPP.PA14", as.character(4:6), sep = "_")
gene_tab <- gene_tab[rownames(gene_enz_tab) %in% unique(gsub("\\_.*", rownames(ccmn_norm_mets_good_old), replacement = "")), ]
gene_tab_filt <- gene_tab_filt[rownames(gene_tab_filt) %in% unique(gsub("\\_.*", rownames(ccmn_norm_mets_good_old), replacement = "")), ]
gene_enz_tab <- gene_enz_tab[rownames(gene_enz_tab) %in% unique(gsub("\\_.*", rownames(ccmn_norm_mets_good_old), replacement = "")), ]
gene_enz_tab_filt <- gene_enz_tab_filt[rownames(gene_enz_tab_filt) %in% unique(gsub("\\_.*", rownames(ccmn_norm_mets_good_old), replacement = "")), ]
#gene_tab_filt <- gene_tab_filt[-grep("M6075", rownames(gene_tab_filt)), ]
#gene_enz_tab <- gene_enz_tab[-grep("M6075", rownames(gene_enz_tab)), ]
#gene_enz_tab_filt <- gene_enz_tab_filt[-grep("M6075", rownames(gene_enz_tab_filt)), ]

metClusts_oldGood <- getMetClusts(ccmn_norm_mets_good_old)
names(metClusts_oldGood)[1] <- "PA14.jbx"
metClusts_oldGood4G <- metClusts_oldGood


# Only 2 major clusters:
metClusts_oldGood <- as.factor(gsub("\\..*", "", metClusts_oldGood))
save(metClusts_oldGood, file = "metClusts_oldGood.RData")

# All genes
contTab_oldGood_allGenes <- getContTab(gene_tab, metClusts_oldGood)

# Removing PA14 and, after that, all across-strain common genes
contTab_oldGood_filt <- getContTab(gene_tab_filt, metClusts_oldGood)

# Only enzymatic genes (the ones with EC number)
contTab_oldGood_enz <- getContTab(gene_enz_tab, metClusts_oldGood)

# Enzimatic genes without PA14 and that are not shared across all strains
contTab_oldGood_filtEnz <- getContTab(gene_enz_tab_filt, metClusts_oldGood)

####################
#
# Fisher test between 4 groups of the tables obtained in previous step
#
####################################################################################################################################

# Load predictors FELLA OPLS-DA results obtained with swarming/metabolome supervised analysis
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis/tab_OPLSDAQuant.RData")



fish_oldGood_allGen <- doFisher(contTab_oldGood_allGenes, metClusts_oldGood)                                # 180 sign genes
fish_oldGood_filt <- doFisher(contTab_oldGood_filt, metClusts_oldGood)         # 139 sign genes
fish_oldGood_enz <- doFisher(contTab_oldGood_enz, metClusts_oldGood)                                        # 12 sign genes
fish_oldGood_filtEnz <- doFisher(contTab_oldGood_filtEnz, metClusts_oldGood)   # 9 sign genes
save(fish_oldGood_allGen, file = "fish_oldGood_allGen.RData")
save(fish_oldGood_filt, file = "fish_oldGood_filt.RData")
save(fish_oldGood_enz, file = "fish_oldGood_enz.RData")
save(fish_oldGood_filtEnz, file = "fish_oldGood_allGen.RData")

fish_oldGood_filtEnzKEGGIDs <- sub("\\).*", "", sub(".*\\(", "", fish_oldGood_filtEnz$sign.genes[, 1]))
fish_oldGood_filtEnzKEGGIDs <- gsub("EC ", replacement = "ec:", fish_oldGood_filtEnzKEGGIDs)
fish_oldGood_filtEnzKEGGIDsPaths <- sapply(fish_oldGood_filtEnzKEGGIDs, keggLink, target = "pathway")
gsub("path:ec|path:map", replacement = "pae",  fish_oldGood_filtEnzKEGGIDsPaths[[2]]) %in% tab_OPLSDAQuant$KEGG.id # --> those genes don't participate in any of the affected pathways

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
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

# Remove ambiguous metabolites from ccmn normalized data 

ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[, !is.na(dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)])]
colnames(ccmn_norm_mets_good_old) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
ccmn_norm_mets_good_oldAmbigRem <- colnames(ccmn_norm_mets_good_old)
names(ccmn_norm_mets_good_oldAmbigRem) <- ccmn_norm_mets_good_oldAmbigRem
ccmn_norm_mets_good_oldAmbigRem <- rmAmbig(ccmn_norm_mets_good_oldAmbigRem)
ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[, colnames(ccmn_norm_mets_good_old) %in% ccmn_norm_mets_good_oldAmbigRem]
colnames(ccmn_norm_mets_good_old) <- dictionary$`KEGG IDs`[match(colnames(ccmn_norm_mets_good_old), dictionary$Consensus)]

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

unique(as.vector(mannWhitPerGeneFilt_oldGood_allGenes$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_allGenes$as_mat))]) # 60/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_filt$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_filt$as_mat))])         # 60/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_enz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_enz$as_mat))])           # 57/86 mets
unique(as.vector(mannWhitPerGeneFilt_oldGood_filtEnz$as_mat)[!is.na(as.vector(mannWhitPerGeneFilt_oldGood_filtEnz$as_mat))])   # 57/86 mets

####################################################################################################################################
### Build a list with genes whose presence/absence individually was found to be correlated with     ################################
### differences in abundance of each one of the metabolites of the 3 sets of differential           ################################
### metabolites. Also with the set of predictors for swarming according to OPLS-DA quant            ################################
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
### to affect that metabolite, including the strains that belong to both of the clusters            ################################
### to whom the set of differental metabolites belong. Include a total of times those genes appear  ################################
### in each strain of the matrix.                                                                   ################################
####################################################################################################################################

geneMatMets_oldGood_allGenes <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_allGenes,
                                                  genePresAbsObjkt = gene_tab,
                                                  metClustObjkt = metClusts_oldGood4G)
geneMatMets_oldGood_filt <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_filt,
                                              genePresAbsObjkt = gene_tab_filt,
                                              metClustObjkt = metClusts_oldGood4G[2:length(metClusts_oldGood4G)])
geneMatMets_oldGood_enz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_enz,
                                             genePresAbsObjkt = gene_enz_tab,
                                             metClustObjkt = metClusts_oldGood4G)
geneMatMets_oldGood_filtEnz <- getGeneMatsPerMet(DGenesPerDMetsObjkt = diffGenesPerDiffMet_oldGood_filtEnz,
                                                 genePresAbsObjkt = gene_enz_tab_filt,
                                                 metClustObjkt = metClusts_oldGood4G[2:length(metClusts_oldGood4G)])

# Load OPLS-DA quant predictors for swarming: we're gonna filter the matrixes to keep only the ones of the mets that 
# predict swarming

load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis/OPLSDAQuantResultKEGGIDs.RData")

# See overlap diffmets and OPLS-DA
overlapOPLSDADiffMets <- list("diff mets C1 & C2" = diffMets_oldGood[, 1],
                              "Swarming Predictors (OPLS-DA)" = OPLSDAQuantResultKEGGIDs)

library(Vennerable)

venndiffMetsOPLSDA <- Venn(overlapOPLSDADiffMets)
tiff(filename = "venndiffMetsOPLSDA.tiff", height = 1400, width = 1800, res = 300)
plot(venndiffMetsOPLSDA, doWeights = T, type = "circles")
dev.off()
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

swarmRelated <- diffGenesPerMet_oldGood_filt$Total$`1&2`[names(diffGenesPerMet_oldGood_filt$Total$`1&2`) %in% OPLSDAQuantResultKEGGIDs]

swarmRelated <- unlist(swarmRelated)[-grep("group", unlist(swarmRelated))]

swarmRelated_enz <- diffGenesPerMet_oldGood_filtEnz$Total$`1&2`[names(diffGenesPerMet_oldGood_filtEnz$Total$`1&2`) %in% OPLSDAQuantResultKEGGIDs]

swarmRelated_enz <- unlist(swarmRelated_enz)[-grep("group", unlist(swarmRelated_enz))]

swarmRel_enz_ECs <- apply(dictEnzymes[dictEnzymes$Gene %in% unlist(swarmRelated), 34:ncol(dictEnzymes)], 1, function(x) unique(x)[!is.na(unique(x))][1])

sum(tab_OPLSDAQuant$KEGG.id %in% swarmRel_enz_ECs)
swarmRel_enz_ECs
# Build heatmap with presence/absence of the genes that affect diffmets between 2 major groups, with strains ordered according 
# to metabolic clusters.


C1_2Genes_old_filtEnz <- unique(unlist(diffGenesPerMet_oldGood_filtEnz$Total$`1&2`))[unique(unlist(diffGenesPerMet_oldGood_filtEnz$Total$`1&2`)) != "Any relationship between the gene and the metabolite"]
presAbsC1_2Genes_old_filtEnz <- gene_enz_tab_filt[, C1_2Genes_old_filtEnz]
clustOrdC1_2_old <- names(sort(metClusts_oldGood4G[1:length(metClusts_oldGood4G)]))
presAbsC1_2Genes_old_filtEnz <- t(presAbsC1_2Genes_old_filtEnz[clustOrdC1_2_old, ])
save(presAbsC1_2Genes_old_filtEnz, file = "presAbsC1_2Genes_old_filtEnz.RData")
save(C1_2Genes_old_filtEnz, file = "C1_2Genes_old_filtEnz.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/presAbsC1_2Genes_old_filtEnz.RData")

tiff("minedDiffGenes_old_filtEnz.tiff", width = 10000, height = 4000, units = "px", pointsize = 100)
heatmap.2(presAbsC1_2Genes_old_filtEnz, Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = c("blue", "red"), ColSideColors = cols[match(clustOrdC1_2_old, rownames(gene_tab))], notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Genes", main = "Presence/absence of diff genes", margins = c(6, 65), 
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

rownames(presAbsC1_2Genes_old_filtEnz)[dictEnzymes$ECnums[match(rownames(presAbsC1_2Genes_old_filtEnz), dictEnzymes$Gene)] %in% tab_OPLSDAQuant$KEGG.id]

dictEnzymes$ECnums[which(dictEnzymes$Gene == "leuB_1")]

rownames(presAbsC1_2Genes_old_filtEnz) %in% dictEnzymes$Gene
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