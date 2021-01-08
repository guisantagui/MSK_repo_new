setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/rnaSeqAnalysis")

if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(airway)) BiocManager::install("airway")
library(airway)
if(!require(DESeq2)) BiocManager::install("DESeq2")
library(DESeq2)
if(!require(gplots)) install.packages("gplots")
if(!require(ropls)) install.packages("ropls")
library(ropls)

rnaSeq <- readxl::read_xlsx("../../rnaSeq/ctsLenCore3.xlsx")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/rhamnMat.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/pcomp1AreaPercCircularity.RData")
prokkAnot <- read.csv("C:/Users/Guillem/Documents/PhD/comput/prokkAnnotStuff/untitled_folder/gene_presence_absence.csv")

rnaSeq

swarmData <- pcomp1Means*-1

genes <- rnaSeq$Gene

rhamnMat$strains <- gsub("PA14A|PA14B", "PA14", rhamnMat$strains)

rnaSeq <- rnaSeq[, 2:ncol(rnaSeq)]

rnaSeq <- as.data.frame(rnaSeq)


# Reorder columns to make PA14 be the first (to use it as reference for obtaining the log2fold changes)
rnaSeq <- rnaSeq[, c(6, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 17, 12, 13, 14, 15, 16, 18, 
                     19, 20, 21, 22, 28, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33)]


rhamn2cats <- rhamnMat$rhamn2cats[match(gsub(".-|(PA14).*", 
                                             colnames(rnaSeq), 
                                             rep = "\\1"),
                                        rhamnMat$strains)]

rhamn3cats <- rhamnMat$rhamn3cats[match(gsub(".-|(PA14).*", 
                                             colnames(rnaSeq), 
                                             rep = "\\1"),
                                        rhamnMat$strains)]

swarm <- swarmData[match(gsub(".-|(PA14).*", 
                              colnames(rnaSeq), 
                              rep = "\\1"),
                         names(swarmData))]

rhamn2cats
rhamnBin <- rep("rhamnPos", ncol(rnaSeq))
rhamnBin[rhamn2cats == 0] <- "rhamnNeg"
rhamnBin <- as.factor(rhamnBin)


# create DESeq object (for running DESeq2 pipeline)
rnaSeqColData <- data.frame(phenotype = rhamnBin)
rownames(rnaSeqColData) <- colnames(rnaSeq)

rnaSeqRounded <- round(rnaSeq)

rownames(rnaSeqRounded) <- genes 

DESeqRNASeq <- DESeqDataSetFromMatrix(rnaSeqRounded, rnaSeqColData, design = ~ phenotype)

# Filter rows with all zeros (it's core so none will be removed)
DESeqRNASeq <- DESeqRNASeq[colSums(counts(DESeqRNASeq)) > 1, ]
nrow(DESeqRNASeq)

# Log2 transform
DESeqRNASeqRlog <- rlog(DESeqRNASeq, blind=TRUE)

plot(log2(counts(DESeqRNASeq, normalized=F)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(DESeqRNASeqRlog)[,1:2],
     pch=16, cex=0.3)

# Dist matrix across strains 

sampleDistsRNASeq <- dist( t( assay(DESeqRNASeqRlog) ) )

if(!require(pheatmap)) install.packages("pheatmap")
library(pheatmap)
if(!require(RColorBrewer)) install.packages("RColorBrewer")
library(RColorBrewer)

sampleDistMatRNASeq <- as.matrix(sampleDistsRNASeq)
rownames(sampleDistMatRNASeq) <- paste(colnames(DESeqRNASeqRlog), DESeqRNASeqRlog$phenotype, sep="-")
colnames(sampleDistMatRNASeq) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatRNASeq,
         clustering_distance_rows=sampleDistsRNASeq,
         clustering_distance_cols=sampleDistsRNASeq,
         col=colors)


PCA_DESeqRNASeq <- plotPCA(DESeqRNASeqRlog, intgroup = "phenotype", returnData = T)

tiff(filename = "PCA_DESeq.tiff")
plotPCA(DESeqRNASeqRlog, intgroup = "phenotype")
dev.off()

diffAnalRNASeq <- DESeq(DESeqRNASeq)

resultsDiffAnalRNASeq <- results(diffAnalRNASeq, alpha = 0.01)

# Consider significative the ones with a p-value below 0.01 and a log2 change above 1 or below -1.
signDESeqResults <- as.data.frame(resultsDiffAnalRNASeq[resultsDiffAnalRNASeq$padj < .01, ])

signDESeqResults <- signDESeqResults[signDESeqResults$log2FoldChange > 1 |signDESeqResults$log2FoldChange < -1,]

signDESeqResults <- signDESeqResults[order(signDESeqResults$log2FoldChange), ]

View(signDESeqResults)

which(rownames(signDESeqResults) == "ibpB")

which(rownames(signDESeqResults) == "group_9868")

# Remove hypothetical proteins 
signDESeqResultsNoHypProt <- signDESeqResults[-grep("group", rownames(signDESeqResults)), ]

signDESeqResultsNoHypProt$Annotation <- prokkAnot$Annotation[match(rownames(signDESeqResultsNoHypProt), prokkAnot$Gene)]

View(signDESeqResultsNoHypProt)

save(signDESeqResultsNoHypProt, file = "signDESeqResultsNoHypProt.RData")
write.csv(signDESeqResultsNoHypProt, "signDESeqResultsNoHypProt.csv")

# PQS genes are quorum sensing genes and are found lower in rhamn non producers

pqsGenes <- cbind.data.frame(rnaSeq_log[, grep("pqs", colnames(rnaSeq_log))], rnaSeq_log$rhamn2cats)

table(resultsDiffAnalRNASeq$padj < .05)


which(rownames(signDESeqResultsNoHypProt) == "nicP_10")

which(rownames(signDESeqResultsNoHypProt) == "group_9868")

# pipeline by hand
rnaSeq <- as.data.frame(t(rnaSeq))

colnames(rnaSeq) <- genes

# There are still some genes that are not core (have zero counts)

rnaSeq <- rnaSeq[, -which(apply(rnaSeq, 2, function(x) sum(x == 0) > 0))]


log2(rnaSeq)
rnaSeq_log <- log2(rnaSeq)


rnaSeq_log$rhamn2cats <- rhamn2cats
rnaSeq_log$rhamn3cats <- rhamn3cats
rnaSeq_log$swarm <- swarm



pCompRnaSeq <- prcomp(rnaSeq_log[, 1:(ncol(rnaSeq_log)-3)])


pdf("PCArnaSeq_rhamn.pdf", height = 12, width = 12)
fviz_pca_ind(pCompRnaSeq, 
             col.ind = rnaSeq_log$rhamn2cats,
             gradient.cols = c("blue", "red"),
             legend.title = "rhamnolipid production")
dev.off()

pdf("PCArnaSeq_rhamn3cats.pdf", height = 12, width = 12)
fviz_pca_ind(pCompRnaSeq, 
             col.ind = rnaSeq_log$rhamn3cats,
             gradient.cols = c("blue", "green", "red"),
             legend.title = "rhamnolipid production")
dev.off()

pdf("PCArnaSeq_swarm.pdf", height = 12, width = 12)
fviz_pca_ind(pCompRnaSeq, 
             col.ind = rnaSeq_log$swarm,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "swarming score")
dev.off()

pdf("biplotRnaSeq.pdf", height = 12, width = 12)
fviz_pca_biplot(pCompRnaSeq)
dev.off()

OPLSDA_RNASeq <- opls(rnaSeq_log[, 1:(ncol(rnaSeq_log) - 3)],
                      as.factor(rnaSeq_log$rhamn2cats), 
                      predI = 1, 
                      orthoI = NA,
                      permI = 50)

head(OPLSDA_RNASeq@loadingMN[order(OPLSDA_RNASeq@loadingMN), ], 50)

OPLSDA_RNASeq@loadingMN[rownames(OPLSDA_RNASeq@loadingMN) == "acrA", ]

pVals <- c()
for(i in 1:(ncol(rnaSeq_log)-3)){
        mWhitTest <- wilcox.test(rnaSeq_log[rnaSeq_log$rhamn2cats == 0, i], 
                                 rnaSeq_log[rnaSeq_log$rhamn2cats == 1, i])$p.value
        pVals <- c(pVals, mWhitTest)
}
names(pVals) <- colnames(rnaSeq_log[1:(ncol(rnaSeq_log)-3)])

pValsAdj <- p.adjust(pVals, method = "BH")

sign <- pValsAdj[pValsAdj < 0.05]

signTab <- data.frame(gene = names(sign),
                      p.val = sign,
                      annotation = prokkAnot$Annotation[match(names(sign), prokkAnot$Gene)])

View(signTab)

wilcox.test(rnaSeq_log[rnaSeq_log$rhamn2cats == 0, grep("lasR", colnames(rnaSeq_log))], 
            rnaSeq_log[rnaSeq_log$rhamn2cats == 1, grep("lasR", colnames(rnaSeq_log))])$p.value

rnaSeqAnnot <- prokkAnot$Annotation[match(colnames(rnaSeq_log)[1:(ncol(rnaSeq_log)-3)], prokkAnot$Gene)]
names(rnaSeqAnnot) <- colnames(rnaSeq_log)[1:(ncol(rnaSeq_log)-3)]

quorumGenesRnaSeq <- rnaSeqAnnot[grep("quorum", tolower(rnaSeqAnnot))]

wilcox.test(rnaSeq_log[rnaSeq_log$rhamn2cats == 0, names(quorumGenesRnaSeq)[1]], 
            rnaSeq_log[rnaSeq_log$rhamn2cats == 1, names(quorumGenesRnaSeq)[1]])$p.value

wilcox.test(rnaSeq_log[rnaSeq_log$rhamn2cats == 0, names(quorumGenesRnaSeq)[2]], 
            rnaSeq_log[rnaSeq_log$rhamn2cats == 1, names(quorumGenesRnaSeq)[2]])$p.value

wilcox.test(rnaSeq_log[rnaSeq_log$rhamn2cats == 0, names(quorumGenesRnaSeq)[3]], 
            rnaSeq_log[rnaSeq_log$rhamn2cats == 1, names(quorumGenesRnaSeq)[3]])$p.value

data("airway")
se <- airway
se

se$dex <- relevel(se$dex, "untrt")
se$dex

round( colSums(assay(se)) / 1e6, 1 )

se$BioSample

countdata <- assay(se)
head(countdata, 3)

coldata <- colData(se)

(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~ cell + dex))

signDESeqResultsPQSPHN <- signDESeqResultsNoHypProt[grep("pqs|phn|las|rhl", tolower(rownames(signDESeqResultsNoHypProt))), ]

View(signDESeqResultsPQSPHN)

View(signDESeqResultsPQSPHN)







