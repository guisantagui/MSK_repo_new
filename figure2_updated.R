####################################################################################################################
##### Figure 2 paper 
####################################################################################################################

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure2")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/pcomp1AreaPercCircularity.RData")


rownames(ccmnNormMets)[grep(pattern = "W70322", x = rownames(ccmnNormMets))] <- paste("W70332", 
                                                                                      gsub(pattern = ".*_", 
                                                                                           replacement = "", 
                                                                                           rownames(ccmnNormMets)[grep(pattern = "W70322", 
                                                                                                                       x = rownames(ccmnNormMets))]), 
                                                                                      sep = "_")


groups <- unique(gsub("_.*", replacement = "", rownames(ccmnNormMets)))

sapply(groups, function(x) sd(ccmnNormMets[grep(x, rownames(ccmnNormMets)), 1]))

SDs <- apply(ccmnNormMets, MARGIN = 2, FUN = function(y) sapply(groups, function(x) sd(y[grep(x, names(y))])))
tail(SDs, 10)
apply(SDs, 2, rank)[9, ]



doAnovas <- function(metMat){
        groups <- gsub("_.*", replacement = "", rownames(metMat))
        pVals <- c()
        aovAllRes <- list()
        for(i in 1:ncol(metMat)){
               data <- cbind.data.frame(metMat[, i], groups)
               colnames(data) <- c("abundance", "groups")
               res.aov <- aov(abundance ~ groups, data = data)
               summ <- summary(res.aov)
               pVal <- summ[[1]]$`Pr(>F)`[1]
               pVals <- c(pVals, pVal)
               aovAllRes[[i]] <- res.aov
        }
        names(pVals) <- colnames(metMat)
        names(aovAllRes) <- colnames(metMat)
        return(list("p.values" = pVals, "allANOVA" = aovAllRes))
}
anovasStrains <- doAnovas(ccmnNormMets)


data <- cbind.data.frame(ccmnNormMets, groups)
colnames(data) <- make.names(colnames(data))

resANOVA <- aov(as.formula(paste(colnames(data)[1],"groups", sep = " ~ ")), data = data)

for(i in 1:length(anovasStrains$allANOVA)){
        plot(anovasStrains$allANOVA[[i]], 1, main = names(anovasStrains$allANOVA[i]))
}
plot(resANOVA, 1)

doGroupAnovas <- function(metMat){
        rep1 <- which(((1:28) + 2) %% 3 == 0)
        rep2 <- which(((1:28) + 1) %% 3 == 0)
        rep3 <- which((1:28) %% 3 == 0)
        groups <- c(rep("rep1", 28), 
                    rep("rep2", 28),
                    rep("rep3", 28))
        pVals <- c()
        for(i in 1:ncol(metMat)){
                abVec <- c(metMat[rep1, i],
                           metMat[rep2, i],
                           metMat[rep3, i])
                data <- cbind.data.frame(abVec, groups)
                colnames(data) <- c("abundance", "groups")
                res.aov <- kruskal.test(x = data$abundance,
                                        g = data$groups)
                #summ <- summary(res.aov)
                pVal <- res.aov$p.value
                pVals <- c(pVals, pVal)
        }
        names(pVals) <- colnames(metMat)
        return(pVals)
}

doGroupAnovas(ccmnNormMets)

# Create colorscale and remove sample 25 (it's very different to the replicates)
if(!require(gplots)) install.packages("gplots")
library(gplots)
if(!require(dendextend)) install.packages("dendextend")
library(dendextend)

source("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

ccmnNormMets_quant <- quantNorm(ccmnNormMets)

cols <- topo.colors(length(groups))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(ccmnNormMets))))
for(ii in 1:nrow(ccmnNormMets)) {
        selected <- which(groups == gsub("\\_.*", "", rownames(ccmnNormMets))[ii])
        Cols[ii] <- cols[selected]
}
ccmnNormMets <- ccmnNormMets[-25, ]
Cols <- Cols[-25]

names(pcomp1Means) <- gsub("W70322", "W70332", names(pcomp1Means))

pcomp1Filt <- pcomp1Means[match(gsub("\\_.*|(PA14).*", 
                                     rownames(ccmnNormMets), 
                                     rep = "\\1"), 
                                names(pcomp1Means))]
pcomp1Filt <- pcomp1Filt * -1



# PCA
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(readxl)) install.packages("readxl")
library(readxl)

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/swarmMeans.RData")

names(swarmMeans) <- gsub(".*_", replacement = "", names(swarmMeans))

swarmMeansFilt <- swarmMeans[match(gsub("\\_.*|(PA14).*", 
                                        rownames(ccmnNormMets), 
                                        rep = "\\1"), 
                                   names(swarmMeans))]
pcaMets <- prcomp(ccmnNormMets)
fviz_eig(pcaMets)


PCA_swarms <- fviz_pca_ind(pcaMets, 
                           col.ind = swarmMeansFilt,
                           gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
                           legend.title = "Log Swarming Area")


tiff("pcaSwarmMets.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = swarmMeansFilt,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Log Swarming Area")
dev.off()

tiff("pcaSwarmMets_pcomp1.tiff", res = 300, height = 5000, width = 5000)
fviz_pca_ind(pcaMets, 
             col.ind = pcomp1Filt,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Swarming (P.Comp. 1 Max. Length and circularity)")
dev.off()

pdf("pcaSwarmMets_pcomp1.pdf", height = 20, width = 20)
fviz_pca_ind(pcaMets, 
             col.ind = pcomp1Filt,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Swarming (P.Comp. 1 Max Length and circularity)")
dev.off()

mets <-read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/newData/metaboliteTable.xlsx")

mets$mutant[grep(pattern = "W70322", x = mets$mutant)] <- paste("W70332", 
                                                                gsub(pattern = ".*_", 
                                                                     replacement = "", 
                                                                     rownames(ccmnNormMets)[grep(pattern = "W70322", 
                                                                                                 x = mets$mutant)]), 
                                                                sep = "_")

batches <- c()
batchNames <- c()
for(i in seq_along(groups)){
        subMat <- mets[grep(groups[i], mets$mutant), ]
        reps <- unique(subMat$Replicate)
        for(j in seq_along(reps)){
                batch <- unique(subMat$batch[grep(reps[j], subMat$Replicate)])
                batches <- c(batches, batch)
                batchNames <- c(batchNames, 
                                paste(groups[i],
                                      as.character(j),
                                      sep = "_"))
        }
}
names(batches) <- batchNames
batches <- batches[-25]

tiff("pcaBatchMets.tiff", res = 300, height = 3000, width = 3000)
fviz_pca_ind(pcaMets,
             col.ind = batches,
             legend.title = "Batches")
dev.off()




# Do heatmaps

heatMapNoRowNorm <- heatmap.2(as.matrix(t(ccmnNormMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                              density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
                              col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
                              ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 10), 
                              cex.main = 10,
                              keysize = 0.7,
                              cexRow = 0.7,
                              cexCol = 1.2,
                              scale = "none",
                              #colCol = colCols,
                              #cellnote = round(as.matrix(t(ccmnNormMets)), 2),
                              #notecex = 0.7,
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

tiff(filename = "noNormHeatmap.tiff", width = 2000, height = 2000)
heatmap.2(as.matrix(t(ccmnNormMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = colorpanel(75, low = "red", mid = "black", high = "green"), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          scale = "none",
          #colCol = colCols,
          #cellnote = round(as.matrix(t(ccmnNormMets)), 2),
          #notecex = 0.7,
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

plot(hclust(dist(ccmnNormMets, method = "euclidean"), method = "ward.D"))



colThick <- set(heatMapNoRowNorm$colDendrogram, "branches_lwd", 5)
rowThick <- set(heatMapNoRowNorm$rowDendrogram, "branches_lwd", 5)

tiff("heatmapCCMN.tiff", width = 10000, height = 6000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(ccmnNormMets)), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          scale = "row",
          Colv = colThick,
          Rowv = rowThick,
          #colCol = colCols,
          #cellnote = round(as.matrix(t(ccmnNormMets)), 2),
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

dev.off()


# Load the solved ambig vector 
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)


ccmnNormMetsNoAmbig <- ccmnNormMets[, -grep("_?", colnames(ccmnNormMets), fixed = T)]
ccmnNormMetsNoAmbigQuant <- quantNorm(ccmnNormMetsNoAmbig)

metNames <- keggList("compound")[match(dictionary$`KEGG IDs`[match(colnames(ccmnNormMetsNoAmbigQuant), 
                                                                   dictionary$Consensus)], 
                                       gsub("cpd:", 
                                            "", 
                                            names(keggList("compound"))))]

metNames <- sapply(metNames, function(x) unlist(strsplit(x, ";"))[1])

metNames[17] <- "AICAR"

colnames(ccmnNormMetsNoAmbigQuant)[!is.na(metNames)] <- metNames[!is.na(metNames)]
colnames(ccmnNormMetsNoAmbig)[!is.na(metNames)] <- metNames[!is.na(metNames)]

# Change names to the ones Kyu Rhee said
KR_metNames <- readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/KRNewMetNames/metabolitesName-Xavier lab-KR.xlsx")

KR_metNames <- as.data.frame(KR_metNames)

colnames(ccmnNormMetsNoAmbig) <- KR_metNames$newName[match(colnames(ccmnNormMetsNoAmbig), KR_metNames$oldName)]
colnames(ccmnNormMetsNoAmbigQuant) <- KR_metNames$newName[match(colnames(ccmnNormMetsNoAmbigQuant), KR_metNames$oldName)]


# Export curated dataset
write.csv(ccmnNormMetsNoAmbig, "metsCCMNNormCurated.csv")

#colOrder <- order.dendrogram(heatMapNoRowNorm$colDendrogram)
#ccmnNormMetsNoAmbigQuant <- ccmnNormMetsNoAmbigQuant[colOrder, ]

colCols <- Cols

heatMapRowNorm <- heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
                            density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
                            col = redgreen(75), breaks = 76, ColSideColors = Cols[match(rownames(ccmnNormMetsNoAmbigQuant),
                                                                                        rownames(ccmnNormMets))], 
                            notecol = NULL, trace = "none", xlab = "Strains", 
                            ylab = "Metabolites", 
                            #main = "CCMN normalized",
                            margins = c(10, 16), 
                            cex.main = 20,
                            keysize = 0.7,
                            cexRow = 0.7,
                            cexCol = 1.2,
                            colCol = colCols[match(rownames(ccmnNormMetsNoAmbigQuant),
                                                   rownames(ccmnNormMets))],
                            cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
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

#heatMapRowNorm <- heatmap.2(as.matrix(t(ccmnNormMetsNoAmbig)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
#                           density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
#                            col = redgreen(75), breaks = 76, ColSideColors = Cols[match(rownames(ccmnNormMetsNoAmbigQuant),
#                                                                                        rownames(ccmnNormMets))], 
#                            notecol = NULL, trace = "none", xlab = "Strains", 
#                            ylab = "Metabolites", 
#                            #main = "CCMN normalized",
#                            margins = c(10, 16), 
#                            cex.main = 20,
#                            keysize = 0.7,
#                            cexRow = 0.7,
#                            cexCol = 1.2,
#                            colCol = colCols[match(rownames(ccmnNormMetsNoAmbigQuant),
#                                                   rownames(ccmnNormMets))],
#                            cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
#                            notecex = 0.7,
#                            key.xtickfun=function() {
#                                    cex <- par("cex")*par("cex.axis")
#                                    side <- 1
#                                    line <- 0
#                                    col <- par("col.axis")
#                                    font <- par("font.axis")
#                                    mtext("low", side=side, at=0, adj=0,
#                                          line=line, cex=cex, col=col, font=font)
#                                    mtext("high", side=side, at=1, adj=1,
#                                          line=line, cex=cex, col=col, font=font)
#                                    return(list(labels=FALSE, tick=FALSE))
#                            })

colThick <- set(heatMapNoRowNorm$colDendrogram, "branches_lwd", 1)
rowThickQuant <- set(heatMapRowNorm$rowDendrogram, "branches_lwd", 1)

tiff("heatmapCCMN_quant.tiff", width = 10000, height = 6000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols#[match(rownames(ccmnNormMetsNoAmbigQuant),
                                                               #      rownames(ccmnNormMets))]
          , 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", 
          #main = "CCMN normalized",
          Colv = colThick,
          Rowv = rowThickQuant,
          margins = c(15, 25), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 1.5,
          cexCol = 2,
          #colCol = colCols#[match(rownames(ccmnNormMetsNoAmbigQuant),
                          #       rownames(ccmnNormMets))]
          #,
          #cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
          notecex = 1.2,
          key = T,
          key.title = "Metabolite Levels",
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(4.5, 2.5, 10, 0)),
          key.xtickfun=function() {
                  cex <- par("cex")*2*par("cex.axis")
                  side <- 1
                  line <- 1
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=F, tick=F))
          }
)

dev.off()

# Create color vectors for rhamnolipid production and swarming

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/rhamnMat.RData")

gradient <- colorRampPalette(colors = c("#003CFF", "#66FF00", "#FF0000"))
swarmCols_PC1 <- gradient(nrow(ccmnNormMetsNoAmbigQuant))[cut(pcomp1Filt + abs(min(pcomp1Filt)), 83)] 
rhamn2CatFilt <- rhamnMat$rhamn2cats[match(gsub("\\_.*", "", rownames(ccmnNormMets)), 
                                           rhamnMat$strains)]
rhamn3CatFilt <- rhamnMat$rhamn3cats[match(gsub("\\_.*", "", rownames(ccmnNormMets)), 
                                           rhamnMat$strains)]
rhamn2CatCols <- rep("red", length(rhamn2CatFilt))
rhamn2CatCols[rhamn2CatFilt == 0] <- "blue"

rhamn3CatCols <- rep("red", length(rhamn3CatFilt))
rhamn3CatCols[rhamn3CatFilt == 0] <- "blue"
rhamn3CatCols[rhamn3CatFilt == 1] <- "green"

pdf("heatmapCCMN_quant.pdf", width = 60, height = 40, pointsize = 25)
heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = rhamn3CatCols#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #      rownames(ccmnNormMets))]
          , 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", 
          cex.axis = 30,
          #main = "CCMN normalized",
          Colv = colThick,
          Rowv = rowThickQuant,
          margins = c(15, 40), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 1.5,
          cexCol = 2,
          colCol = swarmCols_PC1#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #       rownames(ccmnNormMets))]
          ,
          #cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
          notecex = 1.2,
          key = T,
          key.title = "Metabolite Levels",
          key.par=list(mgp=c(1.5, 0.5, 0),
                       mar=c(2.5, 0.5, 8, 2)),
          key.xtickfun=function() {
                  cex <- par("cex")*2*par("cex.axis")
                  side <- 1
                  line <- 1
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=F, tick=F))
          }
)

dev.off()


pdf("heatmapCCMN_quantScaleAlter.pdf", width = 13, height = 11, pointsize = 2.5)
#par(oma = c(0, 0, 0, 0));
heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(12), 
          breaks = quantile(as.vector(t(ccmnNormMetsNoAmbigQuant)),
                            probs = c(0, 0.025, 0.05, 0.1, 0.225, 0.375, 0.5, 0.625, 0.775, 0.9, 0.95, 0.975, 1)), 
          ColSideColors = rhamn3CatCols#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #      rownames(ccmnNormMets))]
          , 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", 
          cex.axis = 300,
          #main = "CCMN normalized",
          Colv = colThick,
          Rowv = rowThickQuant,
          margins = c(13, 42), 
          cex.main = 20,
          keysize = 0.5,
          cexRow = 2.4,
          cexCol = 2.2,
          colCol = swarmCols_PC1#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #       rownames(ccmnNormMets))]
          ,
          #cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
          notecex = 2,
          key = T,
          key.title = "Metabolite Levels",
          key.par=list(mgp=c(2.0, 0.5, 0),
                       mar=c(3.5, 0.5, 8, 2)),
          key.xtickfun=function() {
                  cex <- par("cex")*2*par("cex.axis")
                  side <- 1
                  line <- 1
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=F, tick=F))
          }
)

dev.off()


pdf("heatmapCCMN_quantScaleAlter.pdf", width = 7, height = 7, pointsize = 2.5)
#par(oma = c(0, 0, 0, 0));
heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(12), 
          breaks = quantile(as.vector(t(ccmnNormMetsNoAmbigQuant)),
                            probs = c(0, 0.025, 0.05, 0.1, 0.225, 0.375, 0.5, 0.625, 0.775, 0.9, 0.95, 0.975, 1)), 
          ColSideColors = rhamn3CatCols#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #      rownames(ccmnNormMets))]
          , 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", 
          cex.axis = 300,
          #main = "CCMN normalized",
          Colv = colThick,
          Rowv = rowThickQuant,
          margins = c(1.2, 5.1), 
          cex.main = 20,
          keysize = 0.5,
          cexRow = 2.2,
          cexCol = 2.0,
          colCol = swarmCols_PC1#[match(rownames(ccmnNormMetsNoAmbigQuant),
          #       rownames(ccmnNormMets))]
          ,
          #cellnote = round(as.matrix(t(ccmnNormMets_quant)), 2),
          notecex = 2,
          key = T,
          key.title = "Metabolite Levels",
          key.par=list(mgp=c(2.0, 0.5, 0),
                       mar=c(3.5, 0.5, 8, 2)),
          key.xtickfun=function() {
                  cex <- par("cex")*2*par("cex.axis")
                  side <- 1
                  line <- 1
                  col <- par("col.axis")
                  font <- par("font.axis")
                  mtext("low", side=side, at=0, adj=0,
                        line=line, cex=cex, col=col, font=font)
                  mtext("high", side=side, at=1, adj=1,
                        line=line, cex=cex, col=col, font=font)
                  return(list(labels=F, tick=F))
          }
)

dev.off()

# Heatmaps with correlation as distance measure. 

heatMapNoRowNorm <- heatmap.2(as.matrix(t(ccmnNormMets)), Rowv = T, distfun = function(x) as.dist(cor(t(x))), 
                              density.info = "none", hclust = function(x) hclust(x, method = "complete"), dendrogram = "both", 
                              col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
                              ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
                              cex.main = 20,
                              keysize = 0.7,
                              cexRow = 0.7,
                              cexCol = 1.2,
                              scale = "row",
                              #colCol = colCols,
                              cellnote = round(as.matrix(t(ccmnNormMets)), 2),
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

colThick <- set(heatMapNoRowNorm$colDendrogram, "branches_lwd", 5)
rowThick <- set(heatMapNoRowNorm$rowDendrogram, "branches_lwd", 5)

tiff("heatmapCCMN.tiff", width = 10000, height = 6000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(ccmnNormMets)), distfun = function(x) as.dist(cor(t(x))), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          scale = "row",
          Colv = colThick,
          Rowv = rowThick,
          #colCol = colCols,
          cellnote = round(as.matrix(t(ccmnNormMets)), 2),
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

dev.off()

