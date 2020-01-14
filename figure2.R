####################################################################################################################
##### Figure 2 paper 
####################################################################################################################

setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/figure2")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/dictionary/dictionary.RData")
dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]

ccmnNormMets <- ccmn_norm_mets_good_old
colnames(ccmnNormMets) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
ccmnNormMets <- ccmnNormMets[, !is.na(colnames(ccmnNormMets))]

groups <- unique(gsub("_.*", replacement = "", rownames(ccmnNormMets)))

sapply(groups, function(x) sd(ccmnNormMets[grep(x, rownames(ccmnNormMets)), 1]))

SDs <- apply(ccmnNormMets, MARGIN = 2, FUN = function(y) sapply(groups, function(x) sd(y[grep(x, names(y))])))
tail(SDs, 10)
apply(SDs, 2, rank)[9, ]



doAnovas <- function(metMat){
        groups <- gsub("_.*", replacement = "", rownames(metMat))
        pVals <- c()
        for(i in 1:ncol(metMat)){
               data <- cbind.data.frame(metMat[, i], groups)
               colnames(data) <- c("abundance", "groups")
               res.aov <- aov(abundance ~ groups, data = data)
               summ <- summary(res.aov)
               pVal <- summ[[1]]$`Pr(>F)`[1]
               pVals <- c(pVals, pVal)
        }
        names(pVals) <- colnames(metMat)
        return(pVals)
}
doAnovas(ccmnNormMets)

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



if(!require(gplots)) install.packages("gplots")
library(gplots)
source("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/MSK_repo_new/diffMetAnal_functions.R")

ccmnNormMets_quant <- quantNorm(ccmnNormMets)

cols <- topo.colors(length(groups))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(ccmnNormMets))))
for(ii in 1:nrow(ccmnNormMets)) {
        selected <- which(groups == gsub("\\_.*", "", rownames(ccmnNormMets))[ii])
        Cols[ii] <- cols[selected]
}

heatMapNoRowNorm <- heatmap.2(as.matrix(t(ccmnNormMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
                              density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
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

tiff("heatmapCCMN.tiff", width = 10000, height = 4000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(ccmnNormMets)), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
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

dev.off()

ccmnNormMetsNoAmbig <- ccmnNormMets[, -grep("?", colnames(ccmnNormMets), fixed = T)]
ccmnNormMetsNoAmbigQuant <- quantNorm(ccmnNormMetsNoAmbig)

colOrder <- order.dendrogram(heatMapNoRowNorm$colDendrogram)
ccmnNormMetsNoAmbigQuant <- ccmnNormMetsNoAmbigQuant[colOrder, ]

tiff("heatmapCCMN_quant.tiff", width = 10000, height = 4000, units = "px", pointsize = 50)
heatmap.2(as.matrix(t(ccmnNormMetsNoAmbigQuant)), Rowv = T, Colv = F, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "row", 
          col = redgreen(75), breaks = 76, ColSideColors = Cols[match(rownames(ccmnNormMetsNoAmbigQuant),
                                                                      rownames(ccmnNormMets))], 
          notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 16), 
          cex.main = 20,
          keysize = 0.7,
          cexRow = 0.7,
          cexCol = 1.2,
          #colCol = colCols,
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

dev.off()

