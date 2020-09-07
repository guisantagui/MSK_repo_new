setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/KerryBoylePaper")

metabolomics <- readxl::read_xlsx("metabolomics.xlsx", sheet = 4)

metabolomics <- as.data.frame(metabolomics)

rownames(metabolomics) <- metabolomics[, 1]

metabolomics <- metabolomics[, -c(1:5, ncol(metabolomics))]

View(metabolomics)

library(gplots)

metabolomics[is.na(metabolomics)] <- 0


View(metabolomics)


heatmap.2(as.matrix(metabolomics), Rowv = T, distfun = function(x) dist(x, method = "euclidean"), 
          density.info = "none", hclust = function(x) hclust(x, method = "ward.D"), dendrogram = "both", 
          col = redgreen(75), breaks = 76, notecol = NULL, trace = "none", xlab = "Strains", 
          ylab = "Metabolites", main = "CCMN normalized", margins = c(10, 10), 
          cex.main = 10,
          keysize = 0.7,
          cexRow = 0.5,
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
