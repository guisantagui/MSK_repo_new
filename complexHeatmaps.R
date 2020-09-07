setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/swarmAnalysis")

# Load data

load(file = "geneEnzTabSignOPLSDA.RData")
load("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")
load("pcomp1AreaPercCircularity.RData")
swarmMeans <- pcomp1Means
names(swarmMeans) <- gsub("W70322", replacement = "W70332", names(swarmMeans))
swarmMeans <- swarmMeans[match(rownames(geneEnzTabSignOPLSDA), names(swarmMeans))]

# Build a version of previous heatmap with pathway information
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

if(!require(devtools)) install.packages("devtools")
library(devtools)

if(!require(ComplexHeatmap)) install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

heatmapSignEnzsOPLSDA <- Heatmap(t(as.matrix(geneEnzTabSignOPLSDA)), 
                                 clustering_method_rows = "ward.D", 
                                 cluster_columns = F, 
                                 show_heatmap_legend = F)

OPLSDASignEnzsECnums <- dictEnzymes$ECnums[match(colnames(geneEnzTabSignOPLSDA)[row_order(heatmapSignEnzsOPLSDA)], 
                                                 dictEnzymes$Gene)]
pathsSignEnzs <- sapply(OPLSDASignEnzsECnums, keggLink, target = "pathway")
allPathsSignEnzs <- gsub("map", replacement = "pau", unique(unlist(pathsSignEnzs)))[gsub("map", 
                                                                                         replacement = "pau", 
                                                                                         unique(unlist(pathsSignEnzs))) %in% keggLink(target = "pathway", 
                                                                                                                                      source = "pau")]
annotMat_pathsSignEnzs <- matrix(data = rep("0", length(pathsSignEnzs)*length(allPathsSignEnzs)),
                                 ncol = length(allPathsSignEnzs), 
                                 nrow = length(pathsSignEnzs), 
                                 dimnames = list(paste("ec:", OPLSDASignEnzsECnums, sep = ""),
                                                 allPathsSignEnzs))
for(i in seq_along(pathsSignEnzs)){
        idx <- match(gsub("path:pau", 
                          replacement = "", 
                          allPathsSignEnzs), 
                     gsub("path:map|path:ec", 
                          replacement = "",
                          pathsSignEnzs[[i]]))
        idx <- idx[!is.na(idx)]
        annotMat_pathsSignEnzs[i, idx] <- "1"
}

for(i in seq_along(pathsSignEnzs)){
        idx <- gsub("path:pau", 
                    replacement = "", 
                    allPathsSignEnzs) %in% 
                gsub("path:map|path:ec", 
                     replacement = "",
                     pathsSignEnzs[[i]])
        idx <- as.character(as.numeric(idx))
        annotMat_pathsSignEnzs[i, ] <- idx
}

as.character(annotMat_pathsSignEnzs)
colnames(annotMat_pathsSignEnzs)
annot <- rowAnnotation("Purine metabolism" = annotMat_pathsSignEnzs[, 1],
                       #pau01100 = annotMat_pathsSignEnzs[, 2],
                       "Lipid metabolism" = annotMat_pathsSignEnzs[, 3],
                       "Biotin metabolism" = annotMat_pathsSignEnzs[, 4],
                       "Tyr metabolism" = annotMat_pathsSignEnzs[, 5],
                       "Monobactam biosynthesis" = annotMat_pathsSignEnzs[, 6],
                       "Lys biosynthesis" = annotMat_pathsSignEnzs[, 7],
                       #pau01110 = annotMat_pathsSignEnzs[, 8],
                       #pau01120 = annotMat_pathsSignEnzs[, 9],
                       #pau01130 = annotMat_pathsSignEnzs[, 10],
                       "Arg and Pro metabolism" = annotMat_pathsSignEnzs[, 11],
                       "Gly, Ser and Thr metabolism" = annotMat_pathsSignEnzs[, 12],
                       "Fatty acid degradation" = annotMat_pathsSignEnzs[, 13],
                       "Caprolactam degradation" = annotMat_pathsSignEnzs[, 14],
                       "Glyoxylate and dicarboxylate metabolism" = annotMat_pathsSignEnzs[, 15],
                       "AA-tRNA biosynthesis" = annotMat_pathsSignEnzs[, 16],
                       "Pentose phosphate pathway" = annotMat_pathsSignEnzs[, 17],
                       "Folate biosynthesis" = annotMat_pathsSignEnzs[, 18],
                       "Amino sugar and nucleotide sugar metabolism" = annotMat_pathsSignEnzs[, 19],
                       "Galactose metabolism" = annotMat_pathsSignEnzs[, 20],
                       "Sulfur metabolism" = annotMat_pathsSignEnzs[, 21],
                       "Propanoate metabolism" = annotMat_pathsSignEnzs[, 22],
                       "Cyanoamino acid metabolism" = annotMat_pathsSignEnzs[, 23],
                       "Starch and sucrose metabolism" = annotMat_pathsSignEnzs[, 24],
                       "Benzoate degradation" = annotMat_pathsSignEnzs[, 25],
                       "Vitamin B6 metabolism" = annotMat_pathsSignEnzs[, 26],
                       "Phe metabolism" = annotMat_pathsSignEnzs[, 27],
                       col = list(pau00230 = c("1" = "red", "0" = "blue"),
                                  #pau01100 = c("1" = "red", "0" = "blue"),
                                  pau00061 = c("1" = "red", "0" = "blue"),
                                  pau00780 = c("1" = "red", "0" = "blue"),
                                  pau00350 = c("1" = "red", "0" = "blue"),
                                  pau00261 = c("1" = "red", "0" = "blue"),
                                  pau00300 = c("1" = "red", "0" = "blue"),
                                  #pau01110 = c("1" = "red", "0" = "blue"),
                                  #pau01120 = c("1" = "red", "0" = "blue"),
                                  #pau01130 = c("1" = "red", "0" = "blue"),
                                  pau00330 = c("1" = "red", "0" = "blue"),
                                  pau00260 = c("1" = "red", "0" = "blue"),
                                  pau00071 = c("1" = "red", "0" = "blue"),
                                  pau00930 = c("1" = "red", "0" = "blue"),
                                  pau00630 = c("1" = "red", "0" = "blue"),
                                  pau00970 = c("1" = "red", "0" = "blue"),
                                  pau00030 = c("1" = "red", "0" = "blue"),
                                  pau00790 = c("1" = "red", "0" = "blue"),
                                  pau00520 = c("1" = "red", "0" = "blue"),
                                  pau00052 = c("1" = "red", "0" = "blue"),
                                  pau00920 = c("1" = "red", "0" = "blue"),
                                  pau00640 = c("1" = "red", "0" = "blue"),
                                  pau00460 = c("1" = "red", "0" = "blue"),
                                  pau00500 = c("1" = "red", "0" = "blue"),
                                  pau00362 = c("1" = "red", "0" = "blue"),
                                  pau00750 = c("1" = "red", "0" = "blue"),
                                  pau00360 = c("1" = "red", "0" = "blue")),
                       simple_anno_size = unit(0.12, "cm"))

annot <- rowAnnotation("gene metabolic role" = cbind("Purine metabolism" = as.numeric(annotMat_pathsSignEnzs[, 1]),
                                                     #pau01100 = annotMat_pathsSignEnzs[, 2],
                                                     "Lipid metabolism" = as.numeric(annotMat_pathsSignEnzs[, 3]),
                                                     "Biotin metabolism" = as.numeric(annotMat_pathsSignEnzs[, 4]),
                                                     "Tyr metabolism" = as.numeric(annotMat_pathsSignEnzs[, 5]),
                                                     "Monobactam biosynthesis" = as.numeric(annotMat_pathsSignEnzs[, 6]),
                                                     "Lys biosynthesis" = as.numeric(annotMat_pathsSignEnzs[, 7]),
                                                     #pau01110 = annotMat_pathsSignEnzs[, 8],
                                                     #pau01120 = annotMat_pathsSignEnzs[, 9],
                                                     #pau01130 = annotMat_pathsSignEnzs[, 10],
                                                     "Arg and Pro metabolism" = as.numeric(annotMat_pathsSignEnzs[, 11]),
                                                     "Gly, Ser and Thr metabolism" = as.numeric(annotMat_pathsSignEnzs[, 12]),
                                                     "Fatty acid degradation" = as.numeric(annotMat_pathsSignEnzs[, 13]),
                                                     "Caprolactam degradation" = as.numeric(annotMat_pathsSignEnzs[, 14]),
                                                     "Glyoxylate and dicarboxylate metabolism" = as.numeric(annotMat_pathsSignEnzs[, 15]),
                                                     "AA-tRNA biosynthesis" = as.numeric(annotMat_pathsSignEnzs[, 16]),
                                                     "Pentose phosphate pathway" = as.numeric(annotMat_pathsSignEnzs[, 17]),
                                                     "Folate biosynthesis" = as.numeric(annotMat_pathsSignEnzs[, 18]),
                                                     "Amino sugar and nucleotide sugar metabolism" = as.numeric(annotMat_pathsSignEnzs[, 19]),
                                                     "Galactose metabolism" = as.numeric(annotMat_pathsSignEnzs[, 20]),
                                                     "Sulfur metabolism" = as.numeric(annotMat_pathsSignEnzs[, 21]),
                                                     "Propanoate metabolism" = as.numeric(annotMat_pathsSignEnzs[, 22]),
                                                     "Cyanoamino acid metabolism" = as.numeric(annotMat_pathsSignEnzs[, 23]),
                                                     "Starch and sucrose metabolism" = as.numeric(annotMat_pathsSignEnzs[, 24]),
                                                     "Benzoate degradation" = as.numeric(annotMat_pathsSignEnzs[, 25]),
                                                     "Vitamin B6 metabolism" = as.numeric(annotMat_pathsSignEnzs[, 26]),
                                                     "Phe metabolism" = as.numeric(annotMat_pathsSignEnzs[, 27])),
                       simple_anno_size = unit(0.3, "cm"),
                       show_legend = F)

swarmAnnot <- HeatmapAnnotation("swarming" = swarmMeans*-1)


heatmapSignEnzsOPLSDA <- Heatmap(t(as.matrix(geneEnzTabSignOPLSDA)), 
                                 clustering_method_rows = "ward.D", 
                                 cluster_columns = F, 
                                 show_heatmap_legend = F,
                                 right_annotation = annot,
                                 top_annotation = swarmAnnot)
tiff("geneEnzTabSignPaths.tiff", width = 10, height = 10, units = "in", res = 300,pointsize = 8)
heatmapSignEnzsOPLSDA
dev.off()

ha = rowAnnotation(foo = anno_points(runif(10), ylim = c(0, 1),
                                     width = unit(2, "cm"),
                                     axis_param = list(
                                             side = "bottom",
                                             at = c(0, 0.5, 1), 
                                             labels = c("zero", "half", "one"),
                                             labels_rot = 45
                                     ))
)

anno = anno_points(annotMat_pathsSignEnzs)
draw(anno, test = "anno_points")
anno = anno_points(matrix(runif(20), nc = 2), pch = 1:2)
draw(anno, test = "matrix")
