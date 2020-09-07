setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mithrilAnalysis")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/unNormMetImpData.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/statMetricsSwarm.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/statMetricsRhamn2Cats.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/statMetricsRhamn3Cats.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/swarm_signMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamn2Cats_signMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamn3Cats_signMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat.RData")

dim(ccmnNormMets)

rownames(ccmnNormMets) <- gsub("W70322", "W70332", rownames(ccmnNormMets))
rownames(unNormMetImpData) <- gsub("W70322", "W70332", rownames(unNormMetImpData))

ccmnNormMetsNoAmbig <- ccmnNormMets[, -grep("_?", colnames(ccmnNormMets), fixed = T)]

noNormMetsNoAmbig <- unNormMetImpData[, -grep("_?", colnames(unNormMetImpData), fixed = T)]

defNames <- dictionary$definitiveNames[match(colnames(ccmnNormMetsNoAmbig), dictionary$Consensus)]

defNamesNoNorm <- dictionary$definitiveNames[match(colnames(noNormMetsNoAmbig),
                                                   dictionary$Consensus)]

colnames(ccmnNormMetsNoAmbig)[!is.na(defNames)] <- defNames[!is.na(defNames)]

colnames(noNormMetsNoAmbig)[!is.na(defNamesNoNorm)] <- defNamesNoNorm[!is.na(defNamesNoNorm)]

keggIDs <- dictionary$`KEGG IDs`[match(colnames(ccmnNormMetsNoAmbig),
                                       dictionary$definitiveNames)]

keggIDsNoNorm <- dictionary$`KEGG IDs`[match(colnames(noNormMetsNoAmbig),
                                             dictionary$definitiveNames)]

mithFormat <- as.data.frame(t(ccmnNormMetsNoAmbig))

mithFormat <- cbind.data.frame(rownames(mithFormat),
                               keggIDs, 
                               mithFormat)

mithFormatNoNorm <- as.data.frame(t(noNormMetsNoAmbig))

kruskRes <- apply(mithFormatNoNorm, 1, function(x) kruskal.test(unlist(x), 
                                                                as.factor(gsub("\\_.*",
                                                                               "",
                                                                               colnames(mithFormatNoNorm))))$p.value)


kruskal.test(unlist(mithFormatNoNorm[1, ]), as.factor(gsub("\\_.*", "", colnames(mithFormatNoNorm))))

mithFormatNoNorm

mithFormatNoNorm <- cbind.data.frame(rownames(mithFormatNoNorm),
                                     keggIDsNoNorm, 
                                     mithFormatNoNorm)

colnames(mithFormat)[1:2] <- c("BIOCHEMICAL", "KEGG")

colnames(mithFormatNoNorm)[1:2] <- c("BIOCHEMICAL", "KEGG")

# Remove non changing mets from un normalized data

mithFormatNoNorm <- mithFormatNoNorm[- which(p.adjust(kruskRes, method = "BH") > 0.05), ]

statMetsSwarm <- statMetrics[-grep("_?", statMetrics$metabolites,
                                   fixed = T), ]
statMetsRhamn2 <- statMetricsRhamn2Cats[-grep("_?", statMetricsRhamn2Cats$metabolites,
                                              fixed = T), ]
statMetsRhamn3 <- statMetricsRhamn3Cats[-grep("_?", statMetricsRhamn3Cats$metabolites,
                                              fixed = T), -3]

defNamesStatMats <- dictionary$definitiveNames[match(statMetsSwarm$metabolites, 
                                                     dictionary$Consensus)]

defNamesStatMats[is.na(defNamesStatMats)] <- as.character(statMetsSwarm$metabolites[is.na(defNamesStatMats)])

statMetsSwarm$metabolites <- defNamesStatMats
statMetsRhamn2$metabolites <- defNamesStatMats
statMetsRhamn3$metabolites <- defNamesStatMats




RhamnPosMet <- ccmnNormMetsNoAmbig[rhamnSwarmMat$rhamn2cats[match(gsub("\\_.*|(PA14).*", "\\1", rownames(ccmnNormMetsNoAmbig)), 
                                                                  rhamnSwarmMat$strain)] == 1, ]

RhamnNegMet <- ccmnNormMetsNoAmbig[rhamnSwarmMat$rhamn2cats[match(gsub("\\_.*|(PA14).*", "\\1", rownames(ccmnNormMetsNoAmbig)), 
                                                                  rhamnSwarmMat$strain)] == 0, ]

RhamnPosMetNoNorm <- noNormMetsNoAmbig[rhamnSwarmMat$rhamn2cats[match(gsub("\\_.*|(PA14).*", "\\1", rownames(noNormMetsNoAmbig)), 
                                                                rhamnSwarmMat$strain)] == 1, ]

RhamnNegMetNoNorm <- noNormMetsNoAmbig[rhamnSwarmMat$rhamn2cats[match(gsub("\\_.*|(PA14).*", "\\1", rownames(noNormMetsNoAmbig)), 
                                                                rhamnSwarmMat$strain)] == 0, ]

RhamnPosMeds <- apply(RhamnPosMet, 2, median)
RhamnNegMeds <- apply(RhamnNegMet, 2, median)

RhamnPosMedsNoNorm <- apply(RhamnPosMetNoNorm, 2, median)
RhamnNegMedsNoNorm <- apply(RhamnNegMetNoNorm, 2, median)

logFCRhamn2 <- c()
for(i in 1:length(RhamnPosMeds)){
        logFC <- log((RhamnNegMeds[i]/RhamnPosMeds[i]), base = 2)
        logFCRhamn2 <- c(logFCRhamn2, logFC)
}
names(logFCRhamn2) <- names(RhamnPosMeds)

logFCRhamn2NoNorm <- c()
for(i in 1:length(RhamnPosMedsNoNorm)){
        logFC <- log((RhamnNegMedsNoNorm[i]/RhamnPosMedsNoNorm[i]), base = 2)
        logFCRhamn2NoNorm <- c(logFCRhamn2NoNorm, logFC)
}
names(logFCRhamn2NoNorm) <- names(RhamnPosMedsNoNorm)

mannWhitNoNorm <- c()
for(i in 1:ncol(RhamnPosMetNoNorm)){
        pval <- wilcox.test(RhamnPosMetNoNorm[, i],
                            RhamnNegMetNoNorm[, i])$p.value
        mannWhitNoNorm <- c(mannWhitNoNorm, pval)
}
names(mannWhitNoNorm) <- colnames(RhamnPosMetNoNorm)
mannWhitNoNormPadj <- p.adjust(mannWhitNoNorm, method = "BH")

which(mannWhitNoNormPadj < 0.05)


mithFormat <- cbind.data.frame(mithFormat, statMetsRhamn2$mannWhit, logFCRhamn2)

colnames(mithFormat)[(ncol(mithFormat)-1):ncol(mithFormat)] <- c("p.adj", "log2fc")

mithFormatNoNorm <- cbind.data.frame(mithFormatNoNorm, 
                                     mannWhitNoNormPadj[match(rownames(mithFormatNoNorm),
                                                              names(mannWhitNoNormPadj))], 
                                     logFCRhamn2NoNorm[match(rownames(mithFormatNoNorm),
                                                             names(logFCRhamn2NoNorm))])

colnames(mithFormatNoNorm)[(ncol(mithFormatNoNorm)-1):ncol(mithFormatNoNorm)] <- c("p.adj", "log2fc")


mithFormatRhamPos <- mithFormat[!colnames(mithFormat) %in% rownames(RhamnNegMet)]
mithFormatRhamNeg <- mithFormat[!colnames(mithFormat) %in% rownames(RhamnPosMet)]

mithFormatRhamPosNoNorm <- mithFormatNoNorm[!colnames(mithFormatNoNorm) %in% rownames(RhamnNegMetNoNorm)]
mithFormatRhamNegNoNorm <- mithFormatNoNorm[!colnames(mithFormatNoNorm) %in% rownames(RhamnPosMetNoNorm)]

write.csv(mithFormat, file = "mithFormat.csv")

write.csv(mithFormatNoNorm, file = "mithFormatNoNorm.csv")

write.csv(mithFormatRhamPos, file = "mithFormatRhamPos.csv")
write.csv(mithFormatRhamNeg, file = "mithFormatRhamNeg.csv")
write.csv(mithFormatRhamPosNoNorm, file = "mithFormatRhamPosNoNorm.csv")
write.csv(mithFormatRhamNegNoNorm, file = "mithFormatRhamNegNoNorm.csv")

read.csv("mithFormatRhamPos.csv")
