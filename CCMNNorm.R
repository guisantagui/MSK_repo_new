setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/metabolomicsPackageFormat.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/newData/metabolitePeaks_newData.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure2/solvedAmbigMets.RData")

if(!require(crmn)) install.packages("crmn")
if(!require(plotly)) install.packages("plotly")
if(!require(AUC)) install.packages("AUC")
if(!require(GGally)) install.packages("GGally")
if(!require(rmarkdown)) install.packages("rmarkdown")
if(!require(impute)) BiocManager::install("impute")
if(!require(pcaMethods)) BiocManager::install("pcaMethods")
if(!require(NormalizeMets)) install.packages("C:/Users/Guillem/Documents/PhD/comput/NormalizeMets_0.25.tar.gz",
                                             repos = NULL, 
                                             type = "source")
library(NormalizeMets)
if(!require(cluster)) install.packages("cluster")
library(cluster)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

metvar <- function(valuemat, groups) {
        ns <- nrow(valuemat)
        nmet <- ncol(valuemat)
        unig <- unique(groups)
        ng <- length(unig)
        meanvalues <- matrix(nrow=ng, ncol=nmet)
        row.names(meanvalues) <- unig
        colnames(meanvalues) <- colnames(valuemat)
        meandifvalues <- matrix(nrow=ns, ncol=nmet)
        meancvvalues <- matrix(nrow=ns,ncol=nmet)
        for (i in 1:ng) {
                meanvalues[i, ] <- colMeans(valuemat[groups == unig[i], ], na.rm = T)
                for (j in 1:nmet) {
                        meandifvalues[groups == unig[i], j] <-
                                valuemat[groups == unig[i], j] - meanvalues[i, j]
                        meancvvalues[groups == unig[i], j] <-
                                (valuemat[groups == unig[i], j] - meanvalues[i, j])/meanvalues[i, j]
                }
        }
        
        ssr <- vector()
        sst <- vector()
        ssbet <- vector()
        n1 <- vector()
        n2 <- vector()
        fstat <- vector()
        kruskals <- vector()
        kruskalp <- vector()
        kruskalpadj <- vector()
        metmean <- colMeans(valuemat, na.rm = T)
        metnames <- colnames(valuemat)
        rangefc <- vector()
        
        for (i in 1:nmet) {
                ssr[i] <- sum(meandifvalues[, i]^2, na.rm = T)
                sst[i] <- sum((valuemat[, i] - metmean[i])^2, na.rm = T)
                ssbet[i] <- sst[i] - ssr[i]
                n1[i] <- sum(!is.na(meanvalues[, i])) - 1
                n2[i] <- sum(!is.na(valuemat[, i])) - n1[i]
                fstat[i] <- (ssbet[i]/n1[i])/(ssr[i]/n2[i])
                datavec <- valuemat[, i]
                groupvec <- groups[!is.na(datavec)]
                datavec <- datavec[!is.na(datavec)]
                kresult <- kruskal.test(datavec,groupvec)
                kruskalp[i] <- kresult$p.value
                kruskals[i] <- kresult$statistic
                rangefc[i] <- 10^(max(meanvalues[, i]) - min(meanvalues[, i]))
        }
        kruskalpadj <- p.adjust(kruskalp, method = "BH")
        metvariation <- data.frame(names = metnames, ssr, ssbet, sst,
                                   n1, n2, f = fstat, k = kruskals, p = kruskalp, 
                                   padj = kruskalpadj, mean = metmean, range = rangefc)
        
        
        list(metvar=metvariation, groupmeans=meanvalues, residues=meandifvalues, cv=meancvvalues)
        
}


neibcheck <- function(valuemat) {
        gclass <- c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,
                    11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,
                    19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24,25,25,25,26,26,26,27,27,27)
        straindist <- as.matrix(dist(valuemat, method="euclidean"))
        score <- rep(0, 28)
        rscore <- rep(0, 28)
        silh <- silhouette(gclass, dmatrix = straindist)
        mgclass <- matrix(rep(gclass, 84), nrow=84, ncol=84)
        tmgclass <- t(mgclass)
        diffclass <- (mgclass != tmgclass)
        sameclassdist <- mean(straindist[!diffclass])
        difclassdist <- mean(straindist[diffclass])
        for (i in 1:28) {
                subdist <- straindist[(1+3*(i-1)):(3*i), ]
                points <- 0
                rpoints <- 0
                for (j in 1:3){
                        rankvec <- rank(subdist[j, ])
                        dvec <- subdist[j, ]
                        mvec <- dvec
                        mvec[dvec == 0] <- max(dvec)
                        dvec[dvec == 0] <- min(mvec)
                        dvec <- (dvec - min(mvec))/(max(dvec) - min(mvec))
                        points <- points + sum(rankvec[(1+3*(i-1)):(3*i)]) - 6
                        rpoints <- rpoints + sum(dvec[(1+3*(i-1)):(3*i)])
                }
                score[i] <- points
                rscore[i] <- rpoints
        }
        score <- score/3
        rscore <- rscore/6
        list(score = score, rscore = rscore, silh = silh[,3],
             sdist = sameclassdist, ddist = difclassdist)
}


colnames(metabolomicsPackageFormat)[2:ncol(metabolomicsPackageFormat)] <- dictionary$Consensus[match(colnames(metabolomicsPackageFormat)[2:ncol(metabolomicsPackageFormat)], 
                                                                                                     dictionary$`Old Data Names`)]

metabolomicsPackageFormat <- metabolomicsPackageFormat[, !is.na(colnames(metabolomicsPackageFormat))]

solvedAmbigMets <- solvedAmbigMets[-4]

# Change names of solved ambig mets. 
colnames(metabolomicsPackageFormat)[match(names(solvedAmbigMets), colnames(metabolomicsPackageFormat))] <- solvedAmbigMets

newDatMetNames <- metabolitePeaks_newData$`Compound Name`

metsNew <- t(metabolitePeaks_newData[, 2:(ncol(metabolitePeaks_newData)-1)])

colnames(metsNew) <- newDatMetNames

metsNewAdded <- as.data.frame(metsNew[, dictionary$State[match(colnames(metsNew), dictionary$`New Data Names`)] == "Added"])

metsNewAdded$`hexose diphosphate`[metsNewAdded$`hexose diphosphate` == min(metsNewAdded$`hexose diphosphate`)] <- NA

metabolomicsPackageFormat <- cbind.data.frame(metabolomicsPackageFormat,
                                              metsNewAdded)


# Do scatter plot of metabolite quality according to NA proportion
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

metQual <- data.frame("metabolite" = colnames(metabolomicsPackageFormat[2:ncol(metabolomicsPackageFormat)]),
                      "missingProp" = apply(metabolomicsPackageFormat[, 2:ncol(metabolomicsPackageFormat)], 
                                            2, 
                                            function(x) sum(is.na(x))/nrow(metabolomicsPackageFormat)),
                      "rangePeaks" = apply(metabolomicsPackageFormat[, 2:ncol(metabolomicsPackageFormat)],
                                           2,
                                           function(x) max(x[!is.na(x)]) - min(x[!is.na(x)])))

metqualDefNames <- dictionary$definitiveNames[match(metQual$metabolite, dictionary$Consensus)]
metqualDefNames[is.na(metqualDefNames)] <- as.character(metQual$metabolite[is.na(metqualDefNames)])

#Change names to the ones suggested by Kyu Rhee
KR_metNames <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/KRNewMetNames/metabolitesName-Xavier lab-KR.xlsx"))
KR_metNames$oldName[grep("Hexose phosphates", KR_metNames$oldName)] <- "Hexose phosphate"

metqualDefNames[!is.na(match(metqualDefNames, 
                             KR_metNames$oldName))] <- KR_metNames$newName[match(metqualDefNames, 
                                                                                 KR_metNames$oldName)[!is.na(match(metqualDefNames, 
                                                                                                                   KR_metNames$oldName))]]

metQual$metabolite <- metqualDefNames

ggplot(data = metQual, aes(x = missingProp, y = rangePeaks)) + geom_point() +
        geom_text(label=metQual$metabolite, nudge_x  = 0.04, nudge_y = 1) + 
        ylab("Peak range") + 
        xlab("Missing Value Proportion") + 
        theme_minimal()
ggsave(filename = "metQuality.pdf", width = 6, height = 6)

metQualMaxPeak <- data.frame("metabolite" = colnames(metabolomicsPackageFormat[2:ncol(metabolomicsPackageFormat)]),
                             "missingProp" = apply(metabolomicsPackageFormat[, 2:ncol(metabolomicsPackageFormat)], 
                                                   2, 
                                                   function(x) sum(is.na(x))/nrow(metabolomicsPackageFormat)),
                             "rangePeaks" = apply(metabolomicsPackageFormat[, 2:ncol(metabolomicsPackageFormat)],
                                                  2,
                                                  function(x) max(x[!is.na(x)])))

metQualMaxPeak$metabolite <- metqualDefNames

ggplot(data = metQualMaxPeak, aes(x = missingProp, y = rangePeaks)) + geom_point() +
        geom_text(label=metQual$metabolite, nudge_x  = 0.04, nudge_y = 1) + 
        ylab("Maximum Peak Area") + 
        xlab("Missing Value Proportion") + 
        theme_minimal()
ggsave(filename = "metQualityMaxPeak.pdf", width = 6, height = 6)

# Substitute NAs: the ones that doesn't appear in any replicates by zero, the other ones by replicates' mean.

strains <- gsub("\\_.*|(PA14).*", 
                rownames(metabolomicsPackageFormat), 
                rep = "\\1")

metMat <- metabolomicsPackageFormat

hasNAs <- colnames(metMat)[2:ncol(metMat)][which(apply(metMat[, 2:ncol(metMat)], 
                                                       2, 
                                                       function(x) sum(is.na(x))) > 0)]

# Look at compounds with NAs and annotation to see if makes sense that this compounds are not present 
# in the strains that are NAs. 

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")

hasNAsKEGGIDs <- dictionary$`KEGG IDs`[match(hasNAs, dictionary$Consensus)]

hasNAsKEGGIDs[!is.na(dictionary$`KEGG IDs`[match(hasNAs, 
                                                 dictionary$definitiveNames)])] <- dictionary$`KEGG IDs`[match(hasNAs, 
                                                                                                               dictionary$definitiveNames)][!is.na(dictionary$`KEGG IDs`[match(hasNAs, 
                                                                                                                                                                               dictionary$definitiveNames)])]
hasNAsKEGGIDs <- gsub("*", "", hasNAsKEGGIDs, fixed = T)

hasNAsKEGGIDs <- hasNAsKEGGIDs[!is.na(hasNAsKEGGIDs)] 

relEnzsNACpds <- sapply(hasNAsKEGGIDs, function(x) gsub("ec:", "", keggLink(target = "enzyme", source = x)))

# 2-Acetyl-aminoadipate_?

dictEnzymes[match(relEnzsNACpds[[1]], dictEnzymes$ECnums), ] # Any related enzyme with this compound

# 2,6-Diaminoheptanedioate 1_?

dictEnzymes[match(relEnzsNACpds[[2]], dictEnzymes$ECnums), ] # All the strains have the enzymes. Has many NAs

cbind(metabolomicsPackageFormat$`2,6-Diaminoheptanedioate 1_?`, metabolomicsPackageFormat$Group)

# 4-Amino-4-cyanobutanoic acid_?

dictEnzymes[match(relEnzsNACpds[[3]], dictEnzymes$ECnums), ] # Any strain has the enzymes (EC 3.5.5.1)

cbind(metabolomicsPackageFormat$`4-Amino-4-cyanobutanoic acid_?`, metabolomicsPackageFormat$Group)

# Acetymuramate

dictEnzymes[match(relEnzsNACpds[[4]], dictEnzymes$ECnums), ] # All the strains have the enzyme 

cbind(metabolomicsPackageFormat$Acetylmuramate, metabolomicsPackageFormat$Group)

# Homo-cisaconitate

dictEnzymes[match(relEnzsNACpds[[5]], dictEnzymes$ECnums), ] # Any strain has the enzymes. Has many NAs. Here NAs could be zeros.

cbind(metabolomicsPackageFormat$`Homo-cisaconitate 2`, metabolomicsPackageFormat$Group) 

# d-ala-d-ala

dictEnzymes[match(relEnzsNACpds[[6]], dictEnzymes$ECnums), ] # All the strains have the enzymes. Not many NAs

cbind(metabolomicsPackageFormat$`d-ala-d-ala 2_?`, metabolomicsPackageFormat$Group)

# Deoxyguanosine

dictEnzymes[match(relEnzsNACpds[[7]], dictEnzymes$ECnums), ] # EC3.1.3.5, 3.1.5.1 & 2.4.2.1 are in all strains. EC 2.4.2.4 is not present in 
                                                             # X78812, X9820, T63266, PA14, W16497, T52373, M74707, F5677, F63912, H27930,
                                                             # H5708, M1608, M37351. Not many NAs. The strains that have are others.

cbind(metabolomicsPackageFormat$`deoxyguanosine 1_?`, metabolomicsPackageFormat$Group)

# Dimethylmalate

dictEnzymes[match(relEnzsNACpds[[8]], dictEnzymes$ECnums), ] # Any strain has the enzymes. Many NAs in all replicates 

cbind(metabolomicsPackageFormat$`Dimethylmalate 1_?`, metabolomicsPackageFormat$Group)

# Formyl-methionine

dictEnzymes[match(relEnzsNACpds[[9]], dictEnzymes$ECnums), ] # Any strain has the enzymes. Not many NAs

cbind(metabolomicsPackageFormat$`Formyl-methionine`, metabolomicsPackageFormat$Group)

# Glutamine

dictEnzymes[match(relEnzsNACpds[[10]], dictEnzymes$ECnums), ] # All the strains have the enzymes

# Guanosine

dictEnzymes[match(relEnzsNACpds[[11]], dictEnzymes$ECnums), ] # All the strains have the enzymes

# Homovanillate

dictEnzymes[match(relEnzsNACpds[[12]], dictEnzymes$ECnums), ] # All the strains have the enzymes

# Hypoxantine

dictEnzymes[match(relEnzsNACpds[[13]], dictEnzymes$ECnums), ] # EC 2.4.2.44, 2.4.2.8, 3.5.4.2, 3.2.2.8, 1.17.1.4 & 2.4.2.1 present
                                                              # in all strains. EC 2.4.2.4 is not present in 
                                                              # X78812, X9820, T63266, PA14, W16497, T52373, M74707, F5677, F63912, H27930,
                                                              # H5708, M1608, M37351. Very few NAs

cbind(metabolomicsPackageFormat$hypoxanthine, metabolomicsPackageFormat$Group)

# Norepinephrine

dictEnzymes[match(relEnzsNACpds[[14]], dictEnzymes$ECnums), ] # Any strain has the enzymes. Very few NAs

cbind(metabolomicsPackageFormat$Norepinephrine, metabolomicsPackageFormat$Group)

# Phosphoenolpyruvate

dictEnzymes[match(relEnzsNACpds[[15]], dictEnzymes$ECnums), ] # All the strains have the enzymes, except EC 4.1.1.31, which is not present in 
                                                              # W60856. Many NAs

cbind(metabolomicsPackageFormat$Phosphoenolpyruvate, metabolomicsPackageFormat$Group)

# UDP-N-acetylmuramate

dictEnzymes[match(relEnzsNACpds[[16]], dictEnzymes$ECnums), ] # All the strains have the enzymes. Many NAs


imputed <- metMat
for(i in seq_along(hasNAs)){
        impMet <- metMat[[hasNAs[i]]]
        for(j in seq_along(strains)){
                strainSel <- gsub("\\_.*|(PA14).*",
                                  rownames(metMat),
                                  rep = "\\1") %in% strains[j]
                subVec <- impMet[strainSel]
                if(sum(is.na(subVec)) == length(subVec)){
                        impMet[strainSel] <- rep(NA, length(subVec))
                }else{
                        meanImp <- mean(subVec[!is.na(subVec)])
                        impMet[strainSel] <- rep(meanImp, length(subVec))
                }
        }
        imputed[[hasNAs[i]]] <- impMet 
}

metQualPostImp <- data.frame("metabolite" = colnames(imputed[2:ncol(imputed)]),
                             "missingProp" = apply(imputed[, 2:ncol(imputed)], 
                                                   2, 
                                                   function(x) sum(is.na(x))/nrow(imputed)),
                             "rangePeaks" = apply(imputed[, 2:ncol(imputed)],
                                                  2,
                                                  function(x) max(x[!is.na(x)]) - min(x[!is.na(x)])))

metQualPostImp$metabolite <- metqualDefNames

ggplot(data = metQualPostImp, aes(x = missingProp, y = rangePeaks)) + geom_point() +
        geom_text(label=metQualPostImp$metabolite, nudge_x  = 0.04, nudge_y = 1) + 
        ylab("Peak range") + 
        xlab("Missing Value Proportion") + 
        theme_minimal()
ggsave(filename = "metQualityPostImp.pdf", width = 6, height = 6)

impNAProp <- apply(imputed, 2, function(x) sum(is.na(x))/length(x))

# remove metabolites over the threshold
toRemove <- names(impNAProp[which(impNAProp > 0.00)])

imputed[is.na(imputed)] <- 1
# Remove these 2 compounds because they have more than 80% of NAs

imputed <- imputed[, !colnames(imputed) %in% toRemove]

metvarprenorm <- metvar(imputed[, 2:ncol(imputed)],
                        imputed[, 1])

#detach("package:NormalizeMets", unload = T)
#library(metabolomics)
#metvarprenorm <- metvar(LogTransform(MissingValues(imputed[, 2:ncol(imputed)], column.cutoff = 0.8,
#                                                   group.cutoff = 0.2)$output,
#                                     base = 10)$output,
#                        imputed[, 1])

#detach("package:metabolomics", unload = T)
#library(NormalizeMets)

featuredataMet <- imputed[, 2:ncol(imputed)]
sampledataMet <- data.frame(Group = imputed[, 1], 
                            Species = c(rep("P. aeruginosa", length(imputed[, 1]))))

IS <- c(rep(0, length(metvarprenorm$metvar$padj)))
for(i in 1:length(metvarprenorm$metvar$padj)){
        if(metvarprenorm$metvar$padj[i]>0.05){
                IS[i] <- 1
        }
}

metabolitedataMet <- data.frame(names = colnames(featuredataMet), IS = IS)

alldataMet <- list(featuredata = featuredataMet, sampledata = sampledataMet, 
                   metabolitedata = metabolitedataMet)

unNormMetImpData <- featuredataMet

save(unNormMetImpData, file = "unNormMetImpData.RData")

#LogTransform Data

logdataMet<- LogTransform(alldataMet$featuredata, zerotona = T)

# Do boxplots of samples and metabolites

metnames=colnames(logdataMet$featuredata)
samplename=rownames(logdataMet$featuredata)
logdatatable=data.frame(sample=rep('.',length(samplename)*length(metnames)),
                        strain=rep('.',length(samplename)*length(metnames)),
                        metabolite=rep('.',length(samplename)*length(metnames)),
                        peak=rep(0,length(samplename)*length(metnames)),
                        significant=rep(0,length(samplename)*length(metnames)),stringsAsFactors = F)
line=0
for (i in 1:length(samplename)){
        for (j in seq_along(metnames)){
                line=line+1
                logdatatable[line,1]=samplename[i]
                logdatatable[line,2]=as.character(logdataMet$featuredata[i,1])
                logdatatable[line,3]=metnames[j]
                logdatatable[line,4]=logdataMet$featuredata[i,j]
                if (metvarprenorm$metvar$padj[metvarprenorm$metvar$names==metnames[j]]<0.05){
                        logdatatable[line,5]=1
                }
        }
}

save(logdatatable, file = "logdatatable.RData")

ggplot(data = logdatatable, mapping = aes(x = sample, y = peak)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Abundance",
                x= "Sample"
        )

ggsave("sample_boxplot.pdf",width = 9, height = 5)  


ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=median), y = peak, color=significant)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Abundance",
                x= "Metabolite"
        )

ggsave("metabolite_boxplot.pdf",width = 7, height = 11)  


missingMet <-  MissingValues(logdataMet$featuredata, sampledataMet,metabolitedataMet,
                             feature.cutof=0.8, sample.cutoff=0.2, method="knn")


neibprenorm <- neibcheck(missingMet$featuredata)

# Normalize with different methods
facmat=matrix(0, nrow=84, ncol=28)
for (i in 1:28){
        facmat[(1+3*(i-1)):(3*i),i]=1
}


ccmnNormMets <- NormQcmets(missingMet$featuredata, method = "ccmn", 
                           qcmets = which(missingMet$metabolitedata$IS == 1), factors = sampledataMet$Group,
                           saveoutput = T)
metvarCcmnNorm <- metvar(ccmnNormMets$featuredata, metabolomicsPackageFormat$Group)
neibCcmnNorm <- neibcheck(ccmnNormMets$featuredata)

nomisNormMets <- NormQcmets(missingMet$featuredata, method = "nomis", 
                            qcmets = which(missingMet$metabolitedata$IS == 1))
metvarNomisNorm <- metvar(nomisNormMets$featuredata, metabolomicsPackageFormat$Group)
neibNomisNorm <- neibcheck(nomisNormMets$featuredata)

ruvRandNormMets <- NormQcmets(missingMet$featuredata, method = "ruvrand", k = 2, 
                              qcmets = which(missingMet$metabolitedata$IS == 1),
                              plotk = T)
metvarRuvRandNorm <- metvar(ruvRandNormMets$featuredata, metabolomicsPackageFormat$Group)
neibRuvRandNorm <- neibcheck(ruvRandNormMets$featuredata)

isNormMets <- NormQcmets(missingMet$featuredata, method = "is", 
                         isvec = missingMet$featuredata[, which(missingMet$metabolitedata$IS == 1)[1]])
metvarIsNorm <- metvar(isNormMets$featuredata, metabolomicsPackageFormat$Group)
neibIsNorm <- neibcheck(isNormMets$featuredata)

medianNormMets <- NormScaling(missingMet$featuredata, method = "median")
metvarMedianNorm <- metvar(medianNormMets$featuredata, metabolomicsPackageFormat$Group)
neibMedianNorm <- neibcheck(medianNormMets$featuredata)

meanNormMets <- NormScaling(missingMet$featuredata, method = "mean")
metvarMeanNorm <- metvar(meanNormMets$featuredata, metabolomicsPackageFormat$Group)
neibMeanNorm <- neibcheck(meanNormMets$featuredata)

sumNormMets <- NormScaling(missingMet$featuredata, method = "sum")
metvarSumNorm <- metvar(sumNormMets$featuredata, metabolomicsPackageFormat$Group)
neibSumNorm <- neibcheck(sumNormMets$featuredata)


NormalizeMethod <- c("pre", "nomis", "ccmn", "ruvrandom", "ruvrandclust", "is", "median", "mean", "sum")
NormalizeReport <- data.frame(method=NormalizeMethod, silh=rep(0,9), score=rep(0,9),
                              rscore=rep(0,9), difmet=rep(0,9), sd=rep(0,9), 
                              dd=rep(0,9), stringsAsFactors = F)

NormalizeReport[1,2] = mean(neibprenorm$silh)
NormalizeReport[1,3] = mean(neibprenorm$score)
NormalizeReport[1,4] = mean(neibprenorm$rscore)
NormalizeReport[1,5] = sum(metvarprenorm$metvar$padj<0.01)
NormalizeReport[1,6] = neibprenorm$sdist
NormalizeReport[1,7] = neibprenorm$ddist


NormalizeReport[2,2] = mean(neibNomisNorm$silh)
NormalizeReport[2,3] = mean(neibNomisNorm$score)
NormalizeReport[2,4] = mean(neibNomisNorm$rscore)
NormalizeReport[2,5] = sum(metvarNomisNorm$metvar$padj<0.01)
NormalizeReport[2,6] = neibNomisNorm$sdist
NormalizeReport[2,7] = neibNomisNorm$ddist


NormalizeReport[3,2] = mean(neibCcmnNorm$silh)
NormalizeReport[3,3] = mean(neibCcmnNorm$score)
NormalizeReport[3,4] = mean(neibCcmnNorm$rscore)
NormalizeReport[3,5] = sum(metvarCcmnNorm$metvar$padj<0.01)
NormalizeReport[3,6] = neibCcmnNorm$sdist
NormalizeReport[3,7] = neibCcmnNorm$ddist

NormalizeReport[4,2] = mean(neibRuvRandNorm$silh)
NormalizeReport[4,3] = mean(neibRuvRandNorm$score)
NormalizeReport[4,4] = mean(neibRuvRandNorm$rscore)
NormalizeReport[4,5] = sum(metvarRuvRandNorm$metvar$padj<0.01)
NormalizeReport[4,6] = neibRuvRandNorm$sdist
NormalizeReport[4,7] = neibRuvRandNorm$ddist

NormalizeReport[6,2] = mean(neibIsNorm$silh)
NormalizeReport[6,3] = mean(neibIsNorm$score)
NormalizeReport[6,4] = mean(neibIsNorm$rscore)
NormalizeReport[6,5] = sum(metvarIsNorm$metvar$padj<0.01, na.rm = T)
NormalizeReport[6,6] = neibIsNorm$sdist
NormalizeReport[6,7] = neibIsNorm$ddist

NormalizeReport[7,2] = mean(neibMedianNorm$silh)
NormalizeReport[7,3] = mean(neibMedianNorm$score)
NormalizeReport[7,4] = mean(neibMedianNorm$rscore)
NormalizeReport[7,5] = sum(metvarMedianNorm$metvar$padj<0.01)
NormalizeReport[7,6] = neibMedianNorm$sdist
NormalizeReport[7,7] = neibMedianNorm$ddist

NormalizeReport[8,2] = mean(neibMeanNorm$silh)
NormalizeReport[8,3] = mean(neibMeanNorm$score)
NormalizeReport[8,4] = mean(neibMeanNorm$rscore)
NormalizeReport[8,5] = sum(metvarMeanNorm$metvar$padj<0.01)
NormalizeReport[8,6] = neibMeanNorm$sdist
NormalizeReport[8,7] = neibMeanNorm$ddist

NormalizeReport[9,2] = mean(neibSumNorm$silh)
NormalizeReport[9,3] = mean(neibSumNorm$score)
NormalizeReport[9,4] = mean(neibSumNorm$rscore)
NormalizeReport[9,5] = sum(metvarSumNorm$metvar$padj<0.01)
NormalizeReport[9,6] = neibSumNorm$sdist
NormalizeReport[9,7] = neibSumNorm$ddist

forReport <- NormalizeReport[c(1, 2, 3, 4, 7, 8),]

tiff("normmethod_dist_normalizeMets.tiff", res = 300, height = 1500, width = 1500)
ggplot(data = forReport, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=NormalizeMethod[c(1, 2, 3, 4, 7, 8)]),nudge_y=0.025) +
        labs(
                y= "Mean distance between strains",
                x= "Mean distance between replicates"
        )
dev.off()

# Assess silhouette


silhvarNormalizeMets <- data.frame(method=rep("pre",length(neibprenorm$silh)),
                                   silh=neibprenorm$silh)
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("nomis",length(neibNomisNorm$silh)),
                                         silh=neibNomisNorm$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets, 
                              data.frame(method=rep("ccmn",length(neibCcmnNorm$silh)),
                                         silh=neibCcmnNorm$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("ruvrandom",length(neibRuvRandNorm$silh)),
                                         silh=neibRuvRandNorm$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("median",length(neibMedianNorm$silh)),
                                         silh=neibMedianNorm$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("mean",length(neibMeanNorm$silh)),
                                         silh=neibMeanNorm$silh))

tiff("normmethod_silh_normalizeMets.tiff", res = 300, height = 1500, width = 1500)
ggplot(data = silhvarNormalizeMets, mapping = aes(x = reorder(method, silh, FUN=median), y = silh)) +
        geom_boxplot() +
        #geom_line(data=normreport,mapping=aes(x=method,y=silh,group=1)) +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Silhouete",
                x= "Method"
        )
dev.off()
ggsave("normmethod_silh_normalizeMets.pdf") 


#Within-group RLA plot

RlaPlots(featuredata = missingMet$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Pre-normalization data",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = isNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, IS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = nomisNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, NOMIS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = ccmnNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, CCMN",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = ruvRandNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, RUV Random",
         saveplot = T, savetype = "tiff")
RlaPlots(featuredata = medianNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Median",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = meanNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Mean",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = sumNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Sum",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
#Across groups RLA plots
RlaPlots(featuredata = missingMet$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Pre-normalization Data",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = isNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, IS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = nomisNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, NOMIS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = ccmnNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, CCMN",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = ruvRandNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, RUV Random",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = medianNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Median",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = meanNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Mean",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = sumNormMets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Sum",
         saveplot = T, savetype = "tiff", showlegend = T,
         ylim = c(-2, 2))    

# p value histograms

qplot(metvarprenorm$metvar$padj, geom = "histogram",
      main = "Pre-normalization", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5), 
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_pre_oldData.jpeg", device = "jpeg")

qplot(metvarNomisNorm$metvar$padj, geom = "histogram",
      main = "NOMIS", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_nomis_oldData.jpeg", device = "jpeg")

qplot(metvarCcmnNorm$metvar$padj, geom = "histogram",
      main = "CCMN", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ccmn.jpeg", device = "jpeg")

qplot(metvarRuvRandNorm$metvar$padj, geom = "histogram",
      main = "RUV Random", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ruvrand.jpeg", device = "jpeg")

qplot(metvarMedianNorm$metvar$padj, geom = "histogram",
      main = "Median", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_median.jpeg", device = "jpeg")

qplot(metvarMeanNorm$metvar$padj, geom = "histogram",
      main = "Mean", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_mean.jpeg", device = "jpeg")


ccmnNormMets <- ccmnNormMets$featuredata
colnames(ccmnNormMets) <- colnames(imputed)[match(colnames(ccmnNormMets), make.names(colnames(imputed)))]

meanNormMets <- meanNormMets$featuredata
colnames(meanNormMets) <- colnames(imputed)[match(colnames(meanNormMets), make.names(colnames(imputed)))]

nomisNormMets <- nomisNormMets$featuredata
colnames(nomisNormMets) <- colnames(imputed)[match(colnames(nomisNormMets), make.names(colnames(imputed)))]

write.csv(ccmnNormMets, file = "ccmnNormMets.csv")

save(ccmnNormMets, file = "ccmnNormMets.RData")
save(meanNormMets, file = "meanNormMets.RData")
save(nomisNormMets, file = "nomisNormMets.RData")

