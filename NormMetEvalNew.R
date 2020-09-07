setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/newData")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script prepares the Mass Spectrometry metabolomic data to a format appropriate for "metabolomics" and           #################################
# "NormalizeMets" R packages, for then normalize it through different methods and later evaluate the performance of    #################################
# these methods through different approaches for, finally, exporting the data normalized by CCMN method (the best      #################################
# one) in a CSV file. It differs to NormMethEval.R in that the input it takes is a curated version of the metabolomic  #################################
# data used in past version of the script.                                                                             #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

if(!require(readxl)) install.packages("readxl")
library(readxl)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)
if(!require(UpSetR)) install.packages("UpSetR")
library(UpSetR)

########################################################################################################################################################
#
# Read the data and prepare the format.
#
########################################################################################################################################################

#Load data


mets <-read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/newData/metaboliteTable.xlsx")

strain_names <- unique(mets$mutant)
sampleNames <- gsub(unique(gsub(gsub(mets$sample, pattern = "Negative", replacement = ""), pattern = "Positive", replacement = "")), pattern = "day1", replacement = "_")

positives <- mets[which(mets$mode == "Positive"), ]
posMat <- NULL
posBatches <- NULL
for(i in 1:length(strain_names)){
        for(j in 1:3){
                if(is.null(posMat) == T){
                        posMat <- positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 1]
                        posBatches <- positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 9]
                }else{
                        posMat <- cbind(posMat, positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 1])
                        posBatches <- cbind(posBatches, positives[positives$sample == paste(strain_names[i], "day1PositiveR", j, sep = ""), 9])
                }
        }
}
posMetNames <- unique(positives$metabolite)
colnames(posMat) <- sampleNames
posMat <- cbind.data.frame(posMetNames, as.data.frame(posMat))
colnames(posMat)[1] <- "Compound Name"
posFormulas <- c()
for(i in seq_along(posMetNames)){
        posFormulas<- c(posFormulas, unique(positives$formula[which(positives$metabolite == posMetNames[i])]))
}

posMat <- cbind(posMat, posFormulas)
colnames(posMat)[86] <- "Formula"
colnames(posBatches) <- sampleNames

negatives <- mets[which(mets$mode == "Negative"), ]
negMat <- NULL
negBatches <- NULL
for(i in 1:length(strain_names)){
        for(j in 1:3){
                if(is.null(negMat) == T){
                        negMat <- negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 1]
                        negBatches <- negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 9]
                }else{
                        negMat <- cbind(negMat, negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 1])
                        negBatches <- cbind(negBatches, negatives[negatives$sample == paste(strain_names[i], "day1NegativeR", j, sep = ""), 9])
                }
        }
}
negMetNames <- unique(negatives$metabolite)
colnames(negMat) <- sampleNames
negMat <- cbind.data.frame(negMetNames, as.data.frame(negMat))
colnames(negMat)[1] <- "Compound Name"
negFormulas <- c()
for(i in seq_along(negMetNames)){
        negFormulas<- c(negFormulas, unique(negatives$formula[which(negatives$metabolite == negMetNames[i])]))
}

negMat <- cbind(negMat, negFormulas)
colnames(negMat)[86] <- "Formula"
colnames(negBatches) <- sampleNames

metabolitePeaks <- merge(posMat, negMat, all = T)
metabolitePeaks_newData <- metabolitePeaks
compoundNames <- metabolitePeaks$`Compound Name`
compoundNames_newData <- compoundNames
save(compoundNames_newData, file = "compoundNames_newData.RData")
save(metabolitePeaks_newData, file = "metabolitePeaks_newData.RData")

batches <- rbind(posBatches, negBatches)
batches <- apply(batches, 2, unique)

#Save data as .txt
write.table(metabolitePeaks_newData, 
            sep = "\t", 
            file = "/Users/santamag/Desktop/GUILLEM/wrkng_dirs_newData/normMetAnal_ND/curatePeaks_ND.txt", 
            row.names = F)

if(!require(gplots)) install.packages("gplots")
library(gplots)
if(!require(remotes)) install.packages("remotes")
library(remotes)
if(!require(xml2)) install.packages("xlm2")
if(!require(BiocManager)) install.packages("BiocManager")
if(!require(pcaMethods)) BiocManager::install("pcaMethods")
library(pcaMethods)
if(!require(limma)) BiocManager::install("limma")
if(!require(metabolomics)) install_version("metabolomics", "0.1.4")
library(metabolomics)
if(!require(grDevices)) install.packages("grDevices")
library(grDevices)
if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

#Convert to appropriate format for metabolomics package
sampleNames <- cbind(t(t(sampleNames)), c(1:84), gsub("\\_.*", "", sampleNames))
rownames(sampleNames) <- paste(sampleNames[, 3], sampleNames[, 2], sep = "_")
sampleNames <- sampleNames[, -1]
colnames(sampleNames) <- c("Sample", "Group")
sampleNames <- data.frame(sampleNames)
rownames(metabolitePeaks) <- metabolitePeaks[, 1]

metabolomicsPackageFormat <- data.frame(t(metabolitePeaks[, -c(1, ncol(metabolitePeaks))]))
rownames(metabolomicsPackageFormat) <- rownames(sampleNames)
metabolomicsPackageFormat <- merge(sampleNames, metabolomicsPackageFormat, by = "row.names",
                                   all.x = T, sort = F)
metabolomicsPackageFormat <- metabolomicsPackageFormat[, -c(1, 2)]
row.names(metabolomicsPackageFormat) <- row.names(sampleNames)
colnames(metabolomicsPackageFormat)[2:ncol(metabolomicsPackageFormat)] <- as.character(metabolitePeaks$`Compound Name`)
metabolomicsPackageFormat
write.csv(metabolomicsPackageFormat, file = "metabolomicsPackageFormat_newData.csv")

#Plot mean, median and sd of each metabolite
par(mar = c(1, 1, 1, 1))
meanMetaboliteLeves <- MmsPlot(metabolomicsPackageFormat, variables = "metabolites", main = "Mean")$output
meanMetaboliteLeves$metabolite <- row.names(meanMetaboliteLeves)

#Plot median front of mean
p <- ggplot(meanMetaboliteLeves, aes(Mean, Median))
p + geom_point()

#Plot mean peak area of each metabolite
p <- ggplot(meanMetaboliteLeves, aes(x = reorder(metabolite, -Mean), y = Mean))
p + geom_col() + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("Metabolites") +
        ylab("Mean Peak Area") +
        ggtitle("Mean Metabolite Levels") +
        ggsave("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/newData/meanPeakAreaRaw_ND.jpg")

#Make colormap for each strain
unique.groups <- unique(metabolomicsPackageFormat$Group)
cols <- topo.colors(length(unique.groups))
cols[1] <- '#000000'
cols[2] <- '#555555'

Cols <- c(rep(NA, length(rownames(metabolomicsPackageFormat))))
for(ii in 1:length(metabolomicsPackageFormat[, 1])) {
        selected <- which(unique.groups == metabolomicsPackageFormat[ii, 1])
        Cols[ii] <- cols[selected]
}

#Plot HeatMap of raw data
HeatMap(metabolomicsPackageFormat,
        scale = 'row', 
        ColSideColors = Cols,
        colramp = bluered(75))

#Replace missing values/zeros
metabolomicsPackageFormat[metabolomicsPackageFormat == 0] <- NA
missing <- MissingValues(metabolomicsPackageFormat, group.cutoff = 0.2, column.cutoff = 0.8)$output
length(which(is.na(missing))) / length(which(is.na(metabolomicsPackageFormat)))

#Log transform data & normalize (RUV2 is not really being applied-->Bad code)
log.data <- LogTransform(missing)$output
norm_ruv2 <- Normalise(missing)$output

#Plot Log transformed data 
HeatMap(log.data, 
        scale = 'column',
        ColSideColors = Cols, 
        colramp = bluered(75))

#Plot normalised data
HeatMap(norm_ruv2, 
        scale = 'column', 
        ColSideColors = Cols,
        colramp = bluered(75))


########################################################################################################################################################
#
# Normalize with different normalization methods and compare them with "NormalizeMets" Package
#
########################################################################################################################################################
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
library("KEGGREST")
if(!require(cluster)) install.packages("cluster")
library(cluster)
if(!require(ruv)) install.packages("ruv")
library(ruv)

# This function computes different statistics on the input data
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

#Make function that checks the validity of normalization method
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

########################################################################################################################################################
#
# Normalize with different methods and analyze the performance with "NormalizeMets" Package
#
########################################################################################################################################################

# Determine nonchanging metabolites with a Kruskall-Wallis test through metvar function (after replacing missing values and log transforming data)

metabolomicsPackageFormat[metabolomicsPackageFormat == 0] <- NA
missing <- MissingValues(metabolomicsPackageFormat, group.cutoff = 0.2, column.cutoff = 0.8)$output
length(which(is.na(missing))) / length(which(is.na(metabolomicsPackageFormat)))

log.data <- LogTransform(missing)$output
#Check variance and distance between groups before normalization
metvarprenorm <- metvar(log.data[, 2:ncol(log.data)], log.data[, 1])
neibprenorm <- neibcheck(log.data[, 2:ncol(log.data)])

metnames <- colnames(log.data[,2:ncol(log.data)])
samplename <- row.names(log.data)
logdatatable=data.frame(sample=rep('.',84*(ncol(log.data)-1)),
                        strain=rep('.',84*(ncol(log.data)-1)),
                        metabolite=rep('.',84*(ncol(log.data)-1)),
                        peak=rep(0,84*(ncol(log.data)-1)),
                        significant=rep(0,84*(ncol(log.data)-1)),
                        stringsAsFactors = F)
line=0
for (i in 1:84){
        for (j in 2:ncol(log.data)){
                line=line+1
                logdatatable[line,1]=samplename[i]
                logdatatable[line,2]=as.character(log.data[i,1])
                logdatatable[line,3]=metnames[j-1]
                logdatatable[line,4]=log.data[i,j]
                if (metvarprenorm$metvar$padj[metvarprenorm$metvar$names==metnames[j-1]]<0.05){
                        logdatatable[line,5]=1
                }
        }
}

ggplot(data = logdatatable, mapping = aes(x = sample, y = peak)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Abundance",
                x= "Sample"
        )


ggsave("sample_boxplot_newData.pdf",width = 9, height = 5)  


ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=median), y = peak, color=significant)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Abundance",
                x= "Metabolite"
        )

ggsave("metabolite_boxplot_newData.pdf",width = 7, height = 11)  

ggplot(data = logdatatable[logdatatable$significant==0,], mapping = aes(x = strain, y = peak)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_grid(metabolite ~ ., scales="free") +
        labs(
                y= "Abundance",
                x= "Strain"
        )

ggsave("const_metaboliteds_newData.pdf",width = 7, height = 14)

metvarprenorm <- metvar(log.data[, 2:ncol(log.data)], log.data[, 1])
neibprenorm <- neibcheck(log.data[, 2:ncol(log.data)])

IS <- c(rep(0, ncol(metabolomicsPackageFormat) - 1))
for(i in 1:(ncol(metabolomicsPackageFormat) - 1)){
        if(metvarprenorm$metvar$padj[i]>0.05){
                IS[i] <- 1
        }
}

# The following is the vector of IS (nonchanging mets) that were used by mistake in the past analysis. 

IS2 <- c(rep(0, ncol(metabolomicsPackageFormat) - 1))
for(i in 1:(ncol(metabolomicsPackageFormat) - 1)){
        if(colnames(metabolomicsPackageFormat)[i + 1] == "Acetyl.2..6.diaminopimelate" |
           colnames(metabolomicsPackageFormat)[i + 1] == "glutamate" |
           colnames(metabolomicsPackageFormat)[i + 1] == "Glycerate" |
           colnames(metabolomicsPackageFormat)[i + 1] == "Leucine" |
           colnames(metabolomicsPackageFormat)[i + 1] == "methylmaleate.1" |
           colnames(metabolomicsPackageFormat)[i + 1] == "S.Adenosylmethionine" |
           colnames(metabolomicsPackageFormat)[i + 1] == "Tyrosine"){
                IS2[i] <- 1
        }
}





# Adapt to format
detach("package:metabolomics", unload = T)
if(!require(NormalizeMets)) install.packages("NormalizeMets")
library(NormalizeMets)

featuredataMet <- metabolomicsPackageFormat[, 2:dim(metabolomicsPackageFormat)[2]]

sampledataMet <- data.frame(Group = metabolomicsPackageFormat[, 1], 
                            Batch = batches,
                            Species = c(rep("P. aeruginosa", length(metabolomicsPackageFormat[, 1]))))

metabolitedataMet <- data.frame(names = colnames(featuredataMet), IS = IS)

alldataMet <- list(featuredata = featuredataMet, sampledata = sampledataMet, 
                   metabolitedata = metabolitedataMet)

save(alldataMet, file = "allDataMet.RData")

#LogTransform Data

logdataMet<- LogTransform(alldataMet$featuredata, zerotona = T)
dataview(logdataMet$featuredata)

#Replace Missing

missingMet <-  MissingValues(logdataMet$featuredata, sampledataMet,metabolitedataMet,
                             feature.cutof=0.8, sample.cutoff=0.2, method="knn")
dataview(missingMet$featuredata)

# Normalization by different methods

# Build a factor matrix of 0s and 1s that tells to what strain (group) the sample corresponds (for ruv2 method)
facmat=matrix(0, nrow=84, ncol=28)
for (i in 1:28){
        facmat[(1+3*(i-1)):(3*i),i]=1
}


metvarPreNormNormalizeMets <- metvar(missingMet$featuredata, log.data[, 1])
neibPreNormNormalizeMets <- neibcheck(missingMet$featuredata)

norm_nomis_normalizemets <- NormQcmets(missingMet$featuredata, method = "nomis", 
                                       qcmets = which(missingMet$metabolitedata$IS == 1))
metvarNomisNormalizeMets <- metvar(norm_nomis_normalizemets$featuredata, log.data[, 1])
neibNomisNormalizemets <- neibcheck(norm_nomis_normalizemets$featuredata)

norm_ccmn_normalizemets <- NormQcmets(missingMet$featuredata, method = "ccmn", 
                                      qcmets = which(missingMet$metabolitedata$IS == 1), factors = sampledataMet$Group,
                                      saveoutput = T)
metvarCcmnNormalizeMets <- metvar(norm_ccmn_normalizemets$featuredata, log.data[, 1])
neibCcmnNormalizeMets <- neibcheck(norm_ccmn_normalizemets$featuredata)

norm_ruv2_normalizemets <- NormQcmets(missingMet$featuredata, factormat = facmat, method = "ruv2", k = 2, 
                                      qcmets = which(missingMet$metabolitedata$IS == 1))
metvarRuv2NormalizeMets <- metvar(norm_ruv2_normalizemets$coefficients, log.data[, 1])
neibRuv2NormalizeMets <- neibcheck(norm_ruv2_normalizemets$coefficients)

norm_ruvrandom_normalizemets <- NormQcmets(missingMet$featuredata, method = "ruvrand", k = 2, 
                                           qcmets = which(missingMet$metabolitedata$IS == 1), plotk = T)
metvarRuvRandNormalizeMets <- metvar(norm_ruvrandom_normalizemets$featuredata, log.data[, 1])
neibRuvRandNormalizeMets <- neibcheck(norm_ruvrandom_normalizemets$featuredata)

norm_ruvrandclust_normalizemets <- NormQcmets(missingMet$featuredata, method = "ruvrandclust", k = 2, 
                                              qcmets = which(missingMet$metabolitedata$IS == 1))
metvarRuvRandClustNormalizeMets <- metvar(norm_ruvrandclust_normalizemets$featuredata, log.data[, 1])
neibRuvRandClustNormalizeMets <- neibcheck(norm_ruvrandclust_normalizemets$featuredata)

norm_is_normalizemets <- NormQcmets(missingMet$featuredata, method = "is", 
                                    isvec = missingMet$featuredata[, which(missingMet$metabolitedata$IS == 1)[1]])
metvarIsNormalizeMets <- metvar(norm_is_normalizemets$featuredata, log.data[, 1])
neibIsNormalizeMets <- neibcheck(norm_is_normalizemets$featuredata)

norm_median_normalizemets <- NormScaling(missingMet$featuredata, method = "median")
metvarMedianNormalizeMets <- metvar(norm_median_normalizemets$featuredata, log.data[, 1])
neibMedianNormalizeMets <- neibcheck(norm_median_normalizemets$featuredata)

norm_mean_normalizemets <- NormScaling(missingMet$featuredata, method = "mean")
metvarMeanNormalizeMets <- metvar(norm_mean_normalizemets$featuredata, log.data[, 1])
neibMeanNormalizeMets <- neibcheck(norm_mean_normalizemets$featuredata)

norm_sum_normalizemets <- NormScaling(missingMet$featuredata, method = "sum")
metvarSumNormalizeMets <- metvar(norm_sum_normalizemets$featuredata, log.data[, 1])
neibSumNormalizeMets <- neibcheck(norm_sum_normalizemets$featuredata)


# ASSESS NORMALIZATION METHODS

#Fit normalized results to a linear model

unadjFit <- LinearModelFit(featuredata = missingMet$featuredata, 
                           factormat = facmat, ruv2 = F)
normNomisFit <- LinearModelFit(featuredata = norm_nomis_normalizemets$featuredata,
                               factormat = facmat, ruv2 = F)
normCcmnFit <- LinearModelFit(featuredata = norm_ccmn_normalizemets$featuredata,
                              factormat = facmat, ruv2 = F)
normRuv2Fit <- LinearModelFit(featuredata = missingMet$featuredata, 
                              factormat = facmat, ruv2 = T, k = 2, 
                              qcmets = which(metabolitedataMet$IS == 1))
normRuvRandFit <- LinearModelFit(featuredata = norm_ruvrandom_normalizemets$featuredata, 
                                 factormat = facmat, ruv2 = F)
normRuvRandClustFit <- LinearModelFit(featuredata = norm_ruvrandclust_normalizemets$featuredata,
                                      factormat = facmat, ruv2 = F)
normIsFit <- LinearModelFit(featuredata = norm_is_normalizemets$featuredata, 
                            factormat = facmat, ruv2 = F)
normMedianFit <- LinearModelFit(featuredata = norm_median_normalizemets$featuredata,
                                factormat = facmat, ruv2 = F)
normMeanFit <- LinearModelFit(featuredata = norm_mean_normalizemets$featuredata,
                              factormat = facmat, ruv2 = F)
normSumFit <- LinearModelFit(featuredata = norm_sum_normalizemets$featuredata,
                             factormat = facmat, ruv2 = F)

#Examine residuals with RLA plots

listResidData <- list(unadjusted = unadjFit$residuals, nomis = normNomisFit$residuals,
                      ccmn = normCcmnFit$residuals, ruv2 = normRuv2Fit$residuals, 
                      ruvRand = normRuvRandFit$residuals, ruvRandClust = normRuvRandClustFit$residuals,
                      is = normIsFit$residuals, median = normMedianFit$residuals, 
                      mean = normMeanFit$residuals, sum = normSumFit$residuals)

CompareRlaPlots(listResidData, groupdata = sampledataMet$Group, yrange = c(-3, 3), 
                normmeth = c("unadjusted:", "nomis:", "ccmn:", "ruv2:", "ruvand:", "ruvrandclust", 
                             "is:", "median:", "mean:", "sum:"), cols = rainbow(28), 
                savenoninteractive = TRUE, interactivesavename = "compmethods.png")

factor(as.data.frame(log.data[, 1])[, 1], levels = unique(as.data.frame(log.data[, 1])[, 1]))


unadjFit <- LinearModelFit(log.data[, 2:ncol(log.data)], factormat = facmat, ruv2 = F)
normisFit <- LinearModelFit(norm_is_normalizemets$featuredata, factormat = facmat, ruv2 = F)
normccmnFit <- LinearModelFit(norm_ccmn_normalizemets$featuredata, factormat = facmat, ruv2 = F)
unadjFit$coefficients[, "x1"]
lcoef_x1 <- list(unadjusted=unadjFit$coefficients[,"x1"],
                 is_x1=normisFit$coefficients[,"x1"],
                 ccmn_x1=normccmnFit$coefficients[,"x1"])
lpvals_x1 <- list(unadjusted=unadjFit$p.value[,"x1"],
                  is=normisFit$p.value[,"x1"],
                  ccmn=normccmnFit$p.value[,"x1"])

negcontrols <- colnames(log.data)[(which(metvarprenorm$metvar$padj>0.05)) + 1]

CompareVolcanoPlots(lcoef=lcoef_x1, 
                    lpvals_x1, 
                    normmeth = c(":unadjusted", ":is", ":ccmn"),
                    xlab="Coef",
                    negcontrol=negcontrols)

#Normalize report
NormalizeMethod <- c("pre", "nomis", "ccmn", "ruvrandom", "ruvrandclust", "is", "median", "mean", "sum")
NormalizeReport <- data.frame(method=NormalizeMethod, silh=rep(0,9), score=rep(0,9),
                              rscore=rep(0,9), difmet=rep(0,9), sd=rep(0,9), 
                              dd=rep(0,9), stringsAsFactors = F)

NormalizeReport[1,2] = mean(neibPreNormNormalizeMets$silh)
NormalizeReport[1,3] = mean(neibPreNormNormalizeMets$score)
NormalizeReport[1,4] = mean(neibPreNormNormalizeMets$rscore)
NormalizeReport[1,5] = sum(metvarPreNormNormalizeMets$metvar$padj<0.01)
NormalizeReport[1,6] = neibPreNormNormalizeMets$sdist
NormalizeReport[1,7] = neibPreNormNormalizeMets$ddist


NormalizeReport[2,2] = mean(neibNomisNormalizemets$silh)
NormalizeReport[2,3] = mean(neibNomisNormalizemets$score)
NormalizeReport[2,4] = mean(neibNomisNormalizemets$rscore)
NormalizeReport[2,5] = sum(metvarNomisNormalizeMets$metvar$padj<0.01)
NormalizeReport[2,6] = neibNomisNormalizemets$sdist
NormalizeReport[2,7] = neibNomisNormalizemets$ddist


NormalizeReport[3,2] = mean(neibCcmnNormalizeMets$silh)
NormalizeReport[3,3] = mean(neibCcmnNormalizeMets$score)
NormalizeReport[3,4] = mean(neibCcmnNormalizeMets$rscore)
NormalizeReport[3,5] = sum(metvarCcmnNormalizeMets$metvar$padj<0.01)
NormalizeReport[3,6] = neibCcmnNormalizeMets$sdist
NormalizeReport[3,7] = neibCcmnNormalizeMets$ddist

NormalizeReport[4,2] = mean(neibRuvRandNormalizeMets$silh)
NormalizeReport[4,3] = mean(neibRuvRandNormalizeMets$score)
NormalizeReport[4,4] = mean(neibRuvRandNormalizeMets$rscore)
NormalizeReport[4,5] = sum(metvarRuvRandNormalizeMets$metvar$padj<0.01)
NormalizeReport[4,6] = neibRuvRandNormalizeMets$sdist
NormalizeReport[4,7] = neibRuvRandNormalizeMets$ddist

NormalizeReport[5,2] = mean(neibRuvRandClustNormalizeMets$silh)
NormalizeReport[5,3] = mean(neibRuvRandClustNormalizeMets$score)
NormalizeReport[5,4] = mean(neibRuvRandClustNormalizeMets$rscore)
NormalizeReport[5,5] = sum(metvarRuvRandClustNormalizeMets$metvar$padj<0.01)
NormalizeReport[5,6] = neibRuvRandClustNormalizeMets$sdist
NormalizeReport[5,7] = neibRuvRandClustNormalizeMets$ddist

NormalizeReport[6,2] = mean(neibIsNormalizeMets$silh)
NormalizeReport[6,3] = mean(neibIsNormalizeMets$score)
NormalizeReport[6,4] = mean(neibIsNormalizeMets$rscore)
NormalizeReport[6,5] = sum(metvarIsNormalizeMets$metvar$padj<0.01, na.rm = T)
NormalizeReport[6,6] = neibIsNormalizeMets$sdist
NormalizeReport[6,7] = neibIsNormalizeMets$ddist

NormalizeReport[7,2] = mean(neibMedianNormalizeMets$silh)
NormalizeReport[7,3] = mean(neibMedianNormalizeMets$score)
NormalizeReport[7,4] = mean(neibMedianNormalizeMets$rscore)
NormalizeReport[7,5] = sum(metvarMedianNormalizeMets$metvar$padj<0.01)
NormalizeReport[7,6] = neibMedianNormalizeMets$sdist
NormalizeReport[7,7] = neibMedianNormalizeMets$ddist

NormalizeReport[8,2] = mean(neibMeanNormalizeMets$silh)
NormalizeReport[8,3] = mean(neibMeanNormalizeMets$score)
NormalizeReport[8,4] = mean(neibMeanNormalizeMets$rscore)
NormalizeReport[8,5] = sum(metvarMeanNormalizeMets$metvar$padj<0.01)
NormalizeReport[8,6] = neibMeanNormalizeMets$sdist
NormalizeReport[8,7] = neibMeanNormalizeMets$ddist

NormalizeReport[9,2] = mean(neibSumNormalizeMets$silh)
NormalizeReport[9,3] = mean(neibSumNormalizeMets$score)
NormalizeReport[9,4] = mean(neibSumNormalizeMets$rscore)
NormalizeReport[9,5] = sum(metvarSumNormalizeMets$metvar$padj<0.01)
NormalizeReport[9,6] = neibSumNormalizeMets$sdist
NormalizeReport[9,7] = neibSumNormalizeMets$ddist

forReport <- NormalizeReport[c(1, 2, 3, 4, 5, 7, 8),]

ggplot(data = forReport, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=NormalizeMethod[c(1, 2, 3, 4, 5, 7, 8)]),nudge_y=0.025) +
        labs(
                y= "Mean distance between strains",
                x= "Mean distance between replicates"
        )

ggsave("normmethod_dist_normalizeMets_newData.pdf")

silhvarNormalizeMets <- data.frame(method=rep("pre",length(neibPreNormNormalizeMets$silh)),
                                   silh=neibPreNormNormalizeMets$silh)
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("nomis",length(neibNomisNormalizemets$silh)),
                                         silh=neibNomisNormalizemets$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets, 
                              data.frame(method=rep("ccmn",length(neibCcmnNormalizeMets$silh)),
                                         silh=neibCcmnNormalizeMets$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("ruvrandom",length(neibRuvRandNormalizeMets$silh)),
                                         silh=neibRuvRandNormalizeMets$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("ruvrandclust",length(neibRuvRandClustNormalizeMets$silh)),
                                         silh=neibRuvRandClustNormalizeMets$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("median",length(neibMedianNormalizeMets$silh)),
                                         silh=neibMedianNormalizeMets$silh))
silhvarNormalizeMets <- rbind(silhvarNormalizeMets,
                              data.frame(method=rep("mean",length(neibMeanNormalizeMets$silh)),
                                         silh=neibMeanNormalizeMets$silh))


ggplot(data = silhvarNormalizeMets, mapping = aes(x = reorder(method, silh, FUN=median), y = silh)) +
        geom_boxplot() +
        #geom_line(data=normreport,mapping=aes(x=method,y=silh,group=1)) +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Silhouete",
                x= "Method"
        )

ggsave("normmethod_silh_normalizeMets_newData.pdf") 

#Histograms of pvalues to assess effectivity of normalization
get_hist <- function(p) {
        d <- ggplot_build(p)$data[[1]]
        data.frame(x = d$x, xmin = d$xmin, xmax = d$xmax, y = d$y)
}
set.seed(1)
x = runif(100, 0, 10)
y = cumsum(x)
dataframe <- data.frame(x = sort(x), y = y)

p <- ggplot(data = dataframe, aes(x = x)) + 
        geom_histogram(aes(y = cumsum(..count..)), binwidth = 1, boundary = 0,
                       color = "black", fill = "white")
p <- qplot(metvarPreNormNormalizeMets$metvar$padj, geom = "histogram",
           main = "Pre-normalization", 
           xlab = "p-value",
           ylab = "Frequency",
           fill = I("red"),
           alpha = I(.5))

hist = get_hist(p)
head(hist$x)

head(hist$y)

head(hist$xmax)

head(hist$xmin)


qplot(metvarPreNormNormalizeMets$metvar$padj, geom = "histogram",
      main = "Pre-normalization", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5), 
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_pre_newData.jpeg", device = "jpeg")

qplot(metvarNomisNormalizeMets$metvar$padj, geom = "histogram",
      main = "NOMIS", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_nomis_newData.jpeg", device = "jpeg")

p <- qplot(metvarCcmnNormalizeMets$metvar$padj, geom = "histogram",
           main = "CCMN", 
           xlab = "p-value",
           ylab = "Frequency",
           fill = I("red"),
           alpha = I(.5))

hist = get_hist(p)
head(hist$x)

head(hist$y)

head(hist$xmax)

head(hist$xmin)

qplot(metvarCcmnNormalizeMets$metvar$padj, geom = "histogram",
      main = "CCMN", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ccmn_newData.jpeg", device = "jpeg")

qplot(metvarRuvRandNormalizeMets$metvar$padj, geom = "histogram",
      main = "RUV Random", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ruvrand_newData.jpeg", device = "jpeg")

qplot(metvarRuvRandClustNormalizeMets$metvar$padj, geom = "histogram",
      main = "RUV Random Clust", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ruvrandclust_newData.jpeg", device = "jpeg")

qplot(metvarMedianNormalizeMets$metvar$padj, geom = "histogram",
      main = "Median", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_median_newData.jpeg", device = "jpeg")

qplot(metvarMeanNormalizeMets$metvar$padj, geom = "histogram",
      main = "Mean", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_mean_newData.jpeg", device = "jpeg")


#RLA Plots for assessing the effectiveness of the methods 

#Within-group RLA plot

RlaPlots(featuredata = missingMet$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Pre-normalization data",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_is_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, IS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_nomis_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, NOMIS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ccmn_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, CCMN",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ruv2_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, RUV2",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ruvrandom_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, RUV Random",
         saveplot = T, savetype = "tiff")
RlaPlots(featuredata = norm_ruvrandclust_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, RUV Random Clust",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_median_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Median",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_mean_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "wg",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Within-Group, Mean",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_sum_normalizemets$featuredata, 
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
RlaPlots(featuredata = norm_is_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, IS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_nomis_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, NOMIS",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ccmn_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, CCMN",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ruv2_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, RUV2",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ruvrandom_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, RUV Random",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_ruvrandclust_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, RUV Random Clust",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_median_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Median",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_mean_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Mean",
         saveplot = T, savetype = "tiff",
         ylim = c(-2, 2))
RlaPlots(featuredata = norm_sum_normalizemets$featuredata, 
         groupdata = sampledataMet$Group, type = "ag",
         interactiveplot = F, interactiveonly = F, 
         cols = rainbow(28), minoutlier = 0.5, las = 3,
         plotname = "Across-Group, Sum",
         saveplot = T, savetype = "tiff", showlegend = T,
         ylim = c(-2, 2))    

#Seems that CCMN is the best one

ccmn_norm_mets_newData <- norm_ccmn_normalizemets$featuredata
colnames(ccmn_norm_mets_newData) <- as.character(missingMet$metabolitedata$names[missingMet$metabolitedata$IS != 1])
save(ccmn_norm_mets_newData, file = "ccmn_norm_mets_newData.RData")

write.csv(norm_ccmn_normalizemets$featuredata, file = "normalized_ccmn_ND.csv")
write.csv(norm_nomis_normalizemets$featuredata, file = "normalized_nomis_ND.csv")
write.csv(norm_mean_normalizemets$featuredata, file = "normalized_mean_ND.csv")
write.csv(norm_ruvrandclust_normalizemets$featuredata, file = "normalized_ruvrandclust_ND.csv")