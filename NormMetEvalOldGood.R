setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script prepares the Mass Spectrometry metabolomic data to a format appropriate for "metabolomics" and           #################################
# "NormalizeMets" R packages, for then normalize it through different methods and later evaluate the performance of    #################################
# these methods through different approaches for, finally, exporting the data normalized by CCMN method (the best      #################################
# one) in a CSV file                                                                                                   #################################
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

#Load Names
sampleNames <- read_excel("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/metabolomics_jx/curated_peaks/sampleName1.xlsx")
colnames(sampleNames) <- paste(sampleNames[1, ], sampleNames[2, ], sep = '_')

#Load data
negativePeaks <- read_excel("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/metabolomics_jx/curated_peaks/ClinicalNeg-2.xlsx")
colnames(negativePeaks)[6:89] <- colnames(sampleNames)
negativePeaks$negativePeaks = 1
negativePeaks$mode <- "negative"

positivePeaks1 <- read_excel("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/metabolomics_jx/curated_peaks/ClinicalPos-4.1.xlsx")
colnames(positivePeaks1)[6:89] <- colnames(sampleNames)
positivePeaks1$positivePeaks1 = 1
positivePeaks1$mode <- "positive1"

positivePeaks2 <- read_excel("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/metabolomics_jx/curated_peaks/ClinicalPos-4.2.xlsx")
colnames(positivePeaks2)[6:89] <- colnames(sampleNames)
positivePeaks2$positivePeaks2 = 1
positivePeaks2$mode <- "positive2"

metabolitePeaks <- merge(negativePeaks, positivePeaks1, all = T)
metabolitePeaks <- merge(metabolitePeaks, positivePeaks2, all = T)
metabolitePeaks[is.na(metabolitePeaks)] <- 0

upset(metabolitePeaks)

#Drop variables that we don't want
drops <- c("CAS ID", "negativePeaks", "positivePeaks1", "positivePeaks2")
metabolitePeaks <- metabolitePeaks[, !(names(metabolitePeaks) %in% drops)]
metabolitePeaks_preCur <- metabolitePeaks
save(metabolitePeaks_preCur, file = "metabolitePeaks_oldData.RData")


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
sampleNames <- t(sampleNames)
sampleNames <- sampleNames[, c(1, 2)]
colnames(sampleNames) <- c("Group", "Sample")
sampleNames <- data.frame(sampleNames[, c("Sample", "Group")])

rownames(metabolitePeaks) <- metabolitePeaks[, 1]
metabolomicsPackageFormat <- data.frame(t(metabolitePeaks[, -c(1, 2, 3, 4, 89)]))
metabolomicsPackageFormat <- merge(sampleNames, metabolomicsPackageFormat, by = "row.names",
                                   all.x = T, sort = F)
metabolomicsPackageFormat <- metabolomicsPackageFormat[, -c(1, 2)]
row.names(metabolomicsPackageFormat) <- row.names(sampleNames)
colnames(metabolomicsPackageFormat)[2:ncol(metabolomicsPackageFormat)] <- metabolitePeaks$`Compound Name`
metabolomicsPackageFormat

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
        ggsave("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/normMetAnal/oldDataGood/meanPeakAreaRaw_oldData.jpg")

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

save(metabolomicsPackageFormat, file = "metabolomicsPackageFormat.RData")
########################################################################################################################################################
#
# Normalize with different normalization methods and compare them with "metabolomics" Package
#
########################################################################################################################################################
if(!require(KEGGREST)) BiocManager::install("KEGGREST")
library("KEGGREST")

#log10 transform data 
log10meatabolomics <- log10(metabolomicsPackageFormat[,2:97])

#Calculate means for each group
meanvalues_ag <- aggregate(metabolomicsPackageFormat[,2:97],list(metabolomicsPackageFormat[,1]),mean)

#Calculate the mean,, dif and cv matrixes both for normal data & log10 data
meanvalues <- matrix(nrow=28,ncol=96)
meandifvalues <- matrix(nrow=28*3,ncol=96)
meancvvalues <- matrix(nrow=28*3,ncol=96)

meanlogvalues <- matrix(nrow=28,ncol=96)
meandiflogvalues <- matrix(nrow=28*3,ncol=96)
meancvlogvalues <- matrix(nrow=28*3,ncol=96)

for (i in 1:28){
        meanvalues[i, ] <- colMeans(metabolomicsPackageFormat[(1+3*(i-1)):(3*i), 2:97], na.rm = T)
        meanlogvalues[i, ] <- colMeans(log10meatabolomics[(1+3*(i-1)):(3*i), ], na.rm = T)
        for (j in 1:96) {
                meandifvalues[(1+3*(i-1)):(3*i), j] <- 
                        metabolomicsPackageFormat[(1+3*(i-1)):(3*i), j+1] - meanvalues[i, j]
                meancvvalues[(1+3*(i-1)):(3*i),j] <- 
                        (metabolomicsPackageFormat[(1+3*(i-1)):(3*i), j+1] - meanvalues[i, j])/meanvalues[i, j]
                meandiflogvalues[(1+3*(i-1)):(3*i),j] <- 
                        log10meatabolomics[(1+3*(i-1)):(3*i), j] - meanlogvalues[i, j]
                meancvlogvalues[(1+3*(i-1)):(3*i), j] <- 
                        (log10meatabolomics[(1+3*(i-1)):(3*i), j]-meanlogvalues[i, j])/meanlogvalues[i, j]
        }
}

#Calculate statisics and store in dataframe
ssr <- vector()
sst <- vector()
ssbet <- vector()
n1 <- vector()
n2 <- vector()
fstat <- vector()
kruskals <- vector()
kruskalp <-vector()
metmean <- colMeans(log10meatabolomics, na.rm = T)

for (i in 1:96) {
        ssr[i] <- sum(meandiflogvalues[, i]^2, na.rm=T)
        sst[i] <- sum((log10meatabolomics[, i]-metmean[i])^2, na.rm=T)
        ssbet[i] <- sst[i]-ssr[i]
        n1[i] <- sum(!is.na(meanlogvalues[, i]))-1
        n2[i] <- sum(!is.na(log10meatabolomics[, i]))-n1[i]
        fstat[i] <- (ssbet[i]/n1[i])/(ssr[i]/n2[i])
        datavec <- log10meatabolomics[, i]
        groupvec <- metabolomicsPackageFormat[!is.na(datavec), 1]
        datavec <- datavec[!is.na(datavec)]
        kresult <- kruskal.test(datavec,groupvec)
        kruskalp[i] <- kresult$p.value
        kruskals[i] <- kresult$statistic
}
metvariation <- data.frame(ssr, ssbet, sst, n1, n2, fstat, kruskals, kruskalp)

#Make a function that does the same than the previous loop
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

#substitute NAs and log transform data
missing <- MissingValues(metabolomicsPackageFormat, group.cutoff = 0.2, column.cutoff = 0.8)$output
log.data <- LogTransform(missing,base = 10)$output

if(!require(cluster)) install.packages("cluster")
library(cluster)
if(!require(ruv)) install.packages("ruv")
library(ruv)

#Check variance and distance between groups before normalization
metvarprenorm <- metvar(log.data[, 2:97], log.data[, 1])
neibprenorm <- neibcheck(log.data[, 2:97])

#Test different normalization methods
norm_med <- Normalise(log.data,method="median")$output
metvarmed=metvar(norm_med[,2:97],log.data[,1])
neibmed=neibcheck(norm_med[,2:97])

norm_mean <- Normalise(log.data,method="mean")$output
metvarmean=metvar(norm_mean[,2:97],log.data[,1])
neibmean=neibcheck(norm_mean[,2:97])

norm_sum <- Normalise(log.data,method="sum")$output
metvarsum=metvar(norm_sum[,2:97],log.data[,1])
neibsum=neibcheck(norm_sum[,2:97])

norm_ccmn <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$f<2))$output
metvarccmn=metvar(norm_ccmn[,2:88],log.data[,1])
neibccmn=neibcheck(norm_ccmn[,2:88])

norm_ccmn05 <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$padj>0.05))$output
metvarccmn05=metvar(norm_ccmn05[,2:87],log.data[,1])
neibccmn05=neibcheck(norm_ccmn05[,2:87])

norm_ccmn05_3 <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$padj>0.05 & metvarprenorm$metvar$mean>3))$output
metvarccmn05_3=metvar(norm_ccmn05_3[,2:80],log.data[,1])
neibccmn05_3=neibcheck(norm_ccmn05_3[,2:80])

norm_ccmn01_3 <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$padj>0.01 & metvarprenorm$metvar$mean>3))$output
metvarccmn01_3=metvar(norm_ccmn01_3[,2:80],log.data[,1])
neibccmn01_3=neibcheck(norm_ccmn01_3[,2:80])

norm_ccmn05_2 <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$padj>0.05 & metvarprenorm$metvar$mean>2))$output
metvarccmn05_2=metvar(norm_ccmn05_2[,2:87],log.data[,1])
neibccmn05_2=neibcheck(norm_ccmn05_2[,2:87])

norm_ccmn01_2 <- Normalise(log.data,method="ccmn",nc=(metvarprenorm$metvar$padj>0.01 & metvarprenorm$metvar$mean>2))$output
metvarccmn01_2=metvar(norm_ccmn01_2[,2:80],log.data[,1])
neibccmn01_2=neibcheck(norm_ccmn01_2[,2:80])

norm_nomis05 <- Normalise(log.data,method="nomis",nc=(metvarprenorm$metvar$padj>0.05))$output
metvarnomis05=metvar(norm_nomis05[,2:87],log.data[,1])
neibnomis05=neibcheck(norm_nomis05[,2:87])

norm_nomis05_3 <- Normalise(log.data,method="nomis",nc=(metvarprenorm$metvar$padj>0.05 & metvarprenorm$metvar$mean>3))$output
metvarnomis05_3=metvar(norm_nomis05_3[,2:87],log.data[,1])
neibnomis05_3=neibcheck(norm_nomis05_3[,2:87])

#Builds a factor matrix of 0s and 1s that tells to what strain (group) the sample corresponds
facmat=matrix(0, nrow=84, ncol=28)
for (i in 1:28){
        facmat[(1+3*(i-1)):(3*i),i]=1
}

norm_ruv2_05=LinearModelFit(log.data[,2:97],
                            factormat = facmat, covariatemat = NULL, contrastmat = NULL,
                            ruv2 = TRUE, k = 4, nc = metvarprenorm$metvar$padj>0.05,
                            moderated = FALSE, padjmethod = "BH",
                            saveoutput = FALSE, outputname = "results")

boxplot(t(log.data[,2:97]),ylim=c(0,6.5),xaxt="n")
axis(1,at=seq(1,82,3)+1,labels=log.data[seq(1,82,3),1],pos=-0.27,las=2)

boxplot((log.data[,2:97]),ylim=c(0,6.5),xaxt="n")

metnames=colnames(log.data[,2:97])
samplename=row.names(log.data)
logdatatable=data.frame(sample=rep('.',84*96),strain=rep('.',84*96),metabolite=rep('.',84*96),peak=rep(0,84*96),significant=rep(0,84*96),stringsAsFactors = F)
line=0
for (i in 1:84){
        for (j in 2:97){
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

ggsave("sample_boxplot_oldData.pdf",width = 9, height = 5)  


ggplot(data = logdatatable, mapping = aes(x = reorder(metabolite, peak, FUN=median), y = peak, color=significant)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Abundance",
                x= "Metabolite"
        )

ggsave("metabolite_boxplot_oldData.pdf",width = 7, height = 11)  

save(logdatatable, file = "logdatatable.RData")

ggplot(data = logdatatable[logdatatable$significant==0,], mapping = aes(x = strain, y = peak)) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90)) +
        facet_grid(metabolite ~ ., scales="free") +
        labs(
                y= "Abundance",
                x= "Strain"
        )

ggsave("const_metaboliteds_oldData.pdf",width = 7, height = 14)  

nonvarmet=metnames[metvarprenorm$metvar$padj>0.05]
ncset2=which(is.element(metnames,nonvarmet))
nonvarmet=nonvarmet[c(1,4,5,6,7,8)]
ncset=which(is.element(metnames,nonvarmet))



norm_med <- Normalise(log.data,method="median")$output
metvarmed=metvar(norm_med[,2:97],log.data[,1])
neibmed=neibcheck(norm_med[,2:97])

norm_mean <- Normalise(log.data,method="mean")$output
metvarmean=metvar(norm_mean[,2:97],log.data[,1])
neibmean=neibcheck(norm_mean[,2:97])

norm_sum <- Normalise(log.data,method="sum")$output
metvarsum=metvar(norm_sum[,2:97],log.data[,1])
neibsum=neibcheck(norm_sum[,2:97])

norm_ccmn <- Normalise(log.data,method="ccmn",nc=ncset)$output
metvarccmn=metvar(norm_ccmn[,2:91],log.data[,1])
neibccmn=neibcheck(norm_ccmn[,2:91])

norm_ccmn2 <- Normalise(log.data,method="ccmn",nc=ncset2)$output
metvarccmn2=metvar(norm_ccmn2[,2:87],log.data[,1])
neibccmn2=neibcheck(norm_ccmn2[,2:87])


norm_nomis <- Normalise(log.data,method="nomis",nc=ncset)$output
metvarnomis=metvar(norm_nomis[,2:91],log.data[,1])
neibnomis=neibcheck(norm_nomis[,2:91])

norm_is <- Normalise(log.data,method="is",refvec=ncset[3])$output
metvaris=metvar(norm_is[,2:97],log.data[,1])
neibis=neibcheck(norm_is[,2:97])

normmethod=c("pre","median","mean","sum","ccmn","nomis","is","ccmn2", "ruv2")

normreport=data.frame(method=normmethod, silh=rep(0,9), score=rep(0,9), rscore=rep(0,9),difmet=rep(0,9),sd=rep(0,9),dd=rep(0,9), stringsAsFactors = F)

normreport[1,2]=mean(neibprenorm$silh)
normreport[1,3]=mean(neibprenorm$score)
normreport[1,4]=mean(neibprenorm$rscore)
normreport[1,5]=sum(metvarprenorm$metvar$padj<0.01)
normreport[1,6]=neibprenorm$sdist
normreport[1,7]=neibprenorm$ddist


normreport[2,2]=mean(neibmed$silh)
normreport[2,3]=mean(neibmed$score)
normreport[2,4]=mean(neibmed$rscore)
normreport[2,5]=sum(metvarmed$metvar$padj<0.01)
normreport[2,6]=neibmed$sdist
normreport[2,7]=neibmed$ddist


normreport[3,2]=mean(neibmean$silh)
normreport[3,3]=mean(neibmean$score)
normreport[3,4]=mean(neibmean$rscore)
normreport[3,5]=sum(metvarmean$metvar$padj<0.01)
normreport[3,6]=neibmean$sdist
normreport[3,7]=neibmean$ddist

normreport[4,2]=mean(neibsum$silh)
normreport[4,3]=mean(neibsum$score)
normreport[4,4]=mean(neibsum$rscore)
normreport[4,5]=sum(metvarsum$metvar$padj<0.01)
normreport[4,6]=neibsum$sdist
normreport[4,7]=neibsum$ddist

normreport[5,2]=mean(neibccmn$silh)
normreport[5,3]=mean(neibccmn$score)
normreport[5,4]=mean(neibccmn$rscore)
normreport[5,5]=sum(metvarccmn$metvar$padj<0.01)
normreport[5,6]=neibccmn$sdist
normreport[5,7]=neibccmn$ddist

normreport[6,2]=mean(neibnomis$silh)
normreport[6,3]=mean(neibnomis$score)
normreport[6,4]=mean(neibnomis$rscore)
normreport[6,5]=sum(metvarnomis$metvar$padj<0.01)

normreport[6,6]=neibnomis$sdist
normreport[6,7]=neibnomis$ddist

normreport[7,2]=mean(neibis$silh)
normreport[7,3]=mean(neibis$score)
normreport[7,4]=mean(neibis$rscore)
normreport[7,5]=sum(metvaris$metvar$padj<0.01)
normreport[7,6]=neibis$sdist
normreport[7,7]=neibis$ddist




normreport[8,2]=mean(neibccmn2$silh)
normreport[8,3]=mean(neibccmn2$score)
normreport[8,4]=mean(neibccmn2$rscore)
normreport[8,5]=sum(metvarccmn2$metvar$padj<0.01)
normreport[8,6]=neibccmn2$sdist
normreport[8,7]=neibccmn2$ddist

fnormreport=normreport[c(1,2,3,5,6),]

ggplot(data = fnormreport, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=method),nudge_y=0.025) +
        labs(
                y= "Mean distance between strains",
                x= "Mean distance between replicates"
        )

ggsave("normmethod_dist_oldData.pdf")

silhvar=data.frame(method=rep("pre",length(neibprenorm$silh)),silh=neibprenorm$silh)
silhvar=rbind(silhvar,data.frame(method=rep("median",length(neibmed$silh)),silh=neibmed$silh))
silhvar=rbind(silhvar,data.frame(method=rep("mean",length(neibmean$silh)),silh=neibmean$silh))
silhvar=rbind(silhvar,data.frame(method=rep("sum",length(neibsum$silh)),silh=neibsum$silh))
silhvar=rbind(silhvar,data.frame(method=rep("ccmn",length(neibccmn$silh)),silh=neibccmn$silh))
silhvar=rbind(silhvar,data.frame(method=rep("nomis",length(neibnomis$silh)),silh=neibnomis$silh))
silhvar=rbind(silhvar,data.frame(method=rep("is",length(neibis$silh)),silh=neibis$silh))

ggplot(data = silhvar, mapping = aes(x = reorder(method, silh, FUN=median), y = silh)) +
        geom_boxplot() +
        #geom_line(data=normreport,mapping=aes(x=method,y=silh,group=1)) +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Silhouete",
                x= "Method"
        )

ggsave("normmethod_silh_oldData.pdf") 

scorevar=data.frame(method=rep("pre",length(neibprenorm$score)),score=neibprenorm$score)
scorevar=rbind(scorevar,data.frame(method=rep("median",length(neibmed$score)),score=neibmed$score))
scorevar=rbind(scorevar,data.frame(method=rep("mean",length(neibmean$score)),score=neibmean$score))
scorevar=rbind(scorevar,data.frame(method=rep("sum",length(neibsum$score)),score=neibsum$score))
scorevar=rbind(scorevar,data.frame(method=rep("ccmn",length(neibccmn$score)),score=neibccmn$score))
scorevar=rbind(scorevar,data.frame(method=rep("nomis",length(neibnomis$score)),score=neibnomis$score))
scorevar=rbind(scorevar,data.frame(method=rep("is",length(neibis$score)),score=neibis$score))

ggplot(data = scorevar, mapping = aes(x = method, y = score)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Score",
                x= "Method"
        )

metnames2=colnames(norm_ccmn[,2:91])
samplename2=row.names(norm_ccmn)
ccmndatatable=data.frame(sample=rep('.',84*91),strain=rep('.',84*91),
                         metabolite=rep('.',84*91),peak=rep(0,84*91),significant=rep(0,84*91),stringsAsFactors = F)
line=0
for (i in 1:84){
        for (j in 2:91){
                line=line+1
                ccmndatatable[line,1]=samplename2[i]
                ccmndatatable[line,2]=as.character(norm_ccmn[i,1])
                ccmndatatable[line,3]=metnames2[j-1]
                ccmndatatable[line,4]=norm_ccmn[i,j]
                if (metvarccmn$metvar$padj[metvarccmn$metvar$names==metnames2[j-1]]<0.05){
                        ccmndatatable[line,5]=1
                }
        }
}


ggplot(data = ccmndatatable, mapping = aes(x = sample, y = peak)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90)) +
        #coord_flip() +
        labs(
                y= "Abundance",
                x= "Sample"
        )

ggsave("sample_boxplot_ccmn_oldData.pdf",width = 9, height = 5)  


ggplot(data = ccmndatatable, mapping = aes(x = reorder(metabolite, peak, FUN=median), y = peak, color=significant)) +
        geom_boxplot() +
        coord_flip() +
        labs(
                y= "Abundance",
                x= "Metabolite"
        )

ggsave("metabolite_boxplot_ccmn_oldData.pdf",width = 7, height = 11)  

ccmncor=as.dist(1-cor(norm_ccmn[,2:91]))
ccmnlink=hclust(ccmncor,method="average")
plot(ccmnlink, hang=-1)


#RUVrandom method
if(!require(MetNorm)) install.packages("MetNorm")
library(MetNorm)
if(!require(impute)) BiocManager::install("impute")
library(impute)

randomRUV <- NormalizeRUVRand(Y = data.matrix(log.data[, 2:97]), ctl = metvarprenorm$metvar$f<2, k = "3")
forclust <- NormalizeRUVRandClust(data.matrix(randomRUV$newY), maxIter = 200)
data.matrix(log.data[, 2:97])
length(metvarprenorm$metvar$padj)

if(!require(RUVnormalize)) BiocManager::install("RUVnormalize")
library("RUVnormalize")
naiveRandRUV(Y = data.matrix(log.data[, 2:97]), cIdx = which(metvarprenorm$metvar$f<2), k = "6")

########################################################################################################################################################
#
# Normalize with different methods and analyze the performance with "NormalizeMets" Package
#
########################################################################################################################################################

detach("package:metabolomics", unload = T)
if(!require(NormalizeMets)) install.packages("NormalizeMets")
library(NormalizeMets)

#Adapt to format

featuredataMet <- metabolomicsPackageFormat[, 2:97]
sampledataMet <- data.frame(Group = metabolomicsPackageFormat[, 1], 
                            Species = c(rep("P. aeruginosa", length(metabolomicsPackageFormat[, 1]))))

IS <- c(rep(0, 96))
for(i in 1:96){
        if(metvarprenorm$metvar$padj[i]>0.05){
                IS[i] <- 1
        }
}

metabolitedataMet <- data.frame(names = colnames(featuredataMet), IS = IS)

alldataMet <- list(featuredata = featuredataMet, sampledata = sampledataMet, 
                   metabolitedata = metabolitedataMet)

#LogTransform Data

logdataMet<- LogTransform(alldataMet$featuredata, zerotona = T)
dataview(logdataMet$featuredata)

#Replace Missing

missingMet <-  MissingValues(logdataMet$featuredata, sampledataMet,metabolitedataMet,
                             feature.cutof=0.8, sample.cutoff=0.2, method="knn")
#In this step metabolites 17 and 95 ("X4.Amino.4.cyanobutanoic.acid" & "UDP.N.acetylmuramate") are removed (have a lot of missing values)
dataview(missingMet$featuredata)

#Normalization by different methods

####### On the normalization I was doing before the qcmets was set according to metabolitedataMet object. As the MissingValues function removes 2 metabolites
####### ("X4.Amino.4.cyanobutanoic.acid" & "UDP.N.acetylmuramate") from LogDataMet object, and we're applying normalization to missingMet object, the nonchanging
####### metabolites were shifted, so instead of normalizing to "X.Acetyl.L.glutamate.5.semialdehyde", "Acetolactate", "Glucose", "gly" and "lactate" we were normalizing 
####### to "X.Acetyl.L.glutamate.5.semialdehyde", "X4.Aminobutanoate.3", "Acetyl.2.6.diaminopimelate", "glutamate", "Glycerate", "Leucine", "methylmaleate.1", 
####### "S.Adenosylmethionine" and "Tyrosine". Interestingly, with that last set of metabolites (the wrong one), after normalization all replicates clustered together, while
####### when normalizing with the correct set of nonchanging metabolites, there are 3 replicates that are not clustering together.


# This is the submatrix of the metabolites that were misinterpreted as nonchanging.
View(missingMet$featuredata[, c("X.Acetyl.L.glutamate.5.semialdehyde", 
                                "X4.Aminobutanoate.3",
                                "Acetyl.2.6.diaminopimelate", 
                                "glutamate",
                                "Glycerate", 
                                "Leucine", 
                                "methylmaleate.1", 
                                "S.Adenosylmethionine", 
                                "Tyrosine")])

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
facmat=matrix(0, nrow=84, ncol=28)
for (i in 1:28){
        facmat[(1+3*(i-1)):(3*i),i]=1
}

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


unadjFit <- LinearModelFit(log.data[, 2:97], factormat = facmat, ruv2 = F)
normisFit <- LinearModelFit(norm_is_normalizemets$featuredata, factormat = facmat, ruv2 = F)
normccmnFit <- LinearModelFit(norm_ccmn_normalizemets$featuredata, factormat = facmat, ruv2 = F)
unadjFit$coefficients[, "x1"]
lcoef_x1 <- list(unadjusted=unadjFit$coefficients[,"x1"],
                 is_x1=normisFit$coefficients[,"x1"],
                 ccmn_x1=normccmnFit$coefficients[,"x1"])
lpvals_x1 <- list(unadjusted=unadjFit$p.value[,"x1"],
                  is=normisFit$p.value[,"x1"],
                  ccmn=normccmnFit$p.value[,"x1"])

negcontrols <- colnames(log.data)[(which(metvarprenorm$metvar$f<2)) + 1]

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

tiff("normmethod_dist_normalizeMets_oldData.tiff", res = 300, height = 1500, width = 1500)
ggplot(data = forReport, mapping = aes(x = sd, y = dd)) +
        geom_point() +
        geom_label(aes(label=NormalizeMethod[c(1, 2, 3, 4, 5, 7, 8)]),nudge_y=0.025) +
        labs(
                y= "Mean distance between strains",
                x= "Mean distance between replicates"
        )
dev.off()

ggsave("normmethod_dist_normalizeMets_oldData.pdf")
ggsave("normmethod_dist_normalizeMets_oldData.tiff")

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

tiff("normmethod_silh_normalizeMets_oldData.tiff", res = 300, height = 1500, width = 1500)
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
ggsave("normmethod_silh_normalizeMets_oldData.pdf") 

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

ggsave("pvaluehist_pre_oldData.jpeg", device = "jpeg")

qplot(metvarNomisNormalizeMets$metvar$padj, geom = "histogram",
      main = "NOMIS", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_nomis_oldData.jpeg", device = "jpeg")

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

ggsave("pvaluehist_ccmn_oldData.jpeg", device = "jpeg")

qplot(metvarRuvRandNormalizeMets$metvar$padj, geom = "histogram",
      main = "RUV Random", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ruvrand_oldData.jpeg", device = "jpeg")

qplot(metvarRuvRandClustNormalizeMets$metvar$padj, geom = "histogram",
      main = "RUV Random Clust", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_ruvrandclust_oldData.jpeg", device = "jpeg")

qplot(metvarMedianNormalizeMets$metvar$padj, geom = "histogram",
      main = "Median", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_median_oldData.jpeg", device = "jpeg")

qplot(metvarMeanNormalizeMets$metvar$padj, geom = "histogram",
      main = "Mean", 
      xlab = "p-value",
      ylab = "Frequency",
      fill = I("red"),
      alpha = I(.5),
      xlim = c(-0.0106, 0.6),
      ylim = c(0, 85))

ggsave("pvaluehist_mean_oldData.jpeg", device = "jpeg")

#ranks in function of t-value for nochanging metabolites in each strain, for each one
#of the tested methods.
nonchanging <- c("X.Acetyl.L.glutamate.5.semialdehyde", "4.Amino.4.cyanobutanoic.acid", "Acetolactate",
                 "Glucose", "gly", "lactate", "Methylglyoxal", "Pyruvate", "Tyramine")
ranksAllStrains <- list()
for(strain in 1:28){
        preRanks <- c()
        for(i in rownames(unadjFit$t)){
                for(j in nonchanging){
                        if(i == j){
                                preRanks[i] <- rank(unadjFit$t[, strain])[i]
                        }
                }
        }
        
        nomisRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normNomisFit$t)){
                        if(j == i){
                                nomisRanks[i] <- rank(normNomisFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normNomisFit$t)) == F){nomisRanks[j] <- NA}
                }
        }
        
        ccmnRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normCcmnFit$t)){
                        if(j == i){
                                ccmnRanks[i] <- rank(normCcmnFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normCcmnFit$t)) == F){ccmnRanks[j] <- NA}
                }
        }
        
        ruv2Ranks <- c()
        for(j in nonchanging){
                for(i in rownames(normRuv2Fit$t)){
                        if(i == j){
                                ruv2Ranks[i] <-rank(normRuv2Fit$t[, strain])[i]
                        }else if(is.element(j, rownames(normRuv2Fit$t)) == F){ruv2Ranks[j] <- NA}
                }
        }
        
        ruvRandRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normRuvRandFit$t)){
                        if(i == j){
                                ruvRandRanks[i] <- rank(normRuvRandFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normRuvRandFit$t)) == F){ruvRandRanks[j] <- NA}
                }
        }
        
        ruvRandClustRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normRuvRandClustFit$t)){
                        if(i == j){
                                ruvRandClustRanks[i] <- rank(normRuvRandClustFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normRuvRandClustFit$t)) == F){ruvRandClustRanks[j] <- NA}
                }
        }
        
        isRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normIsFit$t)){
                        if(i == j){
                                isRanks[i] <- rank(normIsFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normisFit$t)) == F){isRanks[j] <- NA}
                }
        }
        
        medianRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normMedianFit$t)){
                        if(i == j){
                                medianRanks[i] <- rank(normMedianFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normMedianFit$t)) == F){medianRanks[j] <- NA}
                }
        }
        
        meanRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normMeanFit$t)){
                        if(i == j){
                                meanRanks[i] <- rank(normMeanFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normMeanFit$t)) == F){meanRanks[j] <- NA}
                }
        }
        
        sumRanks <- c()
        for(j in nonchanging){
                for(i in rownames(normSumFit$t)){
                        if(i == j){
                                sumRanks[i] <- rank(normSumFit$t[, strain])[i]
                        }else if(is.element(j, rownames(normSumFit$t)) == F){sumRanks[j] <- NA}
                }
        }
        rankNonChangMets <- rbind(preRanks, nomisRanks, ccmnRanks, 
                                  ruv2Ranks, ruvRandRanks, ruvRandClustRanks,
                                  isRanks, medianRanks, meanRanks, sumRanks)
        
        rownames(rankNonChangMets) <- c("pre", "NOMIS", "CCMN", "RUV2", 
                                        "RUV Random", "RUV Random Clust", 
                                        "IS", "Median", "Mean","Sum")
        ranksAllStrains[[strain]] <- rankNonChangMets
}
names(ranksAllStrains) <- unique(metabolomicsPackageFormat$Group)
sink("ranksAllStrainsMeths_oldData.csv", type = "output")
invisible(lapply(names(ranksAllStrains),
                 function(x) {print(x) 
                         dput(write.csv(ranksAllStrains[[x]]))}))
sink()
ranksNonChangMean <- apply(simplify2array(ranksAllStrains), 1:2, mean)
ranksNonChangSD <- apply(simplify2array(ranksAllStrains), 1:2, sd)

ncmets <- c()
for(i in nonchanging){
        ncmets <- append(ncmets, rep(i, 10))
}

ranks4hist <- data.frame(NCmetabolites = ncmets, 
                         method = rep(rownames(ranksNonChangMean), 9), 
                         ranks = c(ranksNonChangMean),
                         SD = c(ranksNonChangSD))






ggplot(data = ranks4hist, mapping = aes(x=ranks4hist$method, 
                                        y=ranks4hist$ranks, 
                                        fill=ranks4hist$NCmetabolites)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = ranks4hist$ranks - ranks4hist$SD, 
                          ymax = ranks4hist$ranks + ranks4hist$SD), width = .2,
                      position = position_dodge(.9)) +
        coord_flip() +
        theme(legend.position = "bottom", 
              legend.direction = "vertical", 
              legend.text = element_text(size=8),
              legend.title = element_text(size = 9),
              legend.key.size = unit(0.75, "line"),
              legend.justification = c(-1.2, 1),
              legend.margin = margin(t = 0, unit = 'cm')) +
        guides(linetype = guide_legend(nrow = 3)) +
        scale_x_discrete(limits = c("pre", "IS", "NOMIS", "CCMN", "RUV2", 
                                    "RUV Random", "RUV Random Clust", 
                                    "Median", "Mean", "Sum")) +
        ggtitle("Mean ranks accoring t-statistic \nof nonchanging metabolites \nfor each normalization method") +
        labs(
                x = "Normalization Method",
                y = "Mean Rank",
                fill = "Nonchanging metabolites"
        )
ggsave("meanRanksTstatVert_oldData.pdf")
ggsave("meanRanksaTstatVert_oldData.tiff", device = "tiff")   
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
ccmn_norm_mets_good_old <- norm_ccmn_normalizemets$featuredata
colnames(ccmn_norm_mets_good_old) <- as.character(missingMet$metabolitedata$names[missingMet$metabolitedata$IS != 1])
save(ccmn_norm_mets_good_old, file = "ccmn_norm_mets_good_old.RData")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")

colnames(ccmn_norm_mets_good_old) <- dictionary$Consensus[match(colnames(ccmn_norm_mets_good_old), dictionary$`Old Data Names`)]
ccmn_norm_mets_good_old <- ccmn_norm_mets_good_old[, !is.na(colnames(ccmn_norm_mets_good_old))]

write.csv(ccmn_norm_mets_good_old, file = "normalized_ccmn_good.csv")
write.csv(dictionary, file = "dictionary.csv")

write.csv(norm_nomis_normalizemets$featuredata, file = "normalized_nomis_good.csv")
write.csv(norm_mean_normalizemets$featuredata, file = "normalized_mean_good.csv")
write.csv(norm_ruvrandclust_normalizemets$featuredata, file = "normalized_ruvrandclust_good.csv")
