if(!require(UniProt.ws)) BiocManager::install("UniProt.ws")
library(UniProt.ws)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(factoextra)) install.packages("factoextra")
library(factoextra)
if(!require(rtracklayer)) install.packages("rtracklayer")
library(rtracklayer)
if(!require(ropls)) install.packages("ropls")
library(ropls)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal")

fileNames <- list.files("C:/Users/Guillem/Documents/PhD/comput/data/mtbcLineagesDels")
H37Rv_annot <- readGFF("C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot/genes.gff")

allDels <- data.frame()
for(i in 1:length(fileNames)){
        delsLin <- read.csv(paste("C:/Users/Guillem/Documents/PhD/comput/data/mtbcLineagesDels", fileNames[i], sep = "/"))
        rownames(delsLin) <- delsLin[, 1]
        delsLin <- delsLin[, 2:ncol(delsLin)]
        lin <- gsub("nr", "", strsplit(fileNames[i], "_|\\.")[[1]][2])
        colnames(delsLin) <- paste(colnames(delsLin), lin, sep = "_")
        if(nrow(allDels) == 0){
                allDels <- delsLin
        }else{
                allDels <- cbind.data.frame(allDels, delsLin)
        }
}

allDelsPCA <- prcomp(t(allDels))

fviz_pca_ind(allDelsPCA, habillage = sapply(colnames(allDels), function(x) strsplit(x, "_")[[1]][2]))

annot <- sapply(paste("mtu", rownames(allDels), sep = ":"), function(x){
        a <- keggGet(x)[[1]]$ORTHOLOGY 
        if(!is.null(a)){
                x <- a
        }else{
                x <- NA
        }
        })

getGeneDBIDs <- function(RvCodes){
        ort <- c()
        ncbi_gID <- c()
        ncbi_pID <- c()
        uniProt <- c()
        for(i in 1:length(RvCodes)){
                kGet <- keggGet(paste("mtu", RvCodes[i], sep = ":"))[[1]]
                a <- kGet$ORTHOLOGY[1]
                dbs <- kGet$DBLINKS
                if(!is.null(a)){
                        ort <- c(ort, a)
                }else{
                        ort <- c(ort, NA)
                }
                gID <- dbs[grep("GeneID", dbs)]
                pID <- dbs[grep("ProteinID", dbs)]
                up <- dbs[grep("UniProt", dbs)]
                if(length(gID) != 0){
                        ncbi_gID <- c(ncbi_gID, gsub("NCBI-GeneID: ", "", gID))
                }else{
                        ncbi_gID <- c(ncbi_gID, NA)
                }
                if(length(pID) != 0){
                        ncbi_pID <- c(ncbi_pID, gsub("NCBI-ProteinID: ", "", pID))
                }else{
                        ncbi_pID <- c(ncbi_pID, NA)
                }
                if(length(up) != 0){
                        uniProt <- c(uniProt, gsub("UniProt: ", "", up))
                }else{
                        uniProt <- c(uniProt, NA)
                }
        }
        annotDF <- data.frame(ncbi_geneID = ncbi_gID,
                              ncbi_proteinID = ncbi_pID,
                              uniProt = uniProt,
                              orthology = ort)
        
        rownames(annotDF) <- RvCodes
        return(annotDF)
}

annotDF <- getGeneDBIDs(rownames(allDels))

annotDF$product <- H37Rv_annot$product[match(rownames(annotDF), H37Rv_annot$locus_tag) + 1]

ECnum <- sub("\\].*", "", sub(".*\\[", "", as.character(annotDF$orthology))) 

ECnum[!grepl("EC:", ECnum)] <- NA

#strsplit(ECnum, " ")

annotDF$EC_number <- ECnum 

allDels <- cbind.data.frame(annotDF, allDels)

save(allDels, file = "allDels.RData")

allDelsPCA <- prcomp(t(allDels[, 7:ncol(allDels)]))

fviz_pca_ind(allDelsPCA, habillage = sapply(colnames(allDels)[7:ncol(allDels)], function(x) strsplit(x, "_")[[1]][2]))

save(allDelsPCA, file = "mtbcAllDelsPCA.RData")

enzDels <- allDels[!is.na(allDels$EC_number), ]

enzDelsPCA <- prcomp(t(enzDels[, 7:ncol(enzDels)]))

fviz_pca_ind(enzDelsPCA, habillage = sapply(colnames(enzDels)[7:ncol(enzDels)], function(x) strsplit(x, "_")[[1]][2]))

linGroup <- gsub("[[:digit:]]+", "", sapply(colnames(allDels)[7:ncol(allDels)], function(x) strsplit(x, "_")[[1]][2]))

enzDelsOPLSDA <- opls(t(enzDels[, 7:ncol(enzDels)]), linGroup)
