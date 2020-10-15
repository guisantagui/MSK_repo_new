
########################################################################################
### Create a list of metabolites that we want to see in MTBC strains' LC-MS analysis ###
########################################################################################
if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(RMassBank)) BiocManager::install("RMassBank")
library(RMassBank)
if(!require(R.cache)) install.packages("R.cache")
library(R.cache)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcCPDs4MS")

# Import model dataframes

sMtb_mets <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/models/MTBC/sMtb.xlsx", sheet = 2))
sMtb_mets

m2155_mets <- read.csv("m2155_mets.csv")

mtbH37Rv_mets <- read.csv("mtbH37Rv_mets.csv")

iJO1366_mets <- as.data.frame(readxl::read_xls("C:/Users/Guillem/Documents/PhD/comput/models/Ecoli/iJO1366_info.xls", sheet = 4))

# Import modelSEED compound dataframe

modelSEED <- as.data.frame(fread('C:/Users/Guillem/Documents/PhD/comput/data/ModelSEEDDatabase-master/Biochemistry/compounds.tsv'))

getModelSEEDKEGGIDS <- function(x){
        splitAlias <- strsplit(x, split = "|", fixed = T)[[1]]
        posKEGGID <- grep("KEGG", splitAlias)
        if(length(posKEGGID) > 0){
                KEGGID <- gsub("KEGG: ", "", splitAlias[posKEGGID])
        }else{
                KEGGID <- NA
        }
        return(KEGGID)
}

modelSEED$KEGGIDs <- sapply(modelSEED$aliases, getModelSEEDKEGGIDS)

View(modelSEED)

sMtb_mets$`KeGG ID`[sMtb_mets$`KeGG ID` == "none"] <- NA

KEGGids2SEEDids <- function(KEGGID){
        if(!require(R.cache)) install.packages("R.cache") 
        library(R.cache)
        key = list("seedKEGGIds")
        SEED_KEGGIDs <- loadCache(key = key)
        if(is.null(SEED_KEGGIDs)){
                SEED_KEGGIDs <- sapply(modelSEED$KEGGIDs,
                                       function(x) strsplit(x, 
                                                            split = "; ")[[1]][1])
                saveCache(SEED_KEGGIDs, key = key)
        }
        KEGGID <- strsplit(KEGGID, split = "; ")[[1]][1]
        if(!is.na(KEGGID)){
                SEEDid <- modelSEED$id[match(KEGGID, SEED_KEGGIDs)]
        }else{
                SEEDid <- NA
        }
        return(SEEDid)
}

SEEDids2Mass <- function(SEEDID){
        if(!require(R.cache)) install.packages("R.cache") 
        library(R.cache)
        key = list("seedIds")
        SEEDIds <- loadCache(key = key)
        if(is.null(SEEDIds)){
                SEEDIds <- modelSEED$id
                saveCache(SEEDIds, key = key)
        }
        if(!is.na(SEEDID)){
                SEEDMASS <- modelSEED$mass[match(SEEDID, SEEDIds)]
        }else{
                SEEDMASS <- NA
        }
        return(SEEDMASS)
}

KEGGid2MolWeightNExactMass <- function(KEGGID){
        if(!is.na(KEGGID)){
                get <- keggGet(KEGGID)
                molWeight <- as.numeric(get[[1]]$MOL_WEIGHT)
                exactMass <- as.numeric(get[[1]]$EXACT_MASS)
        }else{
                molWeight = NA
                exactMass = NA
        }
        out <- list(molWeight = molWeight, 
                    exactMass = exactMass)
        return(out)
}

sMtb_mets$SEEDid <- sapply(sMtb_mets$`KeGG ID`, KEGGids2SEEDids)
sMtb_mets$seedMass <- sapply(sMtb_mets$SEEDid, SEEDids2Mass)
sMtb_mets_exMassMW <- t(sapply(sMtb_mets$`KeGG ID`, KEGGid2MolWeightNExactMass))
sMtb_mets$exactMass <- sMtb_mets_exMassMW[, 2]
sMtb_mets$molWeight <- sMtb_mets_exMassMW[, 1]

iJO1366_mets$SEEDid <- sapply(iJO1366_mets$`KEGG ID`, KEGGids2SEEDids)
iJO1366_mets$seedMass <- sapply(iJO1366_mets$SEEDid, SEEDids2Mass)
iJO1366_mets_exMassMW <- t(sapply(iJO1366_mets$`KEGG ID`, KEGGid2MolWeightNExactMass))
iJO1366_mets$exactMass <- iJO1366_mets_exMassMW[, 2]
iJO1366_mets$molWeight <- iJO1366_mets_exMassMW[, 1]



























