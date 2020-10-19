
########################################################################################
### Create a list of metabolites that we want to see in MTBC strains' LC-MS analysis ###
########################################################################################
if(!require(data.table)) install.packages("data.table")
library(data.table)
#if(!require(RMassBank)) BiocManager::install("RMassBank")
#library(RMassBank)
if(!require(R.cache)) install.packages("R.cache")
library(R.cache)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcCPDs4MS")

# Import model dataframes

sMtb_mets <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/models/MTBC/sMtb.xlsx", sheet = 2))

m2155_mets <- read.csv("m2155_mets.csv")

mtbH37Rv_mets <- read.csv("mtbH37Rv_mets.csv")

iJO1366_mets <- as.data.frame(readxl::read_xls("C:/Users/Guillem/Documents/PhD/comput/models/Ecoli/iJO1366_info.xls", sheet = 4))

# Remove duplicated compounds 

mtbH37Rv_mets <- mtbH37Rv_mets[match(unique(mtbH37Rv_mets$seedID), mtbH37Rv_mets$seedID), 2:ncol(mtbH37Rv_mets)]

m2155_mets <- m2155_mets[match(unique(m2155_mets$name), m2155_mets$name), 3:ncol(m2155_mets)]

sMtb_mets <- sMtb_mets[match(unique(sMtb_mets$Name), sMtb_mets$Name), 2:ncol(sMtb_mets)]

iJO1366_mets <- iJO1366_mets[match(unique(iJO1366_mets$`Metabolite Name`), iJO1366_mets$`Metabolite Name`), 2:ncol(iJO1366_mets)]
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
                get <- tryCatch(keggGet(KEGGID), error = function(e) return("none"))
                if(get == "none"){
                        molWeight <- NA
                        exactMass <- NA
                }else{
                        molWeight <- as.numeric(get[[1]]$MOL_WEIGHT)
                        exactMass <- as.numeric(get[[1]]$EXACT_MASS)
                }
        }else{
                molWeight = NA
                exactMass = NA
        }
        out <- list(molWeight = molWeight, 
                    exactMass = exactMass)
        out$kegg <- KEGGID
        print(out)
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


m2155_mets$seedMass <- sapply(m2155_mets$seedID, SEEDids2Mass)
m2155_mets_exactMass <- t(sapply(m2155_mets$keggID, KEGGid2MolWeightNExactMass))
m2155_mets$exactMass <- m2155_mets_exactMass[, 2]
m2155_mets$molWeight <- m2155_mets_exactMass[, 1]

# Curate tables: start by removing the compounds that don't have any identifier assigned:

sMtb_mets <- sMtb_mets[!is.na(sMtb_mets$`KeGG ID`) | sMtb_mets$`PubChem ID` != "none" | sMtb_mets$`ChEBI ID` != "none", ]

# Now curate manuallyç

# (R)-3-Hydroxy-3-methyl-2-oxopentanoate don't correspond either to kegg id or formula in sMtb model

# Change the name 
sMtb_mets$Name[sMtb_mets$Name == "(R)-3-Hydroxy-3-methyl-2-oxopentanoate"] <- "(S)-2-Aceto-2-hydroxybutanoate"

# Alpha mycolate: not in kegg, but in chebi: include monoisotopic mass:
sMtb_mets$exactMass[sMtb_mets$Name == "alpha-mycolate C80 (26)"] <- 1165.20545

# apo-[acetyl-CoA:carbon-dioxide ligase (ADP-forming)]: there is no mass information--> remove 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "apo-[acetyl-CoA:carbon-dioxide ligase (ADP-forming)]"), ]

# Acyl-carrier protein: generic compound, no mass info--> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Acyl-carrier protein"), ]

# Amylose monomer: Can't be free in the cell--> remove it 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Amylose monomer"), ]

# [Acetyl-CoA:carbon-dioxide ligase (ADP-forming)]: Generic compound--> Remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "[Acetyl-CoA:carbon-dioxide ligase (ADP-forming)]"), ]

# 2-Demethylmenaquinone (n=7): in kegg there is no mass because there is no n in the entry, but here is 7 so 
# we can calculate
sMtb_mets$exactMass[which(sMtb_mets$Name == "2-Demethylmenaquinone (n=7)")] <- 702.5376

# decaprenyl phosphate: in kegg this compound doesn't appear--> is a lipid. Add monoisotopic mass from chebi/pubChem
sMtb_mets$exactMass[which(sMtb_mets$Name == "decaprenyl phosphate")] <- 776.5872

# Dihydrolipoylprotein: no mass info in KEGG-->generic compound
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Dihydrolipoylprotein"), ]

# decaprenylphoshoryl-5-phosphoribose: not in kegg.---> Get mass from chebi
sMtb_mets$exactMass[which(sMtb_mets$Name == "decaprenylphoshoryl-5-phosphoribose")] <- 990.6115

# Reduced & Oxidized ferredoxin: are proteins, so there is no mass info--->remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Reduced ferredoxin"), ]
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Oxidized ferredoxin"), ]

# Fe2+ and Fe3+: remove one, because the mass is the same so the MS cannot difference the 2 of them
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Fe2+"), ]

# 7,8-dihydromonapterin: not in KEGG. Get monoisotopic mass from Chebi
sMtb_mets$exactMass[which(sMtb_mets$Name == "7,8-dihydromonapterin")] <- 255.0967

# D-glucan monomer: not found in the cell free--> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "D-glucan monomer"), ]

# hexacosanoyl-CoA: not in KEGG--> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "hexacosanoyl-CoA")] <- 1145.5014

# heptadecanoate: not in KEGG---> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "heptadecanoate")] <- 270.25588

# heptanoyl-CoA: not in KEGG---> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "heptanoyl-CoA")] <- 879.20403

# hexacosanoate: not in KEGG---> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "hexacosanoate")] <- 396.39673

# Hexadecanoyl-[acp]: generic compound---> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Hexadecanoyl-[acp]"), ]

# S-(2-Methylpropanoyl)-dihydrolipoamide: Not in KEGG---> add mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "S-(2-Methylpropanoyl)-dihydrolipoamide")] <- 277.11702

# Lipoylprotein: variable mass, so there's no info about it---> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Lipoylprotein"), ]

# Menaquinone: variable mass, so there's no info about it---> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Menaquinone"), ]

# Octadecanoyl-[acyl-carrier protein]: no info about mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Octadecanoyl-[acyl-carrier protein]"), ]

# Malonyl-[acyl-carrier protein]: variable mass---> remove it 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Malonyl-[acyl-carrier protein]"), ]

# Polyphosphate: multiple sizes---> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Polyphosphate"), ]

# Thioredoxin disulfide: No info about mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Thioredoxin disulfide"), ]

# Peptidoglycan: variable mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Peptidoglycan"), ]

# (2E)-Octadecenoyl-[acp]: multiple sizes---> remove it 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "(2E)-Octadecenoyl-[acp]"), ]

# Phosphatidylethanolamine: variable sizes---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Phosphatidylethanolamine"), ]

# PIM2: not in KEGG---> add mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "PIM2")] <- 1176.67843

# S-Aminomethyldihydrolipoylprotein: has R groups---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "S-Aminomethyldihydrolipoylprotein"), ]

# Ubiquinone (n=7) and Ubiquinol (n=7): variable mass---> remove them
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Ubiquinone (n=7)"), ]
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Ubiquinol (n=7)"), ]

# alpha,alpha-Trehalose-2-sulfate: not in KEGG---> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "alpha,alpha-Trehalose-2-sulfate")] <- 422.07303

# SL659: not in KEGG---> get mass from ChEBI
sMtb_mets$exactMass[which(sMtb_mets$Name == "SL659")] <- 659.29542

length(unique(sMtb_mets$exactMass))

sum(is.na(sMtb_mets$exactMass))
