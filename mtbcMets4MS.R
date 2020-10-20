
########################################################################################
### Create a list of metabolites that we want to see in MTBC strains' LC-MS analysis ###
########################################################################################
if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(R.cache)) install.packages("R.cache")
library(R.cache)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(Rdisop)) BiocManager::install("Rdisop")
library(Rdisop)

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

# From the long string of aliases get KEGG ids and add a column to the dataframe with them
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
        #out$kegg <- KEGGID
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

# Now curate manually

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

print(paste("Proportion of compounds with unique mass:", 
            as.character(length(unique(sMtb_mets$exactMass))/length(sMtb_mets$exactMass))))

sMtb_dupMets <- sMtb_mets[unlist(duplicated(sMtb_mets$exactMass)), ]

View(iJO1366_mets)

# 2-Demethylmenaquinol 8: has kegg ID, but as longitude is variable in kegg there is no mass. n=8, so from chebi, we get mass
iJO1366_mets$`KEGG ID`[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinol 8")] <- "C19847"
iJO1366_mets$exactMass[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinol 8")] <- 704.55323

# 2-Demethylmenaquinone 8: get mass from ChEBI
iJO1366_mets$exactMass[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinone 8")] <- 702.53758

# [2Fe-1S] desulfurated iron-sulfur cluster, [2Fe-2S] iron-sulfur cluster, [3Fe-4S] damaged iron-sulfur cluster,
# [4Fe-4S] iron-sulfur cluster: remove these compounds because are part of proteins
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("[2Fe-1S] desulfurated iron-sulfur cluster", 
                                                                          "[2Fe-2S] iron-sulfur cluster", 
                                                                          "[3Fe-4S] damaged iron-sulfur cluster",
                                                                          "[4Fe-4S] iron-sulfur cluster")), ]

save(iJO1366_mets, file = "iJO1366_mets.RData")
load("iJO1366_mets.RData")
# lipoprotein, applipoprotein: generic--->remove it 
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("lipoprotein", 
                                                                          "applipoprotein")), ]

# Fe2+: same mass as Fe3+--->remove it
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` == "Fe2+"), ]

# Remove compounds with "X" or "R" in the formula, as they have undefined mass
iJO1366_mets <- iJO1366_mets[-grep("X|R", iJO1366_mets$`Neutral Formula`), ]

iJO1366_mets$exMassRdisop <- sapply(iJO1366_mets$`Neutral Formula`, function(x) getMolecule(x)$exactmass)

iJO1366_mets$exactMass <- sapply(iJO1366_mets$exactMass, 
                                 function(x) if(length(x) == 0){x <- NA}else{x <- x})


# See if there are discrepancies between the masses given by function and the ones given based in KEGG ids
View(iJO1366_mets[!is.na(iJO1366_mets$exactMass), ][round(iJO1366_mets$exactMass[!is.na(iJO1366_mets$exactMass)], digits = 4) != round(iJO1366_mets$exMassRdisop, digits = 4)[!is.na(iJO1366_mets$exactMass)], ])

# 1,4-dihydroxy-2-napthoyl-CoA: neutral formula here is wrong-charges are applied--> correct it
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "1,4-dihydroxy-2-napthoyl-CoA"] <- "C32H42N7O19P3S"

# ADP-D-glycero-D-manno-heptose: formula has a typo-0 instead of O---> change it 
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "ADP-D-glycero-D-manno-heptose"] <- "C17H27N5O16P2"
iJO1366_mets$exMassRdisop[iJO1366_mets$`Metabolite Name` == "ADP-D-glycero-D-manno-heptose"] <- getMolecule(iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "ADP-D-glycero-D-manno-heptose"])$exactmass

# Adenosine-GDP-cobinamide: keep the same, Rdisop has less decimals, don't know why

# S-Adenosyl-L-methionine: neutral formula is wrong, is the charged one
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "S-Adenosyl-L-methionine"] <- "C15H22N6O5S"

# S-Adenosylmethioninamine: kegg considers an intermediate ion (C14H23N6O3S)---> keep it

# Arsenite: KEGG considers the ion mass (AsO3)---> keep it

# Betaine aldehyde: KEGG considers the ion (C5H12NO)---> keep it

# Choline: KEGG considers the ion (C5H14NO)---> keep it

# D-Carnitine: KEGG considers charged + ion (C7H16NO3)---> keep it

# D-Fructose 1,6-bisphosphate: Rdisop is 0.0001 higher. Keep it

# Ferrichrome: Rdisop is 0.0001 higher. Keep it

# Fe-enterobactin: KEGG cibsiders a negative charged ion (C30H21N3O15)---> keep it 

# ferroxamine: Rdisop is 0.0001 higher. Keep it 

# sn-Glycero-3-phosphocholine: KEGG considers a positive charged ion (C8H21NO6P)---> keep it 

# gamma-butyrobetaine: KEGG considers a positive charged ion (C7H16NO2)---> keep it

# Glycogen: KEGG entry (C00182) is of fixed length but glycogen not---> remove it
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("glycogen")), ]

# Bicarbonate: KEGG considers negative charged ion (HCO3-)---> keep it

# Heme O: KEGG rounds the decimal to the lower---> keep it 

# Hexanoyl-CoA (n-C6:0CoA): neutral formula has an extra O that shouldn't be there 
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "Hexanoyl-CoA (n-C6:0CoA)"] <- "C27H46N7O17P3S"

# KDO(2)-lipid IV(A) with laurate, KDO(2)-lipid (A): discrepancies and probably we won't see them  in MS--> remove
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("KDO(2)-lipid IV(A) with laurate",
                                                                          "KDO(2)-lipid (A)")), ]

# S-Methyl-L-methionine: KEGG considers a + charged ion (C6H14NO2S)---> keep it 

# Mn2+: KEGG is rounding it---> keep it

# molybdenum cofactor: KEGG is rounding it---> keep it

# molybdopterin: discrepancy-in KEGG doesn't have a Cu atom---> remove it 
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("molybdopterin")), ]

# Nicotinamide adenine dinucleotide phosphate: KEGG considers the + ion (C21H29N7O17P3): keep it

# Ammonium: KEGG considers the + ion (HN4): keep it

# trans-Oct-2-enoyl-CoA: Neutral formula has an extra O--> change it
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "trans-Oct-2-enoyl-CoA"] <- "C29H48N7O17P3S"

# Protoheme: KEGG is rounding it slightly differnt---> keep it

# 7-aminomethyl-7-deazaguanine: formula is very different---> remove it
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("7-aminomethyl-7-deazaguanine")), ]

# Siroheme: KEGG considers positive ion (C42H44FeN4O16)---> keep it 

# Selenite: KEGG considers negative ion SeO3---> keep it

# D-Tagatose 1,6-biphosphate: different rounding---> keep it 

# Thiamin monophosphate: KEGG considers negative ion (C12H17N4O4PS)---> keep it

# Thiosulfate: KEGG considers negative ion HS2O3---> keep it 

# Undecaprenyl-diphospho-N-acetylmuramoyl: large and discrepancies---> remove both 
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`KEGG ID` %in% c("C05898",
                                                                  "C05897")), ]
# (-)-Ureidoglycolate: neutral formula is wrong (anion)---> correct it 
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "(-)-Ureidoglycolate"] <- "C3H6N2O4"

# Zinc: weird, in Rdisop there is 71.90470, but when running the function it's ok
iJO1366_mets$exMassRdisop[iJO1366_mets$`Metabolite Name` == "Zinc"] <- getMolecule(iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "Zinc"])$exactmass

# 