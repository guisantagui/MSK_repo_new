############################################################################################################################
############################################################################################################################
### Create a list of metabolites that we want to see in MTBC strains' LC-MS analysis #######################################
############################################################################################################################
############################################################################################################################

if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(R.cache)) install.packages("R.cache")
library(R.cache)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(enviPat)) install.packages("enviPat")
library(enviPat)
if(!require(Rdisop)) BiocManager::install("Rdisop")
library(Rdisop)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcCPDs4MS")

# Import model dataframes and other data
data("isotopes")
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

# From the long string of aliases get KEGG ids and add a column to modelSEED dataframe with them
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

# in this table where there is no KEGG id there is the string "none". Substitute by NA
sMtb_mets$`KeGG ID`[sMtb_mets$`KeGG ID` == "none"] <- NA

# transforms kegg ids to seed ids. 
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

# Gets mass from modelSEED dataframe
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

# retrieves monoisotopic mass and molecular weight from kegg--> problem: usually masses and weights from KEGG are from ions
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
sMtb_mets_exMassMW <- apply(sMtb_mets_exMassMW, 2, as.numeric)
sMtb_mets$mIsotMass_KEGG <- sMtb_mets_exMassMW[, 2]
sMtb_mets$molWeight <- sMtb_mets_exMassMW[, 1]

iJO1366_mets$SEEDid <- sapply(iJO1366_mets$`KEGG ID`, KEGGids2SEEDids)
iJO1366_mets$seedMass <- sapply(iJO1366_mets$SEEDid, SEEDids2Mass)
iJO1366_mets_exMassMW <- t(sapply(iJO1366_mets$`KEGG ID`, KEGGid2MolWeightNExactMass))
iJO1366_mets_exMassMW <- apply(iJO1366_mets_exMassMW, 2, as.numeric)
iJO1366_mets$mIsotMass_KEGG <- iJO1366_mets_exMassMW[, 2]
iJO1366_mets$molWeight <- iJO1366_mets_exMassMW[, 1]


m2155_mets$seedMass <- sapply(m2155_mets$seedID, SEEDids2Mass)
m2155_mets_exactMass <- t(sapply(m2155_mets$keggID, KEGGid2MolWeightNExactMass))
m2155_mets_exactMass <- apply(m2155_mets_exactMass, 2, as.numeric)
m2155_mets$mIsotMass_KEGG <- m2155_mets_exactMass[, 2]
m2155_mets$molWeight <- m2155_mets_exactMass[, 1]


# few compounds have a couple of kegg ids. Keep just the first, as it's the good one
mtbH37Rv_mets$KEGGID <- sapply(mtbH37Rv_mets$KEGGID, function(x) strsplit(as.character(x), split = "; ")[[1]][1])

mtbH37Rv_mets_exactMass <- t(sapply(mtbH37Rv_mets$KEGGID, KEGGid2MolWeightNExactMass))
mtbH37Rv_mets_exactMass <- apply(mtbH37Rv_mets_exactMass, 2, as.numeric)
mtbH37Rv_mets$mIsotMass_KEGG <- mtbH37Rv_mets_exactMass[, 2]
mtbH37Rv_mets$molWeigh <- mtbH37Rv_mets_exactMass[, 1]

############################################################################################################################
##### TABLE CURATION ####################################################################################################### 
############################################################################################################################

############################################################################################################################
### sMtb metabolites                                                                                                       # 
############################################################################################################################

# Remove the compounds that don't have any identifier assigned:
sMtb_mets <- sMtb_mets[!is.na(sMtb_mets$`KeGG ID`) | sMtb_mets$`PubChem ID` != "none" | sMtb_mets$`ChEBI ID` != "none", ]

# Now curate manually

# (R)-3-Hydroxy-3-methyl-2-oxopentanoate don't correspond either to kegg id or formula in sMtb model

# Change the name 
sMtb_mets$Name[sMtb_mets$Name == "(R)-3-Hydroxy-3-methyl-2-oxopentanoate"] <- "(S)-2-Aceto-2-hydroxybutanoate"

# Alpha mycolate: not in kegg, but in chebi: include monoisotopic mass:
sMtb_mets$mIsotMass_KEGG[sMtb_mets$Name == "alpha-mycolate C80 (26)"] <- 1165.20545

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
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "2-Demethylmenaquinone (n=7)")] <- 702.5376

# decaprenyl phosphate: in kegg this compound doesn't appear--> is a lipid. Add monoisotopic mass from chebi/pubChem
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "decaprenyl phosphate")] <- 776.5872

# Dihydrolipoylprotein: no mass info in KEGG-->generic compound
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Dihydrolipoylprotein"), ]

# decaprenylphoshoryl-5-phosphoribose: not in kegg.---> Get mass from chebi
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "decaprenylphoshoryl-5-phosphoribose")] <- 990.6115

# Reduced & Oxidized ferredoxin: are proteins, so there is no mass info--->remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Reduced ferredoxin"), ]
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Oxidized ferredoxin"), ]

# Fe2+ and Fe3+: remove one, because the mass is the same so the MS cannot difference the 2 of them
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Fe2+"), ]

# 7,8-dihydromonapterin: not in KEGG. Get monoisotopic mass from Chebi
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "7,8-dihydromonapterin")] <- 255.0967

# D-glucan monomer: not found in the cell free--> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "D-glucan monomer"), ]

# hexacosanoyl-CoA: not in KEGG--> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "hexacosanoyl-CoA")] <- 1145.5014

# heptadecanoate: not in KEGG---> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "heptadecanoate")] <- 270.25588

# heptanoyl-CoA: not in KEGG---> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "heptanoyl-CoA")] <- 879.20403

# hexacosanoate: not in KEGG---> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "hexacosanoate")] <- 396.39673

# Hexadecanoyl-[acp]: generic compound---> remove
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Hexadecanoyl-[acp]"), ]

# S-(2-Methylpropanoyl)-dihydrolipoamide: Not in KEGG---> add mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "S-(2-Methylpropanoyl)-dihydrolipoamide")] <- 277.11702

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

# Thioredoxin: no info about mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Thioredoxin"), ]

# Thioredoxin disulfide: No info about mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Thioredoxin disulfide"), ]

# Peptidoglycan: variable mass---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Peptidoglycan"), ]

# (2E)-Octadecenoyl-[acp]: multiple sizes---> remove it 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "(2E)-Octadecenoyl-[acp]"), ]

# Phosphatidylethanolamine: variable sizes---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Phosphatidylethanolamine"), ]

# PIM2: not in KEGG---> add mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "PIM2")] <- 1176.67843

# S-Aminomethyldihydrolipoylprotein: has R groups---> remove it
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "S-Aminomethyldihydrolipoylprotein"), ]

# Ubiquinone (n=7) and Ubiquinol (n=7): variable mass---> remove them
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Ubiquinone (n=7)"), ]
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "Ubiquinol (n=7)"), ]

# alpha,alpha-Trehalose-2-sulfate: not in KEGG---> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "alpha,alpha-Trehalose-2-sulfate")] <- 422.07303

# SL659: not in KEGG---> get mass from ChEBI
sMtb_mets$mIsotMass_KEGG[which(sMtb_mets$Name == "SL659")] <- 659.29542

print(paste("Proportion of compounds with unique mass in sMtb:", 
            as.character(length(unique(sMtb_mets$mIsotMass_KEGG))/length(sMtb_mets$mIsotMass_KEGG))))

# Check what are the compounds with the same mass
sMtb_dupMets <- sMtb_mets[unlist(duplicated(sMtb_mets$mIsotMass_KEGG)), ]

# see if formulas are uncharged:
print(paste("Number of charged formulas in sMtb table:", 
            as.character(sum(sapply(sMtb_mets$Formula, function(x) getMolecule(x)$charge)))))

# Get monoisotopic masses from neutral formulas
sMtb_mets$mIsotMass_fromNeutFormula <- sapply(sMtb_mets$Formula, function(x) check_chemform(isotopes = isotopes, chemforms = x)$monoisotopic_mass)

sMtb_mets_massDisc <- sMtb_mets[!is.na(sMtb_mets$mIsotMass_KEGG), ][round(unlist(sMtb_mets$mIsotMass_KEGG[!is.na(sMtb_mets$mIsotMass_KEGG)]), digits = 4) != round(sMtb_mets$mIsotMass_fromNeutFormula, digits = 4)[!is.na(sMtb_mets$mIsotMass_KEGG)], ]

View(sMtb_mets_massDisc)

# Check each compound with discrepancies to see why this happens

# 1-(5-Phosphoribosyl)-5-amino-4-imidazolecarboxamide: difference in rounding---> keep it

# 7,8-dihydromonapterin: difference in rounding---> keep it

# decaprenyl phosphate: ChEBI shows charged formula---> keep it

# 7,8-Didemethyl-8-hydroxy-5-deazariboflavin: difference in rounding---> keep it

# D-Fructose 1,6-bisphosphate: difference in rounding---> keep it

# Mn: difference in rounding---> keep it

# Nitric oxide: difference in rounding---> keep it

# SL659: ChEBI considers the negative ion---> keep it

# alpha,alpha-Trehalose 6,6-bismycolate: discrepancies in formula---> remove 
sMtb_mets <- sMtb_mets[-which(sMtb_mets$Name == "alpha,alpha-Trehalose 6,6-bismycolate"), ]

# Save objects to avoid wasting the time functions need to run
save(sMtb_mets, file = "sMtb_mets.RData")
save(m2155_mets, file = "m2155_mets.RData")
save(mtbH37Rv_mets, file = "mtbH37Rv_mets.RData")
save(iJO1366_mets, file = "iJO1366_mets.RData")

############################################################################################################################
### iJO13666 metabolites                                                                                                   # 
############################################################################################################################

# 2-Demethylmenaquinol 8: has kegg ID, but as longitude is variable in kegg there is no mass. n=8, so from chebi, we get mass
iJO1366_mets$`KEGG ID`[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinol 8")] <- "C19847"
iJO1366_mets$mIsotMass_KEGG[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinol 8")] <- 704.55323

# 2-Demethylmenaquinone 8: get mass from ChEBI
iJO1366_mets$mIsotMass_KEGG[which(iJO1366_mets$`Metabolite Name` == "2-Demethylmenaquinone 8")] <- 702.53758

# [2Fe-1S] desulfurated iron-sulfur cluster, [2Fe-2S] iron-sulfur cluster, [3Fe-4S] damaged iron-sulfur cluster,
# [4Fe-4S] iron-sulfur cluster: remove these compounds because are part of proteins
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("[2Fe-1S] desulfurated iron-sulfur cluster", 
                                                                          "[2Fe-2S] iron-sulfur cluster", 
                                                                          "[3Fe-4S] damaged iron-sulfur cluster",
                                                                          "[4Fe-4S] iron-sulfur cluster")), ]

# lipoprotein, applipoprotein: generic--->remove it 
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` %in% c("lipoprotein", 
                                                                          "applipoprotein")), ]

# Fe2+: same mass as Fe3+--->remove it
iJO1366_mets <- iJO1366_mets[-which(iJO1366_mets$`Metabolite Name` == "Fe2+"), ]

# Remove compounds with "X" or "R" in the formula, as they have undefined mass
iJO1366_mets <- iJO1366_mets[-grep("X|R", iJO1366_mets$`Neutral Formula`), ]

# Get masses from formulas
iJO1366_mets$mIsotMass_fromNeutFormula <- sapply(iJO1366_mets$`Neutral Formula`, 
                                                 function(x) check_chemform(isotopes = isotopes, 
                                                                            chemforms = x)$monoisotopic_mass)

iJO1366_mets$mIsotMass_KEGG <- sapply(iJO1366_mets$mIsotMass_KEGG, 
                                      function(x) if(length(x) == 0){x <- NA}else{x <- x})


# See if there are discrepancies between the masses given by function and the ones given based in KEGG ids
iJO1366_mets_massDisc <- iJO1366_mets[!is.na(iJO1366_mets$mIsotMass_KEGG), ][round(iJO1366_mets$mIsotMass_KEGG[!is.na(iJO1366_mets$mIsotMass_KEGG)], digits = 4) != round(iJO1366_mets$mIsotMass_fromNeutFormula, digits = 4)[!is.na(iJO1366_mets$mIsotMass_KEGG)], ]
View(iJO1366_mets_massDisc)

# See why the discrepancies happen

# 1,4-dihydroxy-2-napthoyl-CoA: neutral formula here is wrong-charges are applied--> correct it
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "1,4-dihydroxy-2-napthoyl-CoA"] <- "C32H42N7O19P3S"

# Adenosyl cobinamide: KEGG considers positive ion with formula C58H84CoN16O11

# Adenosyl cobinamide phosphate: KEGG considers positive ion with formula C58H85CoN16O14P

# ADP-D-glycero-D-manno-heptose: formula has a typo-0 instead of O---> change it 
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "ADP-D-glycero-D-manno-heptose"] <- "C17H27N5O16P2"

# Adenosine-GDP-cobinamide: KEGG considers negative ion---> keep it

# 5-Amino-1-(5-Phospho-D-ribosyl)imidazole-4-carboxamide: rounding different---> keep it 

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

# N-Acetyl-D-galactosamine: formula is wrong
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "(-)-Ureidoglycolate"] <- "C8H15NO6"

# Fe(III)dicitrate: KEGG considers charged ion (negative)---> keep it

# fusidic acid: here it has one less oxigen atom--->correct formula, charged an uncharged
iJO1366_mets$`Neutral Formula`[iJO1366_mets$`Metabolite Name` == "fusidic acid"] <- "C31H48O6"
iJO1366_mets$`Charged Formula`[iJO1366_mets$`Metabolite Name` == "fusidic acid"] <- "C31H48O6"

iJO1366_mets$mIsotMass_fromNeutFormula <- sapply(iJO1366_mets$`Neutral Formula`, 
                                                 function(x) check_chemform(isotopes = isotopes, 
                                                                            chemforms = x)$monoisotopic_mass)

save(iJO1366_mets, file = "iJO1366_mets.RData")

load("iJO1366_mets.RData")
load("sMtb_mets.RData")
load("mtbH37Rv_mets.RData")
load("m2155_mets.RData")


iJO1366_mets_massDisc <- iJO1366_mets[!is.na(iJO1366_mets$mIsotMass_KEGG), ][round(iJO1366_mets$mIsotMass_KEGG[!is.na(iJO1366_mets$mIsotMass_KEGG)], digits = 4) != round(iJO1366_mets$mIsotMass_fromNeutFormula, digits = 4)[!is.na(iJO1366_mets$mIsotMass_KEGG)], ]
View(iJO1366_mets_massDisc)

# Remove unnecessary columns
iJO1366_mets <- iJO1366_mets[-which(colnames(iJO1366_mets) %in% c("seedMass", "Alternate Names", "Compartment"))]

length(unique(iJO1366_mets$`Metabolite Name`))

print(paste("Proportion of compounds with unique mass in iJO1366:", 
            as.character(length(unique(round(iJO1366_mets$mIsotMass_fromNeutFormula, digits = 4)))/length(iJO1366_mets$mIsotMass_fromNeutFormula))))

############################################################################################################################
### m2155 metabolites                                                                                                      # 
############################################################################################################################

# Turn names into character vector 
m2155_mets$name <- as.character(m2155_mets$name)
# Turn KEGG ids into caracter vector

m2155_mets$keggID <- as.character(m2155_mets$keggID)

# Remove metabolites with no name 
m2155_mets <- m2155_mets[-which(m2155_mets$name == ""), ]

# Remove this compound as it's nothing
m2155_mets <- m2155_mets[-which(m2155_mets$name == "C15815"), ]

m2155_mets$formula <- as.character(m2155_mets$formula)
m2155_mets$formula[m2155_mets$formula %in% c("", ".")] <- NA

m2155_mets <- m2155_mets[-grep("n|R|X", m2155_mets$formula), ]

# m2155_mets <- m2155_mets[!is.na(m2155_mets$formula), ]

# m2155_mets <- m2155_mets[m2155_mets$mIsotMass_KEGG != "numeric(0)", ]


m2155_mets$mIsotMass_fromNeutFormula <- sapply(m2155_mets$formula, 
                                               function(x) if(!is.na(x)) check_chemform(isotopes = isotopes, 
                                                                                        chemforms = x)$monoisotopic_mass else x <- NA)

# get formulas from KEGG
m2155_mets$formula_KEGG <- sapply(m2155_mets$keggID, function(x) if(!is.na(x)) keggGet(x)[[1]]$FORMULA else x <- NA)
m2155_mets$formula_KEGG <- unlist(lapply(m2155_mets$formula_KEGG, function(x) if(is.null(x)) x <- NA else x <- x))

# Remove formulas with ambiguous masses
m2155_mets <- m2155_mets[-grep("n|R|X", m2155_mets$formula_KEGG), ]

# remove compounds that don't have formula from KEGG and from seed
m2155_mets <- m2155_mets[!is.na(m2155_mets$formula) & !is.na(m2155_mets$formula_KEGG), ]

m2155_mets$mIsotMass_fromKEGGFormula <- sapply(m2155_mets$formula_KEGG, 
                                               function(x) check_chemform(isotopes = isotopes, 
                                                                          chemforms = x)$monoisotopic_mass)

m2155_mets_massDisc <- m2155_mets[!is.na(m2155_mets$mIsotMass_fromKEGGFormula), ][round(unlist(m2155_mets$mIsotMass_fromKEGGFormula[!is.na(m2155_mets$mIsotMass_fromKEGGFormula)]), 
                                                                                        digits = 4) != round(m2155_mets$mIsotMass_fromNeutFormula, digits = 4)[!is.na(m2155_mets$mIsotMass_fromKEGGFormula)], ]
# A lot of discrepancies---> seems that formulas here are all ions ?¿?¿?¿


dim(m2155_mets)
dim(m2155_mets_massDisc)

print(paste("Proportion of compounds with unique mass in m2155:",
            as.character(length(round(unique(m2155_mets$mIsotMass_fromKEGGFormula), digits = 4))/length(m2155_mets$mIsotMass_fromKEGGFormula))))


############################################################################################################################
### mtbH37Rv metabolites                                                                                                   # 
############################################################################################################################

View(mtbH37Rv_mets)

# remove compounds that have in formula stuff that are not chemical elements
mtbH37Rv_mets <- mtbH37Rv_mets[-grep("n|R|X", mtbH37Rv_mets$formula), ]

# Convert formula to factor and remove compounds without a formula 
mtbH37Rv_mets$formula <- as.character(mtbH37Rv_mets$formula)
mtbH37Rv_mets <- mtbH37Rv_mets[!mtbH37Rv_mets$formula %in% c("", "."), ]

# Get monoisotopic mass from seed formula
mtbH37Rv_mets$mIsotMass_fromFormula <- sapply(mtbH37Rv_mets$formula, 
                                             function(x) check_chemform(isotopes = isotopes, 
                                                                        chemforms = x)$monoisotopic_mass)

# get formulas from KEGG
mtbH37Rv_mets$formula_KEGG <- sapply(mtbH37Rv_mets$KEGGID, function(x) if(!is.na(x)) keggGet(x)[[1]]$FORMULA else x <- NA)
mtbH37Rv_mets <- mtbH37Rv_mets[-grep("n|R|X", mtbH37Rv_mets$formula_KEGG), ]

# get monoisotopic mass from kegg formulas
mtbH37Rv_mets$mIsotMass_fromKEGGFormula <- sapply(mtbH37Rv_mets$formula_KEGG, 
                                                 function(x) if(!is.na(x))check_chemform(isotopes = isotopes, 
                                                                                         chemforms = x)$monoisotopic_mass else x <- NA)

mtbH37Rv_mets_massDisc <- mtbH37Rv_mets[!is.na(mtbH37Rv_mets$mIsotMass_fromKEGGFormula), ][round(mtbH37Rv_mets$mIsotMass_fromKEGGFormula[!is.na(mtbH37Rv_mets$mIsotMass_fromKEGGFormula)], 
                                                                                                 digits = 4) != round(mtbH37Rv_mets$mIsotMass_fromFormula, digits = 4)[!is.na(mtbH37Rv_mets$mIsotMass_fromKEGGFormula)], ]

save(sMtb_mets, file = "sMtb_mets.RData")
save(m2155_mets, file = "m2155_mets.RData")
save(mtbH37Rv_mets, file = "mtbH37Rv_mets.RData")
save(iJO1366_mets, file = "iJO1366_mets.RData")

load("sMtb_mets.RData")
load("m2155_mets.RData")
load("mtbH37Rv_mets.RData")
load("iJO1366_mets.RData")


# Do mixed dataframe from all model lists

# remove charges from m2155 names 
m2155_mets$name <- gsub("\\(\\+\\)|\\(1\\+\\)|\\(2\\+\\)|\\(3\\+\\)|\\(1\\-\\)|\\(2\\-\\)|\\(3\\-\\)|\\(4\\-\\)|\\(5\\-\\)", "", m2155_mets$name)

length(m2155_mets$name)/length(unique(m2155_mets$name))

# These compounds have weird characters---> change them 
m2155_mets$name[m2155_mets$keggID == "C05399"] <- gsub(";", "", keggGet("C05399")[[1]]$NAME[1])
m2155_mets$name[m2155_mets$keggID == "C11508"] <- gsub(";", "", keggGet("C11508")[[1]]$NAME[1])
m2155_mets$name[m2155_mets$keggID == "C05670"] <- gsub(";", "", keggGet("C05670")[[1]]$NAME[1])

allKEGGIDs <- c(sMtb_mets$`KeGG ID`, 
                m2155_mets$keggID, 
                mtbH37Rv_mets$KEGGID, 
                iJO1366_mets$`KEGG ID`)

sum(is.na(allKEGGIDs))/length(allKEGGIDs)

length(unique(allKEGGIDs))


sMtb_mets <- sMtb_mets[, c("Name", "Formula", "KeGG ID", "PubChem ID", "ChEBI ID", "mIsotMass_fromNeutFormula", "molWeight")]
sMtb_mets$`PubChem ID`[sMtb_mets$`PubChem ID` == "none"] <- NA
sMtb_mets$`ChEBI ID`[sMtb_mets$`ChEBI ID` == "none"] <- NA

iJO1366_mets$`PubChem ID` <- sapply(iJO1366_mets$`KEGG ID`, function(x) if(!is.na(x))gsub("pubchem:", "", keggConv("pubchem", paste("cpd:", sep = "", x))[1]) else x <- NA)
iJO1366_mets$`ChEBI ID` <- sapply(iJO1366_mets$`KEGG ID`, function(x) if(!is.na(x))gsub("chebi:", "", keggConv("chebi", paste("cpd:", sep = "", x))[1]) else x <- NA)
iJO1366_mets <- iJO1366_mets[, c("Metabolite Name", "Neutral Formula", "KEGG ID", "PubChem ID", "ChEBI ID","mIsotMass_fromNeutFormula", "molWeight")]

m2155_mets$`PubChem ID` <- sapply(m2155_mets$keggID, function(x) if(!is.na(x))gsub("pubchem:", "", keggConv("pubchem", paste("cpd:", sep = "", x))[1]) else x <- NA)
m2155_mets$`ChEBI ID` <- sapply(m2155_mets$keggID, function(x) if(!is.na(x))gsub("chebi:", "", keggConv("chebi", paste("cpd:", sep = "", x))[1]) else x <- NA)
m2155_mets <- m2155_mets[, c("name", "formula_KEGG", "keggID", "PubChem ID", "ChEBI ID", "mIsotMass_fromKEGGFormula", "molWeight")]

mtbH37Rv_mets$`PubChem ID` <- sapply(mtbH37Rv_mets$KEGGID, function(x) if(!is.na(x))gsub("pubchem:", "", keggConv("pubchem", paste("cpd:", sep = "", x))[1]) else x <- NA)
mtbH37Rv_mets$`ChEBI ID` <- sapply(mtbH37Rv_mets$KEGGID, function(x) if(!is.na(x))gsub("chebi:", "", keggConv("chebi", paste("cpd:", sep = "", x))[1]) else x <- NA)
mtbH37Rv_mets <- mtbH37Rv_mets[, c("seedName", "formula_KEGG", "KEGGID", "PubChem ID", "ChEBI ID", "mIsotMass_fromKEGGFormula", "molWeigh")]

colnames(sMtb_mets) <- c("name", "formula", "KEGG_ID", "PubChem_ID", "ChEBI ID", "monoIsotopicMass", "molecularWeight")
colnames(iJO1366_mets) <- c("name", "formula", "KEGG_ID", "PubChem_ID", "ChEBI ID", "monoIsotopicMass", "molecularWeight")
colnames(m2155_mets) <- c("name", "formula", "KEGG_ID", "PubChem_ID", "ChEBI ID", "monoIsotopicMass", "molecularWeight")
colnames(mtbH37Rv_mets) <- c("name", "formula", "KEGG_ID", "PubChem_ID", "ChEBI ID", "monoIsotopicMass", "molecularWeight")

allMets <- rbind.data.frame(sMtb_mets,
                            iJO1366_mets,
                            m2155_mets,
                            mtbH37Rv_mets)

allMets_wKeggID <- allMets[match(unique(allMets$KEGG_ID[!is.na(allMets$KEGG_ID)]), allMets$KEGG_ID), ]

allMets_woKeggID <- allMets[is.na(allMets$KEGG_ID), ][apply(allMets_woKeggID, 1, function(x) sum(is.na(x)) != 6), ]

length(unique(allMets_woKeggID$name))/length(allMets_woKeggID$name)


allMets <- rbind.data.frame(allMets_wKeggID, allMets_woKeggID)

save(allMets, file = "allMets.RData")
write.csv(allMets, file = "allMets.csv")
