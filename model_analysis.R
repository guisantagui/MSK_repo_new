setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs/modelAnalysis")

if(!require(SBMLR)) BiocManager::install("SBMLR")
library(SBMLR)
library(glue)

PA14 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/PA14.sbml")
F22031 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/F22031.sbml")
F30658 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/F34365.sbml")
F34365 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/F63912.sbml")
F63912 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/F63912.sbml")
F9670 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/F9670.sbml")
H27930 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/H27930.sbml")
H47921 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/H47921.sbml")
H5708 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/H5708.sbml")
M1608 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/M1608.sbml")
M37351 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/M37351.sbml")
M74707 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/M74707.sbml")
S86968 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/S86968.sbml")
T38079 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/T38079.sbml")
T52373 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/T52373.sbml")
T6313 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/T6313.sbml")
T63266 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/T63266.sbml")
W16407 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W16407.sbml")
W25637 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W25637.sbml")
W36662 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W36662.sbml")
W45909 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W45909.sbml")
W60856 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W60856.sbml")
W70322 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W70322.sbml")
W91453 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/W91453.sbml")
X78812 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/X78812.sbml")
X9820 <- readSBML("/Users/santamag/Desktop/GUILLEM/PATRIC_PA_models/X9820.sbml") 

models <- list(
        "PA14" = PA14,
        "F22031" = F22031,
        "F30658" = F30658,
        "F34365" = F34365,
        "F63912" = F63912,
        "F9670" = F9670,
        "H27930" = H27930,
        "H47921" = H47921,
        "H5708" = H5708,
        "M1608" = M1608,
        "M37351" = M37351,
        "M74707" = M74707,
        "S86968" = S86968,
        "T38079" = T38079,
        "T52373" = T52373,
        "T6313" = T6313,
        "T63266" = T63266,
        "W16407" = W16407,
        "W25637" = W25637,
        "W36662" = W36662,
        "W45909" = W45909,
        "W60856" = W60856,
        "W70322" = W70322,
        "W91453" = W91453,
        "X78812" = X78812,
        "X9820" = X9820
)

# This function removes all reactions and the compounds that participate in those reactions that have been 
# gapfilled, remaining only in the model the reactions that have genetic evidence.

filtMods <- function(mods){
        for(j in seq_along(mods)){
                vecPat <- c()
                maxLoop <- max(which(nchar(mods[[j]]$notes) > nchar("PROTEIN_ASSOCIATION:Unknown")))
                mods[[j]]$notes <- mods[[j]]$notes[1:maxLoop]
                for(i in 1:maxLoop){
                        if(grepl("GENE_ASSOCIATION", mods[[j]]$notes[i]) == TRUE & i %% 2 == 1){
                                vecPat <- c(vecPat, TRUE)
                        }
                        else if(grepl("PROTEIN_ASSOCIATION", mods[[j]]$notes[i]) == TRUE & i %% 2 == 0){
                                vecPat <- c(vecPat, TRUE)
                        }else{vecPat <- c(vecPat, FALSE)}
                }
                vecTrue <- which(vecPat)
                if(!all(vecPat) & length(vecTrue) < maxLoop){
                        h <- 1
                        breaks <- c()
                        while(h < length(vecTrue)){
                                M <- max(c(vecTrue[h]:max(vecTrue))[which(c(vecTrue[h]:max(vecTrue))[1:length(vecTrue)] %in% vecTrue)])
                                breaks <- c(breaks, M)
                                h <- which(vecTrue ==M) + 1
                        }
                        toSub <- c()
                        for(k in 1:length(breaks)){
                                toSub <- c(toSub, (breaks[k] + 1), (vecTrue[breaks[k] + 1] - 1))
                        }
                        toSub <- toSub[!is.na(toSub)]
                        for(l in seq_along(toSub)){
                                mods[[j]]$notes[toSub[l] - 1] <- paste(mods[[j]]$notes[toSub[l] - 1], 
                                                                       mods[[j]]$notes[toSub[l]], 
                                                                       sep = "")
                        }
                        mods[[j]]$notes <- mods[[j]]$notes[-toSub]
                }
                even <- function(x) x%%2 == 0
                mods[[j]]$reactions <- mods[[j]]$reactions[1:(length(mods[[j]]$notes)/2)]
                unKnowns <- which(grepl("GENE_ASSOCIATION:Unknown|PROTEIN_ASSOCIATION:Unknown", 
                                        mods[[j]]$notes))
                mods[[j]]$notes <- mods[[j]]$notes[-unKnowns]
                toRemove <- unKnowns[even(unKnowns)]/2
                mods[[j]]$reactions <- mods[[j]]$reactions[-toRemove]
                cpdVec <- c()
                for(m in seq_along(mods[[j]]$reactions)){
                        cpdVec <- c(cpdVec, mods[[j]]$reactions[[m]]$reactants)
                        cpdVec <- c(cpdVec, mods[[j]]$reactions[[m]]$products)
                }
                cpdVec <- unique(cpdVec)
                mods[[j]]$species <- mods[[j]]$species[cpdVec]
        }
        return(mods)
}

filtered <- filtMods(models)

# This function returns a model that corresponds to the reactions that are found across all models introduced in the 
# input (as a list of models), the accessory model ef each one of the strains, and a matrix of presence/absence of
# accessory reactions for each strain.

getCoreModel <- function(filteredModels){
        filt_mods <- filteredModels
        reacList <- list()
        for(i in seq_along(filt_mods)){
                reacList[[i]] <- names(filt_mods[[i]]$reactions)
        }
        coreReacts <- Reduce(intersect, reacList)
        coreModel <- filt_mods[[1]]
        coreModel$reactions <- coreModel$reactions[coreReacts]
        cpdVec <- c()
        for(j in seq_along(coreModel$reactions)){
                cpdVec <- c(cpdVec, coreModel$reactions[[j]]$reactants)
                cpdVec <- c(cpdVec, coreModel$reactions[[j]]$products)
        }
        cpdVec <- unique(cpdVec)
        coreModel$species <- coreModel$species[cpdVec]
        
        accessMods <- filt_mods
        accessReacts <- c()
        for(k in seq_along(accessMods)){
                accessMods[[k]]$reactions <- accessMods[[k]]$reactions[!names(accessMods[[k]]$reactions) %in% coreReacts]
                accessReacts <- c(accessReacts, names(accessMods[[k]]$reactions)) 
                cpdVecAcc <- c()
                for(l in seq_along(accessMods[[k]]$reactions)){
                        cpdVecAcc <- c(cpdVecAcc, accessMods[[k]]$reactions[[l]]$reactants)
                        cpdVecAcc <- c(cpdVecAcc, accessMods[[k]]$reactions[[l]]$products)
                }
                cpdVecAcc <- unique(cpdVecAcc)
                accessMods[[k]]$species <- accessMods[[k]]$species[cpdVecAcc]
        }
        accessReacts <- unique(accessReacts)
        accessMat <- matrix(nrow = length(accessMods), ncol = length(accessReacts))
        for(m in seq_along(accessMods)){
                for(n in seq_along(accessReacts)){
                        if(accessReacts[n] %in% names(accessMods[[m]]$reactions)){
                                accessMat[m, n] <- 1
                        }else{
                                accessMat[m, n] <- 0
                        }
                }
        }
        rownames(accessMat) <- names(accessMods)
        colnames(accessMat) <- accessReacts
        return(list("Core Model" = coreModel,
                    "Accessory Models" = accessMods,
                    "Accessoy Reactions Matrix" = accessMat))
}

coreAccessModels <- getCoreModel(filtered)

SBMLR::saveSBML(coreAccessModels$`Core Model`, filename = "coreModelPA.sbml")

# Here we loog for groups of genes that have same pattern of presence across groups of strains -> 
# Very similar to cluster the strains according to the accessory matrix
accSubMats <- list()
presPats <- unique(apply(coreAccessModels$`Accessoy Reactions Matrix`, 2, glue_collapse))
for(i in seq_along(presPats)){
        equalPres <- which(apply(coreAccessModels$`Accessoy Reactions Matrix`, 2, glue_collapse) == presPats[i])
        accSubMats[[i]] <- coreAccessModels$`Accessoy Reactions Matrix`[, names(equalPres)]
}

distAcc <- dist(coreAccessModels$`Accessoy Reactions Matrix`, method = "euclidean")
HCAAcc <- hclust(distAcc, method = "ward.D")
plot(HCAAcc)


glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["F63912", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["F34365", ])

# F34365 & F63912 are, in terms of accessory reactions, the same.

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["F9670", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["W70322", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("F9670", "W70322"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["F9670", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["W70322", ])]

# F9670 & W70322 differ in 4 reactions: rxn07849_c0 (F9670 has),
# and rxn01599_c0  rxn03610_c0  rxn00540_c0 (the three present in W70322).

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["H5708", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["W16407", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("H5708", "W16407"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["H5708", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["W16407", ])]

# H5708 & W16407 differ in 9 reactions: rxn00714_c0, rxn05467_c0 and rxn10473_c0 (present in H5708) and 
# rxn05236_c0, rxn00117_c0, rxn07849_c0, rxn01453_c0, rxn06038_c0 and rxn05171_c0 (present in W16407).

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["H47921", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["X78812", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("H47921", "X78812"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["H47921", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["X78812", ])]

# H47921 & X78812 differ in 26 reactions: rxn00292_c0, rxn02975_c0, rxn03208_c0, rxn02377_c0, rxn00297_c0, rxn04016_c0,
# rxn04132_c0, rxn05195_c0, rxn00062_c0, rxn04457_c0, rxn02916_c0, rxn03012_c0, rxn10571_c0, rxn05468_c0, rxn04133_c0,
# rxn05467_c0, rxn04675_c0 (present in H47921) and rxn02969_c0, rxn03618_c0, rxn03908_c0, rxn05363_c0, rxn00714_c0, 
# rxn05209_c0, rxn09345_c0, rxn00117_c0, rxn00946_c0 (present in X78812).

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["T38079", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["W45909", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("T38079", "W45909"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["T38079", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["W45909", ])]

# T38079 & W45909 differ in 3 reactions: rxn10181_c0 (present in T38079), and rxn03208_c0, rxn04016_c0 (present in W45909)

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["M37351", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["W25637", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("M37351", "W25637"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["M37351", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["W25637", ])]

# M37351 & W25637 differ in 3 reactions: rxn03091_c0 (present in M37351), and rxn01453_c0, rxn06038_c0 (present in W25637)

glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["H27930", ]) == 
        glue_collapse(coreAccessModels$`Accessoy Reactions Matrix`["W91453", ])

coreAccessModels$`Accessoy Reactions Matrix`[c("H27930", "W91453"), 
                                             which(!coreAccessModels$`Accessoy Reactions Matrix`["H27930", ] == 
                                                           coreAccessModels$`Accessoy Reactions Matrix`["W91453", ])]
# H27930 & W91453 differ in 5 reactions: rxn02290_c0, rxn03091_c0 (present in H27930), and rxn10571_c0,
# rxn01453_c0, rxn06038_c0 (present in W91453).

# Load P. putida model
iJN746 <- readSBML("/Users/santamag/Desktop/GUILLEM/models/iJN746.xml")

# Lets load the newer Putida model
if(!require(rjson)) install.packages("rjson")
library(rjson)

iJN1411 <- fromJSON(file = "/Users/santamag/Desktop/GUILLEM/models/putida_new/SI3_iJN1411final_flux.json")



iJN1411$genes[[1]]$id




















filt_mods <- filtered
reacList <- list()
for(i in seq_along(filt_mods)){
        reacList[[i]] <- names(filt_mods[[i]]$reactions)
}
coreReacts <- Reduce(intersect, reacList)
coreModel <- filt_mods[[1]]
coreModel$reactions <- coreModel$reactions[coreReacts]
cpdVec <- c()
for(j in seq_along(coreModel$reactions)){
        cpdVec <- c(cpdVec, coreModel$reactions[[j]]$reactants)
        cpdVec <- c(cpdVec, coreModel$reactions[[j]]$products)
}
cpdVec <- unique(cpdVec)
coreModel$species <- coreModel$species[cpdVec]
accessMods <- filt_mods
for(k in seq_along(accessMods)){
        accessMods[[k]]$reactions <- accessMods[[k]]$reactions[!names(accessMods$PA14$reactions) %in% coreReacts]
        cpdVecAcc <- c()
        for(l in seq_along(accessMods[[k]]$reactions)){
                cpdVecAcc <- c(cpdVecAcc, accessMods[[k]]$reactions[[l]]$reactants)
                cpdVecAcc <- c(cpdVecAcc, accessMods[[k]]$reactions[[l]]$products)
        }
        cpdVecAcc <- unique(cpdVecAcc)
        accessMods[[k]]$species <- accessMods[[k]]$species[cpdVecAcc]
}








names(filtered$PA14$reactions)
a <- c(1, 2, 3, 4, 8)
b <- c(3, 6, 8, 1)
c <- c(1, 3, 4, 8)
Reduce(intersect, list(a,b,c))

































j <- 7
mods <- models
vecPat <- c()
maxLoop <- max(which(nchar(mods[[j]]$notes) > nchar("PROTEIN_ASSOCIATION:Unknown")))
mods[[j]]$notes <- mods[[j]]$notes[1:maxLoop]
for(i in 1:maxLoop){
        if(grepl("GENE_ASSOCIATION", mods[[j]]$notes[i]) == TRUE & i %% 2 == 1){
                vecPat <- c(vecPat, TRUE)
        }
        else if(grepl("PROTEIN_ASSOCIATION", mods[[j]]$notes[i]) == TRUE & i %% 2 == 0){
                vecPat <- c(vecPat, TRUE)
        }else{vecPat <- c(vecPat, FALSE)}
}
vecTrue <- which(vecPat)
if(!all(vecPat) & length(vecTrue) < maxLoop){
        h <- 1
        breaks <- c()
        while(h < length(vecTrue)){
                M <- max(c(vecTrue[h]:max(vecTrue))[which(c(vecTrue[h]:max(vecTrue))[1:length(vecTrue)] %in% vecTrue)])
                breaks <- c(breaks, M)
                h <- (which(vecTrue == M) + 1)
        }
        toSub <- c()
        for(k in 1:length(breaks)){
                toSub <- c(toSub, (breaks[k] + 1), (vecTrue[breaks[k] + 1] - 1))
        }
        toSub <- toSub[!is.na(toSub)]
        for(l in seq_along(toSub)){
                mods[[j]]$notes[toSub[l] - 1] <- paste(mods[[j]]$notes[toSub[l] - 1], 
                                                       mods[[j]]$notes[toSub[l]], 
                                                       sep = "")
        }
        mods[[j]]$notes <- mods[[j]]$notes[-toSub]
}
even <- function(x) x%%2 == 0
mods[[j]]$reactions <- mods[[j]]$reactions[1:(length(mods[[j]]$notes)/2)]
unKnowns <- which(grepl("GENE_ASSOCIATION:Unknown|PROTEIN_ASSOCIATION:Unknown", 
                        mods[[j]]$notes))
mods[[j]]$notes <- mods[[j]]$notes[-unKnowns]
toRemove <- unKnowns[even(unKnowns)]/2
mods[[j]]$reactions <- mods[[j]]$reactions[-toRemove]
cpdVec <- c()
for(m in seq_along(mods[[j]]$reactions)){
        cpdVec <- c(cpdVec, mods[[j]]$reactions[[m]]$reactants)
        cpdVec <- c(cpdVec, mods[[j]]$reactions[[m]]$products)
}
cpdVec <- unique(cpdVec)
mods[[j]]$species <- mods[[j]]$species[cpdVec]













odd <- function(x) x%%2 != 0

even <- function(x) x%%2 == 0

evenb <- function(x) !odd(x)

even(unKnowns)


















breaks <- c()
h <- 1
while(h < max(vecTrue)){
        M <- max(c(vecTrue[h]:max(vecTrue))[which(c(vecTrue[h]:max(vecTrue))[1:length(vecTrue)] %in% vecTrue)])
        breaks <- c(breaks, M)
        h <- (M + 1)
}
vecTrue
models$F22031$notes[65]

max(which(c(1:length(which(vecPat))) %in% which(vecPat)))
which(vecPat)[66]







paste(X9820$notes[1273], X9820$notes[1274], sep = "")

X9820$reactions[1:(3058/2)]

which(grepl("Unknown", ))






length(F22031$notes[which(nchar(F22031$notes) > nchar("PROTEIN_ASSOCIATION:Unknown"))])




