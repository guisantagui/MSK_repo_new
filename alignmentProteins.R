library(DECIPHER)
library(caret)


protDict <- read.csv("C:/Users/Guillem/Documents/PhD/comput/scriptOrhologues/Source_code_for_Pseudomonas_Metabolomics_Paper-master/find_protein_orthologues/protein_orthologue_dictionary_ref_PA14.csv")
colnames(protDict)[1] <- "UCBPP.PA14"

allProtFiles <- list.files("C:/Users/Guillem/Documents/PhD/comput/scriptOrhologues/Source_code_for_Pseudomonas_Metabolomics_Paper-master/find_protein_orthologues/prots",
                            pattern = "\\_protein.faa$")

allProtSeqs <- list()
for(i in seq_along(allProtFiles)){
        seq <- readAAStringSet(paste("C:/Users/Guillem/Documents/PhD/comput/scriptOrhologues/Source_code_for_Pseudomonas_Metabolomics_Paper-master/find_protein_orthologues/prots",
                                     allProtFiles[i],
                                     sep = "/"))
        allProtSeqs[[i]] <- seq
}
names(allProtSeqs) <- make.names(gsub(".faa", "", allProtFiles))


# This function asks for a fig code of pa14 and aligns all the proteins with this code
alignProt <- function(ref_fig){
        geneFigs <- protDict[grep(ref_fig, protDict$UCBPP.PA14, fixed = T), ]
        strNames <- colnames(geneFigs)
        setProtSeqs <- c()
        for(i in 1:ncol(geneFigs)){
                print(i)
                strProtSeqs <- allProtSeqs[[grep(strNames[i], names(allProtSeqs))]]
                protSeq <- strProtSeqs[[grep(geneFigs[, i], names(strProtSeqs), fixed = T)]]
                protSeq <- as.character(protSeq)
                setProtSeqs <- c(setProtSeqs, protSeq)
        }
        setProtSeqs <- AAStringSet(setProtSeqs)
        names(setProtSeqs) <- strNames
        alignedProts <- AlignSeqs(setProtSeqs)
        return(alignedProts)
}

# Align gcvP_2
fig287.6770.peg.5751_aligned <-  alignProt("fig|287.6770.peg.5751")

fig287.6770.peg.5751_aligned <- fig287.6770.peg.5751_aligned[-grep("PA7|PAO1|W70332", names(fig287.6770.peg.5751_aligned))]

# Build the matrix with the SNPs
fig287.6770.peg.5751_aligned_SNPsMat <- buildSNPsMat(fig287.6770.peg.5751_aligned, rhamnMat)


doRForest(fig287.6770.peg.5751_aligned_SNPsMat)


fig287.6770.peg.5751_SNPsBin <- binaryCodeSNPs(fig287.6770.peg.5751_aligned_SNPsMat)


library(dplyr)

library(randomForest)
library(caret)
library(e1071)


trControl <- trainControl(method = "cv",
                          number = 60,
                          search = "grid")


set.seed(1234)
# Run the model
rf_default_bin <- train(rhamn2cats ~ .,
                        data = fig287.6770.peg.5751_SNPsBin[, 1:(ncol(fig287.6770.peg.5751_SNPsBin)-1)],
                        method = "rf",
                        metric = "Accuracy",
                        trControl = trControl)

rf_pred_bin <- predict(rf_default_bin, fig287.6770.peg.5751_SNPsBin[, 1:(ncol(fig287.6770.peg.5751_SNPsBin)-1)])

c_bin <- confusionMatrix(rf_pred_bin, fig287.6770.peg.5751_SNPsBin$rhamn2cats)

rf_default <- train(rhamn2cats ~ .,
                    data = fig287.6770.peg.5751_aligned_SNPsMat[, 1:(ncol(fig287.6770.peg.5751_aligned_SNPsMat)-1)],
                    method = "rf",
                    metric = "Accuracy",
                    trControl = trControl)

rf_pred <- predict(rf_default, fig287.6770.peg.5751_aligned_SNPsMat[, 1:(ncol(fig287.6770.peg.5751_aligned_SNPsMat)-1)])


# compare predicted outcome and true outcome
c <- confusionMatrix(rf_pred, fig287.6770.peg.5751_aligned_SNPsMat$rhamn2cats)
# Print the results

print(rf_default)

gcvP_2 <- ParsedSubGraphs_allStrains_named_new$gcvP_2$AAAlignment

names(gcvP_2) <- gsub(".fasta|.fna", "", allSeqs)

gcvP_2 <- gcvP_2[-grep("PA14_jbx|W70332", names(gcvP_2))]

gcvP_2_SNPs <- buildSNPsMat(gcvP_2, rhamnMat)

doRForest(gcvP_2_SNPs)

BrowseSeqs(gcvP_2)
BrowseSeqs(fig287.6770.peg.5751_aligned)

rf_gcvP_2 <- train(rhamn2cats ~ .,
                   data = gcvP_2_SNPs[, 1:(ncol(gcvP_2_SNPs)-1)],
                   method = "rf",
                   metric = "Accuracy",
                   trControl = trControl)

rf_gcvP_2_pred <- predict(rf_gcvP_2, gcvP_2_SNPs[, 1:(ncol(gcvP_2_SNPs)-1)])

c_gcvP_2 <- confusionMatrix(rf_gcvP_2_pred, gcvP_2_SNPs$rhamn2cats)
gsub("-", replacement = "", as.character(gcvP_2$F22031)) == gsub("-|*", replacement = "", as.character(fig287.6770.peg.5751_aligned$F22031))


# My alignment and chen's one are the same, with the exception of W70332, which is longer in chen's one
# In gcvP_2 chen's random forest has an accuracy of 0.9667 and in mine 0.76

commonBothSets <- names(fig287.6770.peg.5751_aligned)[names(fig287.6770.peg.5751_aligned) %in% gsub(".fasta|.fna", "", allSeqs)]
alignSameStrain <- list()
for(i in seq_along(commonBothSets)){
       myAlig <- as.character(gcvP_2[[grep(commonBothSets[i], names(gcvP_2))]])
       myAlig <- gsub("-|*", "", myAlig)
       chenAlig <- as.character(fig287.6770.peg.5751_aligned[[grep(commonBothSets[i], names(fig287.6770.peg.5751_aligned))]])
       chenAlig <- gsub("-|*", "", chenAlig)
       print(myAlig)
       print(chenAlig)
       set <- c(myAlig,
                chenAlig)
       set <- AAStringSet(set)
       print(set)
       set <- AlignSeqs(set)
       alignSameStrain[[i]] <- set
}
names(alignSameStrain) <- commonBothSets

BrowseSeqs(alignSameStrain$UCBPP.PA14)
BrowseSeqs(alignSameStrain$F22031)
BrowseSeqs(alignSameStrain$F23197)
BrowseSeqs(alignSameStrain$F30658)
BrowseSeqs(alignSameStrain$F34365)
BrowseSeqs(alignSameStrain$F5677)
BrowseSeqs(alignSameStrain$F63912)
BrowseSeqs(alignSameStrain$F9670)
BrowseSeqs(alignSameStrain$H27930)
BrowseSeqs(alignSameStrain$H47921)
BrowseSeqs(alignSameStrain$H5708)
BrowseSeqs(alignSameStrain$M1608)
BrowseSeqs(alignSameStrain$M37351)
BrowseSeqs(alignSameStrain$M55212)
BrowseSeqs(alignSameStrain$M74707)
BrowseSeqs(alignSameStrain$S86968)
BrowseSeqs(alignSameStrain$T38079)
BrowseSeqs(alignSameStrain$T52373)
BrowseSeqs(alignSameStrain$T6313)
BrowseSeqs(alignSameStrain$T63266)
BrowseSeqs(alignSameStrain$W16407)
BrowseSeqs(alignSameStrain$W25637)
BrowseSeqs(alignSameStrain$W36662)
BrowseSeqs(alignSameStrain$W45909)
BrowseSeqs(alignSameStrain$W60856)
BrowseSeqs(alignSameStrain$W70332)
BrowseSeqs(alignSameStrain$W91453)
BrowseSeqs(alignSameStrain$X78812)
BrowseSeqs(alignSameStrain$X9820)

for(i in seq_along(alignSameStrain)){
        print(buildSNPsMat(alignSameStrain[[i]], rhamnMat))
}
buildSNPsMat(alignSameStrain$UCBPP.PA14, rhamnMat)


aja <- c(gsub("-", replacement = "", as.character(gcvP_2$F22031)), 
         gsub("-|*", replacement = "", as.character(fig287.6770.peg.5751_aligned$F22031)))
aja <- AAStringSet(aja)

BrowseSeqs(AlignSeqs(aja))

gcvP_2_seqInfo <- gffsAllStrains$W70332[grep("gcvP_2", gffsAllStrains$W70332$Name),]

W70332_gcvP_2_gene <- allGenomes$W70332[[1]][6318839:6315963]
W70332_gcvP_2_gene <- complement(W70332_gcvP_2_gene)
translate(W70332_gcvP_2_gene)

myAlig <- as.character(gcvP_2[[grep(commonBothSets[26], names(gcvP_2))]])
chenAlig <- as.character(fig287.6770.peg.5751_aligned[[grep(commonBothSets[26], names(fig287.6770.peg.5751_aligned))]])

set <- c(myAlig, chenAlig)

set <- AAStringSet(set)

