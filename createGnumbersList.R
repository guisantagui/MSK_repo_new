# This script creates lists of all gnumbers in each lineage from mtbc database retrieved from labkey 
# (with the ones already ran removed)

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/createGnumbersList")

dataBase <- as.data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/data/mtbcDatabase/Genomes_2020-11-10_10-55-41.xlsx"))

dataBase <- dataBase[!is.na(dataBase$`PATH MAPPING OUTPUT`), ]

# Some "L2" are encoded as "2"
dataBase$LINEAGE[dataBase$LINEAGE == "2"] <- "L2"
# L9 are in "SUBLINEAGE COLL" column
dataBase$LINEAGE[dataBase$`SUBLINEAGE COLL` == "L9"] <- "L9"

L1 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L1" & !is.na(dataBase$LINEAGE)])
L2 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L2" & !is.na(dataBase$LINEAGE)])
L3 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L3" & !is.na(dataBase$LINEAGE)])
L4 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L4" & !is.na(dataBase$LINEAGE)])
L5 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L5" & !is.na(dataBase$LINEAGE)])
L6 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L6" & !is.na(dataBase$LINEAGE)])
L7 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L7" & !is.na(dataBase$LINEAGE)])
L8 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L8" & !is.na(dataBase$LINEAGE)])
L9 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "L9" & !is.na(dataBase$LINEAGE)])
A1 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "A1" & !is.na(dataBase$LINEAGE)])
A2 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "A2" & !is.na(dataBase$LINEAGE)])
A3 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "A3" & !is.na(dataBase$LINEAGE)])
A4 <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "A4" & !is.na(dataBase$LINEAGE)])
A <- data.frame(gNumber = dataBase$`G NUMBER`[dataBase$LINEAGE == "A" & !is.na(dataBase$LINEAGE)])

LNA <- data.frame(gNumber = dataBase$`G NUMBER`[is.na(dataBase$LINEAGE)])

newGnums <- list(A1 = A1, 
                 A2 = A2, 
                 A3 = A3, 
                 A4 = A4, 
                 L1 = L1, 
                 L2 = L2, 
                 L3 = L3, 
                 L4 = L4, 
                 L5 = L5,
                 L6 = L6, 
                 L7 = L7, 
                 L8 = L8, 
                 L9 = L9, 
                 LNA = LNA)

oldGnumFiles <- list.files("C:/Users/Guillem/Documents/PhD/comput/data/mtbcGnumbers")

doneGnums <- list()
for(i in 1:length(oldGnumFiles)){
        tab <- read.table(paste("C:/Users/Guillem/Documents/PhD/comput/data/mtbcGnumbers", 
                                oldGnumFiles[i], 
                                sep = "/"),
                          sep = "\t")
        gNums <- as.character(tab$V1)
        doneGnums[[i]] <- gNums
}
names(doneGnums) <- gsub("nr.txt", "", oldGnumFiles)

all(doneGnums$A1 %in% newGnums$A1$gNumber)
all(doneGnums$A2 %in% newGnums$A2$gNumber)
all(doneGnums$A3 %in% newGnums$A3$gNumber)
all(doneGnums$A4 %in% newGnums$A4$gNumber)
doneGnums$A4[which(!doneGnums$A4 %in% newGnums$A4$gNumber)]
doneGnums$A4[which(!doneGnums$A4 %in% newGnums$A4$gNumber)] %in% dataBase$`G NUMBER`

all(doneGnums$L1 %in% newGnums$L1$gNumber)
all(doneGnums$L2 %in% newGnums$L2$gNumber)
all(doneGnums$L3 %in% newGnums$L3$gNumber)
all(doneGnums$L4 %in% newGnums$L4$gNumber)
all(doneGnums$L5 %in% newGnums$L5$gNumber)
all(doneGnums$L6 %in% newGnums$L6$gNumber)
all(doneGnums$L7 %in% newGnums$L7$gNumber)
all(doneGnums$L8 %in% newGnums$L8$gNumber)
all(doneGnums$L9 %in% newGnums$L9$gNumber)
doneGnums$L9[which(!doneGnums$L9 %in% newGnums$L9$gNumber)]
doneGnums$L9[which(!doneGnums$L9 %in% newGnums$L9$gNumber)] %in% dataBase$`G NUMBER`

lengthsDone <- c()
lengthsNew <- c()
for(i in 1:13){
        donG <- length(doneGnums[[i]])
        newG <- length(newGnums[[i]]$gNumber)
        lengthsDone <- c(lengthsDone, donG)
        lengthsNew <- c(lengthsNew, newG)
}

doneNewLengths <- data.frame(done = lengthsDone,
                             new = lengthsNew)
rownames(doneNewLengths) <- names(doneGnums)
doneNewLengths

newGnums$A4$gNumber[!newGnums$A4$gNumber %in% doneGnums$A4]


gNums2Do <- list()
for(i in 1:13){
        g <- newGnums[[i]]$gNumber[!newGnums[[i]]$gNumber %in% doneGnums[[i]]]
        gNums2Do[[i]] <- g
}
names(gNums2Do) <- names(newGnums)[1:13]
gNums2Do$LNA <- newGnums$LNA

gNums2Do <- lapply(gNums2Do, as.data.frame)

gNums2Do <- gNums2Do[sapply(gNums2Do, function(x) nrow(x) > 0)]

gNums2Do$A <- A


for(i in 1:length(gNums2Do)){
        write.table(gNums2Do[[i]], 
                    file = paste(names(gNums2Do)[i],
                                 "-2",
                                 ".txt",
                                 sep = ""), 
                    sep = "\t", 
                    quote = F, 
                    row.names = F, 
                    col.names = F)
}
