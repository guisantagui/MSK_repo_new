setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcPlcCovAnal")
gnumStrain <- data.frame(readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/data/mtbcGnumbers/mtbcGenomesInfo.xlsx"))
allDels <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/mtbcDelAnal/allDels.csv")
annotDir <- "C:/Users/Guillem/Documents/PhD/comput/data/MTBC_ancestralAnnot"
annot1 <- as.data.frame(rtracklayer::readGFF(paste(annotDir, "genes.gff", sep = "/")))
annot2 <- as.data.frame(rtracklayer::readGFF(paste(annotDir, "h37RvNCBIAnnot.gff3", sep = "/")))
gnumDir <- "C:/Users/Guillem/Documents/PhD/comput/data/mtbcAllGNums"


repStrains <- c("N0069",
                "N0072",
                "N0157",
                "N0031",
                "N0052",
                "N0145",
                "N0155",
                "N0004",
                "N0054",
                "N1274",
                "N0136",
                "N1216",
                "N1283",
                "N1176",
                "N1268",
                "N1272",
                "N0091",
                "N1201",
                "N1202",
                "N3913")

gnumStrain$G.NUMBER[match(repStrains, gnumStrain$N.NUMBER)]

repStrainsDF <- data.frame(N_NUMBER = repStrains,
                           G_NUMBER = gnumStrain$G.NUMBER[match(repStrains, gnumStrain$N.NUMBER)], 
                           LIN = gnumStrain$LINEAGE[match(repStrains, gnumStrain$N.NUMBER)])

A1Files <- list.files(gnumDir)[grep("A1", list.files(gnumDir))]
A2Files <- list.files(gnumDir)[grep("A2", list.files(gnumDir))]
A3Files <- list.files(gnumDir)[grep("A3", list.files(gnumDir))]
A4Files <- list.files(gnumDir)[grep("A4", list.files(gnumDir))]
L8Files <- list.files(gnumDir)[grep("L8", list.files(gnumDir))]
L9Files <- list.files(gnumDir)[grep("L9", list.files(gnumDir))]


A1Gnums <- as.character(read.table(paste(gnumDir, A1Files, sep = "/"))$V1)
A2Gnums <- c(as.character(read.table(paste(gnumDir, A2Files[1], sep = "/"))$V1),
             as.character(read.table(paste(gnumDir, A2Files[2], sep = "/"))$V1))
A3Gnums <- as.character(read.table(paste(gnumDir, A3Files, sep = "/"))$V1)
A4Gnums <- c(as.character(read.table(paste(gnumDir, A4Files[1], sep = "/"))$V1),
             as.character(read.table(paste(gnumDir, A4Files[2], sep = "/"))$V1))
L8Gnums <- as.character(read.table(paste(gnumDir, L8Files, sep = "/"))$V1)
L9Gnums <- as.character(read.table(paste(gnumDir, L9Files, sep = "/"))$V1)

set.seed(321)
A1Smp <- sample(A1Gnums, 4)
A2Smp <- sample(A2Gnums, 4)
A3Smp <- sample(A3Gnums, 4)
A4Smp <- sample(A4Gnums, 4)
L8Smp <- L8Gnums
L9Smp <- sample(L9Gnums, 4)

# Add the strains (Chimp) that chantal and carlos tested in the lab
A1Smp <- c(A1Smp, "G00857")

missLinsGNums <- c(L8Smp, L9Smp, A1Smp, A2Smp, A3Smp, A4Smp)

missLinsDF <- data.frame(N_NUMBER = gnumStrain$N.NUMBER[match(missLinsGNums, gnumStrain$G.NUMBER)],
                         G_NUMBER = missLinsGNums, 
                         LIN = as.character(gnumStrain$LINEAGE[match(missLinsGNums, gnumStrain$G.NUMBER)]))


missLins <- as.character(missLinsDF$LIN)
missLins[3:6] <- "L9"
missLins[20] <- "A4"
missLinsDF$LIN <- missLins

repStrainsDF <- rbind.data.frame(repStrainsDF, missLinsDF)

covDir <- "C:/Users/Guillem/Documents/PhD/comput/data/mtbcCoverageFiles"
covFiles <- list.files(covDir)

covList <- list()
for(i in 1:length(repStrainsDF$G_NUMBER)){
        covFileIdx <- grep(repStrainsDF$G_NUMBER[i], covFiles)
        if(length(covFileIdx) != 0){
                cov <- read.table(paste(covDir, covFiles[covFileIdx], sep = "/"))
                cov <- cov$V1
        }else{
                cov <- NA
        }
        covList[[i]] <- cov
}
names(covList) <- repStrainsDF$G_NUMBER

# Remove G00573 because it doesn't have coverage file
covList <- covList[unlist(lapply(covList, function(x) unique(!is.na(x))))]

covDF <- as.data.frame(do.call(cbind, covList))
covDF <- cbind.data.frame(as.numeric(rownames(covDF)), covDF)
colnames(covDF)[1] <- "position"

colnames(covDF)[2:ncol(covDF)] <- paste(colnames(covDF)[2:ncol(covDF)], 
                                        repStrainsDF$LIN[match(colnames(covDF)[2:ncol(covDF)], repStrainsDF$G_NUMBER)], 
                                        sep = "_")

getStartEnd <- function(genes){
        genes1 <- genes[genes %in% annot1$locus_tag]
        genes2 <- genes[!genes %in% annot1$locus_tag]
        gene1 <- paste(genes1, collapse = "|")
        startEnd1 <- annot1[grep(gene1, annot1$locus_tag), c("start", "end")]
        rownames(startEnd1) <- genes1
        startEnd <- startEnd1
        if(length(genes2) != 0){
                gene2 <- paste(genes2, collapse = "|")
                startEnd2 <- annot2[grep(gene2, annot2$ID), c("start", "end")]
                rownames(startEnd2) <- genes2
                startEnd <- rbind.data.frame(startEnd1, startEnd2)
        }
        startEnd <- startEnd[order(rownames(startEnd)), ]
        return(startEnd)
}

plcAndFlankCoord <- getStartEnd(c("Rv1754", "Rv1755", "Rv1756", "Rv2345", "Rv2346c", "Rv2347c","Rv2348c", "Rv2349c", "Rv2350c", "Rv2351c", "Rv2352c", "Rv2353c"))

plcABCCordInt <- min(plcAndFlankCoord[4:nrow(plcAndFlankCoord), ]):max(plcAndFlankCoord[4:nrow(plcAndFlankCoord), ])
plcDCordInt <- min(plcAndFlankCoord[1:3, ]):max(plcAndFlankCoord[1:3, ])

plcCovDF <- covDF[c(plcDCordInt, plcABCCordInt), ]


geneNameCov <- rep(NA, nrow(plccovDF))

for(i in 1:nrow(plcAndFlankCoord)){
        trueVec <- plcCovDF$position <= max(plcAndFlankCoord[i, ]) & plcCovDF$position >= min(plcAndFlankCoord[i, ])
        geneNameCov[trueVec] <- rownames(plcAndFlankCoord)[i]
}

plcCovDF$gene <- geneNameCov

rownames(plcCovDF) <- paste(plcCovDF$gene, plcCovDF$position, sep = "_")

heatmap(as.matrix(plcCovDF[, 2:(ncol(plcCovDF)-1)]), Rowv = NA, Colv = NA, scale = "none")


plcCovDFBin <- apply(plcCovDF[, 2:(ncol(plcCovDF) -1)], 2, function(x){
        x[x < 5] <- 0
        x[x >= 5] <- 1 
        return(x)})

pdf("plcGenesDels.pdf")
heatmap(as.matrix(plcCovDFBin), Rowv = NA, Colv = NA, scale = "none")
dev.off()
