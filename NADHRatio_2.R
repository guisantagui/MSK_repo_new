setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/NADHRatio")

if(!require(DescTools)) install.packages("DescTools")
library(DescTools)
library(ggpubr)

od600Range <- readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/070720GlycerolGrowthCurve4NADH_Maurice.xlsx")
luminiscenceRange <- readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/070820LuminiscenseRangePA14_2h_Hildi.xlsx")

od600Range <- od600Range[40:nrow(od600Range), ]
colnames(od600Range) <- as.character(od600Range[1, ])
od600Range <- as.data.frame(apply(od600Range[2:(which(is.na(od600Range$A1))[1] -1), ], 
                                  2, 
                                  as.numeric))

luminiscenceRange <- luminiscenceRange[36:(nrow(luminiscenceRange)-5), ]
colnames(luminiscenceRange) <- as.character(luminiscenceRange[1, ])
luminiscenceRange <- as.data.frame(apply(luminiscenceRange[-1, ], 
                                         2, 
                                         as.numeric))

plot(luminiscenceRange$`Time [s]`, luminiscenceRange$L9)

# Maximum is about 12500 seconds

standardConc <- c(1000000, 100000, 10000, 1000, 500, 250, 100, 50, 25, 10, 5, 0)

luminiscenceRange$`Time [s]`/60



as.vector(luminiscenceRange[which(luminiscenceRange$`Time [s]`/60 == Closest(luminiscenceRange$`Time [s]`/60, 
                                                                   90)), grep("9", colnames(luminiscenceRange))])[]




timeInts <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 390, 420)

for(i in 1:length(timeInts)){
        plot(standardConc[5:12], 
             as.vector(luminiscenceRange[which(luminiscenceRange$`Time [s]`/60 == Closest(luminiscenceRange$`Time [s]`/60, 
                                                                                          timeInts[i])), grep("9", colnames(luminiscenceRange))])[5:12])
}

for(i in 1:length(timeInts)){
        plot(standardConc[5:12], 
             as.vector(luminiscenceRange[which(luminiscenceRange$`Time [s]`/60 == Closest(luminiscenceRange$`Time [s]`/60, 
                                                                                          timeInts[i])), grep("10", colnames(luminiscenceRange))])[5:12])
}

for(i in 1:length(timeInts)){
        plot(standardConc[5:12], 
             as.vector(luminiscenceRange[which(luminiscenceRange$`Time [s]`/60 == Closest(luminiscenceRange$`Time [s]`/60, 
                                                                                          timeInts[i])), grep("11", colnames(luminiscenceRange))])[5:12])
}

NADH071420 <-readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/071420LuminiscenseNADH_Hildi.xlsx")

NADH071420 <- NADH071420[36:(nrow(NADH071420)-5), ]
colnames(NADH071420) <- as.character(NADH071420[1, ])
NADH071420 <- as.data.frame(apply(NADH071420[-1, ], 
                                  2, 
                                  as.numeric))

NADH071620 <-readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/071620LuminiscenseNADH_Hildi.xlsx")

NADH071620 <- NADH071620[36:(nrow(NADH071420)-5), ]
colnames(NADH071620) <- as.character(NADH071620[1, ])
NADH071620 <- as.data.frame(apply(NADH071620[-1, ], 
                                  2, 
                                  as.numeric))

NADH072120 <-readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/072120LuminiscenseNADH_Hildi.xlsx")

NADH072120 <- NADH072120[36:(nrow(NADH072120)-5), ]
colnames(NADH072120) <- as.character(NADH072120[1, ])
NADH072120 <- as.data.frame(apply(NADH072120[-1, ], 
                                  2, 
                                  as.numeric))

NADH072320 <-readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/072320LuminiscenseNADH_Hildi.xlsx")

NADH072320 <- NADH072320[36:(nrow(NADH072320)-5), ]
colnames(NADH072320) <- as.character(NADH072320[1, ])
NADH072320 <- as.data.frame(apply(NADH072320[-1, ], 
                                  2, 
                                  as.numeric))

NADH072720 <-readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/072720LuminiscenseNADH_Hildi.xlsx")

NADH072720 <- NADH072720[36:(nrow(NADH072720)-5), ]
colnames(NADH072720) <- as.character(NADH072720[1, ])
NADH072720 <- as.data.frame(apply(NADH072720[-1, ], 
                                  2, 
                                  as.numeric))


getNADHDATA <- function(experimentMatrix, minutes, strains, standConc, substractZero = T, standConNAD = NULL){
        standConc <- sort(standConc, decreasing = T)
        if(!is.null(standConNAD)){
                standConNAD <- sort(standConNAD, decreasing = T)
        }
        letters <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P")
        position <- which(experimentMatrix$`Time [s]`/60 == Closest(experimentMatrix$`Time [s]`/60, minutes))
        strainsVec <- experimentMatrix[position, 4:tail(grep("18", colnames(experimentMatrix)), 1)]
        strainMat <- c()
        count <- 0
        countLet <- 1
        while(count < length(strains)){
                for(i in 1:3){
                        if(count < length(strains)){
                                nums <- (i*6-5):(6*i)
                                colPos <- paste(letters[countLet], nums, sep = "")
                                strainMat <- cbind(strainMat, as.numeric(strainsVec[, colPos]))
                                count <- count + 1
                        }else{
                                break()
                        }
                }
                countLet <- countLet + 1
        }
        colnames(strainMat) <- strains
        NAD <- strainMat[1:6 %% 2 == 1, ]
        rownames(NAD) <- c("R1", "R2", "R3")
        NADH <- strainMat[1:6 %% 2 == 0, ]
        rownames(NADH) <- c("R1", "R2", "R3")
        
        standCols <- sapply(c(20, 22, 24), function(x) paste(letters[1:length(standConc)], x, sep = ""))
        standMat <- matrix(nrow = 3, ncol = length(standConc))
        for(i in 1:ncol(standCols)){
                for(j in 1:nrow(standCols)){
                        standMat[i, j] <- experimentMatrix[position, colnames(experimentMatrix) == standCols[j, i]]
                }
        }
        rownames(standMat) <- c("R1", "R2", "R3")
        colnames(standMat) <- as.character(standConc)
        
        if(!is.null(standConNAD)){
                standColsNAD <- sapply(c(19, 21, 23), function(x) paste(letters[1:length(standConNAD)], x, sep = ""))
                standMatNAD <- matrix(nrow = 3, ncol = length(standConNAD))
                for(i in 1:ncol(standColsNAD)){
                        for(j in 1:nrow(standColsNAD)){
                                standMatNAD[i, j] <- experimentMatrix[position, colnames(experimentMatrix) == standColsNAD[j, i]]
                        }
                }
                rownames(standMatNAD) <- c("R1", "R2", "R3")
                colnames(standMatNAD) <- as.character(standConNAD)
        }
        
        if(substractZero == T && is.null(standConNAD) == T){
                zero = mean(standMat[, ncol(standMat)])
                NAD = NAD - zero
                NADH = NADH - zero
                standMat = standMat - zero
        }else if(substractZero == T && is.null(standConNAD) == F){
                zero = mean(standMat[, ncol(standMat)])
                zeroNAD = mean(standMatNAD, ncol(standMatNAD))
                NAD = NAD - zeroNAD
                NADH = NADH - zero
                standMat = standMat - zero
                standMatNAD = standMatNAD - zeroNAD
        }
        
        NADH_NAD_ratio <- NADH/NAD
        if(is.null(standConNAD) == T){
                NAD_NADH <- list(NAD = NAD,
                                 NADH = NADH,
                                 NADH_NAD_ratio = NADH_NAD_ratio,
                                 standMat = standMat)
        }else if(is.null(standConNAD) == F){
                NAD_NADH <- list(NAD = NAD,
                                 NADH = NADH,
                                 NADH_NAD_ratio = NADH_NAD_ratio,
                                 standMat = standMat,
                                 standMatNAD = standMatNAD)
        }
        return(NAD_NADH)
}

NADH071420_data <- getNADHDATA(NADH071420, 
                               minutes = 30, 
                               strains =  c("PA14", "F22031", "F34365", "F63921", "F9670", "H27930", "H5708", "M6075", "S86968", "T6313"),
                               standConc = c(1000, 500, 400, 200, 100, 50, 25, 10, 0))

doScaleReg <- function(NADHData,
                       YposFormula = 320, 
                       YposStat = 250){
        conc <- as.numeric(colnames(NADHData$standMat))
        avLum <- apply(NADHData$standMat, 2, mean)
        scaleDF <- data.frame(conc = conc,
                              avLum = avLum, 
                              R1 = NADHData$standMat[1, ],
                              R2 = NADHData$standMat[2, ],
                              R3 = NADHData$standMat[3, ])
       p <- ggscatter(scaleDF, x = "conc", y = "avLum", add = "reg.line") +
               stat_cor(label.x = 0, label.y = YposStat) +
               stat_regline_equation(label.x = 0, label.y = YposFormula)
       return(p)
}

doScaleReg(NADH071420_data)

timePointList <- list()
timePoints <- c(30, 60, 90, 120, 150, 180, 210, 240)
for(i in 1:length(timePoints)){
        print(i)
        D <- getNADHDATA(NADH071420, 
                         minutes = timePoints[i], 
                         strains =  c("PA14", "F22031", "F34365", "F63921", "F9670", "H27930", "H5708", "M6075", "S86968", "T6313"),
                         standConc = c(1000, 500, 400, 200, 100, 50, 25, 10, 0))
        timePointList[[i]] <- D
}


doScaleReg(timePointList[[1]])
doScaleReg(timePointList[[2]], YposFormula = 500)
doScaleReg(timePointList[[3]], YposFormula = 500)
doScaleReg(timePointList[[4]], YposFormula = 500)
doScaleReg(timePointList[[5]], YposStat = 700) #----> Best fit 
doScaleReg(timePointList[[6]], YposStat = 1000)
doScaleReg(timePointList[[7]])
doScaleReg(timePointList[[8]])

NADH071420_data <- getNADHDATA(NADH071420, 
                               minutes = 180, 
                               strains =  c("PA14", "F22031", "F34365", "F63921", "F9670", "H27930", "H5708", "M6075", "S86968", "T6313"),
                               standConc = c(1000, 500, 400, 200, 100, 50, 25, 10, 0))


NADH071620_data <- getNADHDATA(NADH071620, 
                               minutes = 30, 
                               strains =  c("PA14", "H47912", "W45909", "W60856", "W70332", "X78812", "X9820"),
                               standConc = c(10000, 1000, 500, 400, 200, 100, 50, 25, 10, 0))

timePointList <- list()
for(i in 1:length(timePoints)){
        print(i)
        D <- getNADHDATA(NADH071620, 
                         minutes = timePoints[i], 
                         strains =  c("PA14", "H47912", "W45909", "W60856", "W70332", "X78812", "X9820"),
                         standConc = c(10000, 1000, 500, 400, 200, 100, 50, 25, 10, 0))
        timePointList[[i]] <- D
}

doScaleReg(timePointList[[1]], YposFormula = 1000)
doScaleReg(timePointList[[2]], YposFormula = 1000) #----> Best fit 
doScaleReg(timePointList[[3]], YposFormula = 1000)
doScaleReg(timePointList[[4]], YposFormula = 1000)
doScaleReg(timePointList[[5]], YposStat = 10000) 
doScaleReg(timePointList[[6]], YposStat = 1000)
doScaleReg(timePointList[[7]], YposFormula = 1000)
doScaleReg(timePointList[[8]], YposFormula = 1000)

NADH071620_data <- getNADHDATA(NADH071620, 
                               minutes = 60, 
                               strains =  c("PA14", "H47912", "W45909", "W60856", "W70332", "X78812", "X9820"),
                               standConc = c(10000, 1000, 500, 400, 200, 100, 50, 25, 10, 0))

boxplot(as.numeric(NADH071620_data$NADH_NAD_ratio[, c(1, 2, 3, 6, 7)]),
        NADH071620_data$NADH_NAD_ratio[, 5],
        NADH071620_data$NADH_NAD_ratio[, 4])


timePointList <- list()
for(i in 1:length(timePoints)){
        print(i)
        D <- getNADHDATA(NADH072120, 
                         minutes = timePoints[i], 
                         strains =  c("T38079", 
                                      "H27930", 
                                      "M74707", 
                                      "T52373", 
                                      "W25637", 
                                      "F22031", 
                                      "PA14", 
                                      "F34365", 
                                      "F63912", 
                                      "H5708", 
                                      "S86968", 
                                      "T6313", 
                                      "M1608", 
                                      "F30658", 
                                      "M37351", 
                                      "T63266", 
                                      "W16407", 
                                      "F9670", 
                                      "M6075", 
                                      "W36662"),
                         standConc = c(500, 250, 100, 10, 0))
        timePointList[[i]] <- D
}

doScaleReg(timePointList[[1]], YposFormula = 1000)
doScaleReg(timePointList[[2]], YposFormula = 1000) #----> Best fit 
doScaleReg(timePointList[[3]], YposFormula = 1000)
doScaleReg(timePointList[[4]], YposFormula = 1000)
doScaleReg(timePointList[[5]], YposStat = 10000) 
doScaleReg(timePointList[[6]], YposStat = 1000)
doScaleReg(timePointList[[7]], YposFormula = 1000)
doScaleReg(timePointList[[8]], YposFormula = 1000)

NADH072120_data <- getNADHDATA(NADH072120, 
                               minutes = 60, 
                               strains =  c("T38079", 
                                            "H27930", 
                                            "M74707", 
                                            "T52373", 
                                            "W25637", 
                                            "F22031", 
                                            "PA14", 
                                            "F34365", 
                                            "F63912", 
                                            "H5708", 
                                            "S86968", 
                                            "T6313", 
                                            "M1608", 
                                            "F30658", 
                                            "M37351", 
                                            "T63266", 
                                            "W16407", 
                                            "F9670", 
                                            "M6075", 
                                            "W36662"),
                               standConc = c(500, 250, 100, 10, 0))


timePointList <- list()
for(i in 1:length(timePoints)){
        print(i)
        D <- getNADHDATA(NADH072320, 
                         minutes = timePoints[i],
                         strains =  c("PAO1", 
                                      "T38079", 
                                      "F5677", 
                                      "PA7", 
                                      "W91453", 
                                      "F23197", 
                                      "PA14", 
                                      "M55212"),
                         standConc = c(500, 250, 100, 10, 0),
                         standConNAD = c(1000, 500, 400, 200, 100, 50, 25, 10, 0),
                         substractZero = F)
        timePointList[[i]] <- D
}

doScaleReg(timePointList[[1]], YposFormula = 1000)
doScaleReg(timePointList[[2]], YposFormula = 1000) 
doScaleReg(timePointList[[3]], YposFormula = 1000)
doScaleReg(timePointList[[4]], YposFormula = 1000)
doScaleReg(timePointList[[5]], YposStat = 10000) 
doScaleReg(timePointList[[6]], YposStat = 1000) #--> Best fit
doScaleReg(timePointList[[7]], YposFormula = 1000)
doScaleReg(timePointList[[8]], YposFormula = 1000)


NADH072320_data <- getNADHDATA(NADH072320, 
                               minutes = 180,
                               strains =  c("PAO1", 
                                            "T38079", 
                                            "F5677", 
                                            "PA7", 
                                            "W91453", 
                                            "F23197", 
                                            "PA14", 
                                            "M55212"),
                               standConc = c(500, 250, 100, 10, 0),
                               standConNAD = c(1000, 500, 400, 200, 100, 50, 25, 10, 0),
                               substractZero = F)


timePointList <- list()
for(i in 1:length(timePoints)){
        print(i)
        D <- getNADHDATA(NADH072720, 
                         minutes = timePoints[i],
                         strains =  c("F22031", 
                                      "F34365", 
                                      "F9670", 
                                      "H5708", 
                                      "M6075", 
                                      "S86968", 
                                      "T52373", 
                                      "T63266", 
                                      "W25637", 
                                      "F63912", 
                                      "H27930", 
                                      "M37351", 
                                      "M74707", 
                                      "T6313", 
                                      "PA14", 
                                      "W16407", 
                                      "T38079",
                                      "W91453", 
                                      "M1608", 
                                      "W36662", 
                                      "F30658", 
                                      "M1608 4C", 
                                      "W36662 4C", 
                                      "F30658 4C", 
                                      "T38079 4C", 
                                      "W91453 4C"),
                         standConc = c(500, 250, 100, 10, 0),
                         standConNAD = c(1000, 500, 400, 200, 100, 50, 25, 10, 0), 
                         substractZero = F)
        timePointList[[i]] <- D
}

doScaleReg(timePointList[[1]], YposFormula = 1000)
doScaleReg(timePointList[[2]], YposFormula = 1000) 
doScaleReg(timePointList[[3]], YposFormula = 1000)
doScaleReg(timePointList[[4]], YposFormula = 1000)
doScaleReg(timePointList[[5]], YposStat = 10000) 
doScaleReg(timePointList[[6]], YposStat = 1000)#----> Best fit 
doScaleReg(timePointList[[7]], YposFormula = 1000)
doScaleReg(timePointList[[8]], YposFormula = 1000)


NADH072720_data <- getNADHDATA(NADH072720, 
                               minutes = 180,
                               strains =  c("F22031", 
                                            "F34365", 
                                            "F9670", 
                                            "H5708", 
                                            "M6075", 
                                            "S86968", 
                                            "T52373", 
                                            "T63266", 
                                            "W25637", 
                                            "F63912", 
                                            "H27930", 
                                            "M37351", 
                                            "M74707", 
                                            "T6313", 
                                            "PA14", 
                                            "W16407", 
                                            "T38079",
                                            "W91453", 
                                            "M1608", 
                                            "W36662", 
                                            "F30658", 
                                            "M1608 4C", 
                                            "W36662 4C", 
                                            "F30658 4C", 
                                            "T38079 4C", 
                                            "W91453 4C"),
                               standConc = c(500, 250, 100, 10, 0),
                               standConNAD = c(1000, 500, 400, 200, 100, 50, 25, 10, 0), 
                               substractZero = F)


boxplot(c(NADH072120_data$NADH_NAD_ratio[-2, 12],
         as.numeric(NADH072120_data$NADH_NAD_ratio[, c(1, 3, 4, 5, 6, 7, 8, 18)])),
        c(NADH072120_data$NADH_NAD_ratio[-1, 10],
          as.numeric(NADH072120_data$NADH_NAD_ratio[, c(14, 15, 16, 17, 19)])),
        c(as.numeric(NADH072120_data$NADH_NAD_ratio[, c(2, 9, 11, 13)]),
          NADH072120_data$NADH_NAD_ratio[-3, 20]),
        names = c("+", "+/-", "-"))

# Plot the batches to see if there was time between batches influenced nadh/nad+ concentration

boxplot(as.numeric(NADH072120_data$NADH_NAD_ratio[, 1]),
        as.numeric(NADH072120_data$NADH_NAD_ratio[, 2]),
        as.numeric(NADH072120_data$NADH_NAD_ratio[, c(3, 4, 5, 6)]),
        as.numeric(NADH072120_data$NADH_NAD_ratio[, c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]),
        as.numeric(NADH072120_data$NADH_NAD_ratio[, c(18, 19, 20)]))

boxplot(as.numeric(NADH072120_data$NADH[, 1]),
        as.numeric(NADH072120_data$NADH[, 2]),
        as.numeric(NADH072120_data$NADH[, c(3, 4, 5, 6)]),
        as.numeric(NADH072120_data$NADH[, c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]),
        as.numeric(NADH072120_data$NADH[, c(18, 19, 20)]))

boxplot(as.numeric(NADH072120_data$NAD[, 1]),
        as.numeric(NADH072120_data$NAD[, 2]),
        as.numeric(NADH072120_data$NAD[, c(3, 4, 5, 6)]),
        as.numeric(NADH072120_data$NAD[, c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]),
        as.numeric(NADH072120_data$NAD[, c(18, 19, 20)]))



allExps <- readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves/Data_All_200801_withConc.xlsx")

allExps <- as.data.frame(allExps)

exps <- unique(allExps$date)

exp_71620Standard_nadh <- allExps[allExps$date == exps[1], ][grep("standard", allExps[allExps$date == exps[1], ]$strain), c(1, 2, 4)]

exp_71620Standard_nadh$strain <- as.numeric(gsub("NADHstandard", "", exp_71620Standard_nadh$strain))

exp_72020Standard_nadh <- allExps[allExps$date == exps[2], ][grep("standard", allExps[allExps$date == exps[2], ]$strain), c(1, 2, 4)]

exp_72020Standard_nadh$strain <- as.numeric(gsub("NADHstandard", "", exp_72020Standard_nadh$strain))

exp_72320Standard_nadh <- allExps[allExps$date == exps[3], ][grep("NADHstandard", allExps[allExps$date == exps[3], ]$strain), c(1, 2, 4)]

exp_72320Standard_nadh$strain <- as.numeric(gsub("NADHstandard", "", exp_72320Standard_nadh$strain))

exp_72320Standard_nad <- allExps[allExps$date == exps[3], ][grep("NADstandard", allExps[allExps$date == exps[3], ]$strain), c(9, 2, 4)]

exp_72320Standard_nad$strain <- as.numeric(gsub("NADstandard", "", exp_72320Standard_nad$strain))

exp_72620Standard_nadh <- allExps[allExps$date == exps[4], ][grep("NADHstandard", allExps[allExps$date == exps[4], ]$strain), c(1, 2, 4)]

exp_72620Standard_nadh$strain <- as.numeric(gsub("NADHstandard", "", exp_72620Standard_nadh$strain))

exp_72620Standard_nad <- allExps[allExps$date == exps[4], ][grep("NADstandard", allExps[allExps$date == exps[4], ]$strain), c(9, 2, 4)]

exp_72620Standard_nad$strain <- as.numeric(gsub("NADstandard", "", exp_72620Standard_nad$strain))



nadhStandards <-  rbind.data.frame(exp_71620Standard_nadh,
                                   exp_72020Standard_nadh,
                                   exp_72320Standard_nadh,
                                   exp_72620Standard_nadh)

nadhStandards$date <- as.factor(nadhStandards$date)

colnames(nadhStandards) <- c("luminescence", "date", "nadh_concentration")

nadStandards <- rbind.data.frame(exp_72320Standard_nad,
                                 exp_72620Standard_nad)

nadStandards$date <- as.factor(nadStandards$date)

colnames(nadStandards) <- c("luminescence", "date", "nad_concentration")


if(!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

ggplot(nadhStandards, aes(x = nadh_concentration, y = luminescence, color = date) ) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = c("red", "blue", "green", "purple"))

ggsave("nadhScales.tiff")

ggplot(nadStandards, aes(x = nad_concentration, y = luminescence, color = date) ) +
        geom_point() +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = c("green", "purple"))

ggsave("nadScales.tiff")



exp1 <- allExps[allExps$date == exps[1], ]

zeroExp1 <- mean(nadhStandards[nadhStandards$date == exps[1], ][nadhStandards[nadhStandards$date == exps[1], ]$nadh_concentration == 0, ]$luminescence)

exp1$nadh <- exp1$nadh - zeroExp1
exp1$nad <- exp1$nad - zeroExp1

exp1 <- exp1[-grep("standard", exp1$strain), ]

exp1$nadh <- exp1$nadh/exp1$od
exp1$nad <- exp1$nad/exp1$od

normF_PA14NADHExp1 <- mean(exp1[exp1$strain == "PA14", ]$nadh)
normF_PA14NADExp1 <- mean(exp1[exp1$strain == "PA14", ]$nad)

exp1$nadh <- exp1$nadh/normF_PA14NADHExp1
exp1$nad <- exp1$nad/normF_PA14NADExp1

exp1NADH_NAD_ratio <- exp1$nadh/exp1$nad

exp1 <- cbind.data.frame(exp1, exp1NADH_NAD_ratio)

exp1 <- exp1[, -c(5, 10, 11)]

colnames(exp1)[ncol(exp1)] <- "NADH_NAD_ratio"



exp2 <- allExps[allExps$date == exps[2], ]

zeroExp2 <- mean(nadhStandards[nadhStandards$date == exps[2], ][nadhStandards[nadhStandards$date == exps[2], ]$nadh_concentration == 0, ]$luminescence)

exp2$nadh <- exp2$nadh - zeroExp2
exp2$nad <- exp2$nad - zeroExp2

exp2 <- exp2[-grep("standard", exp2$strain), ]

exp2$nadh <- exp2$nadh/exp2$od
exp2$nad <- exp2$nad/exp2$od

normF_PA14NADHExp2 <- mean(exp2[exp2$strain == "PA14", ]$nadh)
normF_PA14NADExp2 <- mean(exp2[exp2$strain == "PA14", ]$nad)

exp2$nadh <- exp2$nadh/normF_PA14NADHExp2
exp2$nad <- exp2$nad/normF_PA14NADExp2

exp2NADH_NAD_ratio <- exp2$nadh/exp2$nad

exp2 <- cbind.data.frame(exp2, exp2NADH_NAD_ratio)

exp2 <- exp2[, -c(5, 10, 11)]

colnames(exp2)[ncol(exp2)] <- "NADH_NAD_ratio"



exp3 <- allExps[allExps$date == exps[3], ]

zeroExp3 <- mean(nadhStandards[nadhStandards$date == exps[3], ][nadhStandards[nadhStandards$date == exps[3], ]$nadh_concentration == 0, ]$luminescence)

exp3$nadh <- exp3$nadh - zeroExp3
exp3$nad <- exp3$nad - zeroExp3

exp3 <- exp3[-grep("standard", exp3$strain), ]

exp3$nadh <- exp3$nadh/exp3$od
exp3$nad <- exp3$nad/exp3$od

normF_PA14NADHExp3 <- mean(exp3[exp3$strain == "PA14", ]$nadh)
normF_PA14NADExp3 <- mean(exp3[exp3$strain == "PA14", ]$nad)

exp3$nadh <- exp3$nadh/normF_PA14NADHExp3
exp3$nad <- exp3$nad/normF_PA14NADExp3

exp3NADH_NAD_ratio <- exp3$nadh/exp3$nad

exp3 <- cbind.data.frame(exp3, exp3NADH_NAD_ratio)

exp3 <- exp3[, -c(5, 10, 11)]

colnames(exp3)[ncol(exp3)] <- "NADH_NAD_ratio"


exp4 <- allExps[allExps$date == exps[4], ]

zeroExp4 <- mean(nadhStandards[nadhStandards$date == exps[4], ][nadhStandards[nadhStandards$date == exps[4], ]$nadh_concentration == 0, ]$luminescence)

exp4$nadh <- exp4$nadh - zeroExp4
exp4$nad <- exp4$nad - zeroExp4

exp4 <- exp4[-grep("standard", exp4$strain), ]

exp4$nadh <- exp4$nadh/exp4$od
exp4$nad <- exp4$nad/exp4$od

normF_PA14NADHExp4 <- mean(exp4[exp4$strain == "PA14", ]$nadh)
normF_PA14NADExp4 <- mean(exp4[exp4$strain == "PA14", ]$nad)

exp4$nadh <- exp4$nadh/normF_PA14NADHExp4
exp4$nad <- exp4$nad/normF_PA14NADExp4

exp4NADH_NAD_ratio <- exp4$nadh/exp4$nad

exp4 <- cbind.data.frame(exp4, exp4NADH_NAD_ratio)

exp4 <- exp4[, -c(5, 10, 11)]

exp4 <- exp4[!exp4$protocol == "P1", ]

colnames(exp4)[ncol(exp4)] <- "NADH_NAD_ratio"

expsNormRatiosData <- rbind.data.frame(exp1,
                                       exp2,
                                       exp3, 
                                       exp4)

expsNormRatiosData$nadh <- expsNormRatiosData$nadh + round(abs(min(expsNormRatiosData$nadh)), digits = 4)+1
expsNormRatiosData$nad <- expsNormRatiosData$nad + round(abs(min(expsNormRatiosData$nad)), digits = 4)+1

expsNormRatiosData$NADH_NAD_ratio <- expsNormRatiosData$nadh/expsNormRatiosData$nad

ggplot(expsNormRatiosData, aes(x = rhl2cats, y = NADH_NAD_ratio, color = rhl2cats)) + 
        geom_boxplot()

ggsave("NADH_NAD_ratio_boxplot.tiff")

wilcox.test(expsNormRatiosData$NADH_NAD_ratio[expsNormRatiosData$rhl2cats == 0], 
            expsNormRatiosData$NADH_NAD_ratio[expsNormRatiosData$rhl2cats == 1])


ggplot(exp1, aes(x = rhl2cats, y = NADH_NAD_ratio, color = rhl2cats)) + 
        geom_boxplot()

ggplot(exp2, aes(x = rhl2cats, y = NADH_NAD_ratio, color = rhl2cats)) + 
        geom_boxplot()

ggplot(exp3, aes(x = rhl2cats, y = NADH_NAD_ratio, color = rhl2cats)) + 
        geom_boxplot()

ggplot(exp4, aes(x = rhl2cats, y = NADH_NAD_ratio, color = rhl2cats)) + 
        geom_boxplot()
