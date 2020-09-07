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

NADH071620 <- NADH071420[36:(nrow(NADH071420)-5), ]
colnames(NADH071620) <- as.character(NADH071620[1, ])
NADH071620 <- as.data.frame(apply(NADH071620[-1, ], 
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
        if(is.null(standConNAD) == F){
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
                                standMatNAD[i, j] <- experimentMatrix[position, colnames(experimentMatrix) == standCols[j, i]]
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
                                 standMatNADH = standMat,
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
                              R1 = NADHData$standMat[, 1],
                              R2 = NADHData$standMat[, 2],
                              R3 = NADHData$standMat[, 3])
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



boxplot(c(6.9101, 7.6228, 4.4505,
          5.3847, 5.2449, 5.3461,
          3.9892, 5.1938,
          7.1745, 8.1230, 7.2127,
          5.7316, 6.9046, 7.2310,
          9.3150, 8.1775, 9.3508,
          6.5093, 7.8384, 7.6957, 
          3.1830, 5.1488, 5.5099,
          3.3634, 2.7942, 3.1892),
        c(3.3387, 1.4738, 3.0303, 
          4.3091, 4.0473, 3.6863,
          7.5340, 4.9329, 4.6677))

wilcox.test(c(6.9101, 7.6228, 4.4505,
              5.3847, 5.2449, 5.3461,
              3.9892, 5.1938,
              7.1745, 8.1230, 7.2127,
              5.7316, 6.9046, 7.2310,
              9.3150, 8.1775, 9.3508,
              6.5093, 7.8384, 7.6957, 
              3.1830, 5.1488, 5.5099,
              3.3634, 2.7942, 3.1892),
            c(3.3387, 1.4738, 3.0303, 
              4.3091, 4.0473, 3.6863,
              7.5340, 4.9329, 4.6677))

t.test(c(6.9101, 7.6228, 4.4505,
         5.3847, 5.2449, 5.3461,
         3.9892, 5.1938,
         7.1745, 8.1230, 7.2127,
         5.7316, 6.9046, 7.2310,
         9.3150, 8.1775, 9.3508,
         6.5093, 7.8384, 7.6957, 
         3.1830, 5.1488, 5.5099,
         3.3634, 2.7942, 3.1892),
       c(3.3387, 1.4738, 3.0303, 
         4.3091, 4.0473, 3.6863,
         7.5340, 4.9329, 4.6677))

hist(c(6.9101, 7.6228, 4.4505,
       5.3847, 5.2449, 5.3461,
       3.9892, 5.1938,
       7.1745, 8.1230, 7.2127,
       5.7316, 6.9046, 7.2310,
       9.3150, 8.1775, 9.3508,
       6.5093, 7.8384, 7.6957, 
       3.1830, 5.1488, 5.5099,
       3.3634, 2.7942, 3.1892))

shapiro.test(c(6.9101, 7.6228, 4.4505,
               5.3847, 5.2449, 5.3461,
               3.9892, 5.1938,
               7.1745, 8.1230, 7.2127,
               5.7316, 6.9046, 7.2310,
               9.3150, 8.1775, 9.3508,
               6.5093, 7.8384, 7.6957, 
               3.1830, 5.1488, 5.5099,
               3.3634, 2.7942, 3.1892))



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





