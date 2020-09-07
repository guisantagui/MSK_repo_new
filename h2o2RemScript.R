setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/h2o2Removal")

if(!require(ggplot2)) install.packages("ggplot2")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat_2.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat.RData")

h2o2_113018 <- readxl::read_xlsx("181130_glycerol_Amplex.xlsx", sheet = 2)

h2o2_113018 <- as.data.frame(apply(h2o2_113018, 2, as.numeric))

strains_113018 <- c('PA14','W25637','W36662','M1608','W91453','F22031','F9670','M37351','T38079','X9820')

forPlot_113018 <- data.frame(time = c(h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`,
                                      h2o2_113018$`Time [s]`),
                             amplex = c(h2o2_113018$B12 - h2o2_113018$B2,
                                        h2o2_113018$B12 - h2o2_113018$B3,
                                        h2o2_113018$B12 - h2o2_113018$B4,
                                        h2o2_113018$B12 - h2o2_113018$B5,
                                        h2o2_113018$B12 - h2o2_113018$B6,
                                        h2o2_113018$B12 - h2o2_113018$B7,
                                        h2o2_113018$B12 - h2o2_113018$B8,
                                        h2o2_113018$B12 - h2o2_113018$B9,
                                        h2o2_113018$B12 - h2o2_113018$B10,
                                        h2o2_113018$B12 - h2o2_113018$B11),
                             rhamn = c(rep("PA14", length(h2o2_113018$`Time [s]`)),
                                       rep("W25637", length(h2o2_113018$`Time [s]`)),
                                       rep("W36662", length(h2o2_113018$`Time [s]`)),
                                       rep("M1608", length(h2o2_113018$`Time [s]`)),
                                       rep("W91453", length(h2o2_113018$`Time [s]`)),
                                       rep("F22031", length(h2o2_113018$`Time [s]`)),
                                       rep("F9670", length(h2o2_113018$`Time [s]`)),
                                       rep("M37351", length(h2o2_113018$`Time [s]`)),
                                       rep("T38079", length(h2o2_113018$`Time [s]`)),
                                       rep("pos", length(h2o2_113018$`Time [s]`))))

ggplot(forPlot_113018, aes(x= time, y = amplex, colour = rhamn)) +
        geom_line()

h2o2_062420 <- readxl::read_xlsx("062420H2O2Removal_Hildi.xlsx")
h2o2_063020_Hildi <- readxl::read_xlsx("063020H2O2Production_Hildi.xlsx")
h2o2_063020_Maurice <- readxl::read_xlsx("063020H2O2Production_Maurice.xlsx")

h2o2_062420_od600 <- h2o2_062420[59:347, 1:99]
colnames_h2o2_062420 <- as.character(h2o2_062420_od600[1, ])
colnames(h2o2_062420_od600) <- colnames_h2o2_062420
h2o2_062420_od600 <- as.data.frame(apply(h2o2_062420_od600[-1, ], 2, as.numeric))

h2o2_062420_amplex <- as.data.frame(apply(h2o2_062420[642:929, 1:99], 2, as.numeric))
colnames(h2o2_062420_amplex) <- colnames_h2o2_062420

h2o2_063020_Hildi_od600 <- h2o2_063020_Hildi[60:348, 1:75]
colnames_h2o2_063020_Hildi <- as.character(h2o2_063020_Hildi_od600[1, ])
colnames(h2o2_063020_Hildi_od600) <- colnames_h2o2_063020_Hildi
h2o2_063020_Hildi_od600 <- as.data.frame(apply(h2o2_063020_Hildi_od600[-1, ], 2, as.numeric))

h2o2_063020_Hildi_amplex <- as.data.frame(apply(h2o2_063020_Hildi[643:930, 1:75], 2, as.numeric))
colnames(h2o2_063020_Hildi_amplex) <- colnames_h2o2_063020_Hildi

h2o2_063020_Maurice_od600 <- h2o2_063020_Maurice[59:347, 1:83]
colnames_h2o2_063020_Maurice <- as.character(h2o2_063020_Maurice_od600[1, ])
colnames(h2o2_063020_Maurice_od600) <- colnames_h2o2_063020_Maurice
h2o2_063020_Maurice_od600 <- as.data.frame(apply(h2o2_063020_Maurice_od600[-1, ], 2, as.numeric))

h2o2_063020_Maurice_amplex <- as.data.frame(apply(h2o2_063020_Maurice[642:929, 1:83], 2, as.numeric))
colnames(h2o2_063020_Maurice_amplex) <- colnames_h2o2_063020_Maurice




getReplicatesRemov <- function(amplexMat, strainNames){
        mat <- amplexMat[, 4:ncol(amplexMat)]
        strainDFs <- list()
        for(i in 1:(length(strainNames))){
                repPos <- c((i*8-7):(i*8))[2:4]
                strainDF <- data.frame(time = rep(amplexMat$`Time [s]`, 3),
                                       amplex = c(mat[, 2] - mat[, repPos[1]],
                                                  mat[, 3] - mat[, repPos[2]],
                                                  mat[, 4] - mat[, repPos[3]]),
                                       repl = c(rep("R1", length(amplexMat$`Time [s]`)),
                                                rep("R2", length(amplexMat$`Time [s]`)),
                                                rep("R3", length(amplexMat$`Time [s]`))))
                strainDFs[[i]] <- strainDF
        }
        names(strainDFs) <- strainNames
        return(strainDFs)
}

repsRem113018 <- getReplicatesRemov(h2o2_113018, c('PA14',
                                                   'W25637',
                                                   'W36662',
                                                   'M1608',
                                                   'W91453',
                                                   'F22031',
                                                   'F9670',
                                                   'M37351',
                                                   'T38079',
                                                   'X9820'))

repsRem062420 <- getReplicatesRemov(h2o2_062420_amplex, c("PBS",
                                                          "PA14", 
                                                          "F30658", 
                                                          "F34365", 
                                                          "F63912",
                                                          "H27930",
                                                          "H47921",
                                                          "H5708",
                                                          "M6075",
                                                          "M74707",
                                                          "S86968",
                                                          "H2O2"))

repsRem063020_Maurice <- getReplicatesRemov(h2o2_063020_Maurice_amplex, c("PBS_2",
                                                                          "PA14_2",
                                                                          "F23197",
                                                                          "F5677",
                                                                          "M55212",
                                                                          "PA7",
                                                                          "PAO1",
                                                                          "T52373",
                                                                          "T6313",
                                                                          "H2O2_2"))

repsRem063020_Hildi <- getReplicatesRemov(h2o2_063020_Hildi_amplex, c("PBS_3",
                                                                      "PA14_3",
                                                                      "T63266",
                                                                      "W16407",
                                                                      "W45909",
                                                                      "W60856",
                                                                      "W70332",
                                                                      "X78812",
                                                                      "H2O2_3"))

getPlotsReps <- function(listReps){
        plotList <- list()
        for(i in 1:length(listReps)){
                p <- ggplot(listReps[[i]], aes(x = time, y = amplex, colour = repl))
                plotList[[i]] <- p
        }
        names(plotList) <- names(listReps)
        return(plotList)
}

plotsReps062420 <- getPlotsReps(repsRem062420)

plotList_113018 <- list()
for(i in 1:length(repsRem113018)){
        p <- ggplot(repsRem113018[[i]], aes(x= time, y = amplex, colour = repl)) +
                geom_line()
        plotList_113018[[i]] <- p
}
plotList_113018[[6]]

plotList <- list()
for(i in 1:length(repsRem062420)){
        p <- ggplot(repsRem062420[[i]], aes(x= time, y = amplex, colour = repl)) +
                geom_line()
        plotList[[i]] <- p
}
plotList[[2]]

plotList_063020_Maurice <- list()
for(i in 1:length(repsRem063020_Maurice)){
        p <- ggplot(repsRem063020_Maurice[[i]], aes(x= time, y = amplex, colour = repl)) +
                geom_line()
        plotList_063020_Maurice[[i]] <- p
}
plotList_063020_Maurice[[2]]

plotList_063020_Hildi <- list()
for(i in 1:length(repsRem063020_Hildi)){
        p <- ggplot(repsRem063020_Hildi[[i]], aes(x= time, y = amplex, colour = repl)) +
                geom_line()
        plotList_063020_Hildi[[i]] <- p
}
plotList_063020_Hildi[[2]]

getMeanRemoval <- function(listReps){
        strains <- names(listReps)
        allStrainsMeans <- data.frame()
        for(i in 1:length(listReps)){
                repDF <- data.frame(R1 = listReps[[i]]$amplex[listReps[[i]]$repl == "R1"],
                                    R2 = listReps[[i]]$amplex[listReps[[i]]$repl == "R2"],
                                    R3 = listReps[[i]]$amplex[listReps[[i]]$repl == "R3"])
                meanDF <- data.frame(time = listReps[[1]]$time[1:nrow(repDF)],
                                     amplex = apply(repDF, 1, mean),
                                     strain <- rep(strains[i], nrow(repDF)))
                if(nrow(allStrainsMeans) == 0){
                        allStrainsMeans <- meanDF
                }else{
                        allStrainsMeans <- rbind.data.frame(allStrainsMeans, meanDF)
                }
        }
        colnames(allStrainsMeans) <- c("time", "h2o2removal", "strain")
        rhamnProd <- rhamnSwarmMat_2$glyc_rhamn3cats[match(allStrainsMeans$strain, 
                                                           rhamnSwarmMat_2$strain)]
        rhamnProd[is.na(rhamnProd)] <- 1
        rhamnProd <- as.factor(rhamnProd)
        allStrainsMeans <- cbind.data.frame(allStrainsMeans, rhamnProd)
        colnames(allStrainsMeans)[4] <- c("rhamn")
        return(allStrainsMeans)
}

meanRemoval062420 <- getMeanRemoval(repsRem062420)

meanRemoval062420_norm <- meanRemoval062420
for(i in 1:length(unique(meanRemoval062420$strain))){
        meanRemoval062420_norm[(i*288-287):(i*288), ]$h2o2removal <- meanRemoval062420_norm[(i*288-287):i*288, ]$h2o2removal/meanRemoval062420[grep("PA14", meanRemoval062420$strain), ]$h2o2removal
}

meanRemoval062420[grep("PA14", meanRemoval062420$strain), ]$h2o2removal

meanRemoval063020_Maurice <- getMeanRemoval(repsRem063020_Maurice)
meanRemoval063020_Hildi <- getMeanRemoval(repsRem063020_Hildi)

ggplot(meanRemoval062420, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()

ggsave("h2o2Removal062420.pdf")

ggplot(meanRemoval063020_Maurice, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()

ggsave("h2o2Removal063020_maurice.pdf")

ggplot(meanRemoval063020_Hildi, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()

ggsave("h2o2Removal063020_hildi.pdf")

ggplot(meanRemoval062420, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()

ggsave("h2o2Removal062520.pdf")

allStrainsRemoval <- rbind.data.frame(meanRemoval062420,
                                      meanRemoval063020_Maurice,
                                      meanRemoval063020_Hildi)

ggplot(allStrainsRemoval, aes(x= time, y = h2o2removal, colour = rhamn)) +
        geom_line(aes(group = strain))

ggsave("h2o2RemovalAllStrainsRhamn.pdf")



ggplot(allStrainsRemoval, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()

ggsave("h2o2RemovalAllStrains.pdf")

allH2O2Controls <- rbind.data.frame(meanRemoval062420[grep("H2O2", meanRemoval062420$strain), ],
                                    meanRemoval063020_Maurice[grep("H2O2", meanRemoval063020_Maurice$strain), ],
                                    meanRemoval063020_Hildi[grep("H2O2", meanRemoval063020_Hildi$strain), ])

ggplot(allH2O2Controls, aes(x= time, y = h2o2removal, colour = strain)) +
        geom_line()
