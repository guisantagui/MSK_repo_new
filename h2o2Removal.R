#croissance

setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/h202Removal")

h2o2Rem <- readxl::read_xlsx("181130_glycerol_Amplex.xlsx")

plot(h2o2Rem$`Time [s]`, h2o2Rem$A13)

C1 <- h2o2Rem[, grep("1", colnames(h2o2Rem), fixed = T)[2:4]]
PA14 <- h2o2Rem[, grep("2", colnames(h2o2Rem), fixed = T)[2:4]]
W25637 <- h2o2Rem[, grep("3", colnames(h2o2Rem), fixed = T)[2:4]]
W36662 <- h2o2Rem[, grep("4", colnames(h2o2Rem), fixed = T)[2:4]]
M1608 <- h2o2Rem[, grep("5", colnames(h2o2Rem), fixed = T)[2:4]]
W91453 <- h2o2Rem[, grep("6", colnames(h2o2Rem), fixed = T)[2:4]]
F22031 <- h2o2Rem[, grep("7", colnames(h2o2Rem), fixed = T)[2:4]]
F9670 <- h2o2Rem[, grep("8", colnames(h2o2Rem), fixed = T)[2:4]]
M37351 <- h2o2Rem[, grep("9", colnames(h2o2Rem), fixed = T)[2:4]]
T38079 <- h2o2Rem[, grep("10", colnames(h2o2Rem), fixed = T)[2:4]]
X9820 <- h2o2Rem[, grep("11", colnames(h2o2Rem), fixed = T)[2:4]]
C2 <- h2o2Rem[, grep("12", colnames(h2o2Rem), fixed = T)[2:4]]


strains <- list(PA14 = PA14, 
                W25637 = W25637,
                W36662 = W36662,
                M1608 = M1608,
                W91453 = W91453,
                F22031 = F22031,
                F9670 = F9670,
                M37351 = M37351,
                T38079 = T38079, 
                X9820 = X9820)

lapply(strains, function(x) x - C1)

apply(C1, 2, function(x) plot(h2o2Rem$`Time [s]`, x))

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat.RData")

rhamnSwarmMat$rhamn2cats[match(names(strains), rhamnSwarmMat$strain)]

rhamnPosH2O2 <- strains[rhamnSwarmMat$rhamn2cats[match(names(strains), rhamnSwarmMat$strain)] == 1]

rhamnNegH2O2 <- strains[rhamnSwarmMat$rhamn2cats[match(names(strains), rhamnSwarmMat$strain)] == 0]

boxPlotDF <- data.frame("H2O2 removal" = c(unlist(rhamnPosH2O2),
                                           unlist(rhamnNegH2O2)),
                        "rhamnolipid production" = factor(c(rep("rhamnolipid producer", 
                                                                length(unlist(rhamnPosH2O2))),
                                                            rep("rhamnolipid non-producer", 
                                                                length(unlist(rhamnNegH2O2))))))

ggplot(boxPlotDF, aes(x=rhamnolipid.production, y=H2O2.removal)) + 
        geom_boxplot()

wilcox.test(unlist(rhamnPosH2O2),
            unlist(rhamnNegH2O2))$p.value

boxplot(unlist(rhamnPosH2O2),
        unlist(rhamnNegH2O2))

ggsave("h2o2RemVSRhamn.png")

