setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurvesFeats")

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3/rhamnSwarmMat_2.RData")

gCurvFeats <- read.csv("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurvesFeats/tblgcfeatures.csv")

p1GRateMax <- data.frame(strain = gCurvFeats$strain[gCurvFeats$phase == 1],
                         maxGRate = gCurvFeats$specific_growth_rate_max[gCurvFeats$phase == 1],
                         rhamn = rhamnSwarmMat_2$glyc_rhamn2cats[match(gCurvFeats$strain[gCurvFeats$phase == 1],
                                                                rhamnSwarmMat_2$strain)])

p1GRateMax <- p1GRateMax[!is.na(p1GRateMax$rhamn),]

mean(p1GRateMax$maxGRate[p1GRateMax$rhamn == 1])

max(p1GRateMax$maxGRate[p1GRateMax$rhamn == 0])

boxplot()