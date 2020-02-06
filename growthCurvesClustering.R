setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/growthCurvesClustering")
if(!require(kml)) install.packages("kml")
library(kml)

# load curves
if(!require(readxl)) install.packages("readxl")
library(readxl)

hildy_010720 <- read_xlsx("/Users/santamag/010720GlycerolGrowthCurve_Hildi_mean.xlsx")
maurice_010720 <- read_xlsx("/Users/santamag/010720GlycerolGrowthCurve_Maurice_mean.xlsx")
hildy_011520 <- read_xlsx("/Users/santamag/011520GlycerolGrowthCurve_Hildi_mean.xlsx")
maurice_011520 <- read_xlsx("/Users/santamag/011520GlycerolGrowthCurve_Maurice_mean.xlsx")
hildy_012120 <- read_xlsx("/Users/santamag/012120GlycerolGrowthCurve_Hildi_mean.xlsx")
spark_012120 <- read_xlsx("/Users/santamag/012120GlycerolGrowthCurve_Spark_mean.xlsx")

# Sparks output is different to the one of the Hildi's and Maurice's one (Values are usually higher
# we observed that PA14 control doesn't replicate with the other ones), so we obtained a linear 
# transformation vor converting Spark to Hildi values

spark2hildy <- function(spark){
        trans <- function(x) {
                y <- 0.7075*x + 0.0114
                return(y)
        }
        transMat <- apply(spark[, 3:ncol(spark)], 2, trans)
        transMat <- cbind.data.frame(spark[, 1:2],
                                     transMat)
        return(transMat)
}

spark2hildyBlanked <- function(spark){
        trans <- function(x) {
                y <- 0.7038*x + 0.0036
                return(y)
        }
        transMat <- apply(spark[, 3:ncol(spark)], 2, trans)
        transMat <- cbind.data.frame(spark[, 1:2],
                                     transMat)
        return(transMat)
}

spark_012120 <- spark2hildyBlanked(spark_012120)
allGrowthCurves <- cbind.data.frame(hildy_010720, 
                                    maurice_010720, 
                                    hildy_011520, 
                                    maurice_011520, 
                                    hildy_012120, 
                                    spark_012120)

allGrowthCurves <- allGrowthCurves[, -c(1, 13, 14, 17, 18, 22, 23, 34, 35, 39, 40)]
View(allGrowthCurves)
plot(allGrowthCurves$Time, allGrowthCurves$PA14)
plot(allGrowthCurves$Time, allGrowthCurves$PA14.1)
plot(allGrowthCurves$Time, allGrowthCurves$PA14.2)
plot(allGrowthCurves$Time, allGrowthCurves$PA14.3)
plot(allGrowthCurves$Time, allGrowthCurves$PA14.4)
plot(allGrowthCurves$Time, allGrowthCurves$PA14.5)

PA14GrowthCurves <- cbind.data.frame(allGrowthCurves$PA14,
                                     allGrowthCurves$PA14.1,
                                     allGrowthCurves$PA14.2,
                                     allGrowthCurves$PA14.3,
                                     allGrowthCurves$PA14.4,
                                     allGrowthCurves$PA14.5)
PA14meansMean <- apply(PA14GrowthCurves, 1, mean)
plot(allGrowthCurves$Time, PA14meansMean)

allGrowthCurves <- allGrowthCurves[, !colnames(allGrowthCurves) %in% c("PA14", 
                                                                       "PA14.1", 
                                                                       "PA14.2", 
                                                                       "PA14.3", 
                                                                       "PA14.4", 
                                                                       "PA14.5")]

allGrowthCurves <- cbind.data.frame(allGrowthCurves, PA14meansMean)
colnames(allGrowthCurves)[ncol(allGrowthCurves)] <- "PA14"

cldObjkt <- cld(t(allGrowthCurves[, 2:ncol(allGrowthCurves)]), 
                colnames(allGrowthCurves[2:ncol(allGrowthCurves)]),
                allGrowthCurves$Time)

X11(type="Xlib")
kml(cldObjkt, toPlot = "both")

cldObjkt["idAll"]
cldObjkt["varNames"]
cldObjkt["traj"]

X11(type="Xlib")
choice(cldObjkt, typeGraph = "bmp")

part2 <- partition(clusters=rep(1:2,3),cldObjkt)
part3 <- partition(clusters=rep(1:3,2),cldObjkt)




























