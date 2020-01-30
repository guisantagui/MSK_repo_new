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
spark_010720 <- read_xlsx("/Users/santamag/012120GlycerolGrowthCurve_Spark_mean.xlsx")
