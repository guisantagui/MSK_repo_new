setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/SPIA")

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

if(!require(SPIA)) BiocManager::install("SPIA")
library(SPIA)
