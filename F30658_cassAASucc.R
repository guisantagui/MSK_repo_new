setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growthCurves")

casSuccCurves <- readxl::read_xlsx("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/growth_curves_scripts/031320Cas_casSuccGrowthCurve_Hildi_mean.xlsx")

if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggplot)) install.packages("ggplot")

F30658 <- data.frame(time = c(casSuccCurves$Time, casSuccCurves$Time), 
                     OD = c(casSuccCurves$F30658, casSuccCurves$F30658_suc),
                     medium = c(rep("casAA", length(casSuccCurves$F30658)),
                                rep("casAA + succ", length(casSuccCurves$F30658_suc))))



p1 <- ggplot(F30658, aes(x= time, y = OD, colour = medium)) +
        geom_line()

p1

ggsave("F30658CasAA_casAASucc.pdf")

