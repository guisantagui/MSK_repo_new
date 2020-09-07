
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/CCMNNorm/ccmnNormMets.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/rhamnMat.RData")

glnGluRatio <- data.frame(ratio = ccmnNormMets$glutamine/ccmnNormMets$glutamate,
                          rhamn = rhamnMat$rhamn2cats[match(gsub("\\_.*", 
                                                                 rownames(ccmnNormMets), 
                                                                 rep = "\\1"),
                                                            rhamnMat$strains)])

wilcox.test(glnGluRatio$ratio[glnGluRatio$rhamn == 0],
            glnGluRatio$ratio[glnGluRatio$rhamn == 1])

t.test(glnGluRatio$ratio[glnGluRatio$rhamn == 0],
       glnGluRatio$ratio[glnGluRatio$rhamn == 1])

rownames(ccmnNormMets) <- gsub("W70322", "W70332", rownames(ccmnNormMets))

boxplot(glnGluRatio$ratio[glnGluRatio$rhamn == 0],
        glnGluRatio$ratio[glnGluRatio$rhamn == 1])
