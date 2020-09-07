setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/gc_ms")

gcMsPA <- read.csv("C:/Users/Guillem/Documents/PhD/comput/GC_MS/peak_areas_from_test_diff_conc.csv", stringsAsFactors = F)

gcMsPA <- gcMsPA[-1, -6]

gcMsPA[, 7:ncol(gcMsPA)] <- apply(gcMsPA[, 7:ncol(gcMsPA)], 2, as.numeric)

samps <- unique(sapply(gcMsPA$X.2, function(x) strsplit(x, "_")[[1]][2]))
sampID <- sapply(gcMsPA$X.2, function(x) strsplit(x, "_")[[1]][2])


pVals <- c()
for(i in 7:ncol(gcMsPA)){
        sampList <- list()
        for(j in 1:length(samps)){
                sampvec <- gcMsPA[sampID == samps[j], i]
                print(sampvec)
                sampList[[j]] <- sampvec
        }
        kResult <- kruskal.test(sampList)$p.value
        pVals <- c(pVals, kResult)
}
names(pVals) <- colnames(gcMsPA)[7:ncol(gcMsPA)]

p.adjust(pVals, method = "BH")
