setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/FELLA_analysis")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes as input the differential metabolites observed between the 4 major clusters and does a             #################################
# pathway enrichment using FELLA algorithm                                                                             #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)

set.seed(1)

#Load the differential metabolites between clusters 

diffMets <- read.csv("/Users/santamag/Desktop/GUILLEM/wrkng_dirs/diffMetAnalGood/diffMets_allClusts.csv")
load("../diffMetAnal/oldDataGood/diffMets_oldGood.RData")
load("../diffMetAnal/newData/diffMets_newData.RData")
load("../diffMetAnal/hybrid/diffMets_hybrid.RData")

graph <- buildGraphFromKEGGREST(
        organism = "pae",
        filter.path = c("01100", "01200", "01210", "01212", "01230")
)

buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = NULL,
        internalDir = TRUE,
        matrices = "none",
        normality = "diffusion",
        niter = 100)


alias2entrez <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
entrez2ec <- KEGGREST::keggLink("enzyme", "pae")
entrez2path <- KEGGREST::keggLink("pathway", "pae")
fella.data <- loadKEGGdata(
        #databaseDir = "created__2019-03-25;meta__pae_Release_89.0_03_23_Mar_19",
        internalDir = T,
        loadMatrix = "none"
)

fella.data

id.cpd <- getCom(fella.data, level = 5, format = "id") %>% names
id.rx <- getCom(fella.data, level = 4, format = "id") %>% names
id.ec <- getCom(fella.data, level = 3, format = "id") %>% names

cpd.clusts1_2 <- diffMets_oldGood[, 1]

analysis.clusts1_2 <- enrich(
        compounds = cpd.clusts1_2,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.clusts1_2 %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.clusts1_2)

tiff("FELLA_clusts_1_2_old.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plot(
        analysis.clusts1_2,
        method = "diffusion",
        data = fella.data,
        nlimit = 350,
        plotLegend = FALSE)
dev.off()

g_1_2 <- generateResultsGraph(
        object = analysis.clusts1_2,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_1_2

tiff("FELLA_clusts_1_2_old_2.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_1_2
        #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_1_2 <- generateResultsTable(
        object = analysis.clusts1_2,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

cpd.clusts1.1_1.2 <- diffMets_oldGood[!is.na(diffMets_oldGood[, 2]), 2]

analysis.clusts1.1_1.2 <- enrich(
        compounds = cpd.clusts1.1_1.2,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.clusts1.1_1.2 %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.clusts1.1_1.2)

tiff("FELLA_clusts_1.1_1.2_old.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plot(
        analysis.clusts1.1_1.2,
        method = "diffusion",
        data = fella.data,
        nlimit = 350,
        plotLegend = FALSE)
dev.off()

g_1.1_1.2 <- generateResultsGraph(
        object = analysis.clusts1.1_1.2,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_1.1_1.2

tiff("FELLA_clusts_1.1_1.2_old_2.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_1.1_1.2
        #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_1.1_1.2 <- generateResultsTable(
        object = analysis.clusts1.1_1.2,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

cpd.clusts2.1_2.2 <- diffMets_oldGood[!is.na(diffMets_oldGood[, 3]), 3]

analysis.clusts2.1_2.2 <- enrich(
        compounds = cpd.clusts2.1_2.2,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.clusts2.1_2.2 %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.clusts2.1_2.2)

tiff("FELLA_clusts_2.1_2.2_old.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plot(
        analysis.clusts2.1_2.2,
        method = "diffusion",
        data = fella.data,
        nlimit = 350,
        plotLegend = FALSE)
dev.off()

g_2.1_2.2 <- generateResultsGraph(
        object = analysis.clusts2.1_2.2,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_2.1_2.2

tiff("FELLA_clusts_2.1_2.2_old_2.tiff", width = 4000, height = 4000, units = "px", pointsize = 50)
plotGraph(
        g_2.1_2.2
        #vertex.label.cex = vertex.label.cex)
)
dev.off()

tab_2.1_2.2 <- generateResultsTable(
        object = analysis.clusts2.1_2.2,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

pathways1_2_old <- tab_1_2[tab_1_2$Entry.type == "pathway", ]
pathways1_2_old[, 3] <- gsub(" -.*", pathways1_2_old[, 3], replacement = "")

pathways1.1_1.2_old <- tab_1.1_1.2[tab_1.1_1.2$Entry.type == "pathway", ]
pathways1.1_1.2_old[, 3] <- gsub(" -.*", pathways1.1_1.2_old[, 3], replacement = "")

pathways2.1_2.2_old <- tab_2.1_2.2[tab_2.1_2.2$Entry.type == "pathway", ]
pathways2.1_2.2_old[, 3] <- gsub(" -.*", pathways2.1_2.2_old[, 3], replacement = "")

write.csv(pathways1_2_old, file = "pathways1_2_old.csv")
write.csv(pathways1.1_1.2_old, file = "pathways1.1_1.2_old.csv")
write.csv(pathways2.1_2.2_old, file = "pathways2.1_2.2_old.csv")
