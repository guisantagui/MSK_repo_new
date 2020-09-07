setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/sAerous")

fluxes <- readxl::read_xlsx("exchange_fluxes.xlsx")


wilcox.test(as.numeric(unlist(fluxes[2, 3:ncol(fluxes)])),
            as.numeric(unlist(fluxes[4, 3:ncol(fluxes)])))

t.test(as.numeric(unlist(fluxes[2, 3:ncol(fluxes)])),
       as.numeric(unlist(fluxes[4, 3:ncol(fluxes)])))

boxplot(as.numeric(unlist(fluxes[2, 3:ncol(fluxes)])),
        as.numeric(unlist(fluxes[4, 3:ncol(fluxes)])))

if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)


graph <- buildGraphFromKEGGREST(
        organism = "sao",
        filter.path = c("01100", "01200", "01210", "01212", "01230")
)

graphSAureus <- graph

save(graphSAureus, file = "graphSAureus.RData")

buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = NULL,
        internalDir = TRUE,
        matrices = "none",
        normality = "diffusion",
        niter = 100)

fella.data <- loadKEGGdata(
        #databaseDir = "created__2019-03-25;meta__pae_Release_89.0_03_23_Mar_19",
        internalDir = T,
        loadMatrix = "none"
)


fella.data


# Dimethylamine, Pyruvate, ethanol, trehalose, acetoin, acetate, lactate
ndh2a_wt <- c("C00543", "C00022", "C00469", "C01083", "C00810", "C00033", "C00186") 
# Methionine, formate, dimethylamine, pyruvate, ethanol, trehalose, acetoin
ndh2b_wt <- c("C00073", "C00058", "C00543", "C00022", "C00469", "C01083", "C00810") 
# acetate, lactate, acetoin, ethanol, formate, methionine
ndh2a_ndh2b <- c("C00033", "C00186", "C00810", "C00469", "C00058", "C00073")

analysis.ndh2a_wt <- enrich(
        compounds = ndh2a_wt,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.ndh2b_wt <- enrich(
        compounds = ndh2b_wt,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.ndh2a_ndh2b <- enrich(
        compounds = ndh2a_ndh2b,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

g_ndh2a_wt <- generateResultsGraph(
        object = analysis.ndh2a_wt,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)

save(g_ndh2a_wt, file = "g_ndh2a_wt.RData")

pdf("FELLA_ndh2a_wt.pdf", width = 15, height = 15)
plotGraph(
        g_ndh2a_wt
)
dev.off()

g_ndh2b_wt <- generateResultsGraph(
        object = analysis.ndh2b_wt,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)

save(g_ndh2b_wt, file = "g_ndh2b_wt.RData")

pdf("FELLA_ndh2b_wt.pdf", width = 15, height = 15)
plotGraph(
        g_ndh2b_wt
)
dev.off()

g_ndh2a_ndh2b <- generateResultsGraph(
        object = analysis.ndh2a_ndh2b,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
save(g_ndh2a_ndh2b, file = "g_ndh2a_ndh2b.RData")

pdf("FELLA_ndh2a_ndh2b.pdf", width = 15, height = 15)
plotGraph(
        g_ndh2a_ndh2b
)
dev.off()

tab_ndh2a_wt <- generateResultsTable(
        object = analysis.ndh2a_wt,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_ndh2a_wt, file = "tab_ndh2a_wt.csv")

tab_ndh2b_wt <- generateResultsTable(
        object = analysis.ndh2b_wt,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_ndh2b_wt, file = "tab_ndh2b_wt.csv")

tab_ndh2a_ndh2b <- generateResultsTable(
        object = analysis.ndh2a_ndh2b,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_ndh2a_ndh2b, file = "tab_ndh2a_ndh2b.csv")
