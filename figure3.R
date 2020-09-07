setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure3")



if(!require(ropls)) BiocManager::install("ropls")
library(ropls)

if(!require(factoextra)) install.packages("factoextra")
library(factoextra)

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/swarmAnalysis/pcomp1AreaPercCircularity.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/normMetAnal/oldDataGood/ccmn_norm_mets_good_old.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/dictionary/dictionary.RData")

strains <- unique(gsub("\\_.*|(PA14).*", 
                       rownames(OPLSDAScores), 
                       rep = "\\1"))

swarmMeans <- pcomp1Means*-1

rownames(ccmn_norm_mets_good_old) <- gsub("W70322", replacement = "W70332", rownames(ccmn_norm_mets_good_old))
names(swarmMeans) <- gsub("W70322", replacement = "W70332", names(swarmMeans))

swarmMeansFilt <- swarmMeans
swarmMeansFilt <- swarmMeansFilt[match(gsub("\\_.*|(PA14).*", 
                                            rownames(ccmn_norm_mets_good_old), 
                                            rep = "\\1"), names(swarmMeans))]
length(swarmMeansFilt)

ccmn_norm_mets_good_old$swarmData <- swarmMeansFilt


# Add binarized swarming data
swarmBin <- rep(NA, nrow(ccmn_norm_mets_good_old))
swarmBin[which(swarmMeansFilt > 0.19)] <- 1
swarmBin[is.na(swarmBin)] <- 0
swarmBin <- as.factor(swarmBin)
ccmn_norm_mets_good_old$swarmDatBin <- swarmBin

# Add rhamnolipid production data

glycRhamn <- c(2, 2, 2, 2, 2, 0, 2, 0, 2, 1, 0, 1, 1, 2, 0, 2, 2, 2, 1, 1, 2, 0, 2, 0, 1, 2, 2, 2)
glycRhamn <- c(1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1)

names(glycRhamn) <- strains

glycRhamn <- glycRhamn[match(gsub("\\_.*|(PA14).*", 
                                  rownames(OPLSDAScores), 
                                  rep = "\\1"), 
                             names(glycRhamn))]

ccmn_norm_mets_good_old$rhamn <- glycRhamn

# Do PCA with swarming color code (figure S3)

PCAMetSwarms <- prcomp(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old)-3)])

pdf("PCAMetSwarms.pdf", height = 12, width = 12)
fviz_pca_ind(PCAMetSwarms, 
             col.ind = ccmn_norm_mets_good_old$swarmData,
             gradient.cols = c("#003CFF", "#66FF00", "#FF0000"),
             legend.title = "Swarming Score")
dev.off()


swarmData

# OPLS-DA

# Qualitative and with partition

OPLSDA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 3)],
               ccmn_norm_mets_good_old$swarmDatBin, 
               predI = 1, 
               orthoI = NA,
               subset = "odd"
)

trainVi <- getSubsetVi(OPLSDA)
table(ccmn_norm_mets_good_old$swarmDatBin[trainVi], fitted(OPLSDA))

confusion <- table(ccmn_norm_mets_good_old$swarmDatBin[-trainVi],
                   predict(OPLSDA, ccmn_norm_mets_good_old[-trainVi, 1:(ncol(ccmn_norm_mets_good_old)-3)]))

(confusion[1, 1] + confusion[2, 2])/sum(confusion)*100 # --> 95% of accuracy

# Quantitative without partition
pdf(file = "OPLSDA_swarm.pdf", width = 10, height = 10)
OPLSDA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 3)],
               ccmn_norm_mets_good_old$swarmData, 
               predI = 1, 
               orthoI = NA
)
dev.off()

# Obtain barplot

load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/figure2/solvedAmbigMets.RData")

solvedAmbigMets <- solvedAmbigMets[-4]

OPLSDALoads <- getLoadingMN(OPLSDA)
rownames(OPLSDALoads) <- dictionary$Consensus[match(make.names(rownames(OPLSDALoads)), make.names(dictionary$`Old Data Names`))]
OPLSDALoads <- OPLSDALoads[!is.na(rownames(OPLSDALoads)), ]

# Change names of solved ambig mets. 
names(OPLSDALoads)[match(names(solvedAmbigMets), names(OPLSDALoads))] <- solvedAmbigMets


# Remove ambiguous mets
OPLSDALoads <- OPLSDALoads[-grep("_?", names(OPLSDALoads), fixed = T)]
defNames <- dictionary$definitiveNames[match(names(OPLSDALoads), dictionary$Consensus)]

names(OPLSDALoads)[!is.na(defNames)] <- defNames[!is.na(defNames)]

OPLSDALoads <- cbind.data.frame(names(OPLSDALoads), OPLSDALoads)
colnames(OPLSDALoads) <- c("metabolite", "loading")

ggplot(data = OPLSDALoads,
       aes(x = reorder(metabolite, loading), y = loading)) +
        labs(x = "Metabolites", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity") +
        coord_flip()
ggsave(filename = "OPLSDALoadsBarplot.pdf", limitsize = F)

# Let's do a fancier version of the OPLS-DA plot

OPLSDAScores <- cbind(getScoreMN(OPLSDA), getScoreMN(OPLSDA, orthoL = T)[, 1])

colnames(OPLSDAScores)[2] <- "o1"

OPLSDAScores <- as.data.frame(OPLSDAScores)

OPLSDAScores$swarmDat <- swarmMeansFilt


#low = "#003CFF", 
#mid = "#66FF00",
#high = "#FF0000",
#space = "Lab")

# get means of scores for plotting purposes

OPLSDAScoresMeans <- c()
for(i in seq_along(strains)){
        strMean <- apply(OPLSDAScores[gsub("\\_.*|(PA14).*", 
                                           rownames(OPLSDAScores), 
                                           rep = "\\1") %in% strains[i], ], 2, mean)
        OPLSDAScoresMeans <- rbind(OPLSDAScoresMeans, strMean)
}
rownames(OPLSDAScoresMeans) <- strains
OPLSDAScoresMeans <- as.data.frame(OPLSDAScoresMeans)


OPLSDAScorePlot<-ggplot(OPLSDAScores, aes(x=p1, y=o1, color=swarmDat)) + geom_point()


mid<-(max(OPLSDAScores$swarmDat) - min(OPLSDAScores$swarmDat))/2 + min(OPLSDAScores$swarmDat)
OPLSDAScorePlot+scale_color_gradient2(midpoint=mid, 
                          low="#003CFF",
                          mid="#66FF00",
                          high="#FF0000",
                          space ="Lab",
                          limits = c(min(OPLSDAScores$swarmDat),
                                     max(OPLSDAScores$swarmDat))) +
        geom_text(label=rownames(OPLSDAScores), nudge_y = 0.3) + 
        labs(col="Swarming Score") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")

ggsave("OPLSDAFancy_swarm.pdf", height = 10, width = 10)

OPLSDAScoreMeansPlot<-ggplot(OPLSDAScoresMeans, aes(x=p1, y=o1, color=swarmDat)) + 
        geom_point() +
        scale_color_gradient2(midpoint=mid, 
                              low="#003CFF",
                              mid="#66FF00",
                              high="#FF0000",
                              space ="Lab",
                              limits = c(min(OPLSDAScores$swarmDat),
                                         max(OPLSDAScores$swarmDat))) +
        geom_text(label=rownames(OPLSDAScoresMeans), nudge_y = 0.3) + 
        labs(col="Swarming Score") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("t1 (6%)") +
        ylab("to1")

OPLSDAScoreMeansPlot

ggsave("OPLSDAFancy_swarm_means.pdf", height = 10, width = 10)

OPLSDA

# OPLS-DA using rhamnolipid production as response variable.
pdf(file = "OPLSDA_rhamn.pdf", width = 10, height = 10)
rhamnOPLSDA <- opls(ccmn_norm_mets_good_old[, 1:(ncol(ccmn_norm_mets_good_old) - 3)],
                    ccmn_norm_mets_good_old$rhamn, 
                    predI = 1, 
                    orthoI = NA)
dev.off()

rhamnOPLSDALoads <- getLoadingMN(rhamnOPLSDA)

rownames(rhamnOPLSDALoads) <- dictionary$Consensus[match(make.names(rownames(rhamnOPLSDALoads)), 
                                                         make.names(dictionary$`Old Data Names`))]
rhamnOPLSDALoads <- rhamnOPLSDALoads[!is.na(rownames(rhamnOPLSDALoads)), ]

# Change names of solved ambig mets. 
names(rhamnOPLSDALoads)[match(names(solvedAmbigMets), names(rhamnOPLSDALoads))] <- solvedAmbigMets


# Remove ambiguous mets
rhamnOPLSDALoads <- rhamnOPLSDALoads[-grep("_?", names(rhamnOPLSDALoads), fixed = T)]
defNames <- dictionary$definitiveNames[match(names(rhamnOPLSDALoads), dictionary$Consensus)]

names(rhamnOPLSDALoads)[!is.na(defNames)] <- defNames[!is.na(defNames)]

rhamnOPLSDALoads <- cbind.data.frame(names(rhamnOPLSDALoads), rhamnOPLSDALoads)
colnames(rhamnOPLSDALoads) <- c("metabolite", "loading")

ggplot(data = rhamnOPLSDALoads,
       aes(x = reorder(metabolite, loading), y = loading)) +
        labs(x = "Metabolites", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity") +
        coord_flip()
ggsave(filename = "rhamnOPLSDALoadsBarplot.pdf", limitsize = F)




rhamnOPLSDAScores <- cbind(getScoreMN(rhamnOPLSDA), getScoreMN(rhamnOPLSDA, orthoL = T)[, 1])

colnames(rhamnOPLSDAScores)[2] <- "o1"

rhamnOPLSDAScores <- as.data.frame(rhamnOPLSDAScores)

rhamnOPLSDAScores$rhamn <- ccmn_norm_mets_good_old$rhamn


#low = "#003CFF", 
#mid = "#66FF00",
#high = "#FF0000",
#space = "Lab")

# get means of scores for plotting purposes

rhamnOPLSDAScoresMeans <- c()
for(i in seq_along(strains)){
        strMean <- apply(rhamnOPLSDAScores[gsub("\\_.*|(PA14).*", 
                                                rownames(rhamnOPLSDAScores), 
                                                rep = "\\1") %in% strains[i], ], 2, mean)
        rhamnOPLSDAScoresMeans <- rbind(rhamnOPLSDAScoresMeans, strMean)
}
rownames(rhamnOPLSDAScoresMeans) <- strains
rhamnOPLSDAScoresMeans <- as.data.frame(rhamnOPLSDAScoresMeans)


rhamnOPLSDAScorePlot <- ggplot(rhamnOPLSDAScores, aes(x=p1, y=o1, color=rhamn)) + geom_point()


rhamnOPLSDAScorePlot+scale_color_gradient2(midpoint=0.5, 
                                      low="blue",
                                      mid="white",
                                      high="red",
                                      space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScores), nudge_y = 0.3) + 
        labs(col="Rhamnolipid Production") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed")

ggsave("OPLSDAFancy_rhamn.pdf", height = 10, width = 10)

rhamnOPLSDAScoreMeansPlot <- ggplot(rhamnOPLSDAScoresMeans, aes(x=p1, y=o1, color=rhamn)) + 
        geom_point() +
        scale_color_gradient2(midpoint=0.5, 
                              low="blue",
                              mid="white",
                              high="red",
                              space ="Lab") +
        geom_text(label=rownames(rhamnOPLSDAScoresMeans), nudge_y = 0.3) + 
        labs(col="Rhamnolipid Production") + theme_minimal() +
        geom_vline(xintercept = 0, linetype = "dashed") +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        xlab("t1 (3%)") +
        ylab("to1")

rhamnOPLSDAScoreMeansPlot

ggsave("OPLSDAFancy_rhamn_means.pdf", height = 10, width = 10)

# Do FELLA for mets with highest loading:
if(!require(FELLA)) BiocManager::install("FELLA")
library(FELLA)
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)
if(!require(igraph)) install.packages("igraph")
library(igraph)
if(!require(magrittr)) install.packages("magrittr")
library(magrittr)

# Get top swarming predictors compounds
n <- 8
extremeLoadComps <- rbind.data.frame(head(OPLSDALoads[order(OPLSDALoads$loading), ], n),
                                     tail(OPLSDALoads[order(OPLSDALoads$loading), ], n))

KEGGIDs <- dictionary$`KEGG IDs`[match(rownames(extremeLoadComps), dictionary$definitiveNames)]

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

loadKEGGdata()

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

analysis.rfeRes<- enrich(
        compounds = KEGGIDs,
        data = fella.data,
        method = "diffusion",
        approx = "normality")

analysis.rfeRes %>%
        getInput %>%
        getName(data = fella.data)
getExcluded(analysis.rfeRes)

g_rfeSwarm <- generateResultsGraph(
        object = analysis.rfeRes,
        method = "diffusion",
        nlimit = 350,
        data = fella.data)
g_rfeSwarm

pdf("FELLA_clusts_rfeRes_swarm.pdf", width = 15, height = 15)
plotGraph(
        g_rfeSwarm
        #vertex.label.cex = vertex.label.cex)
)

dev.off()

tab_rfeSwarm <- generateResultsTable(
        object = analysis.rfeRes,
        data = fella.data,
        method = "diffusion",
        nlimit = 500)

write.csv(tab_rfeSwarm, file = "tab_rfeSwarm.csv")
save(tab_rfeSwarm, file = "tab_rfeSwarm.RData")

KEGGIDs_all <- dictionary$`KEGG IDs`[match(rownames(OPLSDALoads), dictionary$definitiveNames)]
KEGGIDs_all <- KEGGIDs_all[order(OPLSDALoads$loading)]
KEGGIDs_all <- KEGGIDs_all[!is.na(KEGGIDs_all)]

pathsPerMet <- list()
for(i in seq_along(KEGGIDs_all)){
        print(i)
        paths <- keggLink("pathway", 
                          KEGGIDs_all[i])[gsub("path:map", 
                                               "pae", 
                                               keggLink("pathway", 
                                                        KEGGIDs_all[i])) %in% tab_rfeSwarm$KEGG.id]
        pathsPerMet[[i]] <- paths
}
names(pathsPerMet) <- KEGGIDs_all
pathsPerMet

allFELLAPaths <- unique(unlist(pathsPerMet))

metsInPaths <- list()

for(i in seq_along(allFELLAPaths)){
        allMetsinPath <- keggLink("compound", allFELLAPaths[i])
        allMetsinPath <- gsub("cpd:", "", allMetsinPath)
        mets <- allMetsinPath[allMetsinPath %in% KEGGIDs_all]
        metsInPaths[[i]] <- mets
}
names(metsInPaths) <- allFELLAPaths
metsInPaths$anyRelatedPath <- KEGGIDs_all[!KEGGIDs_all %in% unique(unlist(metsInPaths))]

loadsPaths <- data.frame()
for(i in seq_along(names(metsInPaths))){
        mets <- metsInPaths[[i]]
        metsDefNames <- dictionary$definitiveNames[match(mets, dictionary$`KEGG IDs`)]
        loadsSubMat <- OPLSDALoads[match(metsDefNames, OPLSDALoads$metabolite), ]
        loadsSubMat <- cbind.data.frame(loadsSubMat, rep(names(metsInPaths)[i], nrow(loadsSubMat)))
        loadsPaths <- rbind.data.frame(loadsPaths, loadsSubMat)
}

colnames(loadsPaths)[3] <- "Pathway"

allPathways <- keggList("pathway")

loadsPaths$Pathway <- allPathways[match(loadsPaths$Pathway, names(allPathways))]
loadsPaths$Pathway[is.na(loadsPaths$Pathway)] <- "Any related pathway"

loadsPaths$metabolite <- make.unique(as.character(loadsPaths$metabolite))
loadsPaths$metabolite <- factor(loadsPaths$metabolite, levels = loadsPaths$metabolite)


ggplot(data = loadsPaths,
       aes(x = metabolite, y = loading, fill = Pathway)) +
        labs(x = "Metabolite", y = "OPLS-DA Loading") +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip()

ggsave(filename = "loadsOPLSDAPaths.pdf", height = 20, width = 12)

specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

# Grouped
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
        geom_bar(position="dodge", stat="identity")

