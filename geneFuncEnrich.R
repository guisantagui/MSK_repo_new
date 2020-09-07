setwd("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/rnaSeqAnalysis")

if(!requireNamespace("BiocManager", quietly=TRUE))
        install.packages("BiocManager")
if(!require(FGNet)) BiocManager::install("FGNet")
library(FGNet)

if(!require(RGtk2)) BiocManager::install("RGtk2")
if(!require(RCurl)) BiocManager::install("RCurl")
#if(!require(RDAVIDWebService)) BiocManager::install("RDAVIDWebService")
if(!require(gage)) BiocManager::install("gage")
if(!require(topGO)) BiocManager::install("topGO")
if(!require(GO.db)) BiocManager::install("GO.db")
if(!require(reactome.db)) BiocManager::install("reactome.db")
if(!require(org.Sc.sgd.db)) BiocManager::install("org.Sc.sgd.db")
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")
library(clusterProfiler)
if(!require(AnnotationHub)) BiocManager::install("AnnotationHub")
library(AnnotationHub)
if(!require(MeSH.Pae.PAO1.eg.db)) BiocManager::install("MeSH.Pae.PAO1.eg.db")
if(!require(KEGGREST)) install.packages("KEGGREST")
library(KEGGREST)

load("signDESeqResultsNoHypProt.RData")
load("C:/Users/Guillem/Documents/PhD/comput/wrkng_dirs_clean/genePresAbs/dictEnzymes.RData")

read.csv("C:/Users/Guillem/Documents/PhD/comput/eggNOGMapps/job_MM_lqj2rvp1_annotations.tsv")
prokkAnot <- read.csv("C:/Users/Guillem/Documents/PhD/comput/prokkAnnotStuff/untitled_folder/gene_presence_absence.csv")


geneExpr <- rownames(signDESeqResultsNoHypProt)
geneExpr <- setNames(c(rep(-1, max(which(signDESeqResultsNoHypProt$log2FoldChange < 0))),
                       rep(1, length(geneExpr) - max(which(signDESeqResultsNoHypProt$log2FoldChange < 0)))), 
                     geneExpr)

hub <- AnnotationHub()


query(hub, "coli")

query(hub, c("OrgDb", "Pseudomonas"))

PA <- hub[["AH10565"]]

PP <- hub[["AH77379"]]

EColi <- hub[["AH10493"]]

Cgriseus <- hub[["AH48061"]]

s <- sample(keys(org.Sc.sgd.db, keytype = "ENTREZID"), 100, )

underGenes <- gsub("\\_.*", "", names(geneExpr[geneExpr == -1]))

enrichGO(underGenes, 
         OrgDb = PA, 
         pvalueCutoff = 1,
         qvalueCutoff = 1, 
         ont = "ALL", 
         readable = T,
         keyType = "PSEUDOMONAS_AERUGINOSA")

enrichGO(s, 
         OrgDb = org.Sc.sgd.db, 
         pvalueCutoff = 1,
         qvalueCutoff = 1)


query(AnnotationHub(), "Pseudomonas aeruginosa")
cls <- columns(MeSH.Pae.PAO1.eg.db)
cls
kts <- keytypes(MeSH.Pae.PAO1.eg.db)
kt <- kts[2]
kts
ks <- head(keys(MeSH.Pae.PAO1.eg.db, keytype=kts[2]))
ks
res <- select(MeSH.Pae.PAO1.eg.db, keys=ks, columns=cls, keytype=kt)
head(res)
dbconn(MeSH.Pae.PAO1.eg.db)
dbfile(MeSH.Pae.PAO1.eg.db)
dbschema(MeSH.Pae.PAO1.eg.db)
dbInfo(MeSH.Pae.PAO1.eg.db)
species(MeSH.Pae.PAO1.eg.db)

colsPA <- columns(PA)
colsPA
ktsPA <- keytypes(PA)
ktPA <- ktsPA[2]
ksPA <- head(keys(PA, keytype = ktsPA[27]))

select(PA, keys = ksPA, columns = c("PAN_TROGLODYTES", "PSEUDOMONAS_AERUGINOSA"), keytype = "PSEUDOMONAS_AERUGINOSA")


ail <- query(hub, c('inparanoid8', 'ailuropoda'))

library(org.Sc.sgd.db)
geneLabels <- unlist(as.list(org.Sc.sgdGENENAME))

?setNames

c(rep(1, 19), rep(-1, 5))

which(signDESeqResultsNoHypProt$log2FoldChange < 0)


TCAEnzs <- gsub("ec:", "", keggLink("enzyme", "map00020"))


genesInTCACyc <- dictEnzymes$Gene[dictEnzymes$ECnums %in% TCAEnzs]

signDESeqResultsNoHypProt[which(rownames(signDESeqResultsNoHypProt) %in% genesInTCACyc),]
