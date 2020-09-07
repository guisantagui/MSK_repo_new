setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/putidaPA14_overlap")
genes_iJN1411 <- read.csv("genes_in_iJN1411.csv")
if(!require(readxl)) install.packages("readxl")
library(readxl)
genomeCompPA14_put <- read_xls("genome_comparison.xls")

colnames(genomeCompPA14_put) <- genomeCompPA14_put[1, ]
genomeCompPA14_put <- genomeCompPA14_put[-1, ]

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

if(!require(SBMLR)) BiocManager::install("SBMLR")
library(SBMLR)

PA14Model <- readSBML("../../PATRIC_PA_models/PA14.sbml")
genesPA14Mod <- PA14Model$notes[grep("GENE_ASSOCIATION", PA14Model$notes)]
genesPA14Mod <- gsub("GENE_ASSOCIATION:", "", genesPA14Mod)
genesPA14Mod <- gsub("Unknown", "", genesPA14Mod)
genesPA14Mod <- gsub("and", "", genesPA14Mod)
genesPA14Mod <- gsub("or", "", genesPA14Mod)
genesPA14Mod <- gsub("(", "", genesPA14Mod, fixed = T)
genesPA14Mod <- gsub(")", "", genesPA14Mod, fixed = T)

gPA14Mod <- c()
for(i in seq_along(genesPA14Mod)){
        g <- strsplit(genesPA14Mod[i], " ")
        gPA14Mod <- c(gPA14Mod, g)
}
gPA14Mod <- unlist(gPA14Mod)
gPA14Mod <- gPA14Mod[gPA14Mod != ""]
gPA14Mod <- unique(gPA14Mod)



subComparison <- genomeCompPA14_put[genomeCompPA14_put$comp_genome_1_patric_id %in% gPA14Mod, ]

inBoth <- subComparison$ref_genome_locus_tag[subComparison$comp_genome_1_hit == "bi (<->)"]
inBoth <- inBoth[!is.na(inBoth)] 

subCompBi <- subComparison[subComparison$comp_genome_1_hit == "bi (<->)", ]

gPA14Mod_subst <- gPA14Mod
for(i in seq_along(gPA14Mod_subst)){
        if(gPA14Mod_subst[i] %in% subCompBi$comp_genome_1_patric_id){
                eqName <-  subCompBi$ref_genome_locus_tag[subCompBi$comp_genome_1_patric_id == gPA14Mod_subst[i]]
                if(!is.na(eqName)){
                        gPA14Mod_subst[i] <- eqName
                }
        }
}
 
length(grep("PP", gPA14Mod_subst))           # ---> 759 genes out of 976 genes that appear in PA14 model have 
                                             #      a match to P. putida genome according to proteome comparison
                                             #      (77.8% of the genes in PA14 model).
sum(genes_iJN1411[, 1] %in% gPA14Mod_subst)  #      603 of those 759 genes appear in iJN1411 model. That means 
                                             #      that the 42.7% of the genes that appear in iJN1411 also 
                                             #      appear in PA14 automaticcally generated model. 
                                   
sum(genes_iJN1411[, 1] %in% inBoth)/nrow(genes_iJN1411) # ---> 42.7% of overlap
sum(genes_iJN1411[, 1] %in% gPA14Mod_subst)/nrow(genes_iJN1411) # ---> 42.7% of overlap: SAME

if(!require(graph)) BiocManager::install("graph")
if(!require(RBGL)) BiocManager::install("RBGL")
if(!require(Vennerable)) install.packages("Vennerable", repos="http://R-Forge.R-project.org")
library(Vennerable)

compGenes <- list(PA14 = gPA14Mod_subst, iJN1411 = genes_iJN1411[, 1])
vennD <- Venn(Sets = compGenes)
tiff(filename = "PA14_iJN1411_overlap.tiff", height = 1300, width = 1800, res = 300)
plot(vennD, doWeights = TRUE, type = "circles")
dev.off()
