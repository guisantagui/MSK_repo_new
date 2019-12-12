setwd("/Users/santamag/Desktop/GUILLEM/wrkng_dirs_clean/sequence_analysis")

########################################################################################################################################################
########################################################################################################################################################
#                                                                                                                      #################################
# This script takes as input FASTA genomes of the strains, and GFFs files with annotation information,                 #################################
# and performs a Multiple Sequence Alignment, by looking for syntenic blocks between strains for later parsing them    #################################
# to the annotation data, obtaining the ParsedSubGraphs_allStrains_named.RData object, which contains the syntenic     #################################
# blocks alineated between all the strains, each one which a gene assigned.                                            #################################
#                                                                                                                      #################################
########################################################################################################################################################
########################################################################################################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

if(!require(DECIPHER)) BiocManager::install("DECIPHER")
library(DECIPHER)
if(!require(devtools)) install.packages('devtools')
library(devtools)
if(!require(Heron)) install_github("npcooley/Heron")
library(Heron)
if(!require(igraph)) install.packages('igraph')
library(igraph)
if(!require(ape)) install.packages('ape')
library(ape)

########################################################################################################################################################
#
# Create a SQLite database with the FASTA files of the 26 strains (all except PA14, because the FASTA is 
# structured in contigs).---> ONLY RUN THIS IF IT'S THE FIRST TIME RUNNING THE SCRIPT
#
########################################################################################################################################################
db_full_genomes <- paste(paste(getwd(), "/", sep = ""), "PA_strain_full_genomes.sqlite", sep = "")    #-->Only run this if it's the first time running the script

full_seqs <- c()    #-->Only run this if it's the first time running the script
for(i in 1:(length(list.files(path = "../../Sequences/FASTA_full_genomes/my_strains")) - 1)){     #-->Only run this if it's the first time running the script
        full_seqs <- c(full_seqs, paste("../../Sequences/FASTA_full_genomes/my_strains/", 
                                            list.files(path = "../../Sequences/FASTA_full_genomes/my_strains")[i + 1],
                                            sep = ""))
}
names(full_seqs) <- gsub(".fasta", "", list.files(path = "../../Sequences/FASTA_full_genomes/my_strains/Seqs")[1:26])    #-->Only run this if it's the first time running the script



for (i in seq_along(full_seqs)) {   #-->Only run this if it's the first time running the script
        Seqs2DB(full_seqs[i], "FASTA", db_full_genomes, names(full_seqs[i]))
}
#######################################################################################################################################################
#
# Find syntenic blocks between pairs of strains
#
#######################################################################################################################################################
synteny_full_seqs <- FindSynteny(db_full_genomes)#-->In this is the first time running this script run those 2 lines
save(synteny_full_seqs, file = "synteny_full_seqs.RData")
load("synteny_full_seqs.RData")

#######################################################################################################################################################
#
# Create GeneCalls object from GFFs
#
#######################################################################################################################################################
geneCallsAllStrains4Heron <- list()
for(i in 1:length(list.files(path = "../../Sequences/FASTA_full_genomes/gff_files"))){
        thing <- readGFF(filepath = paste("../../Sequences/FASTA_full_genomes/gff_files/",
                                          list.files(path = "../../Sequences/FASTA_full_genomes/gff_files")[i],
                                          sep = ""))[, c(4, 5, 7, 10)]
        thing2 <- data.frame("Index" = rep(1, length(thing[, 1])),
                             "Start" = thing[, 1],
                             "Stop" = thing[, 2],
                             "Strand" = thing[, 3], 
                             "Annotation" = thing[, 4])
        strandVec <- c()
        for(h in 1:length(thing[, 3])){
                if(thing[h, 3] == "+"){strandVec[h] <- 0}
                if(thing[h, 3] == "-"){strandVec[h] <- 1}
        }
        thing2[, 4] <- strandVec
        thing3 <- thing2[order(thing2$Start), ]
        gff <- data.frame()
        for(k in 1:5){
                for(j in 1:dim(thing3)[1]){
                        gff[j, k] <- thing3[j, k]
                }
        }
        colnames(gff) <- colnames(thing2)
        geneCallsAllStrains4Heron[[i]] <- gff
}
names(geneCallsAllStrains4Heron) <- gsub(".gff", "", list.files(path = "../../Sequences/FASTA_full_genomes/gff_files"))
save(geneCallsAllStrains4Heron, file = "geneCallsAllStrains4Heron.RData")
load("geneCallsAllStrains4Heron.RData")

#######################################################################################################################################################
#
# Start pipeline
#
#######################################################################################################################################################
overlap_all_strains <- NucleotideOverlap(SyntenyObject = synteny_full_seqs,
                                         GeneCalls = geneCallsAllStrains4Heron,
                                         Verbose = T)
save(overlap_all_strains, file = "overlap_all_strains.RData")

load("overlap_all_strains.RData")

orthologList_allStrains <- GetOrthologSummary(OrthologsObject = overlap_all_strains,
                                              GeneCalls = geneCallsAllStrains4Heron,
                                              DBPath = "PA_strain_full_genomes.sqlite",
                                              SimilarityScores = FALSE,
                                              Type = "AAStringSet",
                                              Verbose = TRUE)
save(orthologList_allStrains, file = "orthologList_allStrains.RData")
load("orthologList_allStrains.RData")

ResolvedOrthologs_allStrains <- ResolveConflicts(SummaryObject = orthologList_allStrains,
                                                 ResolveBy = "TotalCoverage",
                                                 Verbose = TRUE)

save(ResolvedOrthologs_allStrains, file = "ResolvedOrthologs_allStrains.RData")
load("ResolvedOrthologs_allStrains.RData")

Labels_allStrains <- do.call(rbind,
                             sapply(rownames(ResolvedOrthologs_allStrains),
                                    function(x) strsplit(x,
                                                         split = " ")))
Graph_allStrains <- graph_from_edgelist(el = Labels_allStrains,
                                        directed = FALSE)

SubGraphs_allStrains <- groups(components(Graph_allStrains))

print("SubGraphs Generated")

ParsedSubGraphs_allStrains <- ParseSubGraphs(SubGraphs = SubGraphs_allStrains,
                                  DBPATH = "PA_strain_full_genomes.sqlite",
                                  GeneCalls = geneCallsAllStrains4Heron,
                                  CopyType = "Single",
                                  InputType = "Groups",
                                  OriginGraph = Graph_allStrains,
                                  Verbose = TRUE)
load("ParsedSubGraphs_allStrains.RData")

print("SubGraphs Parsed")

if(!require(zoo)) install.packages("zoo")
if(!require(sandwich)) install.packages("sandwich")
if(!require(bfast)) install.packages("bfast")
if(!require(phenopix)) install.packages("phenopix")
if(!require(Kendall)) install.packages("Kendall")
if(!require(ncdf4)) install.packages("ncdf4")
if(!require(quantreg)) install.packages("quantreg")
if(!require(spam)) install.packages("spam")
if(!require(fields)) install.packages("fields")
if(!require(greenbrown)) install.packages("greenbrown", repos="http://R-Forge.R-project.org")
library(greenbrown)

matchAnnot <- function(parsed, calls){
        library(greenbrown)
        geneNames <- c()
        pb = txtProgressBar(min = 0, max = length(parsed), initial = 0)
        for(i in seq_along(parsed)){
                names_g <- c()
                for(j in 1:ncol(parsed[[i]]$DNADist)){
                        names_g <- c(names_g, as.character(calls[[j]][as.integer(gsub(".*_", 
                                                                                      replacement = "", 
                                                                                      colnames(parsed[[i]]$DNADist)[j])), 5]))
                }
                setTxtProgressBar(pb, i)
                if(AllEqual(names_g)){
                        geneNames <- c(geneNames, names_g[1])
                }
                if(!AllEqual(names_g)){
                        genes_pasted <- names_g[1]
                        for(k in 1:(length(unique(names_g)) - 1)){
                                genes_pasted <- paste(genes_pasted, unique(names_g)[k + 1], sep = " /*/ ")
                        }
                        geneNames <- c(geneNames, genes_pasted)
                }
        }
        print("Done!")
        return(geneNames)
}

namesGenes <- matchAnnot(ParsedSubGraphs_allStrains, geneCallsAllStrains4Heron)
namesGenes <- make.unique(namesGenes)

ParsedSubGraphs_allStrains_named <- ParsedSubGraphs_allStrains
names(ParsedSubGraphs_allStrains_named) <- namesGenes

save(ParsedSubGraphs_allStrains_named, file = "ParsedSubGraphs_allStrains_named.RData")
load("ParsedSubGraphs_allStrains_named.RData")

