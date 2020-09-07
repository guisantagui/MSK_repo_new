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
db_full_genomes_new <- paste(paste(getwd(), "/", sep = ""), "PA_strain_full_genomes_new.sqlite", sep = "")    #-->Only run this if it's the first time running the script

full_seqs_new <- c()    #-->Only run this if it's the first time running the script
for(i in 1:(length(list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive")))){     #-->Only run this if it's the first time running the script
        full_seqs_new <- c(full_seqs_new, paste("/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive/", 
                                                list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive")[i],
                                                sep = ""))
}
names(full_seqs_new) <- gsub(".fasta|.fna", "", list.files(path = "/Users/santamag/Desktop/GUILLEM/Sequences/newSeqs/Archive"))    #-->Only run this if it's the first time running the script



for (i in seq_along(full_seqs_new)) {   #-->Only run this if it's the first time running the script
        Seqs2DB(full_seqs_new[i], "FASTA", db_full_genomes_new, names(full_seqs_new[i]))
}
#######################################################################################################################################################
#
# Find syntenic blocks between pairs of strains
#
#######################################################################################################################################################
synteny_full_seqs_new <- FindSynteny(db_full_genomes_new)#-->In this is the first time running this script run those 2 lines
save(synteny_full_seqs_new, file = "synteny_full_seqs_new.RData")
load("synteny_full_seqs_new.RData")

#######################################################################################################################################################
#
# Create GeneCalls object from GFFs
#
#######################################################################################################################################################
prokkAnot <- read.csv("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/untitled_folder/gene_presence_absence.csv")
namesCalls <- make.names(gsub(".gff", "", list.files(path = "/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3")))
geneCallsAllStrains4HeronProkka <- list()
for(i in 1:length(list.files(path = "/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3"))){
        print(i)
        thing <- readGFF(filepath = paste("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3/",
                                          list.files(path = "/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3")[i],
                                          sep = ""))[, c(4, 5, 7, 9)]
        thing2 <- data.frame("Index" = rep(1, length(thing[, 1])),
                             "Start" = thing[, 1],
                             "Stop" = thing[, 2],
                             "Strand" = thing[, 3], 
                             "Annotation" = thing[, 4])
        strandVec <- gsub("+", replacement = "0", thing[, 3], fixed = T)
        strandVec <- gsub("-", replacement = "1", strandVec, fixed = T)
        # many crispr regions are coded with "*" in the GFF file
        strandVec <- gsub("*", replacement = "2", strandVec, fixed = T)
        strandVec <- as.numeric(strandVec)
        print(strandVec)
        #for(h in 1:length(thing[, 3])){
        #        if(thing[h, 3] == "+"){strandVec[h] <- 0}
        #        if(thing[h, 3] == "-"){strandVec[h] <- 1}
        #}
        thing2[, 4] <- strandVec
        thing3 <- thing2[order(thing2$Start), ]
        gff <- thing3
        #gff <- data.frame()
        #for(k in 1:5){
        #        for(j in 1:dim(thing3)[1]){
        #                gff[j, k] <- thing3[j, k]
        #        }
        #}
        #colnames(gff) <- colnames(thing2)
        geneNames <- prokkAnot$Gene[match(gff$Annotation, prokkAnot[[namesCalls[i]]])]
        #geneNames <- geneNames[!is.na(geneNames)]
        gff$Annotation <- geneNames
        gff <- gff[!gff$Strand == 2, ]
        gff <- gff[!is.na(gff$Annotation), ]
        rownames(gff) <- 1:nrow(gff)
        geneCallsAllStrains4HeronProkka[[i]] <- gff
}
names(geneCallsAllStrains4HeronProkka) <- namesCalls
save(geneCallsAllStrains4HeronProkka, file = "geneCallsAllStrains4HeronProkka.RData")
load("geneCallsAllStrains4HeronProkka.RData")


numWeird <- which(unlist(lapply(geneCallsAllStrains4HeronProkka, function(x) anyNA(x$Strand))))
whichWeird <- lapply(geneCallsAllStrains4HeronProkka, function(x) which(is.na(x$Strand)))[unlist(lapply(geneCallsAllStrains4HeronProkka, function(x) anyNA(x$Strand)))]

weirds <- list()
for(i in seq_along(numWeird)){
        gff1 <- readGFF(filepath = paste("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3/",
                                         list.files(path = "/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3")[numWeird[i]], sep = ""))
        weirds[[i]] <- gff1[whichWeird[[i]], ]
        
}
gff1 <- readGFF(filepath = paste("/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3/",
                         list.files(path = "/Users/santamag/Desktop/GUILLEM/prokkAnnotStuff/gff3")[29], sep = ""))

gff1[c(2771, 5300, 5313), ] 

#######################################################################################################################################################
#
# Start pipeline
#
#######################################################################################################################################################
overlap_all_strains <- NucleotideOverlap(SyntenyObject = synteny_full_seqs_new,
                                         GeneCalls = geneCallsAllStrains4HeronProkka,
                                         Verbose = T)
overlap_all_strains_new <- overlap_all_strains
save(overlap_all_strains_new, file = "overlap_all_strains_new.RData")

load("overlap_all_strains_new.RData")

orthologList_allStrains <- GetOrthologSummary(OrthologsObject = overlap_all_strains_new,
                                              GeneCalls = geneCallsAllStrains4HeronProkka,
                                              DBPath = "PA_strain_full_genomes_new.sqlite",
                                              SimilarityScores = FALSE,
                                              Type = "AAStringSet",
                                              Verbose = TRUE)
orthologList_allStrains_new <- orthologList_allStrains
save(orthologList_allStrains_new, file = "orthologList_allStrains_new.RData")
load("orthologList_allStrains_new.RData")

ResolvedOrthologs_allStrains_new <- ResolveConflicts(SummaryObject = orthologList_allStrains_new,
                                                     ResolveBy = "TotalCoverage",
                                                     Verbose = TRUE)

save(ResolvedOrthologs_allStrains_new, file = "ResolvedOrthologs_allStrains_new.RData")
load("ResolvedOrthologs_allStrains_new.RData")

Labels_allStrains <- do.call(rbind,
                             sapply(rownames(ResolvedOrthologs_allStrains_new),
                                    function(x) strsplit(x,
                                                         split = " ")))
Graph_allStrains <- graph_from_edgelist(el = Labels_allStrains,
                                        directed = FALSE)

SubGraphs_allStrains <- groups(components(Graph_allStrains))

print("SubGraphs Generated")

ParsedSubGraphs_allStrains_new <- ParseSubGraphs(SubGraphs = SubGraphs_allStrains,
                                  DBPATH = "PA_strain_full_genomes_new.sqlite",
                                  GeneCalls = geneCallsAllStrains4HeronProkka,
                                  CopyType = "Single",
                                  InputType = "Groups",
                                  OriginGraph = Graph_allStrains,
                                  Verbose = TRUE)
save(ParsedSubGraphs_allStrains_new, file = "ParsedSubGraphs_allStrains_new.RData")
load("ParsedSubGraphs_allStrains_new.RData")

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
                strsWgene <- as.integer(gsub("_.*", replacement = "", colnames(parsed[[i]]$DNADist)))
                for(j in 1:length(strsWgene)){
                        names_g <- c(names_g, as.character(calls[[strsWgene[j]]][as.integer(gsub(".*_", 
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



namesGenes <- matchAnnot(ParsedSubGraphs_allStrains_new, geneCallsAllStrains4HeronProkka)
namesGenes <- make.unique(namesGenes)

ParsedSubGraphs_allStrains_named_new <- ParsedSubGraphs_allStrains_new
names(ParsedSubGraphs_allStrains_named_new) <- namesGenes

save(ParsedSubGraphs_allStrains_named_new, file = "ParsedSubGraphs_allStrains_named_new.RData")
load("ParsedSubGraphs_allStrains_named_new.RData")

