library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")

otu_mat<- read.csv("../data_sets/HMP_V13_OTU_counts.csv")
tax_mat<- read.csv("../data_sets/HMP_V13_taxonomy_fix.csv")
samples_df <- read_excel("../data_sets/HMP_V13_participant_data.csv")