library("phyloseq")
library("microbiome")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

getwd() # checking file path

otu_mat<- read.csv("data_sets/HMP_V13_OTU_counts.csv", header = TRUE)
tax_mat<- read.csv("data_sets/HMP_V13_taxonomy_fix.csv", header = TRUE)
samples_df <- read.csv("data_sets/HMP_V13_participant_data.csv", header = TRUE)

otu_mat_2 <- otu_mat %>%
  tibble::column_to_rownames("PSN") %>%
  t() %>% # transpose for formatting
  as.data.frame()

tax_mat_2 <- tax_mat %>%
  as.data.frame() %>%
  mutate(OTU = OTU_ID) %>%
  tibble::column_to_rownames("OTU_ID")

rownames(samples_df) <- NULL
samples_df_2 <- samples_df %>%
  tibble::column_to_rownames("PSN")

otu_mat_3 <- as.matrix(otu_mat_2)
tax_mat_3 <- as.matrix(tax_mat_2)

OTU = otu_table(otu_mat_3, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat_3)
samples = sample_data(samples_df_2)

carbom <- phyloseq(OTU, TAX, samples)
carbom

carbom.compositional <- microbiome::transform(carbom, "compositional")

list.taxa.passing <- core_members(carbom.compositional, detection = 0/100, prevalence = 10/100)

filtered_physeq <- subset_taxa(carbom.compositional, OTU %in% list.taxa.passing)


