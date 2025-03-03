library("phyloseq")
library("microbiome")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names

getwd() # checking file path

# Load data
otu_mat <- read.csv("data_sets/HMP_V13_OTU_counts.csv", header = TRUE)
tax_mat <- read.csv("data_sets/HMP_V13_taxonomy_fix.csv", header = TRUE)
samples_df <- read.csv("data_sets/HMP_V13_participant_data.csv", header = TRUE)

# Format OTU matrix
otu_mat_2 <- otu_mat %>%
  tibble::column_to_rownames("PSN") %>%
  t() %>%
  as.data.frame()

# Format taxonomy matrix
tax_mat_2 <- tax_mat %>%
  as.data.frame() %>%
  mutate(OTU = OTU_ID) %>%
  tibble::column_to_rownames("OTU_ID")

# Format sample data
rownames(samples_df) <- NULL
samples_df_2 <- samples_df %>%
  tibble::column_to_rownames("PSN")

# Convert to phyloseq objects
otu_mat_3 <- as.matrix(otu_mat_2)
tax_mat_3 <- as.matrix(tax_mat_2)

OTU = otu_table(otu_mat_3, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat_3)
samples = sample_data(samples_df_2)

carbom <- phyloseq(OTU, TAX, samples)

# Define preprocessing combinations
prevalence_options <- c(TRUE, FALSE)
transformations <- c("none", "log10", "compositional", "rarefy")

for (prev in prevalence_options) {
  if (prev) {
    carbom.compositional <- microbiome::transform(carbom, "compositional")
    list.taxa.passing <- core_members(carbom.compositional, detection = 0/100, prevalence = 0/100)
    filtered_physeq <- subset_taxa(carbom.compositional, OTU %in% list.taxa.passing)
  } else {
    filtered_physeq <- carbom
  }
  
  for (trans in transformations) {
    if (trans == "none") {
      transformed_physeq <- filtered_physeq
    } else if (trans == "rarefy") {
      # Apply rarefaction when 'rarefy' is selected
      transformed_physeq <- rarefy_even_depth(filtered_physeq)
    } else {
      # Apply other transformations
      transformed_physeq <- microbiome::transform(filtered_physeq, trans)
    }
    
    # Create file name based on combinations (00, 10, 01, 11, etc.)
    prevalence_code <- ifelse(prev, "1", "0")
    transformation_code <- ifelse(trans == "none", "0", ifelse(trans == "log10", "1", ifelse(trans == "compositional", "2", ifelse(trans == "rarefy", "3", "4"))))
    
    filename <- paste0("data_sets/preprocessed_table_", prevalence_code, transformation_code, ".csv")
    otu_table_transformed <- t(abundances(transformed_physeq))
    write.csv(otu_table_transformed, filename, row.names = TRUE)
  }
}