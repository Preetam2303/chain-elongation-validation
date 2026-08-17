# ==============================================================================
# 1. WORKSPACE SETUP & DEPENDENCIES
# ==============================================================================
# Install missing packages if necessary
#install.packages(c("dplyr", "tibble", "readxl"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("dada2")

library(dada2)
library(dplyr)
library(tibble)

# Set working directory directly to the Grand Merge folder
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
workdir <- "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"
setwd(workdir)
message("Working directory set to: ", getwd())

# ==============================================================================
# 2. DOWNLOAD SILVA v138.1 REFERENCE DATABASES
# ==============================================================================
# Increase timeout limit for large downloads on Windows
options(timeout = max(600, getOption("timeout")))

# URLs for DADA2 formatted SILVA 138.1 databases (Zenodo canonical links)
silva_train_url <- "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
silva_species_url <- "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz"

# File destinations in your working directory
train_file <- "silva_nr99_v138.1_train_set.fa.gz"
species_file <- "silva_species_assignment_v138.1.fa.gz"

# Download files (mode = "wb" is critical for Windows to prevent archive corruption)
if(!file.exists(train_file)) {
  message("Downloading SILVA training set...")
  download.file(silva_train_url, destfile = train_file, mode = "wb")
}

if(!file.exists(species_file)) {
  message("Downloading SILVA species assignment...")
  download.file(silva_species_url, destfile = species_file, mode = "wb")
}

# ==============================================================================
# 3. LOAD DATA AND ASSIGN TAXONOMY
# ==============================================================================
message("Loading master matrix from local NVMe SSD...")
seqtab.nochim <- readRDS("BIOTWIN_Global_Master_Matrix_NoChim.rds")

message("Assigning taxonomy...")
# Taking full advantage of the Intel Ultra 5 and 32GB RAM. 
# RcppParallel will utilize multiple cores here to process the 18,870 ASVs.
taxa <- assignTaxonomy(seqtab.nochim, train_file, multithread = TRUE)

message("Adding species level assignments...")
taxa <- addSpecies(taxa, species_file)

# ==============================================================================
# 4. EXTRACT AND FORMAT TABLES FOR METADATA MERGE
# ==============================================================================
# A. Create a mapping of DNA sequences to clean ASV IDs
raw_sequences <- colnames(seqtab.nochim)
clean_asv_ids <- paste0("ASV_", seq_along(raw_sequences))

# B. Format the ASV Count Matrix
asv_counts <- as.data.frame(seqtab.nochim)
# Rename the columns from raw DNA to clean ASV IDs
colnames(asv_counts) <- clean_asv_ids
# Move the rownames (currently sample names) into a proper column to match your Excel sheet
asv_counts <- asv_counts %>%
  rownames_to_column(var = "Sample_ID")

# C. Format the Taxonomy Table
taxa_df <- as.data.frame(taxa)
# Bind the clean IDs and the raw DNA sequences to the taxonomy table
taxa_df <- taxa_df %>%
  mutate(
    ASV_ID = clean_asv_ids,
    Raw_Sequence = raw_sequences
  ) %>%
  select(ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, Raw_Sequence)
rownames(taxa_df) <- NULL

# ==============================================================================
# 5. EXPORT PREPARED ASSETS
# ==============================================================================
write.csv(asv_counts, "ASV_Count_Matrix_Clean.csv", row.names = FALSE)
write.csv(taxa_df, "ASV_Taxonomy_Master.csv", row.names = FALSE)

message("Pipeline complete. Clean ASV CSVs are saved in BIOTWIN_Grand_Merge.")