# ==============================================================================
# WINDOWS LOCAL RUN: THE GRAND MERGE & GLOBAL CHIMERA REMOVAL (6 Datasets)
# Goal: Unify 121 samples into the final BIOTWIN ML Matrix
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# Create this folder in your Documents manually and paste the 6 .rds files into it
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
workdir <- "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"
if(!dir.exists(workdir)) dir.create(workdir)
setwd(workdir)

message("Working directory set to: ", getwd())

# --- 2. Load the 6 Matrices ---
message("Loading the 6 Master Sequence Tables from the hard drive...")
tab1 <- readRDS("seqtab_Brodowski_PRJNA715197_Cleaned.rds") # Using the one without the Mock!
tab2 <- readRDS("seqtab_Brodowski_2025_Master.rds")
tab3 <- readRDS("seqtab_Brodowski_2022_Master_open_culture.rds")
tab4 <- readRDS("seqtab_Duber_2025_16S_Master.rds")
tab5 <- readRDS("seqtab_Duber_2022_Master.rds")
tab6 <- readRDS("seqtab_Duber_2020_Master.rds")

# --- 3. The Grand Merge ---
message("Merging sequence tables... (Aligning all columns mathematically)")
master_merged <- mergeSequenceTables(tab1, tab2, tab3, tab4, tab5, tab6)

message("Pre-Chimera Matrix Dimensions:")
print(dim(master_merged))

# --- 4. Global Chimera Removal (The Heavy Compute) ---
message("Engaging Global Chimera Removal on 121 samples. CPU maxed...")
# Method 'consensus' is crucial here across multiple merged runs
master_nochim <- removeBimeraDenovo(master_merged, method="consensus", multithread=FALSE, verbose=TRUE)

# --- 5. Save the Final Asset ---
message("Saving the Final Cleaned BIOTWIN Master Matrix...")
saveRDS(master_nochim, "BIOTWIN_Global_Master_Matrix_NoChim.rds")

message("==========================================================")
message("GRAND MERGE COMPLETE. The data is officially ready for ML.")
message("Final Matrix Dimensions:")
print(dim(master_nochim))
message("==========================================================")