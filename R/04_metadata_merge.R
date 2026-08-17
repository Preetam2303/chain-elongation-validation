# ==============================================================================
# FINAL STAGE: MERGE ASV MATRIX WITH EXCEL METADATA (OPERATIONAL_DATA)
# ==============================================================================

library(readxl)
library(dplyr)

# 1. Set Workspace
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
workdir <- "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"
setwd(workdir)

# 2. Load the Clean ASV Matrix
message("Loading ASV Count Matrix...")
asv_matrix <- read.csv("ASV_Count_Matrix_Clean.csv", stringsAsFactors = FALSE)

# 3. Load the Excel Metadata
metadata_file <- "XGBoost_Matrix WITHOUT COMMUNITY DATA.xlsx" 
message("Loading Excel Metadata...")
metadata <- read_excel(metadata_file)

# ------------------------------------------------------------------------------
# THE SCRUBBER: Prevent join errors by stripping hidden characters (\r, \n, spaces)
# ------------------------------------------------------------------------------
message("Scrubbing invisible characters from Sample_IDs...")
metadata$Sample_ID <- gsub("[^a-zA-Z0-9]", "", metadata$Sample_ID)
asv_matrix$Sample_ID <- gsub("[^a-zA-Z0-9]", "", asv_matrix$Sample_ID)

# 4. Perform the Merge
message("Executing Inner Join on Sample_ID...")
final_ml_matrix <- inner_join(metadata, asv_matrix, by = "Sample_ID")

# 5. Validation Checks
message("--- Merge Validation ---")
message("Expected Rows: 122 (Due to the B1/B2 Inoculum Split)")
message("Actual Rows in Final Matrix: ", nrow(final_ml_matrix))
message("Total Columns (Metadata + ASVs): ", ncol(final_ml_matrix))

if(nrow(final_ml_matrix) != 122) {
  warning("Row count mismatch! Check your Excel Sample_IDs for typos.")
} else {
  message("Merge successful. Row count matches perfectly at 122.")
}

# 6. Export the Final Asset
message("Exporting final BIOTWIN ML Matrix...")
write.csv(final_ml_matrix, "BIOTWIN_FINAL_ML_MATRIX.csv", row.names = FALSE)

message("==========================================================")
message("PIPELINE COMPLETE. Data is locked and ready for modeling.")
message("==========================================================")