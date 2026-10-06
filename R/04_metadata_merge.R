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

# Read from, and write into, the primer-trimmed run folder made by 02_grand_merge.R and 03_taxonomic_assignment.R.
# The metadata Excel file and the old BIOTWIN_FINAL_ML_MATRIX.csv stay in workdir itself.
run_dir <- Sys.glob(file.path(workdir, "primer_trimmed_*"))
if(length(run_dir) != 1) stop("Expected exactly one primer_trimmed_* folder in ", workdir, ", found ", length(run_dir))
if(file.exists(file.path(run_dir, "BIOTWIN_FINAL_ML_MATRIX.csv"))) stop("Metadata already merged in ", run_dir, ". Rename or move that folder before rerunning.")
message("Run folder: ", run_dir)

# 2. Load the Clean ASV Matrix
message("Loading ASV Count Matrix...")
asv_matrix <- read.csv(file.path(run_dir, "ASV_Count_Matrix_Clean.csv"), stringsAsFactors = FALSE)

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
message("Expected Rows: 130 (122 raw 16S samples + 8 added by the B1/B2 inoculum split)")
message("Actual Rows in Final Matrix: ", nrow(final_ml_matrix))
message("Total Columns (Metadata + ASVs): ", ncol(final_ml_matrix))

if(nrow(final_ml_matrix) != 130) {
  warning("Row count mismatch! Expected 130. Check your Excel Sample_IDs for typos.")
} else {
  message("Merge successful. Row count matches perfectly at 130.")
}

# Stage 3 check: the metadata part (rows, Sample_IDs, caproate) must be identical to the old matrix,
# because only the ASV columns should change in the primer-trimmed rerun.
old_file <- file.path(workdir, "BIOTWIN_FINAL_ML_MATRIX.csv")
if(file.exists(old_file)) {
  old <- read.csv(old_file, stringsAsFactors = FALSE, check.names = FALSE)
  old$Sample_ID <- gsub("[^a-zA-Z0-9]", "", old$Sample_ID)
  message("--- Comparison with the old matrix ---")
  message("Rows: new ", nrow(final_ml_matrix), ", old ", nrow(old))
  message("Same Sample_IDs in the same order: ", identical(as.character(final_ml_matrix$Sample_ID), as.character(old$Sample_ID)))
  message("Caproate identical row by row: ", isTRUE(all.equal(as.numeric(final_ml_matrix$Caproate), as.numeric(old$Caproate))))
  meta_cols <- intersect(names(metadata), names(old))
  message("All ", length(meta_cols), " metadata columns identical: ",
          isTRUE(all.equal(as.data.frame(lapply(final_ml_matrix[meta_cols], as.character)),
                           as.data.frame(lapply(old[meta_cols], as.character)), check.attributes = FALSE)))
  message("ASV columns: new ", sum(startsWith(names(final_ml_matrix), "ASV_")), ", old ", sum(startsWith(names(old), "ASV_")))
}

# 6. Export the Final Asset
message("Exporting final BIOTWIN ML Matrix...")
write.csv(final_ml_matrix, file.path(run_dir, "BIOTWIN_FINAL_ML_MATRIX.csv"), row.names = FALSE)

message("==========================================================")
message("PIPELINE COMPLETE. Data is locked and ready for modeling.")
message("==========================================================")