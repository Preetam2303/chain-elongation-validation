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

# Read from the primer-trimmed run folder made by 02_grand_merge.R and 03_taxonomic_assignment.R, and write
# the corrected run into a new dated folder inside it (corrected_<date>), so earlier outputs stay untouched.
# The metadata Excel file and the old BIOTWIN_FINAL_ML_MATRIX.csv stay in workdir itself.
run_dir <- Sys.glob(file.path(workdir, "primer_trimmed_*"))
if(length(run_dir) != 1) stop("Expected exactly one primer_trimmed_* folder in ", workdir, ", found ", length(run_dir))
existing <- Sys.glob(file.path(run_dir, "corrected_????-??-??"))
if(length(existing) > 0) stop("A corrected-run folder already exists: ", paste(existing, collapse = ", "), ". Rename or move it (for example add _try1 to its name) before rerunning.")
out_dir <- file.path(run_dir, paste0("corrected_", format(Sys.Date(), "%Y-%m-%d")))
message("Run folder: ", run_dir)
message("Will write to: ", out_dir)

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

# ------------------------------------------------------------------------------
# DATA FIXES: excluded samples and value corrections, kept in the repository's data/corrections/
# (one row per change, with its source and reason). Hanna_2025 (Nanopore) rows are not in this
# matrix; their corrections are applied when the Nanopore rows are merged (Stage 5).
# ------------------------------------------------------------------------------
# EDIT THIS if your repository clone is somewhere else.
repo_dir <- "C:/Users/IISiIS-ZWW-233/Documents/chain-elongation-validation"
excluded <- read.csv(file.path(repo_dir, "data", "corrections", "excluded_samples.csv"), stringsAsFactors = FALSE)
corrections <- read.csv(file.path(repo_dir, "data", "corrections", "metadata_corrections.csv"), stringsAsFactors = FALSE)
corrections <- corrections[corrections$Paper_ID != "Hanna_2025", ]
excluded$Sample_ID <- gsub("[^a-zA-Z0-9]", "", excluded$Sample_ID)
corrections$Sample_ID <- gsub("[^a-zA-Z0-9]", "", corrections$Sample_ID)

# Applies the corrections to a metadata table, after checking each cell still holds the old value.
# An old value of NA means the cell must be blank (NA, empty or the text "NA").
is_blank <- function(x) is.na(x) || as.character(x) %in% c("", "NA")
is_blank_vec <- function(x) is.na(x) | as.character(x) %in% c("", "NA")
apply_corrections <- function(tab, corrections) {
  for(i in seq_len(nrow(corrections))) {
    k <- corrections[i, ]
    row <- which(tab$Sample_ID == k$Sample_ID)
    if(length(row) != 1) stop("Correction ", i, ": expected one row for ", k$Sample_ID, ", found ", length(row))
    current <- tab[[k$column]][row]
    num_old <- suppressWarnings(as.numeric(k$old_value))
    same <- if(is_blank(k$old_value)) is_blank(current)
            else if(is.na(num_old)) identical(as.character(current), as.character(k$old_value))
            else isTRUE(abs(suppressWarnings(as.numeric(current)) - num_old) < 0.01)
    if(!same) stop("Correction ", i, ": ", k$Sample_ID, " ", k$column, " is ", current, ", expected old value ", k$old_value)
    num_new <- suppressWarnings(as.numeric(k$new_value))
    if(is.character(tab[[k$column]]) && !is.na(num_new) && all(is_blank_vec(tab[[k$column]]) | !is.na(suppressWarnings(as.numeric(tab[[k$column]]))))) {
      tab[[k$column]] <- suppressWarnings(as.numeric(tab[[k$column]]))  # a numeric column read as text because of "NA" cells
    }
    tab[[k$column]][row] <- if(is.character(tab[[k$column]])) as.character(k$new_value) else num_new
  }
  tab
}

message("Excluding ", nrow(excluded), " sample(s): ", paste(excluded$Sample_ID, collapse = ", "))
n_excluded_rows <- sum(metadata$Sample_ID %in% excluded$Sample_ID)
metadata <- metadata[!(metadata$Sample_ID %in% excluded$Sample_ID), ]
message("Applying ", nrow(corrections), " value corrections...")
metadata <- apply_corrections(metadata, corrections)

# 4. Perform the Merge
message("Executing Inner Join on Sample_ID...")
final_ml_matrix <- inner_join(metadata, asv_matrix, by = "Sample_ID")

# 5. Validation Checks
# All rows stay in this matrix, including the 17 inoculum rows; 00_prevalence_filter_genus.py leaves the
# inoculum rows out of the modelling matrices (data/corrections/descriptive_only_samples.csv).
expected_rows <- 130 - n_excluded_rows
message("--- Merge Validation ---")
message("Expected Rows: ", expected_rows, " (122 raw 16S samples + 8 added by the B1/B2 inoculum split, minus ", n_excluded_rows, " excluded)")
message("Actual Rows in Final Matrix: ", nrow(final_ml_matrix))
message("Total Columns (Metadata + ASVs): ", ncol(final_ml_matrix))

if(nrow(final_ml_matrix) != expected_rows) {
  warning("Row count mismatch! Expected ", expected_rows, ". Check your Excel Sample_IDs for typos.")
} else {
  message("Merge successful. Row count matches perfectly at ", expected_rows, ".")
}

# Stage 3 check: the metadata part (rows, Sample_IDs, caproate) must be identical to the old matrix
# once the same exclusions and corrections are applied to it, because only the ASV columns and the
# listed fixes should change in the primer-trimmed rerun.
old_file <- file.path(workdir, "BIOTWIN_FINAL_ML_MATRIX.csv")
if(file.exists(old_file)) {
  old <- read.csv(old_file, stringsAsFactors = FALSE, check.names = FALSE)
  old$Sample_ID <- gsub("[^a-zA-Z0-9]", "", old$Sample_ID)
  old <- apply_corrections(old[!(old$Sample_ID %in% excluded$Sample_ID), ], corrections)
  message("--- Comparison with the old matrix (same exclusions and corrections applied) ---")
  message("Rows: new ", nrow(final_ml_matrix), ", old ", nrow(old))
  message("Same Sample_IDs in the same order: ", identical(as.character(final_ml_matrix$Sample_ID), as.character(old$Sample_ID)))
  message("Caproate identical row by row: ", isTRUE(all.equal(as.numeric(final_ml_matrix$Caproate), as.numeric(old$Caproate))))
  meta_cols <- intersect(names(metadata), names(old))
  # Compares numbers as numbers (so 70.74 equals 70.739999999999995) and treats "NA" text, NA and blank as missing.
  same_col <- function(a, b) {
    a <- as.character(a); b <- as.character(b)
    a[is.na(a) | a %in% c("", "NA")] <- NA; b[is.na(b) | b %in% c("", "NA")] <- NA
    na <- suppressWarnings(as.numeric(a)); nb <- suppressWarnings(as.numeric(b))
    both_num <- !is.na(na) & !is.na(nb)
    all(ifelse(both_num, abs(na - nb) < 1e-6, (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)))
  }
  diff_cols <- meta_cols[!mapply(same_col, final_ml_matrix[meta_cols], old[meta_cols])]
  message("All ", length(meta_cols), " metadata columns identical: ", length(diff_cols) == 0,
          if(length(diff_cols) > 0) paste0(" (differ: ", paste(diff_cols, collapse = ", "), ")") else "")
  message("ASV columns: new ", sum(startsWith(names(final_ml_matrix), "ASV_")), ", old ", sum(startsWith(names(old), "ASV_")))
}

# 6. Export the Final Asset
message("Exporting final BIOTWIN ML Matrix...")
# The folder is created only now, so a run that stops on a failed correction leaves nothing behind.
dir.create(out_dir)
write.csv(final_ml_matrix, file.path(out_dir, "BIOTWIN_FINAL_ML_MATRIX.csv"), row.names = FALSE)
message("Saved: ", file.path(out_dir, "BIOTWIN_FINAL_ML_MATRIX.csv"))

message("==========================================================")
message("PIPELINE COMPLETE. Data is locked and ready for modeling.")
message("==========================================================")