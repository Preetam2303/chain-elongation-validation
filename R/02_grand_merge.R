# ==============================================================================
# WINDOWS LOCAL RUN: THE GRAND MERGE & GLOBAL CHIMERA REMOVAL (6 Datasets)
# Goal: Unify the six short-read studies into the final BIOTWIN ML Matrix.
# NOTE: 121/122 below refer to the raw 16S sample count at this stage; the
# final analytical matrix is 130 rows after metadata linkage (Table 1).
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Grand_Merge")
if(!dir.exists(workdir)) dir.create(workdir)
setwd(workdir)

message("Working directory set to: ", getwd())

# Every output of this run goes into a new dated folder, so earlier results are never overwritten.
run_tag <- paste0("primer_trimmed_", format(Sys.Date(), "%Y-%m-%d"))
run_dir <- file.path(workdir, run_tag)
if(file.exists(file.path(run_dir, "BIOTWIN_Global_Master_Matrix_NoChim.rds"))) stop("This run already finished: ", run_dir, ". Rename or move that folder before rerunning.")
if(!dir.exists(run_dir)) dir.create(run_dir)

# --- 2. Load the 6 Matrices ---
# The primer-trimmed tables are read straight from each study's dated run folder (made by scripts 01a-01f).
# Each study folder must contain exactly one primer_trimmed_* run folder holding the table.
find_seqtab <- function(study_folder, file) {
  hits <- Sys.glob(file.path(base_dir, study_folder, "primer_trimmed_*", file))
  if(length(hits) != 1) stop("Expected exactly one ", file, " in ", study_folder, "/primer_trimmed_*, found ", length(hits))
  hits
}
studies <- list(
  Brodowski_PRJNA715197 = c("BIOTWIN_Brodowski_Full", "seqtab_Brodowski_PRJNA715197_Cleaned.rds"), # Using the one without the Mock!
  Brodowski_2025        = c("BIOTWIN_Brodowski_2025", "seqtab_Brodowski_2025_Master.rds"),
  Brodowski_2022        = c("BIOTWIN_Brodowski_2022_open_culture", "seqtab_Brodowski_2022_Master_open_culture.rds"),
  Duber_2025            = c("BIOTWIN_Duber_2025_16S_Only", "seqtab_Duber_2025_16S_Master.rds"),
  Duber_2022            = c("BIOTWIN_Duber_2022_Lactose_insight", "seqtab_Duber_2022_Master.rds"),
  Duber_2020            = c("BIOTWIN_Duber_2020", "seqtab_Duber_2020_Master.rds")
)
message("Loading the 6 primer-trimmed Master Sequence Tables...")
tab_files <- sapply(studies, function(x) find_seqtab(x[1], x[2]))
print(tab_files)
tabs <- lapply(tab_files, readRDS)

# --- 3. The Grand Merge ---
message("Merging sequence tables... (Aligning all columns mathematically)")
master_merged <- do.call(mergeSequenceTables, unname(tabs))

message("Pre-Chimera Matrix Dimensions:")
print(dim(master_merged))

# --- 4. Global Chimera Removal (The Heavy Compute) ---
message("Engaging Global Chimera Removal on ", nrow(master_merged), " samples. CPU maxed...")
# Method 'consensus' is crucial here across multiple merged runs
master_nochim <- removeBimeraDenovo(master_merged, method="consensus", multithread=TRUE, verbose=TRUE)

# --- 5. Save the Final Asset ---
message("Saving the Final Cleaned BIOTWIN Master Matrix...")
saveRDS(master_nochim, file.path(run_dir, "BIOTWIN_Global_Master_Matrix_NoChim.rds"))

message("==========================================================")
message("GRAND MERGE COMPLETE. The data is officially ready for ML.")
message("Final Matrix Dimensions:")
print(dim(master_nochim))
message("==========================================================")

# --- 6. Stage 2 Checks: ASV counts, chimera losses and lengths, new run vs old run ---
# Answer key: about 6,000 ASVs after chimera removal (the old primers-kept run had 18,870).
message("ASVs after chimera removal: ", ncol(master_nochim), "  (answer key: about 6,000; old run 18,870)")

# Per study: reads before/after chimera removal, for this run and for the old run.
# The old run's files are the six .rds tables and the NoChim matrix in BIOTWIN_Grand_Merge itself.
study_of <- unlist(lapply(names(tabs), function(s) setNames(rep(s, nrow(tabs[[s]])), rownames(tabs[[s]]))))
chim <- data.frame(study = names(tabs),
                   samples = sapply(tabs, nrow),
                   reads_before = sapply(names(tabs), function(s) sum(master_merged[names(study_of)[study_of == s], ])),
                   reads_after = sapply(names(tabs), function(s) sum(master_nochim[names(study_of)[study_of == s], ])))
old_nochim_file <- file.path(workdir, "BIOTWIN_Global_Master_Matrix_NoChim.rds")
old_tab_files <- file.path(workdir, sapply(studies, `[`, 2))
if(file.exists(old_nochim_file) && all(file.exists(old_tab_files))) {
  old_tabs <- lapply(old_tab_files, readRDS)
  old_nochim <- readRDS(old_nochim_file)
  chim$old_reads_before <- sapply(old_tabs, sum)
  chim$old_reads_after <- sapply(old_tabs, function(t) sum(old_nochim[intersect(rownames(t), rownames(old_nochim)), ]))
}
chim <- rbind(chim, data.frame(lapply(chim, function(x) if(is.numeric(x)) sum(x) else "All"), check.names = FALSE))
chim$pct_reads_kept <- round(100 * chim$reads_after / chim$reads_before, 1)
if(!is.null(chim$old_reads_before)) chim$old_pct_reads_kept <- round(100 * chim$old_reads_after / chim$old_reads_before, 1)
asv_counts <- data.frame(metric = c("ASVs_before_chimera_removal", "ASVs_after_chimera_removal"),
                         new = c(ncol(master_merged), ncol(master_nochim)))
if(exists("old_nochim")) asv_counts$old <- c(length(unique(unlist(lapply(old_tabs, colnames)))), ncol(old_nochim))
write.csv(chim, file.path(run_dir, "chimera_removal_by_study.csv"), row.names = FALSE)
write.csv(asv_counts, file.path(run_dir, "asv_counts_new_vs_old.csv"), row.names = FALSE)
print(asv_counts)
print(chim, row.names = FALSE)

# ASV length distribution (amplicon without primers: expected about 400-430 bases).
asv_len <- nchar(colnames(master_nochim))
len_tab <- data.frame(length = sort(unique(asv_len)))
len_tab$n_ASVs <- sapply(len_tab$length, function(L) sum(asv_len == L))
len_tab$n_reads <- sapply(len_tab$length, function(L) sum(master_nochim[, asv_len == L]))
write.csv(len_tab, file.path(run_dir, "asv_length_distribution.csv"), row.names = FALSE)
outside <- asv_len < 400 | asv_len > 431
message("ASVs outside 400-431 bases: ", sum(outside), " ASVs, ",
        round(100 * sum(master_nochim[, outside]) / sum(master_nochim), 2), "% of reads")
