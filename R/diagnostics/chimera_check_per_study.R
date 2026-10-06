# ==============================================================================
# DIAGNOSTIC: ARE THE HIGH CHIMERA RATES REAL? (run after 02_grand_merge.R)
# Three studies still lose 37-42% of reads at chimera removal after primer trimming.
# 02_grand_merge.R runs removeBimeraDenovo(method = "consensus") once, on the merged
# table of all six studies, so a sequence's verdict is a vote across every study it
# appears in. This script answers two questions:
#   1. Does running consensus chimera removal per study give different losses?
#   2. In one Brodowski 2022 sample, what are the flagged sequences? The table it
#      writes (every sequence in that sample, its abundance and its flags) is the
#      input for finding each chimera's parents and BLASTing the top ones.
# Read-only for earlier results: all outputs go into a new dated chimera_check_* folder.
# ==============================================================================

library(dada2)

# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Grand_Merge")

# The primer-trimmed Grand Merge run folder made by 02_grand_merge.R.
run_dir <- Sys.glob(file.path(workdir, "primer_trimmed_*"))
if(length(run_dir) != 1) stop("Expected exactly one primer_trimmed_* folder in ", workdir, ", found ", length(run_dir))
out_dir <- file.path(run_dir, paste0("chimera_check_", format(Sys.Date(), "%Y-%m-%d")))
if(dir.exists(out_dir)) stop(out_dir, " already exists. Rename or move it before rerunning.")
dir.create(out_dir)
message("Writing to: ", out_dir)

# Same six pre-chimera tables as 02_grand_merge.R.
find_seqtab <- function(study_folder, file) {
  hits <- Sys.glob(file.path(base_dir, study_folder, "primer_trimmed_*", file))
  if(length(hits) != 1) stop("Expected exactly one ", file, " in ", study_folder, "/primer_trimmed_*, found ", length(hits))
  hits
}
studies <- list(
  Brodowski_PRJNA715197 = c("BIOTWIN_Brodowski_Full", "seqtab_Brodowski_PRJNA715197_Cleaned.rds"),
  Brodowski_2025        = c("BIOTWIN_Brodowski_2025", "seqtab_Brodowski_2025_Master.rds"),
  Brodowski_2022        = c("BIOTWIN_Brodowski_2022_open_culture", "seqtab_Brodowski_2022_Master_open_culture.rds"),
  Duber_2025            = c("BIOTWIN_Duber_2025_16S_Only", "seqtab_Duber_2025_16S_Master.rds"),
  Duber_2022            = c("BIOTWIN_Duber_2022_Lactose_insight", "seqtab_Duber_2022_Master.rds"),
  Duber_2020            = c("BIOTWIN_Duber_2020", "seqtab_Duber_2020_Master.rds")
)
tabs <- lapply(studies, function(x) readRDS(find_seqtab(x[1], x[2])))
merged_nochim <- readRDS(file.path(run_dir, "BIOTWIN_Global_Master_Matrix_NoChim.rds"))

# --- Question 1: consensus per study vs consensus on the merged table ---
message("Running consensus chimera removal per study (several minutes per study)...")
per_study <- lapply(tabs, function(t) removeBimeraDenovo(t, method = "consensus", multithread = TRUE, verbose = TRUE))
q1 <- data.frame(study = names(tabs),
                 samples = sapply(tabs, nrow),
                 reads_before = sapply(tabs, sum),
                 reads_kept_merged_consensus = sapply(tabs, function(t) sum(merged_nochim[rownames(t), ])),
                 reads_kept_per_study_consensus = sapply(per_study, sum))
q1$pct_removed_merged <- round(100 * (1 - q1$reads_kept_merged_consensus / q1$reads_before), 1)
q1$pct_removed_per_study <- round(100 * (1 - q1$reads_kept_per_study_consensus / q1$reads_before), 1)
write.csv(q1, file.path(out_dir, "chimera_removed_per_study_vs_merged.csv"), row.names = FALSE)
print(q1, row.names = FALSE)

# --- Question 2: the flagged sequences in one Brodowski 2022 sample ---
# Pick the Brodowski 2022 sample that loses the most reads to chimera removal on the merged table.
b22 <- tabs[["Brodowski_2022"]]
removed_merged <- !(colnames(b22) %in% colnames(merged_nochim))
lost <- rowSums(b22[, removed_merged, drop = FALSE])
sample_id <- names(which.max(lost))
message("Brodowski 2022 sample examined: ", sample_id)

abund <- b22[sample_id, ]
abund <- sort(abund[abund > 0], decreasing = TRUE)
abund <- setNames(as.integer(abund), names(abund))
# Per-sample verdict, with the same parent rule consensus uses inside each sample (parents >= 1.5x as abundant).
flag_sample <- isBimeraDenovo(abund, minFoldParentOverAbundance = 1.5, multithread = FALSE)  # one sample is quick; multithread here uses mclapply, which fails on Windows
q2 <- data.frame(rank = seq_along(abund),
                 abundance = as.integer(abund),
                 pct_of_sample = round(100 * as.integer(abund) / sum(abund), 3),
                 flagged_in_this_sample = as.logical(flag_sample),
                 removed_merged_consensus = !(names(abund) %in% colnames(merged_nochim)),
                 removed_per_study_consensus = !(names(abund) %in% colnames(per_study[["Brodowski_2022"]])),
                 length = nchar(names(abund)),
                 sequence = names(abund))
write.csv(q2, file.path(out_dir, paste0("brodowski2022_", sample_id, "_sequences_and_flags.csv")), row.names = FALSE)
message("Sequences in sample: ", nrow(q2), "; flagged in this sample: ", sum(q2$flagged_in_this_sample),
        " (", round(100 * sum(q2$abundance[q2$flagged_in_this_sample]) / sum(q2$abundance), 1), "% of reads)")
message("Top 10 flagged sequences by abundance:")
print(head(q2[q2$flagged_in_this_sample, c("rank", "abundance", "pct_of_sample", "removed_merged_consensus", "removed_per_study_consensus")], 10), row.names = FALSE)
message("Done. Send the two CSV files in ", out_dir)
