# DIAGNOSTIC (read-only for earlier results): how much does chimera removal depend on
# minFoldParentOverAbundance? dada2's consensus default is 1.5 (a parent need only be 1.5x
# as abundant as the "chimera"). Runs consensus per study at 1.5, 4 and 8 and writes one CSV
# into a new dated folder. Download this file, then Source it in RStudio.
library(dada2)
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
run_dir <- Sys.glob(file.path(base_dir, "BIOTWIN_Grand_Merge", "primer_trimmed_*"))
if(length(run_dir) != 1) stop("Expected exactly one primer_trimmed_* folder, found ", length(run_dir))
out_dir <- file.path(run_dir, paste0("chimera_fold_check_", format(Sys.Date(), "%Y-%m-%d")))
if(dir.exists(out_dir)) stop(out_dir, " already exists. Rename or move it before rerunning.")
dir.create(out_dir)
studies <- list(
  Brodowski_PRJNA715197 = c("BIOTWIN_Brodowski_Full", "seqtab_Brodowski_PRJNA715197_Cleaned.rds"),
  Brodowski_2025        = c("BIOTWIN_Brodowski_2025", "seqtab_Brodowski_2025_Master.rds"),
  Brodowski_2022        = c("BIOTWIN_Brodowski_2022_open_culture", "seqtab_Brodowski_2022_Master_open_culture.rds"),
  Duber_2025            = c("BIOTWIN_Duber_2025_16S_Only", "seqtab_Duber_2025_16S_Master.rds"),
  Duber_2022            = c("BIOTWIN_Duber_2022_Lactose_insight", "seqtab_Duber_2022_Master.rds"),
  Duber_2020            = c("BIOTWIN_Duber_2020", "seqtab_Duber_2020_Master.rds"))
res <- list()
for(s in names(studies)) {
  f <- Sys.glob(file.path(base_dir, studies[[s]][1], "primer_trimmed_*", studies[[s]][2]))
  if(length(f) != 1) stop("Expected one table for ", s)
  tab <- readRDS(f)
  for(fold in c(1.5, 4, 8)) {
    message(s, ", minFoldParentOverAbundance = ", fold)
    nochim <- removeBimeraDenovo(tab, method = "consensus", minFoldParentOverAbundance = fold, multithread = TRUE)
    res[[length(res) + 1]] <- data.frame(study = s, min_fold_parent = fold, reads_before = sum(tab),
      reads_kept = sum(nochim), pct_removed = round(100 * (1 - sum(nochim) / sum(tab)), 1),
      asvs_before = ncol(tab), asvs_kept = ncol(nochim))
  }
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "chimera_removed_by_min_fold_parent.csv"), row.names = FALSE)
print(res, row.names = FALSE)
message("Done. Send ", file.path(out_dir, "chimera_removed_by_min_fold_parent.csv"))
