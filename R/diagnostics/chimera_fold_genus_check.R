# ==============================================================================
# DIAGNOSTIC (read-only for earlier results): does the chimera setting change the genus picture?
# For the three studies with high chimera losses, runs consensus chimera removal per study at
# minFoldParentOverAbundance = 1.5 (used) and 8, assigns SILVA genera, and compares the mean
# relative abundance of the top genera (and Caproiciproducens) between the two settings.
# Sequences already in ASV_Taxonomy_Master.csv reuse that assignment; only the extra sequences
# kept at 8x are classified here (same SILVA file and seed as 03_taxonomic_assignment.R).
# Output: one CSV in a new dated chimera_genus_check_* folder inside the run folder.
# ==============================================================================

library(dada2)

# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Grand_Merge")
train_file <- file.path(workdir, "silva_nr99_v138.1_train_set.fa.gz")

run_dir <- Sys.glob(file.path(workdir, "primer_trimmed_*"))
if(length(run_dir) != 1) stop("Expected exactly one primer_trimmed_* folder in ", workdir, ", found ", length(run_dir))
if(!file.exists(train_file)) stop("SILVA training set not found: ", train_file)
out_dir <- file.path(run_dir, paste0("chimera_genus_check_", format(Sys.Date(), "%Y-%m-%d")))
if(dir.exists(out_dir)) stop(out_dir, " already exists. Rename or move it before rerunning.")
dir.create(out_dir)
message("Writing to: ", out_dir)

studies <- list(
  Brodowski_PRJNA715197 = c("BIOTWIN_Brodowski_Full", "seqtab_Brodowski_PRJNA715197_Cleaned.rds"),
  Brodowski_2022        = c("BIOTWIN_Brodowski_2022_open_culture", "seqtab_Brodowski_2022_Master_open_culture.rds"),
  Duber_2025            = c("BIOTWIN_Duber_2025_16S_Only", "seqtab_Duber_2025_16S_Master.rds"))
folds <- c(1.5, 8)
n_top <- 10

# --- Chimera removal per study at both settings ---
nochim <- list()
for(s in names(studies)) {
  f <- Sys.glob(file.path(base_dir, studies[[s]][1], "primer_trimmed_*", studies[[s]][2]))
  if(length(f) != 1) stop("Expected one table for ", s)
  tab <- readRDS(f)
  for(fold in folds) {
    message(s, ", minFoldParentOverAbundance = ", fold)
    nochim[[paste(s, fold)]] <- removeBimeraDenovo(tab, method = "consensus", minFoldParentOverAbundance = fold, multithread = TRUE)
  }
}

# --- Genus for every kept sequence ---
master <- read.csv(file.path(run_dir, "ASV_Taxonomy_Master.csv"), stringsAsFactors = FALSE)
label <- function(genus, family) ifelse(is.na(genus) | genus == "",
                                        paste0("Unclassified (", ifelse(is.na(family) | family == "", "unknown family", family), ")"),
                                        genus)
genus_of <- setNames(label(master$Genus, master$Family), master$Raw_Sequence)
all_seqs <- unique(unlist(lapply(nochim, colnames)))
new_seqs <- setdiff(all_seqs, names(genus_of))
message("Sequences to classify that are not in ASV_Taxonomy_Master.csv: ", length(new_seqs))
if(length(new_seqs) > 0) {
  set.seed(100)
  tx <- assignTaxonomy(new_seqs, train_file, multithread = TRUE)
  genus_of <- c(genus_of, setNames(label(tx[, "Genus"], tx[, "Family"]), rownames(tx)))
}

# --- Mean per-sample relative abundance (%) of each genus ---
genus_pct <- function(t) {
  rel <- t / pmax(rowSums(t), 1)
  g <- rowsum(t(rel), genus_of[colnames(t)])  # genera x samples
  sort(100 * rowMeans(g), decreasing = TRUE)
}
res <- list()
for(s in names(studies)) {
  p15 <- genus_pct(nochim[[paste(s, 1.5)]])
  p8  <- genus_pct(nochim[[paste(s, 8)]])
  top <- unique(c(names(p15)[1:min(n_top, length(p15))], names(p8)[1:min(n_top, length(p8))], "Caproiciproducens"))
  get <- function(p, g) ifelse(g %in% names(p), p[g], 0)
  res[[s]] <- data.frame(study = s, genus = top,
                         mean_pct_fold_1.5 = round(get(p15, top), 2),
                         mean_pct_fold_8 = round(get(p8, top), 2),
                         rank_fold_1.5 = match(top, names(p15)),
                         rank_fold_8 = match(top, names(p8)))
  res[[s]]$difference_pct_points <- res[[s]]$mean_pct_fold_8 - res[[s]]$mean_pct_fold_1.5
  res[[s]] <- res[[s]][order(res[[s]]$rank_fold_1.5), ]
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_dir, "top_genera_fold_1.5_vs_8.csv"), row.names = FALSE)
print(res, row.names = FALSE)
message("Done. Send ", file.path(out_dir, "top_genera_fold_1.5_vs_8.csv"))
