# ==============================================================================
# WINDOWS LOCAL RUN: DATASET 6 (Duber et al. 2025 - 53 True 16S Samples)
# Goal: ENA Download -> Filter -> Denoise -> Merge -> SAVE (NO CHIMERA REMOVAL)
# Note: Excluding Shotgun and ITS datasets from the run.
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Duber_2025_16S_Only")
if(!dir.exists(workdir)) dir.create(workdir)
setwd(workdir)

message("Working directory set to: ", getwd())

# --- 2. The ENA Direct Download Function ---
download_ena_fastq <- function(srr) {
  base_url <- "https://ftp.sra.ebi.ac.uk/vol1/fastq"
  dir1 <- substr(srr, 1, 6) 
  dir2 <- paste0("0", substr(srr, 10, 11)) 
  
  url_fwd <- sprintf("%s/%s/%s/%s/%s_1.fastq.gz", base_url, dir1, dir2, srr, srr)
  url_rev <- sprintf("%s/%s/%s/%s/%s_2.fastq.gz", base_url, dir1, dir2, srr, srr)
  
  dest_fwd <- file.path(workdir, paste0(srr, "_1.fastq.gz"))
  dest_rev <- file.path(workdir, paste0(srr, "_2.fastq.gz"))
  
  if(!file.exists(dest_fwd)) {
    message("Downloading FWD read for ", srr, "...")
    download.file(url_fwd, dest_fwd, mode="wb", quiet=TRUE)
  }
  if(!file.exists(dest_rev)) {
    message("Downloading REV read for ", srr, "...")
    download.file(url_rev, dest_rev, mode="wb", quiet=TRUE)
  }
}

# --- 3. Execute Download for 53 16S Samples ---
message("Initiating targeted download sequence for Duber 2025 (16S only)...")

# Dynamically generating ONLY the 53 16S SRR numbers (169 through 221)
duber25_16s_accs <- sprintf("SRR%d", 28047169:28047221)

for(acc in duber25_16s_accs) {
  download_ena_fastq(acc)
}
message("All 53 downloads complete.")

# --- 4. File Routing ---
fnFs <- sort(list.files(workdir, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(workdir, pattern="_2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filt_path <- file.path(workdir, "filtered")
if(!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# --- 5. Filtering & Trimming ---
message("Filtering reads... (Overnight compute engaged)")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

# --- 6. Learning Errors ---
message("Building Error Models across 53 samples...")
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# --- 7. Denoising & Merging ---
message("Denoising and Merging... see you in the morning.")
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

# --- 8. Generate Sequence Table & Save ---
seqtab.master <- makeSequenceTable(mergers)

message("Saving Duber 2025 Master Sequence Table...")
saveRDS(seqtab.master, "seqtab_Duber_2025_16S_Master.rds")

message("==========================================================")
message("PIPELINE COMPLETE. 53-Sample Master table safely saved.")
print(dim(seqtab.master))
message("==========================================================")