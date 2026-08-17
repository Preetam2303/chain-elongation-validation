# ==============================================================================
# WINDOWS LOCAL RUN: DATASET 5 (Brodowski et al. 2022 - 10 Samples)
# Goal: ENA Download -> Filter -> Denoise -> Merge -> SAVE (NO CHIMERA REMOVAL)
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# Clean, dedicated folder for the 2022 dataset
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Brodowski_2022_open_culture")
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

# --- 3. Execute Download for 10 Samples ---
message("Initiating download sequence for Brodowski 2022 samples...")

# SRR numbers natively plugged in:
brodowski22_accs <- c(
  "SRR18962127", "SRR18962128", "SRR18962129", "SRR18962130", "SRR18962131",
  "SRR18962132", "SRR18962133", "SRR18962134", "SRR18962135", "SRR18962136"
)

for(acc in brodowski22_accs) {
  download_ena_fastq(acc)
}
message("Downloads complete.")

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
message("Filtering reads...")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

# --- 6. Learning Errors ---
message("Building Error Models across 10 samples...")
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# --- 7. Denoising & Merging ---
message("Denoising and Merging...")
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

# --- 8. Generate Sequence Table & Save ---
seqtab.master <- makeSequenceTable(mergers)

message("Saving Brodowski 2022 Master Sequence Table...")
saveRDS(seqtab.master, "seqtab_Brodowski_2022_Master_open_culture.rds")

message("==========================================================")
message("PIPELINE COMPLETE. Master table safely saved.")
print(dim(seqtab.master))
message("==========================================================")