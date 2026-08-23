# ==============================================================================
# WINDOWS LOCAL RUN: DATASET 8 (Duber 2020 - 9 Samples)
# Goal: ENA Download (ERR format) -> Filter -> Denoise -> Merge -> SAVE
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# Clean, dedicated folder for Caprobiome II
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
base_dir <- "C:/Users/IISiIS-ZWW-233/Documents"
workdir <- file.path(base_dir, "BIOTWIN_Duber_2020")
if(!dir.exists(workdir)) dir.create(workdir)
setwd(workdir)

message("Working directory set to: ", getwd())

# --- 2. The MODIFIED ENA Direct Download Function (For 10-char ERR) ---
download_ena_fastq <- function(err) {
  base_url <- "https://ftp.sra.ebi.ac.uk/vol1/fastq"
  dir1 <- substr(err, 1, 6) # e.g., "ERR320"
  
  # Crucial fix for 10-character ERR accessions
  dir2 <- paste0("00", substr(err, 10, 10)) 
  
  url_fwd <- sprintf("%s/%s/%s/%s/%s_1.fastq.gz", base_url, dir1, dir2, err, err)
  url_rev <- sprintf("%s/%s/%s/%s/%s_2.fastq.gz", base_url, dir1, dir2, err, err)
  
  dest_fwd <- file.path(workdir, paste0(err, "_1.fastq.gz"))
  dest_rev <- file.path(workdir, paste0(err, "_2.fastq.gz"))
  
  if(!file.exists(dest_fwd)) {
    message("Downloading FWD read for ", err, "...")
    download.file(url_fwd, dest_fwd, mode="wb", quiet=TRUE)
  }
  if(!file.exists(dest_rev)) {
    message("Downloading REV read for ", err, "...")
    download.file(url_rev, dest_rev, mode="wb", quiet=TRUE)
  }
}

# --- 3. Execute Download for 9 Samples ---
message("Initiating download sequence for Caprobiome II...")

# Dynamically generating the 9 ERR numbers
capro2_accs <- sprintf("ERR%d", 3200162:3200170)

for(acc in capro2_accs) {
  download_ena_fastq(acc)
}
message("All 9 downloads complete.")

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
message("Filtering reads... ")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

# --- 6. Learning Errors ---
message("Building Error Models across 9 samples...")
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

message("Saving Caprobiome II Master Sequence Table...")
saveRDS(seqtab.master, "seqtab_Duber_2020_Master.rds")

message("==========================================================")
message("PIPELINE COMPLETE. 9-Sample Master table safely saved.")
print(dim(seqtab.master))
message("==========================================================")