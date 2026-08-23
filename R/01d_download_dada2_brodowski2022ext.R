# ==============================================================================
# WINDOWS LOCAL RUN: FULL 16-SAMPLE PIPELINE (Brodowski PRJNA715197)
# Goal: ENA Download -> Filter -> Denoise -> Merge -> SAVE 
# NOTE: Chimera removal is skipped here per the Master Matrix Protocol.
# ==============================================================================

library(dada2)

# --- 1. Workspace Setup ---
# Creating a dedicated folder for the full Brodowski dataset
# ============================================================================
# EDIT THIS: set to your own local working directory before running.
# This script manages multi-stage SRA download / DADA2 output and is not
# intended as a one-command rerun -- see repo README for context.
# ============================================================================
workdir <- file.path(getwd(), "BIOTWIN_Brodowski_Full")
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

# --- 3. Execute Download for ALL 16 Samples ---
message("Initiating download sequence for 16 samples...")
full_accs <- c(
  "SRR13989639", "SRR13989640", "SRR13989641", "SRR13989642", 
  "SRR13989643", "SRR13989644", "SRR13989645", "SRR13989646",
  "SRR13989647", "SRR13989648", "SRR13989649", "SRR13989650", 
  "SRR13989651", "SRR13989652", "SRR13989653", "SRR13989654"
)

for(acc in full_accs) {
  download_ena_fastq(acc)
}
message("All 16 downloads complete.")

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
message("Filtering reads... (This will take a moment)")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # CRITICAL for Windows

# --- 6. The Heavy Grind: Learning Errors ---
message("Building Error Models across 16 samples. Go hit the gym...")
errF <- learnErrors(filtFs, multithread=FALSE)
errR <- learnErrors(filtRs, multithread=FALSE)

# --- 7. Denoising & Merging ---
message("Denoising and Merging... almost done.")
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

dadaFs <- dada(derepFs, err=errF, multithread=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

# --- 8. Generate Sequence Table & Save ---
seqtab.master <- makeSequenceTable(mergers)

message("Saving Master Sequence Table...")
saveRDS(seqtab.master, "seqtab_Brodowski_PRJNA715197_Master.rds")

message("==========================================================")
message("PIPELINE COMPLETE. Master table safely saved as .rds file.")
message("Dimensions of final matrix:")
print(dim(seqtab.master))
message("==========================================================")
#########################
# ==============================================================================
# POST-PIPELINE CLEANING: REMOVING THE MOCK CONTROL
# ==============================================================================

# 1. Make sure R is looking in your dedicated directory
# If you closed RStudio, run your setwd() line first:
setwd("C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Brodowski_Full")

# 2. Load the master sequence table back into R
seqtab.master <- readRDS("seqtab_Brodowski_PRJNA715197_Master.rds")

# Check the original size before modifications (Should print: 16  [number of ASVs])
message("Original matrix dimensions:")
print(dim(seqtab.master))

# 3. Use a logical filter to slice out the Mock row
# This keeps every row whose name is NOT equal to "SRR13989654"
seqtab.cleaned <- seqtab.master[rownames(seqtab.master) != "SRR13989654", ]

# 4. Safety Check: Verify the new dimensions
# The first number (rows) must be exactly 15 now!
message("New matrix dimensions after removing Mock:")
print(dim(seqtab.cleaned))

# 5. Save the cleaned matrix to a brand new file
# Pro-tip: Saving it under a new name keeps your original 16-sample file safe as a backup
saveRDS(seqtab.cleaned, "seqtab_Brodowski_PRJNA715197_Cleaned.rds")

message("Operation successful! Cleaned matrix saved as seqtab_Brodowski_PRJNA715197_Cleaned.rds")