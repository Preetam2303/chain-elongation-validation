# ==============================================================================
# WINDOWS LOCAL RUN: DATASET 6 (Duber et al. 2025 - 53 True 16S Samples)
# Goal: ENA Download -> Filter -> Denoise -> Merge -> SAVE (NO CHIMERA REMOVAL)
# Note: Excluding Shotgun and ITS datasets from the run.
# ==============================================================================

library(dada2)

# --- Primer removal settings (primer-trimmed rerun) ---
# EDIT THIS: path to cutadapt (5.2, installed in its own conda environment "cutadapt").
cutadapt <- "C:/Users/IISiIS-ZWW-233/anaconda3/envs/cutadapt/Scripts/cutadapt.exe"
primer_fwd <- "CCTACGGGNGGCWGCAG"      # 341F, 17 bases
primer_rev <- "GACTACHVGGGTATCTAATCC"  # 805R, 21 bases
# Every output of this run goes into a new dated folder, so earlier results are never overwritten.
run_tag <- paste0("primer_trimmed_", format(Sys.Date(), "%Y-%m-%d"))

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

# All outputs of this run go into the dated run folder; the old filtered/ folder and .rds are left untouched.
run_dir <- file.path(workdir, run_tag)
if(file.exists(file.path(run_dir, "seqtab_Duber_2025_16S_Master.rds"))) stop("This run already finished: ", run_dir, ". Rename or move that folder before rerunning.")
if(!dir.exists(run_dir)) dir.create(run_dir)

# --- 5. Primer Removal (cutadapt) ---
# The 341F/805R primers are cut off before DADA2. A fixed trimLeft cut is not used because some primer
# copies are 1-2 bases shorter than 17/21 bases; cutadapt finds where the primer ends in each read.
# ^ anchors each primer at the start of the read, -e 2 allows up to 2 mismatches or indels, and
# read pairs missing either primer are discarded (--discard-untrimmed with --pair-filter=any).
# The ~464 bp amplicon is longer than the reads, so the opposite primer never appears at the read ends.
cut_path <- file.path(run_dir, "cutadapt")
if(!dir.exists(cut_path)) dir.create(cut_path)
cutFs <- file.path(cut_path, basename(fnFs))
cutRs <- file.path(cut_path, basename(fnRs))

message("Removing primers with cutadapt...")
cut_reports <- lapply(seq_along(fnFs), function(i) {
  report <- system2(cutadapt, args = c(
    "-g", paste0("^", primer_fwd), "-G", paste0("^", primer_rev),
    "-e", "2", "--discard-untrimmed", "--pair-filter=any", "--report=minimal",
    "-o", shQuote(cutFs[i]), "-p", shQuote(cutRs[i]),
    shQuote(fnFs[i]), shQuote(fnRs[i])), stdout = TRUE)
  if(!is.null(attr(report, "status"))) stop("cutadapt failed on ", basename(fnFs[i]))
  read.delim(text = report)
})

filt_path <- file.path(run_dir, "filtered")
if(!dir.exists(filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# --- 6. Filtering & Trimming ---
message("Filtering reads... (Overnight compute engaged)")
# truncLen is the old 250 minus the primer lengths (17, 21), so reads cover the same 16S stretch as before.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(233,229),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

# --- 7. Learning Errors ---
message("Building Error Models across 53 samples...")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# --- 8. Denoising & Merging ---
message("Denoising and Merging... see you in the morning.")
derepFs <- derepFastq(filtFs, verbose=FALSE)
derepRs <- derepFastq(filtRs, verbose=FALSE)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=FALSE)

# --- 9. Generate Sequence Table & Save ---
seqtab.master <- makeSequenceTable(mergers)

message("Saving Duber 2025 Master Sequence Table...")
saveRDS(seqtab.master, file.path(run_dir, "seqtab_Duber_2025_16S_Master.rds"))

message("==========================================================")
message("PIPELINE COMPLETE. 53-Sample Master table safely saved.")
print(dim(seqtab.master))
message("==========================================================")

# --- 10. Stage 1 Checks & Read Tracking ---
# No ASV should still start with the 341F primer or end with the reverse complement of 805R.
asvs <- colnames(seqtab.master)
message("ASVs starting with 341F (CCTACGGG): ", sum(grepl("^CCTACGGG", asvs)), "  (should be 0)")
message("ASVs ending with 805R reverse complement: ",
        sum(grepl("GGATTAGATACCC[CGT][AGT]GTAGTC$", asvs)), "  (should be 0)")

# Reads kept at each step, per sample. merged_old_run is the same sample in the old (primers kept) run.
getN <- function(x) sum(getUniques(x))
track <- data.frame(sample = sample.names,
                    raw = sapply(cut_reports, `[[`, "in_reads"),
                    primers_removed = sapply(cut_reports, `[[`, "out_reads"),
                    filtered = out[, "reads.out"],
                    denoisedF = sapply(dadaFs, getN),
                    denoisedR = sapply(dadaRs, getN),
                    merged = sapply(mergers, getN))
track$pct_with_primers <- round(100 * track$primers_removed / track$raw, 1)
track$pct_merged_of_filtered <- round(100 * track$merged / track$filtered, 1)
old_file <- file.path(workdir, "seqtab_Duber_2025_16S_Master.rds")
if(file.exists(old_file)) {
  track$merged_old_run <- rowSums(readRDS(old_file))[track$sample]
}
write.csv(track, file.path(run_dir, "read_tracking_Duber_2025_16S_Master.csv"), row.names = FALSE)
print(track)
