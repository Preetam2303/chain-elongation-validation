# 00b_nanopore_merge.py
# Stage 5 of the primer-trimmed rerun: add the Nanopore samples to the new
# Illumina genus matrix.
#
# The Nanopore study (Prusak et al., 2026; Paper_ID Hanna_2025) was not
# reprocessed: it is full-length 16S, never went through the Illumina primer
# and DADA2 steps, and its 124 genera were already bridged to SILVA names
# (synonym dictionary) in the published matrix. So its 54 rows, metadata and
# genus values are taken unchanged from data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv.
# Only the Hanna_2025 rows of data/corrections/metadata_corrections.csv are
# applied here: sample_C2 and sample_C3 are reactor B1, not B2, and on the
# phase-boundary sampling days pH, temperature and HRT take the values of the
# phase that was ending (decision 2026-10-10).
#
# The two platforms are joined as before: an outer join on genus columns, so a
# genus found on only one platform is 0 in the other platform's rows.
#
# Input:
#   BIOTWIN_GENUS_ML_MATRIX.csv in the run folder's corrected_<date> folder (from 00)
#   data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv (published; Nanopore rows only)
#   data/corrections/metadata_corrections.csv (Hanna_2025 rows only)
# Output (written into the same corrected_<date> folder, never over an existing file):
#   BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv  Illumina + Nanopore rows, all genera
#   stage5_genus_summary.csv                  one row per genus: platform, % of reads
#
# Rows: 112 modelled Illumina samples + 54 Nanopore samples = 166 (the 17
# inoculum rows and ERR3200169 are not modelled).

import glob
import os
from pathlib import Path

import pandas as pd

# ============================================================================
# EDIT THIS: the folder that holds your primer_trimmed_* run folder.
# ============================================================================
GRAND_MERGE_DIR = "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"

NANOPORE_ID = "Hanna_2025"
REPO = Path(__file__).resolve().parents[2]
PUBLISHED_PATH = REPO / "data" / "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
CORRECTIONS_PATH = REPO / "data" / "corrections" / "metadata_corrections.csv"

run_dirs = glob.glob(os.path.join(GRAND_MERGE_DIR, "primer_trimmed_*"))
if len(run_dirs) != 1:
    raise SystemExit(f"Expected exactly one primer_trimmed_* folder in {GRAND_MERGE_DIR}, found {len(run_dirs)}")
RUN_DIR = run_dirs[0]
out_dirs = glob.glob(os.path.join(RUN_DIR, "corrected_????-??-??"))
if len(out_dirs) != 1:
    raise SystemExit(f"Expected exactly one corrected_<date> folder in {RUN_DIR}, found {len(out_dirs)}")
OUT_DIR = out_dirs[0]
print(f"Run folder: {RUN_DIR}")
print(f"Corrected-run folder: {OUT_DIR}")

OUT_PATH = os.path.join(OUT_DIR, "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "stage5_genus_summary.csv")
for p in (OUT_PATH, SUMMARY_PATH):
    if os.path.exists(p):
        raise SystemExit(f"{p} already exists. Rename or move it before rerunning.")

illumina = pd.read_csv(os.path.join(OUT_DIR, "BIOTWIN_GENUS_ML_MATRIX.csv"), low_memory=False)
published = pd.read_csv(PUBLISHED_PATH, low_memory=False)
nanopore = published[published["Paper_ID"] == NANOPORE_ID].reset_index(drop=True)

meta_cols = [c for c in illumina.columns if not c.startswith("g__")]
nano_meta = [c for c in nanopore.columns if not c.startswith("g__")]
if meta_cols != nano_meta:
    raise SystemExit(f"Metadata columns differ.\nIllumina: {meta_cols}\nNanopore: {nano_meta}")
if (illumina["Paper_ID"] == NANOPORE_ID).any():
    raise SystemExit("The Illumina matrix already contains Nanopore rows.")

# Keep only genera actually detected in the Nanopore rows (the published matrix
# also carries the old Illumina genera as all-zero columns for these rows).
nano_genera = [c for c in nanopore.columns if c.startswith("g__") and nanopore[c].fillna(0).sum() > 0]
nanopore = nanopore[nano_meta + nano_genera]

# Hanna_2025 corrections, each checked against the value it replaces.
corrections = pd.read_csv(CORRECTIONS_PATH, dtype=str)
corrections = corrections[corrections["Paper_ID"] == NANOPORE_ID]
for _, k in corrections.iterrows():
    rows = nanopore.index[nanopore["Sample_ID"].astype(str).str.strip() == k["Sample_ID"]]
    if len(rows) != 1:
        raise SystemExit(f"Correction for {k['Sample_ID']}: expected one row, found {len(rows)}")
    current = nanopore.at[rows[0], k["column"]]
    cur_num, old_num = pd.to_numeric(current, errors="coerce"), pd.to_numeric(k["old_value"], errors="coerce")
    same = (abs(cur_num - old_num) < 1e-6) if pd.notna(cur_num) and pd.notna(old_num) else str(current).strip() == k["old_value"]
    if not same:
        raise SystemExit(f"Correction for {k['Sample_ID']} {k['column']}: value is {current}, expected {k['old_value']}")
    new_num = pd.to_numeric(k["new_value"], errors="coerce")
    nanopore.at[rows[0], k["column"]] = new_num if pd.notna(new_num) else k["new_value"]
print(f"Applied {len(corrections)} Nanopore correction(s)")

# Outer join on genus columns: a genus absent from one platform is 0 there.
merged = pd.concat([illumina, nanopore], axis=0, ignore_index=True, sort=False)
genus_cols = sorted(c for c in merged.columns if c.startswith("g__"))
merged[genus_cols] = merged[genus_cols].fillna(0)
merged = merged[meta_cols + genus_cols]

illumina_genera = {c for c in illumina.columns if c.startswith("g__")}
is_nano = merged["Paper_ID"] == NANOPORE_ID
ill_total = merged.loc[~is_nano, genus_cols].values.sum()
nano_total = merged.loc[is_nano, genus_cols].values.sum()
summary = pd.DataFrame({
    "genus": genus_cols,
    "in_illumina": [g in illumina_genera for g in genus_cols],
    "in_nanopore": [g in set(nano_genera) for g in genus_cols],
    "pct_of_illumina_reads": [round(merged.loc[~is_nano, g].sum() / ill_total * 100, 4) for g in genus_cols],
    "pct_of_nanopore_reads": [round(merged.loc[is_nano, g].sum() / nano_total * 100, 4) for g in genus_cols],
})
shared = summary[summary["in_illumina"] & summary["in_nanopore"]]

print("--- Stage 5 results ---")
print(f"Rows: {len(merged)} ({(~is_nano).sum()} Illumina + {is_nano.sum()} Nanopore)")
print(f"Genera: {len(genus_cols)}")
print(f"Illumina genera: {len(illumina_genera)}; Nanopore genera: {len(nano_genera)}")
print(f"Shared genera: {len(shared)}, of which named: {(shared['genus'] != 'g__Unclassified').sum()}")
print("Shared genera, % of reads (Illumina / Nanopore):")
for _, r in shared.sort_values("pct_of_illumina_reads", ascending=False).iterrows():
    print(f"  {r['genus']}: {r['pct_of_illumina_reads']:.2f} / {r['pct_of_nanopore_reads']:.2f}")

merged.to_csv(OUT_PATH, index=False)
summary.to_csv(SUMMARY_PATH, index=False)
print(f"Saved: {OUT_PATH}")
print(f"Saved: {SUMMARY_PATH}")
