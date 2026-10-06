# 00_prevalence_filter_genus.py
# Stage 4 of the primer-trimmed rerun: prevalence filter, then genus collapse.
#
# Keeps every ASV present (count > 0) in at least 5% of the Illumina rows
# (7 of 130), then sums the kept ASVs by SILVA genus. The filter runs on ASVs
# BEFORE genus summing, exactly as in the published pipeline.
#
# Input (in the primer_trimmed_* run folder made by the R scripts):
#   BIOTWIN_FINAL_ML_MATRIX.csv   from 04_metadata_merge.R
#   ASV_Taxonomy_Master.csv       from 03_taxonomic_assignment.R
# Output (written into the same run folder, never over an existing file):
#   BIOTWIN_PRUNED_ML_MATRIX.csv  metadata + ASVs passing the filter
#   BIOTWIN_GENUS_ML_MATRIX.csv   metadata + genus sums of those ASVs
#   stage4_genus_summary.csv      one row per genus: ASVs, samples, % of reads,
#                                 and whether it was in the old 65-genus matrix
#
# Answer key (handoff): about 92 genera.

import glob
import os
from pathlib import Path

import pandas as pd

# ============================================================================
# EDIT THIS: the folder that holds your primer_trimmed_* run folder.
# ============================================================================
GRAND_MERGE_DIR = "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"

PREVALENCE = 0.05  # 5% of rows -> 6.5 -> present in >= 7 of 130 samples

run_dirs = glob.glob(os.path.join(GRAND_MERGE_DIR, "primer_trimmed_*"))
if len(run_dirs) != 1:
    raise SystemExit(f"Expected exactly one primer_trimmed_* folder in {GRAND_MERGE_DIR}, found {len(run_dirs)}")
RUN_DIR = run_dirs[0]
print(f"Run folder: {RUN_DIR}")

PRUNED_PATH = os.path.join(RUN_DIR, "BIOTWIN_PRUNED_ML_MATRIX.csv")
GENUS_PATH = os.path.join(RUN_DIR, "BIOTWIN_GENUS_ML_MATRIX.csv")
SUMMARY_PATH = os.path.join(RUN_DIR, "stage4_genus_summary.csv")
for p in (PRUNED_PATH, GENUS_PATH, SUMMARY_PATH):
    if os.path.exists(p):
        raise SystemExit(f"{p} already exists. Rename or move it before rerunning.")

# The old (untrimmed) 65-genus matrix, for comparison only.
OLD_GENUS_PATH = Path(__file__).resolve().parents[2] / "data" / "historical" / "BIOTWIN_GENUS_ML_MATRIX.csv"

print("Loading ML matrix and taxonomy...")
ml = pd.read_csv(os.path.join(RUN_DIR, "BIOTWIN_FINAL_ML_MATRIX.csv"), low_memory=False)
tax = pd.read_csv(os.path.join(RUN_DIR, "ASV_Taxonomy_Master.csv")).set_index("ASV_ID")

meta_cols = [c for c in ml.columns if not c.startswith("ASV_")]
asv_cols = [c for c in ml.columns if c.startswith("ASV_")]
X = ml[asv_cols].apply(pd.to_numeric, errors="coerce").fillna(0)


def genus_name(asv_id):
    # Unclassified genus calls stay as their own bin rather than being dropped.
    g = str(tax.loc[asv_id, "Genus"]).strip()
    if g.lower() in ("nan", "none", ""):
        return "g__Unclassified"
    return "g__" + g.replace("[", "").replace("]", "").replace(" ", "_")


threshold = len(ml) * PREVALENCE
prevalence = (X > 0).sum()
keep = [a for a in asv_cols if prevalence[a] >= threshold]
genus_of = pd.Series({a: genus_name(a) for a in keep})

G = X[keep].T.groupby(genus_of.values).sum().T

pruned = pd.concat([ml[meta_cols], X[keep]], axis=1)
genus = pd.concat([ml[meta_cols], G], axis=1)

# Per-genus summary, so every number below comes from a saved file.
total_reads = X.values.sum()
summary = pd.DataFrame({
    "genus": G.columns,
    "n_asvs_kept": [int((genus_of == g).sum()) for g in G.columns],
    "samples_present": [int((G[g] > 0).sum()) for g in G.columns],
    "pct_of_all_reads": [round(G[g].sum() / total_reads * 100, 4) for g in G.columns],
})
old_genera = set()
if OLD_GENUS_PATH.exists():
    old = pd.read_csv(OLD_GENUS_PATH)
    old_genera = {c for c in old.columns if c.startswith("g__")}
    summary["in_old_65_genus_matrix"] = summary["genus"].isin(old_genera)
summary = summary.sort_values("pct_of_all_reads", ascending=False)

print("--- Stage 4 results ---")
print(f"Rows: {len(ml)}")
print(f"Prevalence threshold: present in >= {threshold:.1f} samples")
print(f"ASVs in matrix: {len(asv_cols):,}")
print(f"ASVs passing the filter: {len(keep):,} (old run: 4,027)")
print(f"Reads in kept ASVs: {X[keep].values.sum() / total_reads * 100:.2f}% of all reads")
print(f"Genera after collapse: {G.shape[1]} (answer key: about 92; old run: 65)")
if old_genera:
    new_genera = set(G.columns)
    print(f"Genera in both old and new: {len(new_genera & old_genera)}")
    print(f"Old genera missing now: {sorted(old_genera - new_genera)}")
    print(f"New genera not in the old 65: {len(new_genera - old_genera)}")

pruned.to_csv(PRUNED_PATH, index=False)
genus.to_csv(GENUS_PATH, index=False)
summary.to_csv(SUMMARY_PATH, index=False)
print(f"Saved: {PRUNED_PATH}")
print(f"Saved: {GENUS_PATH}")
print(f"Saved: {SUMMARY_PATH}")
