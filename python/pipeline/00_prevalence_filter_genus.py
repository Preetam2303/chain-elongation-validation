# 00_prevalence_filter_genus.py
# Stage 4 of the primer-trimmed rerun: prevalence filter, then genus collapse,
# on the modelled samples only.
#
# 1. Inoculum samples (data/corrections/descriptive_only_samples.csv) are left
#    out of the modelling matrices. They are kept for descriptive use only:
#    their genus profiles are written to a separate file (decision 2026-10-10).
# 2. The ZymoBIOMICS spike-in organisms of the Brodowski 2022 batch runs
#    (genera Imtechella and Allobacillus) are removed before the filter.
# 3. Every ASV present (count > 0) in at least 5% of the modelled Illumina
#    samples is kept; each modelled sample counts once. The filter runs on
#    ASVs BEFORE genus summing, as in the published pipeline.
# 4. The kept ASVs are summed by SILVA genus.
#
# Input:
#   corrected_<date>/BIOTWIN_FINAL_ML_MATRIX.csv  from 04_metadata_merge.R
#   ASV_Taxonomy_Master.csv                       from 03_taxonomic_assignment.R
#                                                 (in the primer_trimmed_* run folder)
# Output (written into the same corrected_<date> folder, never over an existing file):
#   BIOTWIN_PRUNED_ML_MATRIX.csv    metadata + ASVs passing the filter (modelled samples)
#   BIOTWIN_GENUS_ML_MATRIX.csv     metadata + genus sums of those ASVs (modelled samples)
#   stage4_genus_summary.csv        one row per genus: ASVs, samples, % of reads,
#                                   and whether it was in the old 65-genus matrix
#   inocula_genus_profiles.csv      the inoculum samples' genus sums, same genera
#                                   (descriptive use only)

import glob
import os
from pathlib import Path

import pandas as pd

# ============================================================================
# EDIT THIS: the folder that holds your primer_trimmed_* run folder.
# ============================================================================
GRAND_MERGE_DIR = "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"

PREVALENCE = 0.05  # 5% of the modelled samples (112 -> present in >= 6)
SPIKE_IN_GENERA = ['Imtechella', 'Allobacillus']
REPO = Path(__file__).resolve().parents[2]
DESCRIPTIVE_ONLY_PATH = REPO / "data" / "corrections" / "descriptive_only_samples.csv"

run_dirs = glob.glob(os.path.join(GRAND_MERGE_DIR, "primer_trimmed_*"))
if len(run_dirs) != 1:
    raise SystemExit(f"Expected exactly one primer_trimmed_* folder in {GRAND_MERGE_DIR}, found {len(run_dirs)}")
RUN_DIR = run_dirs[0]
out_dirs = glob.glob(os.path.join(RUN_DIR, "corrected_????-??-??"))
if len(out_dirs) != 1:
    raise SystemExit(f"Expected exactly one corrected_<date> folder in {RUN_DIR} (made by 04_metadata_merge.R), "
                     f"found {len(out_dirs)}")
OUT_DIR = out_dirs[0]
print(f"Run folder: {RUN_DIR}")
print(f"Corrected-run folder: {OUT_DIR}")

PRUNED_PATH = os.path.join(OUT_DIR, "BIOTWIN_PRUNED_ML_MATRIX.csv")
GENUS_PATH = os.path.join(OUT_DIR, "BIOTWIN_GENUS_ML_MATRIX.csv")
SUMMARY_PATH = os.path.join(OUT_DIR, "stage4_genus_summary.csv")
INOCULA_PATH = os.path.join(OUT_DIR, "inocula_genus_profiles.csv")
for p in (PRUNED_PATH, GENUS_PATH, SUMMARY_PATH, INOCULA_PATH):
    if os.path.exists(p):
        raise SystemExit(f"{p} already exists. Rename or move it before rerunning.")

# The old (untrimmed) 65-genus matrix, for comparison only.
OLD_GENUS_PATH = REPO / "data" / "historical" / "BIOTWIN_GENUS_ML_MATRIX.csv"

print("Loading ML matrix and taxonomy...")
ml = pd.read_csv(os.path.join(OUT_DIR, "BIOTWIN_FINAL_ML_MATRIX.csv"), low_memory=False)
tax = pd.read_csv(os.path.join(RUN_DIR, "ASV_Taxonomy_Master.csv")).set_index("ASV_ID")

meta_cols = [c for c in ml.columns if not c.startswith("ASV_")]
asv_cols = [c for c in ml.columns if c.startswith("ASV_")]
X = ml[asv_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

# --- 1. Inoculum samples: descriptive use only ---
descriptive = pd.read_csv(DESCRIPTIVE_ONLY_PATH)
sample_ids = ml["Sample_ID"].astype(str).str.strip()
is_inoculum = sample_ids.isin(descriptive["Sample_ID"].astype(str).str.strip())
expected = int(descriptive["matrix_rows"].sum())
if is_inoculum.sum() != expected:
    raise SystemExit(f"Found {is_inoculum.sum()} inoculum rows, expected {expected} "
                     f"(data/corrections/descriptive_only_samples.csv)")
missing = set(descriptive["Sample_ID"]) - set(sample_ids)
if missing:
    raise SystemExit(f"Inoculum samples not in the matrix: {sorted(missing)}")
modelled = ~is_inoculum
print(f"Inoculum rows left out of the models: {is_inoculum.sum()} "
      f"({descriptive.shape[0]} sequencing runs); modelled rows: {modelled.sum()}")


def genus_name(asv_id):
    # Unclassified genus calls stay as their own bin rather than being dropped.
    g = str(tax.loc[asv_id, "Genus"]).strip()
    if g.lower() in ("nan", "none", ""):
        return "g__Unclassified"
    return "g__" + g.replace("[", "").replace("]", "").replace(" ", "_")


# --- 2. Spike-in organisms ---
spike = [a for a in asv_cols if str(tax.loc[a, "Genus"]).strip() in SPIKE_IN_GENERA]
carriers = ml.loc[(X[spike] > 0).any(axis=1), ["Paper_ID", "Sample_ID"]]
print(f"Spike-in ASVs removed ({', '.join(SPIKE_IN_GENERA)}): {len(spike)}; "
      f"found in {len(carriers)} rows, studies: {sorted(carriers['Paper_ID'].astype(str).str.strip().unique())}")
asv_cols = [a for a in asv_cols if a not in spike]

# --- 3. Prevalence filter on the modelled samples, each counted once ---
Xm = X.loc[modelled, asv_cols]
threshold = modelled.sum() * PREVALENCE
prevalence = (Xm > 0).sum()
keep = [a for a in asv_cols if prevalence[a] >= threshold]
genus_of = pd.Series({a: genus_name(a) for a in keep})

# --- 4. Genus collapse ---
G = X[keep].T.groupby(genus_of.values).sum().T

pruned = pd.concat([ml[meta_cols], X[keep]], axis=1)[modelled]
genus = pd.concat([ml[meta_cols], G], axis=1)
inocula = genus[is_inoculum]
genus = genus[modelled]

# Per-genus summary over the modelled samples, so every number below comes from a saved file.
Gm = G[modelled]
total_reads = Xm.values.sum()
summary = pd.DataFrame({
    "genus": G.columns,
    "n_asvs_kept": [int((genus_of == g).sum()) for g in G.columns],
    "samples_present": [int((Gm[g] > 0).sum()) for g in G.columns],
    "pct_of_all_reads": [round(Gm[g].sum() / total_reads * 100, 4) for g in G.columns],
})
old_genera = set()
if OLD_GENUS_PATH.exists():
    old = pd.read_csv(OLD_GENUS_PATH)
    old_genera = {c for c in old.columns if c.startswith("g__")}
    summary["in_old_65_genus_matrix"] = summary["genus"].isin(old_genera)
summary = summary.sort_values("pct_of_all_reads", ascending=False)

print("--- Stage 4 results ---")
print(f"Modelled rows: {modelled.sum()} (inoculum rows kept aside: {is_inoculum.sum()})")
print(f"Prevalence threshold: present in >= {threshold:.1f} modelled samples")
print(f"ASVs in matrix: {len(asv_cols) + len(spike):,} (spike-in removed: {len(spike)})")
print(f"ASVs passing the filter: {len(keep):,} (old run: 4,027)")
print(f"Reads in kept ASVs: {Xm[keep].values.sum() / total_reads * 100:.2f}% of the modelled samples' reads")
print(f"Genera after collapse: {G.shape[1]} (old run: 65)")
if old_genera:
    new_genera = set(G.columns)
    print(f"Genera in both old and new: {len(new_genera & old_genera)}")
    print(f"Old genera missing now: {sorted(old_genera - new_genera)}")
    print(f"New genera not in the old 65: {len(new_genera - old_genera)}")

pruned.to_csv(PRUNED_PATH, index=False)
genus.to_csv(GENUS_PATH, index=False)
summary.to_csv(SUMMARY_PATH, index=False)
inocula.to_csv(INOCULA_PATH, index=False)
for p in (PRUNED_PATH, GENUS_PATH, SUMMARY_PATH, INOCULA_PATH):
    print(f"Saved: {p}")
