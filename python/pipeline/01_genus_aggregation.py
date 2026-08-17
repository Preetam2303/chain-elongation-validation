# 01_genus_aggregation.py
# Collapses ASV-level 16S/Nanopore abundance tables to genus level using the
# SILVA taxonomy assignment, resolving the cross-study, cross-platform ASV
# incompatibility problem described in Methods Section 2.4. ASV identity is a
# per-run bioinformatic artifact (denoising pipeline, reference version,
# sequencing run); the same organism receives different ASV IDs in different
# labs, so any cross-study comparison must operate at genus level or above.
#
# Input: the six short-read (Illumina) sequence tables merged in 03_grand_merge.R
#        and 04_metadata_merge.R (BIOTWIN_FINAL_ML_MATRIX.csv), plus the
#        long-read (Nanopore) contribution bridged to SILVA nomenclature
#        (see nanopore/ notebook -- taxonomy synonym dictionary resolving
#        conflicts such as Anaerotignum <-> Anaerostignum,
#        Escherichia <-> Escherichia-Shigella).
# Output: BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv, the 184-sample x 170-genus
#         matrix used by every downstream script in python/pipeline/.

import pandas as pd

MATRIX_PATH = "../../data/historical/BIOTWIN_FINAL_ML_MATRIX.csv"
TAX_PATH = "../../data/historical/ASV_Taxonomy_Master.csv"
OUTPUT_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"

print("Loading Grand Merge Matrix and Taxonomy...")
df = pd.read_csv(MATRIX_PATH)
tax_df = pd.read_csv(TAX_PATH)

# ASV -> genus lookup from the SILVA taxonomy table. Unclassified genus-level
# calls are retained as a distinct "g__Unclassified" bin rather than dropped,
# consistent with the NA-vs-zero handling used throughout this pipeline:
# "not classifiable" is a different fact from "absent."
asv_col = tax_df.columns[0]
genus_col = next((c for c in tax_df.columns if "genus" in str(c).lower()), tax_df.columns[6])

asv_to_genus = {}
for _, row in tax_df.iterrows():
    asv_id = str(row[asv_col]).strip()
    genus_val = str(row[genus_col]).strip()
    if genus_val.lower() in ["nan", "none", ""]:
        asv_to_genus[asv_id] = "g__Unclassified"
    else:
        clean_genus = genus_val.replace("[", "").replace("]", "").replace(" ", "_")
        asv_to_genus[asv_id] = f"g__{clean_genus}"

op_cols = [c for c in df.columns if not c.startswith("ASV_")]
df_op = df[op_cols]
df_asv = df[[c for c in df.columns if c.startswith("ASV_")]]

print("Aggregating ASVs into functional genera...")
df_genus = df_asv.rename(columns=asv_to_genus).groupby(level=0, axis=1).sum()

df_final = pd.concat([df_op, df_genus], axis=1)
df_final.to_csv(OUTPUT_PATH, index=False)
print(f"Collapsed to {df_genus.shape[1]} unique genera across {len(df_final)} samples.")
print(f"Saved: {OUTPUT_PATH}")

# NOTE: the Nanopore (Hanna Prusak, 2026) contribution is merged into this
# same matrix via an outer join (not inner) to preserve the full biological
# and chemical feature space -- an initial inner join collapsed the feature
# space from ~139 to 32 columns and was rejected (Methods 2.4). See the
# nanopore notebook for the taxonomy-bridging and outer-join step specific
# to integrating the long-read contribution.
