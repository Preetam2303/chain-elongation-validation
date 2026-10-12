# 12_table_s2_shap_ranking.py
# Stage 7 of the primer-trimmed rerun (handoff: "Rebuild Table S2"), steps 1
# and 2: rank the ASVs by mean |SHAP| from a saved, dated model run, and write
# the primer-free sequences to BLAST to a FASTA with the date in its name.
# 12b_table_s2_from_blast.py then builds Table S2 from the BLAST text file
# (step 4). The published Table S2 cannot be traced to any saved ranking, so it
# is rebuilt from scratch. Decided by Preetam on 2026-10-12.
#
# Model: one XGBoost fit on all modelled Illumina samples (BIOTWIN_PRUNED), with
# the core set (pH, temperature, HRT, lactate, acetate, ethanol, butyrate,
# propionate, the 2 categoricals) and the ASVs, the shared settings in
# pipeline_settings.py (no subsampling) and TreeSHAP on the same samples. This
# is how Figure 5 ranks genera: a descriptive ranking of the sequences the
# model relies on, with no split, so it is the same every time. The ASVs enter
# as relative abundance (each sample's counts divided by its total over the
# kept ASVs), as the genera do in the headline. Numeric predictors are scaled;
# this does not change the ranks.
#
# The FASTA holds:
#   - the top 20 ASVs by mean |SHAP|, and
#   - the Caproiciproducens ASVs that together hold 90% of that genus's reads
#     (all of its ASVs in the modelled samples, before the prevalence filter;
#     BIOTWIN_ASV_ALL), even outside the top 20, so the strain-level evidence
#     rests on the whole genus and not only on the sequences that rank highest.
#
# For every ASV the table also gives the number of studies it is found in
# (count > 0 in any sample), as Figure 6 does for genera: a top ASV found in
# one study can mark that study rather than the guild. As a stability check, it
# gives the ASV's rank in each leave-one-study-out fit (fitted on the other
# studies' samples, SHAP on those samples, as in 09). An ASV found in one
# study only gets no attribution when that study is held out.
#
# Ties are ranked deterministically (average rank, as in 09); the top 20 are
# the 20 best ranks, ties broken by ASV number.
#
# Inputs: BIOTWIN_PRUNED (00), BIOTWIN_ASV_ALL (04_metadata_merge.R's
# BIOTWIN_FINAL_ML_MATRIX.csv, every ASV) and BIOTWIN_TAXONOMY (the run's
# ASV_Taxonomy_Master.csv, from 03_taxonomic_assignment.R, with the sequences).
# Output (in BIOTWIN_RESULTS_DIR): 12_asv_shap_ranking.csv (every modelled ASV),
# 12_table_s2_asvs.csv (the ASVs in the FASTA) and Table_S2_ASVs_<date>.fasta.

import datetime
import os
import re
import warnings

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from scipy.stats import rankdata
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler

from pipeline_settings import CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, drop_unmodelled, present, save_table

warnings.filterwarnings('ignore')

PRUNED_PATH = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
ASV_ALL_PATH = os.environ.get("BIOTWIN_ASV_ALL", "../../data/historical/BIOTWIN_FINAL_ML_MATRIX.csv")
TAXONOMY_PATH = os.environ.get("BIOTWIN_TAXONOMY", "../../data/historical/ASV_Taxonomy_Master.csv")
N_TOP = 20
WHOLE_GENUS = 'Caproiciproducens'
GENUS_READ_SHARE = 0.90
# Stage 1's check: no ASV may start with the 341F primer or end with the
# reverse complement of 805R (GACTACHVGGGTATCTAATCC).
PRIMER_341F_START = 'CCTACGGG'
PRIMER_805R_RC_END = re.compile(r'GGATTAGATACCC[CGT][AGT]GTAGTC$')


def read_asv_table(path):
    df = pd.read_csv(path, low_memory=False)
    for c in ['Paper_ID', 'Sample_ID', 'BIOREACTOR']:
        df[c] = df[c].astype(str).str.strip()
    return drop_unmodelled(df)


def load_asv_matrix(path):
    """The ASV-level matrix as in 10_table3a_naive_splits.py, plus the ASVs as
    relative abundance."""
    df = read_asv_table(path)
    for c in OPS + CORE_CHEM + [TARGET]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    for c in CAT_COLS:
        df[c] = df[c].astype('category')
    if df[TARGET].isna().any():
        raise SystemExit("Caproate is missing in the ASV-level matrix; every modelled row needs it.")
    asv_cols = [c for c in df.columns if c.startswith('ASV_')]
    counts = df[asv_cols].fillna(0)
    rel = df.copy()
    rel[asv_cols] = counts.div(counts.sum(axis=1).replace(0, 1), axis=0)
    return rel, asv_cols


def asv_ranks(X, y, asv_cols):
    """Fit on these rows, TreeSHAP on the same rows; mean |SHAP| and rank of each ASV."""
    Xs = X.copy()
    num = [c for c in Xs.columns if c not in CAT_COLS]
    Xs[num] = StandardScaler().fit_transform(Xs[num])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(Xs, y)
    mean_abs = pd.Series(np.abs(shap.TreeExplainer(model).shap_values(Xs)).mean(axis=0), index=Xs.columns)[asv_cols]
    return mean_abs, pd.Series(rankdata(-mean_abs.values, method='average'), index=asv_cols)


def asv_number(a):
    return int(a.split('_')[1])


def lowest_taxon(row):
    """The SILVA genus, or the lowest rank assigned when there is no genus."""
    for level in ['Genus', 'Family', 'Order', 'Class', 'Phylum']:
        if isinstance(row[level], str):
            return row[level] if level == 'Genus' else f"{level.lower()} {row[level]}"
    return 'unassigned'


df, asv_cols = load_asv_matrix(PRUNED_PATH)
features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df))) + asv_cols
studies = sorted(df['Paper_ID'].unique())

# Every ASV of the modelled samples, for read shares and study counts
every = read_asv_table(ASV_ALL_PATH).set_index('Sample_ID')
if len(every) != len(df) or set(every.index) != set(df['Sample_ID']):
    raise SystemExit(f"{ASV_ALL_PATH} and {PRUNED_PATH} do not hold the same modelled samples.")
all_asvs = [c for c in every.columns if c.startswith('ASV_')]
all_counts = every[all_asvs].fillna(0)
tax = pd.read_csv(TAXONOMY_PATH).set_index('ASV_ID')
missing = [a for a in all_asvs if a not in tax.index]
if missing:
    raise SystemExit(f"ASVs missing from {TAXONOMY_PATH}: {missing[:10]}")
genus = tax['Genus'].reindex(all_asvs)
asv_reads = all_counts.sum()
genus_reads = asv_reads.groupby(genus.values).sum()  # ASVs without a genus are left out
print(f"ASV-level matrix: {len(df)} samples, {len(studies)} studies, {len(asv_cols)} ASVs after the prevalence "
      f"filter ({len(all_asvs)} before), {len(features) - len(asv_cols)} other predictors")

# ==============================================================================
# 1. THE RANKING: ONE FIT ON ALL SAMPLES
# ==============================================================================
shap_all, rank_all = asv_ranks(df[features], df[TARGET], asv_cols)

# ==============================================================================
# 2. STABILITY: THE SAME RANKING IN EACH LEAVE-ONE-STUDY-OUT FOLD
# ==============================================================================
fold_ranks = {}
for tr, te in LeaveOneGroupOut().split(df, groups=df['Paper_ID']):
    held = df['Paper_ID'].iloc[te].iloc[0]
    fold_ranks[held] = asv_ranks(df[features].iloc[tr], df[TARGET].iloc[tr], asv_cols)[1]
fold_ranks = pd.DataFrame(fold_ranks)[studies]


def describe(asvs):
    """Taxonomy, read share and where each ASV is found."""
    paper = every['Paper_ID']
    return pd.DataFrame({
        'ASV_ID': asvs,
        'silva_taxon': [lowest_taxon(tax.loc[a]) for a in asvs],
        'genus': genus.reindex(asvs).values,
        'n_samples': [int((all_counts[a] > 0).sum()) for a in asvs],
        'n_studies': [paper[all_counts[a] > 0].nunique() for a in asvs],
        'pct_of_genus_reads': [100 * asv_reads[a] / genus_reads[g] if isinstance(g, str) else np.nan
                               for a, g in zip(asvs, genus.reindex(asvs).values)],
    })


# ==============================================================================
# 3. THE RANKING TABLE (EVERY MODELLED ASV)
# ==============================================================================
table = describe(asv_cols)
table.insert(1, 'rank', rank_all.values)
table.insert(2, 'mean_abs_shap', shap_all.values)
table['mean_rel_abundance_pct'] = 100 * df[asv_cols].mean().values
table['loso_folds_in_top20'] = (fold_ranks <= N_TOP).sum(axis=1).values
table['loso_median_rank'] = fold_ranks.median(axis=1).values
for s in studies:
    table[f'rank_without_{s}'] = fold_ranks[s].values
table = table.assign(_n=table['ASV_ID'].map(asv_number)).sort_values(['rank', '_n']).drop(columns='_n')
table = table.reset_index(drop=True)
save_table(table, "12_asv_shap_ranking.csv")

# ==============================================================================
# 4. THE ASVs TO BLAST: TOP 20, AND THE CAPROICIPRODUCENS ASVs WITH 90% OF ITS READS
# ==============================================================================
top = list(table['ASV_ID'].head(N_TOP))
g_reads = asv_reads[(genus == WHOLE_GENUS).values]
g_reads = g_reads[g_reads > 0]
order = sorted(g_reads.index, key=lambda a: (-g_reads[a], asv_number(a)))
cum = g_reads[order].cumsum() / g_reads.sum()
n_genus = int((cum < GENUS_READ_SHARE).sum()) + 1
genus_set = order[:n_genus]

chosen = top + [a for a in genus_set if a not in top]
blast = describe(chosen)
ranked = table.set_index('ASV_ID')


def shap_rank_text(a):
    if a not in ranked.index:
        return 'not kept'
    return f"{ranked.loc[a, 'rank']:g}" if ranked.loc[a, 'mean_abs_shap'] > 0 else 'none'


in_top, in_genus = f"SHAP top {N_TOP}", f"{WHOLE_GENUS} {100 * GENUS_READ_SHARE:.0f}% of reads"
blast.insert(1, 'included_as', [f"{in_top} and {in_genus}" if a in top and a in genus_set
                                else in_top if a in top else in_genus for a in chosen])
blast.insert(2, 'rank', [ranked.loc[a, 'rank'] if a in ranked.index else np.nan for a in chosen])
blast.insert(3, 'mean_abs_shap', [ranked.loc[a, 'mean_abs_shap'] if a in ranked.index else np.nan for a in chosen])
blast['cum_pct_of_genus_reads'] = [100 * cum[a] if a in genus_set else np.nan for a in chosen]
blast['loso_folds_in_top20'] = [ranked.loc[a, 'loso_folds_in_top20'] if a in ranked.index else np.nan
                                for a in chosen]
blast['sequence'] = tax['Raw_Sequence'].reindex(chosen).values
blast['length'] = blast['sequence'].str.len()
bad = [a for a, s in zip(blast['ASV_ID'], blast['sequence'])
       if s.startswith(PRIMER_341F_START) or PRIMER_805R_RC_END.search(s)]
if bad:
    raise SystemExit(f"Primer sequence found in {bad}; the FASTA must hold primer-free sequences.")

print(f"\n--- TOP {N_TOP} ASVs BY MEAN |SHAP| (one fit on all {len(df)} samples; caproate units, mM C) ---")
print(f"studies: in how many of the {len(studies)} studies the ASV is found. LOSO: its rank in each of the "
      f"{len(studies)} leave-one-study-out fits (in the top {N_TOP} in how many; median rank)")
print("-" * 126)
print(f"{'rank':>4}  {'ASV':<9} {'genus (SILVA)':<38} {'mean|SHAP|':>10} {'samples':>8} {'studies':>8} "
      f"{'% of genus':>10} {'LOSO top 20':>12} {'median':>7}")
for _, r in table.head(N_TOP).iterrows():
    g = r['silva_taxon'] if isinstance(r['genus'], str) else f"(no genus; {r['silva_taxon']})"
    pct = f"{r['pct_of_genus_reads']:>9.1f}%" if isinstance(r['genus'], str) else f"{'-':>10}"
    print(f"{r['rank']:>4.0f}  {r['ASV_ID']:<9} {g:<38} {r['mean_abs_shap']:>10.3f} {r['n_samples']:>8} "
          f"{r['n_studies']:>8} {pct} {r['loso_folds_in_top20']:>7} of {len(studies)} "
          f"{r['loso_median_rank']:>7.1f}")
print("-" * 126)
print(f"ASVs with any attribution: {int((table['mean_abs_shap'] > 0).sum())} of {len(table)}.")

print(f"\n--- {WHOLE_GENUS.upper()}: THE {n_genus} ASVs THAT TOGETHER HOLD {100 * GENUS_READ_SHARE:.0f}% OF ITS READS "
      f"({len(g_reads)} ASVs with reads, {int(g_reads.sum())} reads) ---")
print(f"{'ASV':<9} {'% of genus':>10} {'cumulative':>10} {'SHAP rank':>10} {'studies':>8}")
for a in genus_set:
    print(f"{a:<9} {100 * g_reads[a] / g_reads.sum():>9.1f}% {100 * cum[a]:>9.1f}% "
          f"{shap_rank_text(a):>10} {every.loc[all_counts[a] > 0, 'Paper_ID'].nunique():>8}")
print("SHAP rank 'none': the ASV gets no attribution in the fit; 'not kept': removed by the prevalence filter.")
print(f"In the top {N_TOP} as well: {sum(a in top for a in genus_set)}. The FASTA holds {len(chosen)} ASVs.")
save_table(blast, "12_table_s2_asvs.csv")

out_dir = os.environ.get("BIOTWIN_RESULTS_DIR")
fasta_name = f"Table_S2_ASVs_{datetime.date.today():%Y-%m-%d}.fasta"
if out_dir:
    fasta_path = os.path.join(out_dir, fasta_name)
    if os.path.exists(fasta_path):
        raise SystemExit(f"{fasta_path} already exists.")
    with open(fasta_path, 'w', encoding='utf-8', newline='\n') as f:
        for _, r in blast.iterrows():
            # One word per title, as in the FASTA BLASTed on 2026-10-04, e.g.
            # ASV_143_rank_4_Caproiciproducens or ASV_9_Caproiciproducens_6.5pct_of_genus_reads
            taxon = re.sub(r'[^A-Za-z0-9.-]+', '_', r['silva_taxon']).strip('_')
            title = (f"{r['ASV_ID']}_rank_{r['rank']:g}_{taxon}" if r['ASV_ID'] in top
                     else f"{r['ASV_ID']}_{taxon}_{r['pct_of_genus_reads']:.1f}pct_of_genus_reads")
            f.write(f">{title}\n{r['sequence']}\n")
    print(f"  Saved: {fasta_path}")
    print("\nNext (by hand): BLAST this FASTA at blast.ncbi.nlm.nih.gov (Nucleotide BLAST, database "
          "'16S ribosomal RNA sequences (Bacteria and Archaea)', program megablast), then Download All -> Text, "
          "and run 12b_table_s2_from_blast.py on that file.")
else:
    print(f"  (BIOTWIN_RESULTS_DIR not set: {fasta_name} not saved)")
